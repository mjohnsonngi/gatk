package org.broadinstitute.hellbender.tools.walkers.vqsr.scalable;

import com.fasterxml.jackson.annotation.JsonProperty;
import com.fasterxml.jackson.databind.ObjectMapper;
import com.google.common.primitives.Doubles;
import org.apache.commons.lang3.tuple.Pair;
import org.apache.commons.math3.stat.descriptive.moment.Mean;
import org.apache.commons.math3.stat.descriptive.moment.Variance;
import org.apache.commons.math3.stat.descriptive.rank.Median;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hdf5.HDF5File;
import org.broadinstitute.hdf5.HDF5LibException;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.clustering.BayesianGaussianMixtureModeller;
import org.broadinstitute.hellbender.utils.io.IOUtils;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.PrintStream;
import java.io.Serializable;
import java.util.List;
import java.util.stream.IntStream;

final class VariantAnnotationUtils {

    private static final Logger logger = LogManager.getLogger(VariantAnnotationUtils.class);

    static final class Scorer implements Serializable {
        private static final long serialVersionUID = 1L;

        private final Preprocesser preprocesser;
        private final BayesianGaussianMixtureModeller bgmm;

        Scorer(final Preprocesser preprocesser,
               final BayesianGaussianMixtureModeller bgmm) {
            this.preprocesser = preprocesser;
            this.bgmm = bgmm;
        }

        Pair<double[][], double[]> preprocessAndScoreSamples(final double[][] data) {
            final double[][] preprocessedData = preprocesser.transform(data);
            final double[] scores = bgmm.scoreSamples(preprocessedData);
            return Pair.of(preprocessedData, scores);
        }
    }

    /**
     * TODO standardize and median impute; maybe add Z-score truncation?
     */
    static final class Preprocesser implements Serializable {
        private static final long serialVersionUID = 1L;

        private double[] meansForStandardization;
        private double[] standardDeviationsForStandardization;
        private double[] standardizedMediansForImputation;

        Preprocesser() {}

        private static double standardize(final double value,
                                          final double mean,
                                          final double standardDeviation) {
            // if standard deviation is zero, no standardization is needed, but we must guard against divide-by-zero
            return (value - mean) / (standardDeviation < Double.MIN_VALUE ? 1. : standardDeviation);
        }

        double[][] fitTransform(final double[][] data) {
            // TODO validation
            final int nSamples = data.length;
            final int nFeatures = data[0].length;
            meansForStandardization = IntStream.range(0, nFeatures)
                    .mapToDouble(j -> new Mean().evaluate(
                            IntStream.range(0, nSamples).boxed()
                                    .mapToDouble(i -> data[i][j])
                                    .filter(Double::isFinite)
                                    .toArray()))
                    .toArray();
            standardDeviationsForStandardization = IntStream.range(0, nFeatures)
                    .mapToDouble(j -> new Variance().evaluate(
                            IntStream.range(0, nSamples).boxed()
                                    .mapToDouble(i -> data[i][j])
                                    .filter(Double::isFinite)
                                    .toArray()))
                    .map(Math::sqrt)
                    .toArray();

            standardizedMediansForImputation = IntStream.range(0, nFeatures)
                    .mapToDouble(j -> new Median().evaluate(
                            IntStream.range(0, nSamples).boxed()
                                    .mapToDouble(i -> data[i][j])
                                    .filter(Double::isFinite)
                                    .map(x -> standardize(x, meansForStandardization[j], standardDeviationsForStandardization[j]))
                                    .toArray()))
                    .toArray();

            return transform(data);
        }

        double[][] transform(final double[][] data) {
            // TODO validation
            final int nSamples = data.length;
            final int nFeatures = data[0].length;
            double[][] preprocessedData = new double[nSamples][nFeatures];
            for (int i = 0; i < nSamples; i++) {
                for (int j = 0; j < nFeatures; j++) {
                    final double value = data[i][j];
                    final double preprocessedValue = Double.isNaN(value)
                            ? standardizedMediansForImputation[j]
                            : standardize(value, meansForStandardization[j], standardDeviationsForStandardization[j]);
                    if (!Double.isFinite(preprocessedValue)) {
                        throw new GATKException("Encountered non-finite value during preprocessing of annotations.");
                    }
                    preprocessedData[i][j] = preprocessedValue;
                }
            }
            return preprocessedData;
        }
    }

    // TODO perhaps just put this in the BGMM classes?
    // there is the complication of specifying mean and covariance priors differently here,
    // as well as the possible desire to have different defaults;
    // would also require slightly more code to invoke a BGMM if we use a builder for these instead
    static final class Hyperparameters {
        @JsonProperty("n_components")
        int nComponents = 6;
        @JsonProperty("tol")
        double tol = 1.;
        @JsonProperty("reg_covar")
        double regCovar = 1E-6;
        @JsonProperty("max_iter")
        int maxIter = 150;
        @JsonProperty("n_init")
        int nInit = 1;
        @JsonProperty("init_params")
        String initMethod = "kmeans++";
        @JsonProperty("weight_concentration_prior")
        Double weightConcentrationPrior = null;
        @JsonProperty("mean_precision_prior")
        double meanPrecisionPrior = 1.;
        @JsonProperty("mean_prior")
        Double meanPrior = null;
        @JsonProperty("degrees_of_freedom_prior")
        Double degreesOfFreedomPrior = null;
        @JsonProperty("covariance_prior")
        Double covariancePrior = null;
        @JsonProperty("random_state")
        int seed = 0;
        @JsonProperty("warm_start")
        boolean warmStart = false;
        @JsonProperty("verbose_interval")
        int verboseInterval = 5;

        Hyperparameters() {}

        static Hyperparameters readHyperparameters(final File hyperparametersJSONFile) {
            IOUtils.canReadFile(hyperparametersJSONFile);
            try {
                final ObjectMapper mapper = new ObjectMapper();
                 return mapper.readValue(hyperparametersJSONFile, Hyperparameters.class);
            } catch (final Exception e) {
                throw new UserException.BadInput("Could not read hyperparameters JSON.", e);
            }
        }

        BayesianGaussianMixtureModeller.InitMethod getInitMethod() {
            switch (initMethod) {
                case "kmeans++":
                    return BayesianGaussianMixtureModeller.InitMethod.K_MEANS_PLUS_PLUS;
                case "random":
                    return BayesianGaussianMixtureModeller.InitMethod.RANDOM;
                default:
                    throw new UserException.BadInput("Unknown initialization method specified.");
            }
        }

        @Override
        public String toString() {
            return "Hyperparameters{" +
                    "nComponents=" + nComponents +
                    ", tol=" + tol +
                    ", regCovar=" + regCovar +
                    ", maxIter=" + maxIter +
                    ", nInit=" + nInit +
                    ", initMethod='" + initMethod + '\'' +
                    ", weightConcentrationPrior=" + weightConcentrationPrior +
                    ", meanPrecisionPrior=" + meanPrecisionPrior +
                    ", meanPrior=" + meanPrior +
                    ", degreesOfFreedomPrior=" + degreesOfFreedomPrior +
                    ", covariancePrior=" + covariancePrior +
                    ", seed=" + seed +
                    ", warmStart=" + warmStart +
                    ", verboseInterval=" + verboseInterval +
                    '}';
        }
    }

    static <T> void serialize(final File outputFile, 
                              final T object) {
        try (final FileOutputStream fileOutputStream = new FileOutputStream(outputFile);
             final ObjectOutputStream objectOutputStream = new ObjectOutputStream(fileOutputStream)) {
            objectOutputStream.writeObject(object);
        } catch (final IOException e) {
            throw new GATKException(String.format("Exception encountered during serialization of %s to %s: %s",
                    object.getClass(), outputFile.getAbsolutePath(), e));
        }
    }

    static <T> T deserialize(final File inputFile,
                             final Class<T> clazz) {
        try (final FileInputStream fileInputStream = new FileInputStream(inputFile);
             final ObjectInputStream objectInputStream = new ObjectInputStream(fileInputStream)) {
            return clazz.cast(objectInputStream.readObject());
        } catch (final IOException | ClassNotFoundException e) {
            throw new GATKException(String.format("Exception encountered during deserialization from %s: %s",
                    inputFile.getAbsolutePath(), e));
        }
    }

    static double[] readScores(final File inputFile) {
        try (final HDF5File inputHDF5File = new HDF5File(inputFile, HDF5File.OpenMode.READ_ONLY)) {
            IOUtils.canReadFile(inputHDF5File.getFile());
            return inputHDF5File.readDoubleArray("/data/scores");
        } catch (final HDF5LibException exception) {
            throw new GATKException(String.format("Exception encountered during reading of scores from %s: %s",
                    inputFile.getAbsolutePath(), exception));
        }
    }
    
    static void writeScores(final File outputFile,
                            final double[] scores) {
        try (final HDF5File outputHDF5File = new HDF5File(outputFile, HDF5File.OpenMode.CREATE)) {
            outputHDF5File.makeDoubleArray("/data/scores", scores);
        } catch (final HDF5LibException exception) {
            throw new GATKException(String.format("Exception encountered during writing of scores (%s). Output file at %s may be in a bad state.",
                    exception, outputFile.getAbsolutePath()));
        }
        logger.info(String.format("Scores written to %s.", outputFile.getAbsolutePath()));
    }

    // TODO we could just get all this stuff from the VariantDataManager, but let's clean that up later
    static void writeTruthSensitivityTranches(final File outputFile,
                                              final File scoresFile,
                                              final File annotationsFile,
                                              final List<Double> truthSensitivityTranches,
                                              final VariantTypeMode mode) {
        // TODO validate
        final List<Double> scores = Doubles.asList(readScores(scoresFile));
        final List<Boolean> isTruth = VariantDataCollection.readLabel(annotationsFile, VariantDataCollection.TRUTH_LABEL);
        try {
            final PrintStream tranchesStream = new PrintStream(outputFile);

            // Find the score cutoff values which correspond to the various tranches of calls requested by the user
            final int nCallsAtTruth = TruthSensitivityTranche.countCallsAtTruth(scores, isTruth, Double.NEGATIVE_INFINITY);
            final TruthSensitivityTranche.TruthSensitivityMetric metric = new TruthSensitivityTranche.TruthSensitivityMetric(nCallsAtTruth);
            final List<TruthSensitivityTranche> tranches = TruthSensitivityTranche.findTranches(scores, isTruth, truthSensitivityTranches, metric, mode);
            tranchesStream.print(TruthSensitivityTranche.printHeader());
            tranchesStream.print(TruthSensitivityTranche.tranchesString(tranches));
        } catch (final FileNotFoundException e) {
            throw new UserException.CouldNotCreateOutputFile(outputFile, e);
        }
        logger.info(String.format("Truth-sensitivity tranches written to %s.", outputFile.getAbsolutePath()));
    }
}