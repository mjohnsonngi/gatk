package org.broadinstitute.hellbender.tools.walkers.sv;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.*;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.*;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.StructuralVariantDiscoveryProgramGroup;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.tools.sv.SVCallRecord;
import org.broadinstitute.hellbender.tools.sv.SVCallRecordUtils;
import org.broadinstitute.hellbender.tools.sv.cluster.CanonicalSVCollapser;
import org.broadinstitute.hellbender.tools.sv.cluster.CanonicalSVLinkage;
import org.broadinstitute.hellbender.tools.sv.cluster.SVClusterEngineArgumentsCollection;
import org.broadinstitute.hellbender.tools.sv.cluster.SVConcordanceClusterEngine;
import org.broadinstitute.hellbender.tools.walkers.validation.Concordance;
import org.broadinstitute.hellbender.tools.walkers.validation.ConcordanceState;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.variant.VariantContextGetters;
import picard.vcf.*;

import java.util.*;
import java.util.function.Function;
import java.util.stream.Collectors;

/**
 * <p>Clusters structural variants based on coordinates, event type, and supporting algorithms. Primary use cases include:</p>
 * <ul>
 *     <li>
 *         Clustering SVs produced by multiple callers, based on interval overlap, breakpoint proximity, and sample overlap.
 *     </li>
 *     <li>
 *         Merging multiple SV VCFs with disjoint sets of samples and/or variants.
 *     </li>
 *     <li>
 *         Defragmentation of copy number variants produced with depth-based callers.
 *     </li>
 * </ul>
 *
 * <p>The tool generates a new VCF with clusters collapsed into single representative records. By default, a MEMBERS
 * field is generated that lists the input variant IDs contained in that record's cluster.</p>
 *
 * <h3>Inputs</h3>
 *
 * <ul>
 *     <li>
 *         One or more SV VCFs
 *     </li>
 * </ul>
 *
 * <h3>Output</h3>
 *
 * <ul>
 *     <li>
 *         Clustered VCF
 *     </li>
 * </ul>
 *
 * <h3>Usage example</h3>
 *
 * <pre>
 *     gatk SVConcordance \
 *       -V variants1.vcf.gz \
 *       -V variants2.vcf.gz \
 *       --eval-vcf eval.vcf.gz \
 *       -O annotated.vcf.gz
 * </pre>
 *
 * @author Mark Walker &lt;markw@broadinstitute.org&gt;
 */
@CommandLineProgramProperties(
        summary = "Calculates structural variant genotype concordance with a reference call set",
        oneLineSummary = "Calculates structural variant genotype concordance with a reference call set",
        programGroup = StructuralVariantDiscoveryProgramGroup.class
)
@BetaFeature
@DocumentedFeature
public final class SVConcordance extends MultiVariantWalker {

    public static final String EVAL_VCF_LONG_NAME = "eval-vcf";

    @Argument(
            doc = "Output VCF",
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME
    )
    private GATKPath outputFile;

    @Argument(
            doc = "Evaluation VCF to annotate using reference genotypes provided with one or more -"
                    + StandardArgumentDefinitions.VARIANT_SHORT_NAME + " arguments",
            fullName = EVAL_VCF_LONG_NAME
    )
    private GATKPath evalVcf;

    @ArgumentCollection
    private final SVClusterEngineArgumentsCollection clusterParameterArgs = new SVClusterEngineArgumentsCollection();

    private SAMSequenceDictionary dictionary;
    private VariantContextWriter writer;
    private CanonicalSVLinkage<SVCallRecord> linkage;
    private Iterator<VariantContext> evalVariantSource;
    private VariantContext currentEvalVariant;
    private Comparator<VariantContext> variantComparator;
    private SVConcordanceClusterEngine engine;
    private Long nextItemId = 0L;
    private String currentContig = null;
    private GenotypeConcordanceScheme scheme;

    @Override
    public void onTraversalStart() {
        dictionary = getBestAvailableSequenceDictionary();
        if (dictionary == null) {
            throw new UserException("Reference sequence dictionary required");
        }

        linkage = new CanonicalSVLinkage<>(dictionary, false);
        linkage.setDepthOnlyParams(clusterParameterArgs.getDepthParameters());
        linkage.setMixedParams(clusterParameterArgs.getMixedParameters());
        linkage.setEvidenceParams(clusterParameterArgs.getPESRParameters());
        engine = new SVConcordanceClusterEngine(linkage, this::collapse, dictionary);

        variantComparator = IntervalUtils.getDictionaryOrderComparator(dictionary);
        final FeatureDataSource<VariantContext> source = new FeatureDataSource<>(evalVcf.toString());
        final Object headerObject = source.getHeader();
        if ( ! (headerObject instanceof VCFHeader) ) {
            throw new GATKException("Header for " + source.getName() + " is not in VCF header format");
        }
        final VCFHeader header = (VCFHeader) headerObject;
        evalVariantSource = new FeatureDataSource<VariantContext>(evalVcf.toString()).iterator();
        scheme = new GA4GHSchemeWithMissingAsHomRef();

        writer = createVCFWriter(outputFile);
        writer.writeHeader(createHeader(header));
    }

    @Override
    public Object onTraversalSuccess() {
        flushEvalVariants(null);
        flushClusters(true);
        return super.onTraversalSuccess();
    }

    @Override
    public void closeTool() {
        super.closeTool();
        if (writer != null) {
            writer.close();
        }
    }

    @Override
    public void apply(final VariantContext variant, final ReadsContext readsContext,
                      final ReferenceContext referenceContext, final FeatureContext featureContext) {
        flushEvalVariants(variant);
        if (currentContig == null) {
            currentContig = variant.getContig();
        } else if (!currentContig.equals(variant.getContig())) {
            currentContig = variant.getContig();
            flushClusters(true);
        }
        engine.add(nextItemId++, SVCallRecordUtils.create(variant), true);
        flushClusters(false);
    }

    private void flushEvalVariants(final VariantContext variant) {
        while (evalVariantSource.hasNext()) {
            if (currentEvalVariant == null) {
                currentEvalVariant = evalVariantSource.next();
            }
            if (variantComparator.compare(currentEvalVariant, variant) <= 0) {
                if (!currentEvalVariant.getContig().equals(currentContig)) {
                    if (currentContig != null) {
                        flushClusters(true);
                    }
                    currentContig = currentEvalVariant.getContig();
                }
                engine.add(nextItemId++, SVCallRecordUtils.create(currentEvalVariant), false);
                currentEvalVariant = null;
            } else {
                break;
            }
        }
    }

    private void flushClusters(final boolean force) {
        engine.flush(force).stream()
                .map(SVCallRecordUtils::getVariantBuilder)
                .map(VariantContextBuilder::make)
                .forEach(writer::add);
    }

    private VCFHeader createHeader(final VCFHeader header) {
        header.addMetaDataLine(new VCFFormatHeaderLine(GenotypeConcordance.CONTINGENCY_STATE_TAG, 1, VCFHeaderLineType.String, "The genotype concordance contingency state"));
        header.addMetaDataLine(Concordance.TRUTH_STATUS_HEADER_LINE);
        header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.GENOTYPE_CONCORDANCE_INFO, 1, VCFHeaderLineType.Float, "Genotype concordance"));
        header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.NON_REF_GENOTYPE_CONCORDANCE_INFO, 1, VCFHeaderLineType.Float, "Non-ref genotype concordance"));
        header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.HET_PPV_INFO, 1, VCFHeaderLineType.Float, "Heterozygous genotype positive predictive value"));
        header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.HET_SENSITIVITY_INFO, 1, VCFHeaderLineType.Float, "Heterozygous genotype sensitivity"));
        header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.HET_SPECIFICITY_INFO, 1, VCFHeaderLineType.Float, "Heterozygous genotype specificity"));
        header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.HOMVAR_PPV_INFO, 1, VCFHeaderLineType.Float, "Homozygous genotype positive predictive value"));
        header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.HOMVAR_SENSITIVITY_INFO, 1, VCFHeaderLineType.Float, "Homozygous genotype sensitivity"));
        header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.HOMVAR_SPECIFICITY_INFO, 1, VCFHeaderLineType.Float, "Homozygous genotype specificity"));
        header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.VAR_PPV_INFO, 1, VCFHeaderLineType.Float, "Non-ref genotype positive predictive value"));
        header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.VAR_SENSITIVITY_INFO, 1, VCFHeaderLineType.Float, "Non-ref genotype sensitivity"));
        header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.VAR_SPECIFICITY_INFO, 1, VCFHeaderLineType.Float, "Non-ref genotype specificity"));
        header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.CONCORDANT_MEMBERS_INFO, 1, VCFHeaderLineType.Float, "Matching truth set variant ids"));
        return header;
    }

    public SVCallRecord collapse(final SVConcordanceClusterEngine.CrossRefOutputCluster cluster) {
        final SVCallRecord evalRecord = cluster.getEvalItem();
        final Map<String, List<Genotype>> truthGenotypesMap = cluster.getItems().stream()
                .map(SVCallRecord::getGenotypes)
                .flatMap(Collection::stream)
                .collect(Collectors.groupingBy(Genotype::getSampleName));
        final GenotypesContext evalGenotypes = evalRecord.getGenotypes();

        final ArrayList<Genotype> newGenotypes = new ArrayList<>(evalGenotypes.size());
        final GenotypeConcordanceCounts counts = new GenotypeConcordanceCounts();
        for (final Genotype g : evalGenotypes) {
            final GenotypeBuilder builder = new GenotypeBuilder(g);
            final String sample = g.getSampleName();
            final GenotypeConcordanceStates.TruthState truthState = getTruthState(truthGenotypesMap.get(sample), evalRecord);
            final GenotypeConcordanceStates.CallState evalState = getEvalState(g, evalRecord);
            final GenotypeConcordanceStates.TruthAndCallStates state = new GenotypeConcordanceStates.TruthAndCallStates(truthState, evalState);
            counts.increment(state);
            builder.attribute(GenotypeConcordance.CONTINGENCY_STATE_TAG, scheme.getContingencyStateString(truthState, evalState));
            newGenotypes.add(builder.make());
        }
        final SVCallRecord recordWithGenotypes = SVCallRecordUtils.copyCallWithNewGenotypes(evalRecord, GenotypesContext.create(newGenotypes));
        final Map<String, Object> attributes = new HashMap<>(recordWithGenotypes.getAttributes());
        final List<String> members = cluster.getItems().stream().map(SVCallRecord::getId).collect(Collectors.toList());
        final ConcordanceState variantStatus = cluster.getItems().isEmpty() ? ConcordanceState.FALSE_POSITIVE : ConcordanceState.TRUE_POSITIVE;
        final GenotypeConcordanceSummaryMetrics metrics = new GenotypeConcordanceSummaryMetrics(VariantContext.Type.SYMBOLIC, counts, "truth", "eval", true);

        attributes.put(Concordance.TRUTH_STATUS_VCF_ATTRIBUTE, variantStatus.getAbbreviation());
        attributes.put(GATKSVVCFConstants.GENOTYPE_CONCORDANCE_INFO, metrics.GENOTYPE_CONCORDANCE);
        attributes.put(GATKSVVCFConstants.NON_REF_GENOTYPE_CONCORDANCE_INFO, metrics.NON_REF_GENOTYPE_CONCORDANCE);
        attributes.put(GATKSVVCFConstants.HET_PPV_INFO, metrics.HET_PPV);
        attributes.put(GATKSVVCFConstants.HET_SENSITIVITY_INFO, metrics.HET_SENSITIVITY);
        attributes.put(GATKSVVCFConstants.HET_SPECIFICITY_INFO, metrics.HET_SPECIFICITY);
        attributes.put(GATKSVVCFConstants.HOMVAR_PPV_INFO, metrics.HOMVAR_PPV);
        attributes.put(GATKSVVCFConstants.HOMVAR_SENSITIVITY_INFO, metrics.HOMVAR_SENSITIVITY);
        attributes.put(GATKSVVCFConstants.HOMVAR_SPECIFICITY_INFO, metrics.HOMVAR_SPECIFICITY);
        attributes.put(GATKSVVCFConstants.VAR_PPV_INFO, metrics.VAR_PPV);
        attributes.put(GATKSVVCFConstants.VAR_SENSITIVITY_INFO, metrics.VAR_SENSITIVITY);
        attributes.put(GATKSVVCFConstants.VAR_SPECIFICITY_INFO, metrics.VAR_SPECIFICITY);
        attributes.put(GATKSVVCFConstants.CONCORDANT_MEMBERS_INFO, members);

        return SVCallRecordUtils.copyCallWithNewAttributes(recordWithGenotypes, attributes);
    }

    protected GenotypeConcordanceStates.CallState getEvalState(final Genotype genotype, final SVCallRecord record) {
        final StructuralVariantType svtype = record.getType();
        if (svtype == StructuralVariantType.CNV || svtype == StructuralVariantType.DUP || svtype == StructuralVariantType.DEL) {
            return getCNVEvalState(genotype, record);
        } else {
            return getNonCNVEvalState(genotype);
        }
    }

    private GenotypeConcordanceStates.CallState getNonCNVEvalState(final Genotype g) {
        if (g.isHomRef()) {
            return GenotypeConcordanceStates.CallState.HOM_REF;
        } else if (g.isHet()) {
            return GenotypeConcordanceStates.CallState.HET_REF_VAR1;
        } else if (g.isHomVar()) {
            return GenotypeConcordanceStates.CallState.HOM_VAR1;
        } else if (g.isNoCall() || g.getPloidy() == 0) {
            return GenotypeConcordanceStates.CallState.NO_CALL;
        } else if (g.isFiltered()) {
            return GenotypeConcordanceStates.CallState.GT_FILTERED;
        } else {
            throw new IllegalArgumentException("Unsupported eval set genotype: " + g.toBriefString());
        }
    }

    protected GenotypeConcordanceStates.TruthState getTruthState(final List<Genotype> genotypes, final SVCallRecord record) {
        if (genotypes == null || genotypes.isEmpty()) {
            return GenotypeConcordanceStates.TruthState.HOM_REF;
        }
        final StructuralVariantType svtype = record.getType();
        final boolean isCNV = svtype == StructuralVariantType.CNV || svtype == StructuralVariantType.DUP || svtype == StructuralVariantType.DEL;
        final Function<Genotype, GenotypeConcordanceStates.TruthState> mapper = isCNV ? x -> getCNVTruthState(x, record) : this::getNonCNVTruthState;
        final Set<GenotypeConcordanceStates.TruthState> states = genotypes.stream().map(mapper).collect(Collectors.toSet());
        final int numHetRefVar1 = states.contains(GenotypeConcordanceStates.TruthState.HET_REF_VAR1) ? 1 : 0;
        final int numHetVar1Var2 = states.contains(GenotypeConcordanceStates.TruthState.HET_VAR1_VAR2) ? 1 : 0;
        final int numHomVar = states.contains(GenotypeConcordanceStates.TruthState.HOM_VAR1) ? 1 : 0;
        if (numHetRefVar1 + numHetVar1Var2 + numHomVar > 1) {
            return GenotypeConcordanceStates.TruthState.IS_MIXED;
        } else if (numHetRefVar1 == 1) {
            return GenotypeConcordanceStates.TruthState.HET_REF_VAR1;
        } else if (numHetVar1Var2 == 1) {
            return GenotypeConcordanceStates.TruthState.HET_VAR1_VAR2;
        } else if (numHomVar == 1) {
            return GenotypeConcordanceStates.TruthState.HOM_VAR1;
        } else if (states.contains(GenotypeConcordanceStates.TruthState.HOM_REF)) {
            return GenotypeConcordanceStates.TruthState.HOM_REF;
        } else if (states.contains(GenotypeConcordanceStates.TruthState.NO_CALL)) {
            return GenotypeConcordanceStates.TruthState.NO_CALL;
        } else {
            throw new IllegalArgumentException("Unsupported TruthStates in set: " +
                    String.join(", ", states.stream().map(Object::toString).collect(Collectors.toList())));
        }
    }

    private GenotypeConcordanceStates.TruthState getNonCNVTruthState(final Genotype g) {
        if (g.isHomRef()) {
            return GenotypeConcordanceStates.TruthState.HOM_REF;
        } else if (g.isHet()) {
            return GenotypeConcordanceStates.TruthState.HET_REF_VAR1;
        } else if (g.isHomVar()) {
            return GenotypeConcordanceStates.TruthState.HOM_VAR1;
        } else if (g.isNoCall() || g.getPloidy() == 0) {
            return GenotypeConcordanceStates.TruthState.NO_CALL;
        } else if (g.isFiltered()) {
            return GenotypeConcordanceStates.TruthState.GT_FILTERED;
        } else {
            throw new IllegalArgumentException("Unsupported truth set genotype: " + g.toBriefString());
        }
    }

    private CNVAlleleCounts getCNVAlleleCounts(final Genotype g, final SVCallRecord record) {
        final int ploidy = VariantContextGetters.getAttributeAsInt(g, GATKSVVCFConstants.EXPECTED_COPY_NUMBER_FORMAT, g.getPloidy());
        Utils.validate(ploidy <= 2, "Ploidy can be at most 2 but found: " + ploidy);
        Utils.validateArg(g.hasExtendedAttribute(GATKSVVCFConstants.COPY_NUMBER_FORMAT), "CNV genotype has undefined copy number: " + g.toBriefString());
        final int copyNumber = VariantContextGetters.getAttributeAsInt(g, GATKSVVCFConstants.COPY_NUMBER_FORMAT, 0);
        final List<Allele> determinableAlleles = CanonicalSVCollapser.getCNVGenotypeAllelesFromCopyNumber(record.getAltAlleles(), record.getRefAllele(), ploidy, copyNumber);
        int numVarDel = 0;
        int numVarDup = 0;
        int numRef = 0;
        for (int i = 0; i < determinableAlleles.size(); i++) {
            final Allele a = determinableAlleles.get(i);
            if (a.isNonReference()) {
                if (a.equals(Allele.SV_SIMPLE_DEL)) {
                    numVarDel++;
                } else if (a.equals(Allele.SV_SIMPLE_DUP)) {
                    numVarDup++;
                } else if (!a.equals(Allele.NO_CALL)) {
                    throw new IllegalArgumentException("Unsupported CNV alt allele: " + a);
                }
            } else if (a.isReference()) {
                numRef++;
            }
        }
        return new CNVAlleleCounts(ploidy, numVarDel, numVarDup, numRef);
    }

    private GenotypeConcordanceStates.CallState getCNVEvalState(final Genotype g, final SVCallRecord record) {
        final CNVAlleleCounts counts = getCNVAlleleCounts(g, record);
        if (counts.numDel > 0) {
            if (counts.numDup > 0) {
                return GenotypeConcordanceStates.CallState.HET_VAR1_VAR2;
            } else if (counts.numRef > 0) {
                return GenotypeConcordanceStates.CallState.HET_REF_VAR1;
            } else {
                return GenotypeConcordanceStates.CallState.HOM_VAR1;
            }
        } else if (counts.numDup > 0) {
            if (counts.numRef > 0) {
                return GenotypeConcordanceStates.CallState.HET_REF_VAR1;
            } else {
                return GenotypeConcordanceStates.CallState.HOM_VAR1;
            }
        } else if (counts.numRef > 0) {
            return GenotypeConcordanceStates.CallState.HOM_REF;
        } else {
            return GenotypeConcordanceStates.CallState.NO_CALL;
        }
    }

    private GenotypeConcordanceStates.TruthState getCNVTruthState(final Genotype g, final SVCallRecord record) {
        final CNVAlleleCounts counts = getCNVAlleleCounts(g, record);
        if (counts.numDel > 0) {
            if (counts.numDup > 0) {
                return GenotypeConcordanceStates.TruthState.HET_VAR1_VAR2;
            } else if (counts.numRef > 0) {
                return GenotypeConcordanceStates.TruthState.HET_REF_VAR1;
            } else {
                return GenotypeConcordanceStates.TruthState.HOM_VAR1;
            }
        } else if (counts.numDup > 0) {
            if (counts.numRef > 0) {
                return GenotypeConcordanceStates.TruthState.HET_REF_VAR1;
            } else {
                return GenotypeConcordanceStates.TruthState.HOM_VAR1;
            }
        } else if (counts.numRef > 0) {
            return GenotypeConcordanceStates.TruthState.HOM_REF;
        } else {
            return GenotypeConcordanceStates.TruthState.NO_CALL;
        }
    }

    private static class CNVAlleleCounts {
        public final int ploidy;
        public final int numDel;
        public final int numDup;
        public final int numRef;

        public CNVAlleleCounts(final int ploidy, final int numDel, final int numDup, final int numRef) {
            this.ploidy = ploidy;
            this.numDel = numDel;
            this.numDup = numDup;
            this.numRef = numRef;
        }
    }

}
