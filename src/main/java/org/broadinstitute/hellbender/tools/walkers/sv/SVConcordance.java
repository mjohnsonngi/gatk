package org.broadinstitute.hellbender.tools.walkers.sv;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFFormatHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
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
import org.broadinstitute.hellbender.tools.sv.cluster.CanonicalSVLinkage;
import org.broadinstitute.hellbender.tools.sv.cluster.SVClusterEngineArgumentsCollection;
import org.broadinstitute.hellbender.tools.sv.concordance.SVConcordanceClusterEngine;
import org.broadinstitute.hellbender.tools.sv.concordance.SVConcordanceCollapser;
import org.broadinstitute.hellbender.tools.walkers.validation.Concordance;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import picard.vcf.GenotypeConcordance;

import java.util.Comparator;
import java.util.Iterator;

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
    public static final String STRICT_CNV_INTERPRETATION_LONG_NAME = "strict-cnv-interpretation";

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

    /**
     * If used, CNVs diploid genotypes with copy number equal to the ploidy {@link GATKSVVCFConstants#COPY_NUMBER_FORMAT} ==
     * {@link GATKSVVCFConstants#EXPECTED_COPY_NUMBER_FORMAT} will be interpreted as a "mixed" call, resulting in
     * empty concordance status for the genotype. By default, these are treated as homozygous reference if the ploidy
     * is greater than 0.
     */
    @Argument(
            doc = "Enable strict CNV genotype copy state interpretation. See documentation for details.",
            fullName = STRICT_CNV_INTERPRETATION_LONG_NAME,
            optional = true
    )
    private boolean strictCnvInterpretation = false;

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

        final SVConcordanceCollapser collapser = new SVConcordanceCollapser(strictCnvInterpretation);
        engine = new SVConcordanceClusterEngine(linkage, collapser::collapse, dictionary);

        variantComparator = IntervalUtils.getDictionaryOrderComparator(dictionary);
        final FeatureDataSource<VariantContext> source = new FeatureDataSource<>(evalVcf.toString());
        final Object headerObject = source.getHeader();
        if ( ! (headerObject instanceof VCFHeader) ) {
            throw new GATKException("Header for " + source.getName() + " is not in VCF header format");
        }
        final VCFHeader header = (VCFHeader) headerObject;
        evalVariantSource = new FeatureDataSource<VariantContext>(evalVcf.toString()).iterator();

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


}
