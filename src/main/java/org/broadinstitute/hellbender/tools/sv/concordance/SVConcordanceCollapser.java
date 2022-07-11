package org.broadinstitute.hellbender.tools.sv.concordance;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.variant.variantcontext.*;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.tools.sv.SVCallRecord;
import org.broadinstitute.hellbender.tools.sv.SVCallRecordUtils;
import org.broadinstitute.hellbender.tools.sv.cluster.CanonicalSVCollapser;
import org.broadinstitute.hellbender.tools.walkers.validation.Concordance;
import org.broadinstitute.hellbender.tools.walkers.validation.ConcordanceState;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.variant.VariantContextGetters;
import picard.vcf.*;

import java.util.*;
import java.util.stream.Collectors;

public class SVConcordanceCollapser {

    private final GenotypeConcordanceScheme scheme;
    private final boolean strictCnvInterpretation;

    public SVConcordanceCollapser(final boolean strictCnvInterpretation) {
        this.scheme = new SVGenotypeConcordanceScheme();
        this.strictCnvInterpretation = strictCnvInterpretation;
    }

    public SVCallRecord collapse(final SVConcordanceClusterEngine.ConcordanceOutputCluster cluster) {
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
            final GenotypeState truthGenotypeState = collapseGenotypes(truthGenotypesMap.get(sample), evalRecord);
            final GenotypeConcordanceStates.TruthState truthState = toTruthState(truthGenotypeState);
            final GenotypeConcordanceStates.CallState evalState = toCallState(collapseGenotypes(Collections.singletonList(g), evalRecord), truthGenotypeState);
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

    /**
     *
     * @param state
     * @param truthState
     * @return
     */
    private GenotypeConcordanceStates.CallState toCallState(final GenotypeState state, final GenotypeState truthState) {
        if (state == GenotypeState.HOM_REF) {
            return GenotypeConcordanceStates.CallState.HOM_REF;
        } else if (state == GenotypeState.HET_REF_VAR1) {
            return GenotypeConcordanceStates.CallState.HET_REF_VAR1;
        } else if (state == GenotypeState.HET_REF_VAR2) {
            // Shift VAR2->VAR1 if needed
            if (truthState == GenotypeState.HET_REF_VAR2 || truthState == GenotypeState.HOM_VAR2) {
                return GenotypeConcordanceStates.CallState.HET_REF_VAR1;
            } else {
                return GenotypeConcordanceStates.CallState.HET_REF_VAR2;
            }
        } else if (state == GenotypeState.HOM_VAR1) {
            return GenotypeConcordanceStates.CallState.HOM_VAR1;
        } else if (state == GenotypeState.HOM_VAR2) {
            // Shift VAR2->VAR1 if needed
            if (truthState == GenotypeState.HET_REF_VAR2 || truthState == GenotypeState.HOM_VAR2) {
                return GenotypeConcordanceStates.CallState.HOM_VAR1;
            } else {
                return GenotypeConcordanceStates.CallState.HOM_VAR2;
            }
        } else if (state == GenotypeState.HET_VAR1_VAR2) {
            return GenotypeConcordanceStates.CallState.HET_VAR1_VAR2;
        } else if (state == GenotypeState.UNKNOWN) {
            return GenotypeConcordanceStates.CallState.IS_MIXED;
        } else if (state == GenotypeState.NO_CALL) {
            return GenotypeConcordanceStates.CallState.NO_CALL;
        } else {
            throw new IllegalArgumentException("Unsupported genotype state for eval set: " + state.name());
        }
    }

    private GenotypeConcordanceStates.TruthState toTruthState(final GenotypeState state) {
        if (state == GenotypeState.HOM_REF) {
            return GenotypeConcordanceStates.TruthState.HOM_REF;
        } else if (state == GenotypeState.HET_REF_VAR1) {
            return GenotypeConcordanceStates.TruthState.HET_REF_VAR1;
        } else if (state == GenotypeState.HET_REF_VAR2) {
            return GenotypeConcordanceStates.TruthState.HET_REF_VAR1;
        } else if (state == GenotypeState.HOM_VAR1) {
            return GenotypeConcordanceStates.TruthState.HOM_VAR1;
        } else if (state == GenotypeState.HOM_VAR2) {
            return GenotypeConcordanceStates.TruthState.HOM_VAR1;
        } else if (state == GenotypeState.HET_VAR1_VAR2) {
            return GenotypeConcordanceStates.TruthState.HET_VAR1_VAR2;
        } else if (state == GenotypeState.UNKNOWN) {
            return GenotypeConcordanceStates.TruthState.IS_MIXED;
        } else if (state == GenotypeState.NO_CALL) {
            return GenotypeConcordanceStates.TruthState.NO_CALL;
        } else {
            throw new IllegalArgumentException("Unsupported genotype state for truth set: " + state.name());
        }
    }

    @VisibleForTesting
    protected GenotypeState collapseGenotypes(final List<Genotype> genotypes, final SVCallRecord record) {
        if (genotypes == null || genotypes.isEmpty()) {
            // Note Picard GenotypeConcordanceSummaryMetrics does not calculate concordance for MISSING states
            // correctly when using missingSitesFlag is true
            return GenotypeState.HOM_REF;
        }
        final StructuralVariantType svtype = record.getType();
        final boolean isCNV = svtype == StructuralVariantType.CNV || svtype == StructuralVariantType.DUP || svtype == StructuralVariantType.DEL;
        final Set<GenotypeState> states = genotypes.stream().map(g -> isCNV ? getCnvGenotypeState(g, record) : getNonCNVGenotypeState(g)).collect(Collectors.toSet());
        final GenotypeState[] stateValues = GenotypeState.values();
        final boolean[] stateMap = new boolean[GenotypeState.values().length];
        for (int i = 0; i < stateMap.length; i++) {
            if (states.contains(stateValues[i])) {
                stateMap[i] = true;
            }
        }

        final boolean hasVar1 = stateMap[GenotypeState.HOM_VAR1.ordinal()] || stateMap[GenotypeState.HET_VAR1_VAR2.ordinal()] || stateMap[GenotypeState.HET_REF_VAR1.ordinal()];
        final boolean hasVar2 = stateMap[GenotypeState.HOM_VAR2.ordinal()] || stateMap[GenotypeState.HET_VAR1_VAR2.ordinal()] || stateMap[GenotypeState.HET_REF_VAR2.ordinal()];

        if (hasVar1 && hasVar2) {
            return GenotypeState.HET_VAR1_VAR2;
        } else if (hasVar1) {
            if (stateMap[GenotypeState.HET_REF_VAR1.ordinal()]) {
                return GenotypeState.HET_REF_VAR1;
            } else if (stateMap[GenotypeState.HOM_VAR1.ordinal()]) {
                return GenotypeState.HOM_VAR1;
            } else {
                throw new GATKException.ShouldNeverReachHereException("Should be het or hom var1");
            }
        } else if (hasVar2) {
            if (stateMap[GenotypeState.HET_REF_VAR2.ordinal()]) {
                return GenotypeState.HET_REF_VAR2;
            } else if (stateMap[GenotypeState.HOM_VAR2.ordinal()]) {
                return GenotypeState.HOM_VAR2;
            } else {
                throw new GATKException.ShouldNeverReachHereException("Should be het or hom var1");
            }
        } else if (stateMap[GenotypeState.UNKNOWN.ordinal()]) {
            return GenotypeState.UNKNOWN;
        } else if (stateMap[GenotypeState.HOM_REF.ordinal()]) {
                return GenotypeState.HOM_REF;
        } else if (stateMap[GenotypeState.NO_CALL.ordinal()]) {
            return GenotypeState.NO_CALL;
        } else {
            throw new GATKException.ShouldNeverReachHereException("Empty state map");
        }
    }

    @VisibleForTesting
    protected GenotypeState getNonCNVGenotypeState(final Genotype g) {
        if (g.isHomRef()) {
            return GenotypeState.HOM_REF;
        } else if (g.isHet()) {
            return GenotypeState.HET_REF_VAR1;
        } else if (g.isHomVar()) {
            return GenotypeState.HOM_VAR1;
        } else if (g.isNoCall() || g.getPloidy() == 0) {
            return GenotypeState.NO_CALL;
        } else {
            // Something complicated going on, try filtering out missing/null alleles
            final List<Allele> filteredAlleles = new ArrayList<>(g.getAlleles().size());
            boolean foundNull = false;
            for (final Allele a : g.getAlleles()) {
                if (a != null && a != Allele.NO_CALL) {
                    filteredAlleles.add(a);
                } else {
                    foundNull = true;
                }
            }
            if (foundNull) {
                return getNonCNVGenotypeState(new GenotypeBuilder(g).alleles(filteredAlleles).make());
            } else {
                throw new IllegalArgumentException("Could not determine truth state for genotype: " + g);
            }
        }
    }

    @VisibleForTesting
    protected GenotypeState getCnvGenotypeState(final Genotype g, final SVCallRecord record) {
        Utils.validateArg(g.getExtendedAttribute(GATKSVVCFConstants.EXPECTED_COPY_NUMBER_FORMAT) != null,
                "Encountered missing " + GATKSVVCFConstants.EXPECTED_COPY_NUMBER_FORMAT + " for genotype " +
                        g + " in record " + record.getId());
        Utils.validateArg(g.getExtendedAttribute(GATKSVVCFConstants.COPY_NUMBER_FORMAT) != null,
                "Encountered missing " + GATKSVVCFConstants.COPY_NUMBER_FORMAT + " for genotype " +
                        g + " in record " + record.getId());
        final int ploidy = VariantContextGetters.getAttributeAsInt(g, GATKSVVCFConstants.EXPECTED_COPY_NUMBER_FORMAT, 0);
        final int copyNumber = VariantContextGetters.getAttributeAsInt(g, GATKSVVCFConstants.COPY_NUMBER_FORMAT, 0);

        // Short circuit the trivial case
        if (ploidy == 0) {
            return GenotypeState.NO_CALL;
        }
        // Special cases where the alleles cannot technically be unambiguously determined, but in practice we
        // usually want to call them something sensible. Note we know ploidy > 0 due to the check above.
        if (!strictCnvInterpretation) {
            if (copyNumber == ploidy) {
                return GenotypeState.HOM_REF;
            } else if (copyNumber < ploidy) {
                if (copyNumber == 0) {
                    return GenotypeState.HOM_VAR1;
                } else {
                    return GenotypeState.HET_REF_VAR1;
                }
            } else {
                // copyNumber > ploidy case
                if (ploidy == 1) {
                    // Haploid case
                    return GenotypeState.HOM_VAR2;
                } else if (copyNumber - ploidy == 1) {
                    // Diploid case, 1 copy
                    return GenotypeState.HET_REF_VAR2;
                } else {
                    // Diploid case, 2 or more copies
                    return GenotypeState.HOM_VAR2;
                }
            }
        }
        final List<Allele> alleles = CanonicalSVCollapser.getCNVGenotypeAllelesFromCopyNumber(record.getAltAlleles(), record.getRefAllele(), ploidy, copyNumber);
        int numDup = 0;
        int numDel = 0;
        int numRef = 0;
        int numNoCall = 0;
        for (final Allele a : alleles) {
            if (a.isReference()) {
                numRef++;
            } else if (a == Allele.SV_SIMPLE_DEL) {
                numDel++;
            } else if (a == Allele.SV_SIMPLE_DUP) {
                numDup++;
            } else if (a == Allele.NO_CALL) {
                numNoCall++;
            } else {
                throw new IllegalArgumentException("Unexpected allele " + a + " in genotype " + g + " in record " + record.getId());
            }
        }
        if (numDel > 0) {
            if (numDup > 0) {
                return GenotypeState.HET_VAR1_VAR2;
            } else if (numRef > 0) {
                return GenotypeState.HET_REF_VAR1;
            } else {
                return GenotypeState.HOM_VAR1;
            }
        } else if (numDup > 0) {
            if (numRef > 0) {
                return GenotypeState.HET_REF_VAR2;
            } else {
                return GenotypeState.HOM_VAR2;
            }
        } else if (numRef > 0) {
            return GenotypeState.HOM_REF;
        } else if (numNoCall > 0) {
            return GenotypeState.UNKNOWN;
        } else {
            throw new IllegalArgumentException("Unable to determine CNV call type from alleles: [" +
                    String.join(", ", alleles.stream().map(Allele::getDisplayString).collect(Collectors.toList())) + "]");
        }
    }

    enum GenotypeState {
        HOM_REF,
        HOM_VAR1,
        HOM_VAR2,
        HET_REF_VAR1,
        HET_REF_VAR2,
        HET_VAR1_VAR2,
        UNKNOWN,
        NO_CALL
    }

    /**
     * Based on {@link GA4GHSchemeWithMissingAsHomRef} but using VAR1/VAR2 to mean DEL/DUP for multiallelics.
     * Unused rows have been removed for simplicity.
     */
    private class SVGenotypeConcordanceScheme extends GenotypeConcordanceScheme {

        @Override
        protected void initiateScheme() {
            /** ROW STATE                                               MISSING       HOM_REF       HET_REF_VAR1       HET_VAR1_VAR2        HOM_VAR1        NO_CALL        LOW_GQ        LOW_DP        VC_FILTERED   GT_FILTERED   IS_MIXED    **/
            addRow(GenotypeConcordanceStates.CallState.MISSING,         TN_ONLY,      TN_ONLY,      TN_FN,             FN_ONLY,             FN_ONLY,        EMPTY,         NA,        NA,        NA,        NA,        NA);
            addRow(GenotypeConcordanceStates.CallState.HOM_REF,         TN_ONLY,      TN_ONLY,      TN_FN,             FN_ONLY,             FN_ONLY,        EMPTY,         NA,        NA,        NA,        NA,        NA);
            addRow(GenotypeConcordanceStates.CallState.HET_REF_VAR1,    FP_TN,        FP_TN,        TP_TN,             TP_FN,               TP_FN,          EMPTY,         NA,        NA,        NA,        NA,        NA);
            addRow(GenotypeConcordanceStates.CallState.HET_REF_VAR2,    FP_TN,        FP_TN,        FP_TN_FN,          TP_FP_FN,            FP_FN,          EMPTY,         NA,        NA,        NA,        NA,        NA);
            addRow(GenotypeConcordanceStates.CallState.HET_VAR1_VAR2,   FP_ONLY,      FP_ONLY,      TP_FP,             TP_ONLY,             TP_FP_FN,       EMPTY,         NA,        NA,        NA,        NA,        NA);
            addRow(GenotypeConcordanceStates.CallState.HOM_VAR1,        FP_ONLY,      FP_ONLY,      TP_FP,             TP_FN,               TP_ONLY,        EMPTY,         NA,        NA,        NA,        NA,        NA);
            addRow(GenotypeConcordanceStates.CallState.HOM_VAR2,        FP_ONLY,      FP_ONLY,      FP_FN,             TP_FN,               FP_FN,          EMPTY,         NA,        NA,        NA,        NA,        NA);
            addRow(GenotypeConcordanceStates.CallState.NO_CALL,         EMPTY,        EMPTY,        EMPTY,             EMPTY,               EMPTY,          EMPTY,         NA,        NA,        NA,        NA,        NA);
            addRow(GenotypeConcordanceStates.CallState.IS_MIXED,        EMPTY,        EMPTY,        EMPTY,             EMPTY,               EMPTY,          EMPTY,         NA,        NA,        NA,        NA,        NA);
        }
    }
}
