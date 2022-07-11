package org.broadinstitute.hellbender.tools.sv.concordance;

import com.google.common.collect.Lists;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.StructuralVariantType;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.tools.sv.SVCallRecord;
import org.broadinstitute.hellbender.tools.sv.SVTestUtils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.Arrays;

public class SVConcordanceCollapserTest {

    @Test
    public void testCollapse() {
    }

    @DataProvider(name = "testGetTruthStateNonCNVData")
    public Object[][] testGetTruthStateNonCNVData() {
        return new Object[][]{

                ////////////////////////////////////////
                // Empty edge cases
                ////////////////////////////////////////

                // No alleles case
                {
                        new Allele[]{ },
                        new Allele[]{ },
                        SVConcordanceCollapser.GenotypeState.NO_CALL
                },
                // no-call allele
                {
                        new Allele[]{ Allele.NO_CALL },
                        new Allele[]{ },
                        SVConcordanceCollapser.GenotypeState.NO_CALL
                },
                // mixed no-call / ref allele
                {
                        new Allele[]{ Allele.NO_CALL, Allele.REF_N },
                        new Allele[]{ },
                        SVConcordanceCollapser.GenotypeState.HOM_REF
                },
                // One genotype with mixed no-call / var allele
                {
                        new Allele[]{ Allele.NO_CALL, Allele.SV_SIMPLE_INS },
                        new Allele[]{ },
                        SVConcordanceCollapser.GenotypeState.HOM_VAR1
                },

                ////////////////////////////////////////
                // Haploid
                ////////////////////////////////////////

                // hom ref
                {
                        new Allele[]{ Allele.REF_N },
                        new Allele[]{ Allele.REF_N, Allele.SV_SIMPLE_INS },
                        SVConcordanceCollapser.GenotypeState.HOM_REF
                },
                // hom var
                {
                        new Allele[]{ Allele.SV_SIMPLE_INS },
                        new Allele[]{ Allele.REF_N, Allele.SV_SIMPLE_INS },
                        SVConcordanceCollapser.GenotypeState.HOM_VAR1
                },

                ////////////////////////////////////////
                // Diploid
                ////////////////////////////////////////

                // hom ref
                {
                        new Allele[]{ Allele.REF_N, Allele.REF_N },
                        new Allele[]{ Allele.REF_N, Allele.SV_SIMPLE_INS },
                        SVConcordanceCollapser.GenotypeState.HOM_REF
                },
                // hom var
                {
                        new Allele[]{ Allele.SV_SIMPLE_INS, Allele.SV_SIMPLE_INS },
                        new Allele[]{ Allele.REF_N, Allele.SV_SIMPLE_INS },
                        SVConcordanceCollapser.GenotypeState.HOM_VAR1
                },

                ////////////////////////////////////////
                // Edge cases
                ////////////////////////////////////////

                // hom ref with a no-call
                {
                        new Allele[]{ Allele.NO_CALL, Allele.REF_N },
                        new Allele[]{ Allele.REF_N, Allele.SV_SIMPLE_INS },
                        SVConcordanceCollapser.GenotypeState.HOM_REF
                },
                // hom var with a no-call
                {
                        new Allele[]{ Allele.NO_CALL, Allele.SV_SIMPLE_INS },
                        new Allele[]{ Allele.REF_N, Allele.SV_SIMPLE_INS },
                        SVConcordanceCollapser.GenotypeState.HOM_VAR1
                },
                // switch allele order
                {
                        new Allele[]{ Allele.SV_SIMPLE_INS, Allele.NO_CALL },
                        new Allele[]{ Allele.REF_N, Allele.SV_SIMPLE_INS },
                        SVConcordanceCollapser.GenotypeState.HOM_VAR1
                },
        };
    }

    @Test(dataProvider= "testGetTruthStateNonCNVData")
    public void testGetTruthStateNonCNV(final Allele[] allelesArr,
                                        final Allele[] variantAlleles,
                                        final SVConcordanceCollapser.GenotypeState expected) {
        final SVCallRecord record = SVTestUtils.newCallRecordWithAlleles(
                Arrays.asList(Allele.REF_N, Allele.REF_N),
                Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_INS),
                StructuralVariantType.INS,
                null,
                null
        );
        final Genotype g = alleleArrayToGenotype(allelesArr, null, null);

        final SVConcordanceCollapser collapser = new SVConcordanceCollapser(false);
        final SVConcordanceCollapser.GenotypeState actual = collapser.getNonCNVGenotypeState(g);
        Assert.assertEquals(actual, expected);
    }


    @DataProvider(name = "testGetTruthStateCNVData")
    public Object[][] testGetTruthStateCNVData() {
        return new Object[][]{

                ////////////////////////////////////////
                // Empty edge cases
                ////////////////////////////////////////

                // ploidy 0
                {
                        0,
                        0,
                        false,
                        new Allele[]{ Allele.REF_N, Allele.SV_SIMPLE_DEL },
                        StructuralVariantType.DEL,
                        SVConcordanceCollapser.GenotypeState.NO_CALL
                },
                {
                        0,
                        0,
                        false,
                        new Allele[]{ Allele.REF_N, Allele.SV_SIMPLE_DUP },
                        StructuralVariantType.DUP,
                        SVConcordanceCollapser.GenotypeState.NO_CALL
                },
                {
                        0,
                        0,
                        false,
                        new Allele[]{ Allele.REF_N, Allele.SV_SIMPLE_DEL, Allele.SV_SIMPLE_DUP },
                        StructuralVariantType.CNV,
                        SVConcordanceCollapser.GenotypeState.NO_CALL
                },

                ////////////////////////////////////////
                // Haploid
                ////////////////////////////////////////

                // hom ref
                {
                        1,
                        1,
                        false,
                        new Allele[]{ Allele.REF_N, Allele.SV_SIMPLE_DEL },
                        StructuralVariantType.DEL,
                        SVConcordanceCollapser.GenotypeState.HOM_REF
                },
                {
                        1,
                        1,
                        false,
                        new Allele[]{ Allele.REF_N, Allele.SV_SIMPLE_DUP },
                        StructuralVariantType.DUP,
                        SVConcordanceCollapser.GenotypeState.HOM_REF
                },
                {
                        1,
                        1,
                        false,
                        new Allele[]{ Allele.REF_N, Allele.SV_SIMPLE_DEL, Allele.SV_SIMPLE_DUP },
                        StructuralVariantType.CNV,
                        SVConcordanceCollapser.GenotypeState.HOM_REF
                },
                // hom var
                {
                        1,
                        0,
                        false,
                        new Allele[]{ Allele.REF_N, Allele.SV_SIMPLE_DEL },
                        StructuralVariantType.DEL,
                        SVConcordanceCollapser.GenotypeState.HOM_VAR1
                },
                {
                        1,
                        2,
                        false,
                        new Allele[]{ Allele.REF_N, Allele.SV_SIMPLE_DUP },
                        StructuralVariantType.DUP,
                        SVConcordanceCollapser.GenotypeState.HOM_VAR2
                },
                // CNV - hom DEL
                {
                        1,
                        0,
                        false,
                        new Allele[]{ Allele.REF_N, Allele.SV_SIMPLE_DEL, Allele.SV_SIMPLE_DUP },
                        StructuralVariantType.CNV,
                        SVConcordanceCollapser.GenotypeState.HOM_VAR1
                },
                // CNV - hom DUP
                {
                        1,
                        2,
                        false,
                        new Allele[]{ Allele.REF_N, Allele.SV_SIMPLE_DEL, Allele.SV_SIMPLE_DUP },
                        StructuralVariantType.CNV,
                        SVConcordanceCollapser.GenotypeState.HOM_VAR2
                },

                ////////////////////////////////////////
                // Diploid
                ////////////////////////////////////////

                // hom ref
                {
                        2,
                        2,
                        false,
                        new Allele[]{ Allele.REF_N, Allele.SV_SIMPLE_DEL },
                        StructuralVariantType.DEL,
                        SVConcordanceCollapser.GenotypeState.HOM_REF
                },
                {
                        2,
                        2,
                        false,
                        new Allele[]{ Allele.REF_N, Allele.SV_SIMPLE_DUP },
                        StructuralVariantType.DUP,
                        SVConcordanceCollapser.GenotypeState.HOM_REF
                },
                {
                        2,
                        2,
                        false,
                        new Allele[]{ Allele.REF_N, Allele.SV_SIMPLE_DEL, Allele.SV_SIMPLE_DUP },
                        StructuralVariantType.CNV,
                        SVConcordanceCollapser.GenotypeState.HOM_REF
                },

                // hom var
                {
                        2,
                        0,
                        false,
                        new Allele[]{ Allele.REF_N, Allele.SV_SIMPLE_DEL },
                        StructuralVariantType.DEL,
                        SVConcordanceCollapser.GenotypeState.HOM_VAR1
                },
                // Note strict CNV interpretation is off here, so we interpret copy states 4 and 5 as hom dup
                {
                        2,
                        4,
                        false,
                        new Allele[]{ Allele.REF_N, Allele.SV_SIMPLE_DUP },
                        StructuralVariantType.DUP,
                        SVConcordanceCollapser.GenotypeState.HOM_VAR2
                },
                {
                        2,
                        5,
                        false,
                        new Allele[]{ Allele.REF_N, Allele.SV_SIMPLE_DEL, Allele.SV_SIMPLE_DUP },
                        StructuralVariantType.CNV,
                        SVConcordanceCollapser.GenotypeState.HOM_VAR2
                },

                // Het - strict CNV interpretation off
                {
                        2,
                        1,
                        false,
                        new Allele[]{ Allele.REF_N, Allele.SV_SIMPLE_DEL, Allele.SV_SIMPLE_DUP },
                        StructuralVariantType.CNV,
                        SVConcordanceCollapser.GenotypeState.HET_REF_VAR1
                },
                {
                        2,
                        3,
                        false,
                        new Allele[]{ Allele.REF_N, Allele.SV_SIMPLE_DUP },
                        StructuralVariantType.DUP,
                        SVConcordanceCollapser.GenotypeState.HET_REF_VAR2
                },

                ////////////////////////////////////////
                // Strict CNV interpretation mode
                ////////////////////////////////////////

                // Note we cannot determine the true alleles unambiguously, so expect IS_MIXED
                {
                        2,
                        2,
                        true,
                        new Allele[]{ Allele.REF_N, Allele.SV_SIMPLE_DEL, Allele.SV_SIMPLE_DUP },
                        StructuralVariantType.CNV,
                        SVConcordanceCollapser.GenotypeState.UNKNOWN
                },
                {
                        2,
                        3,
                        true,
                        new Allele[]{ Allele.REF_N, Allele.SV_SIMPLE_DEL, Allele.SV_SIMPLE_DUP },
                        StructuralVariantType.CNV,
                        SVConcordanceCollapser.GenotypeState.UNKNOWN
                },
                {
                        2,
                        4,
                        true,
                        new Allele[]{ Allele.REF_N, Allele.SV_SIMPLE_DEL, Allele.SV_SIMPLE_DUP },
                        StructuralVariantType.CNV,
                        SVConcordanceCollapser.GenotypeState.UNKNOWN
                },
                {
                        2,
                        5,
                        true,
                        new Allele[]{ Allele.REF_N, Allele.SV_SIMPLE_DEL, Allele.SV_SIMPLE_DUP },
                        StructuralVariantType.CNV,
                        SVConcordanceCollapser.GenotypeState.UNKNOWN
                },

                // Special cases that are unambiguous
                {
                        2,
                        0,
                        true,
                        new Allele[]{ Allele.REF_N, Allele.SV_SIMPLE_DEL, Allele.SV_SIMPLE_DUP },
                        StructuralVariantType.CNV,
                        SVConcordanceCollapser.GenotypeState.HOM_VAR1
                },
                {
                        2,
                        1,
                        true,
                        new Allele[]{ Allele.REF_N, Allele.SV_SIMPLE_DEL, Allele.SV_SIMPLE_DUP },
                        StructuralVariantType.CNV,
                        SVConcordanceCollapser.GenotypeState.HET_REF_VAR1
                },
                {
                        1,
                        1,
                        true,
                        new Allele[]{ Allele.REF_N, Allele.SV_SIMPLE_DEL, Allele.SV_SIMPLE_DUP },
                        StructuralVariantType.CNV,
                        SVConcordanceCollapser.GenotypeState.HOM_REF
                },
                {
                        1,
                        2,
                        true,
                        new Allele[]{ Allele.REF_N, Allele.SV_SIMPLE_DEL, Allele.SV_SIMPLE_DUP },
                        StructuralVariantType.CNV,
                        SVConcordanceCollapser.GenotypeState.HOM_VAR2
                },

        };
    }

    @Test(dataProvider= "testGetTruthStateCNVData")
    public void testGetTruthStateCNV(final Integer truthExpectedCopyNumberArr,
                                     final Integer truthCopyNumberArr,
                                     final boolean strictCnvInterpretation,
                                     final Allele[] variantAlleles,
                                     final StructuralVariantType svtype,
                                     final SVConcordanceCollapser.GenotypeState expected) {
        final SVCallRecord record = SVTestUtils.newCallRecordWithAlleles(
                Arrays.asList(Allele.NO_CALL, Allele.NO_CALL),
                Arrays.asList(variantAlleles),
                svtype,
                2,
                2
        );
        final Allele[] truthAllelesArr = new Allele[truthExpectedCopyNumberArr];
        Arrays.fill(truthAllelesArr, Allele.NO_CALL);
        final Genotype g = alleleArrayToGenotype(truthAllelesArr, truthExpectedCopyNumberArr, truthCopyNumberArr);

        final SVConcordanceCollapser collapser = new SVConcordanceCollapser(strictCnvInterpretation);
        final SVConcordanceCollapser.GenotypeState actual = collapser.getCnvGenotypeState(g, record);
        Assert.assertEquals(actual, expected);
    }

    private Genotype alleleArrayToGenotype(final Allele[] allelesArr,
                                           final Integer expectedCopyNumber,
                                           final Integer copyNumber) {
        final GenotypeBuilder builder = new GenotypeBuilder("");
        if (expectedCopyNumber != null) {
            builder.attribute(GATKSVVCFConstants.EXPECTED_COPY_NUMBER_FORMAT, expectedCopyNumber);
        }
        if (copyNumber != null) {
            builder.attribute(GATKSVVCFConstants.COPY_NUMBER_FORMAT, copyNumber);
        }
        builder.alleles(Arrays.asList(allelesArr));
        return builder.make();
    }
}
