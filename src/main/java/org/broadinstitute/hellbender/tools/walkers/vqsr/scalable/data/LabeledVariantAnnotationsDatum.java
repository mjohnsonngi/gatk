package org.broadinstitute.hellbender.tools.walkers.vqsr.scalable.data;

import com.google.common.collect.ImmutableList;
import com.google.common.collect.ImmutableSet;
import htsjdk.samtools.util.Locatable;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;

import java.util.List;
import java.util.Set;

public final class LabeledVariantAnnotationsDatum implements Locatable {
    final SimpleInterval interval;
    final Allele refAllele;
    final ImmutableList<Allele> altAlleles;
    final VariantType variantType;
    final ImmutableSet<String> labels;
    final double[] annotations;

    public LabeledVariantAnnotationsDatum(final VariantContext vc,
                                          final List<Allele> altAlleles,
                                          final VariantType variantType,
                                          final Set<String> labels,
                                          final List<String> sortedAnnotationNames,
                                          final boolean useASAnnotations) {
        Utils.validate(!useASAnnotations || altAlleles.size() == 1,
                "Datum should only be associated with one alt allele in allele-specific mode.");
        this.interval = new SimpleInterval(vc);
        this.refAllele = vc.getReference();
        this.altAlleles = ImmutableList.copyOf(altAlleles);
        this.variantType = variantType;
        this.labels = ImmutableSet.copyOf(labels);
        this.annotations = sortedAnnotationNames.stream()
                .mapToDouble(a -> decodeAnnotation(vc, altAlleles, a, useASAnnotations))
                .toArray();
    }

    @Override
    public String getContig() {
        return interval.getContig();
    }

    @Override
    public int getStart() {
        return interval.getStart();
    }

    @Override
    public int getEnd() {
        return interval.getEnd();
    }

    public VariantType getVariantType() {
        return variantType;
    }

    public ImmutableSet<String> getLabels() {
        return labels;
    }

    public double[] getAnnotations() {
        return annotations;
    }

    // code retained from VQSR
    private static double decodeAnnotation(final VariantContext vc,
                                           final List<Allele> altAlleles,
                                           final String annotationName,
                                           final boolean useASAnnotations) {
        double value;
        try {
            //if we're in allele-specific mode and an allele-specific annotation has been requested, parse the appropriate value from the list
            if (useASAnnotations && annotationName.startsWith(GATKVCFConstants.ALLELE_SPECIFIC_PREFIX)) {
                final List<Object> valueList = vc.getAttributeAsList(annotationName);
                final Allele altAllele = altAlleles.get(0);
                //FIXME: we need to look at the ref allele here too
                if (vc.hasAllele(altAllele)) {
                    final int altIndex = vc.getAlleleIndex(altAllele) - 1; //- 1 is to convert the index from all alleles (including reference) to just alternate alleles
                    value = Double.parseDouble((String) valueList.get(altIndex));
                } else {
                    //if somehow our alleles got mixed up
                    throw new IllegalStateException("ExtractAnnotationsVariantDatum allele " + altAllele + " is not contained in the input VariantContext.");
                }
            } else {
                value = vc.getAttributeAsDouble(annotationName, Double.NaN);
            }
            if (Double.isInfinite(value)) {
                value = Double.NaN;
            }
        } catch (final NumberFormatException e) {
            value = Double.NaN;
        }
        return value;
    }
}