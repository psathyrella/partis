package io.github.cmccoy.sam;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import net.sf.samtools.SAMRecord;
import net.sf.picard.util.IntervalTreeMap;
import net.sf.picard.util.Interval;
import static com.google.common.base.Preconditions.*;

public final class VDJAnnotator {

    private final SAMRecord read;
    private final IntervalTreeMap<String> tree;
    private final Map<String, RegionAnnotation> annotations;

    public VDJAnnotator(final SAMRecord read, final IntervalTreeMap<String> tree) {
        checkNotNull(read, "null SAM record");
        checkNotNull(tree, "null interval tree");
        this.read = read;
        this.tree = tree;
        this.annotations = new HashMap<String, RegionAnnotation>();
    }

    public Map<String, RegionAnnotation> annotateVDJ() {
        final String referenceName = read.getReferenceName();
        final int nm = read.getIntegerAttribute("NM");
        final int as = read.getIntegerAttribute("AS");

        final List<AlignedPair> alignedPairs = SAMUtils.getAlignedPairs(read);
        for(final AlignedPair ap : alignedPairs) {
            if(ap.isIndel() || ap.isUnknown()) continue;
            final int qpos = ap.getQueryPosition();
            final int rpos = ap.getReferencePosition();

            final List<String> overlapping = new ArrayList<String>(
                tree.getOverlapping(new Interval(referenceName, rpos + 1, rpos + 1)));
            overlapping.add(new String(referenceName.substring(3, 4)));
            final boolean isMutation = ap.isMutation();

            for(final String region : overlapping) {
                final RegionAnnotation annot;
                if(this.annotations.containsKey(region)) {
                    annot = this.annotations.get(region);
                } else {
                    annot = new RegionAnnotation(referenceName, ap.getQueryPosition(), as, nm);
                    this.annotations.put(region, annot);
                }

                annot.aligned++;
                annot.minqpos = Math.min(qpos, annot.minqpos);
                annot.maxqpos = Math.max(qpos, annot.minqpos);
                if(isMutation) {
                    annot.mismatch++;
                }
            }
        }

        return this.annotations;
    }

    public static class RegionAnnotation {
      public final String name;
      public int aligned;
      public int mismatch;
      public int minqpos;
      public int maxqpos;
      public int alignmentScore;
      public int nm;

      public RegionAnnotation(final String name, final int qpos, final int alignmentScore, final int nm) {
          this.name = name;
          this.mismatch = 0;
          this.aligned = 0;
          this.minqpos = qpos;
          this.maxqpos = qpos;
          this.alignmentScore = alignmentScore;
          this.nm = nm;
      }
    }
}
