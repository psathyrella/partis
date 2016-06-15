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
        final int startPos = read.getAlignmentStart() - 1;
        final int endPos = read.getAlignmentEnd();
        checkState(startPos >= 0);
        checkState(endPos != 0);

        final int refLength = read.getHeader().getSequence(read.getReferenceIndex()).getSequenceLength();


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
                    annot = new RegionAnnotation(referenceName, startPos, refLength - endPos,
                                                 ap.getQueryPosition(), as, nm);
                    this.annotations.put(region, annot);
                }

                annot.aligned++;
                annot.qstart = Math.min(qpos, annot.qstart);
                annot.qend = Math.max(qpos, annot.qstart);
                if(isMutation) {
                    annot.mismatch++;
                }
            }
        }

        return this.annotations;
    }

    public static class RegionAnnotation {
      public final String name;
      public final int eroded5P;
      public final int eroded3P;

      public int aligned;
      public int mismatch;
      public int qstart;
      public int qend;
      public int alignmentScore;
      public int nm;

      public RegionAnnotation(final String name, final int eroded5P, final int eroded3P,
                              final int qpos, final int alignmentScore, final int nm) {
          this.name = name;
          this.eroded5P = eroded5P;
          this.eroded3P = eroded3P;
          this.mismatch = 0;
          this.aligned = 0;
          this.qstart = qpos;
          this.qend = qpos;
          this.alignmentScore = alignmentScore;
          this.nm = nm;
      }

      @Override
      public int hashCode() {
        return java.util.Objects.hash(
            name,
            eroded5P,
            eroded3P,
            mismatch,
            aligned,
            qstart,
            qend,
            alignmentScore,
            nm);
      }
    }
}
