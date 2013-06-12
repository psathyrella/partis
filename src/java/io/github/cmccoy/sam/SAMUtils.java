package io.github.cmccoy.sam;

import java.util.List;
import java.util.ArrayList;

import net.sf.samtools.SAMRecord;
import net.sf.samtools.CigarElement;

import com.google.common.base.Objects;
import com.google.common.collect.ComparisonChain;

public class SAMUtils {
  public static class Mutation implements Comparable<Mutation> {
    private final char type;
    //  public final String reference;
    private final int position;
    private final char wt;
    private final char mut;
    public Mutation(char type, int position, char wt, char mut) {
      this.type = type;
      //  this.reference = reference;
      this.position = position;
      this.wt = wt;
      this.mut = mut;
    }

    @Override
    public String toString() {
      return Objects.toStringHelper(this)
        .add("type", type)
        .add("position", position)
        .add("wt", wt)
        .add("mut", mut)
        .toString();
    }

    @Override
    public int hashCode() {
      return Objects.hashCode(type, position, wt, mut);
    }

    @Override
    public int compareTo(Mutation other) {
      return ComparisonChain.start()
        .compare(this.type, other.type)
        .compare(this.position, other.position)
        .compare(this.wt, other.wt)
        .compare(this.mut, other.mut)
        .result();
    }
    @Override
    public boolean equals(Object obj) {
	if (obj == null) return false;
	if (obj.getClass() != getClass()) return false;
	final Mutation other = (Mutation) obj;

	return other.type == this.type &&
	    other.position == this.position &&
	    other.wt == this.wt &&
	    other.mut == this.mut;
    }

    public char getType() { return type; }
    public int getPosition() { return position; }
    public char getWildType() { return wt; }
    public char getMutant() { return mut; }
  }

  public static List<Mutation> enumerateMutations(final SAMRecord record,
      final byte[] referenceBases) {
    // Always be noting 1-based indexing!
    final int start = record.getAlignmentStart() - 1;
    final byte[] qBases = record.getReadBases();
    final List<Mutation> result = new ArrayList<Mutation>();

    int qi = 0, ri = start;

    for(final CigarElement e : record.getCigar().getCigarElements()) {
      switch (e.getOperator()) {
        case H:
        case S:
          break;
        case I:
          result.add(new Mutation('I', ri, '-', (char)qBases[qi]));
          break;
        case D:
          result.add(new Mutation('D', ri, (char)referenceBases[ri], '-'));
          break;
        case M:
          for(int i = 0; i < e.getLength(); ++i) {
            if(qBases[qi+i] != referenceBases[ri+i]) {
              result.add(new Mutation('M', ri + i,
                    (char)referenceBases[ri+i],
                    (char)qBases[qi+i]));
            }
          }
          break;
        default:
          throw new IllegalArgumentException("unknown " + e.toString());
      }
      final int l = e.getLength();
      if(e.getOperator().consumesReadBases())
        qi += l;
      if(e.getOperator().consumesReferenceBases())
        ri += l;
    }
    return result;
  }
}
