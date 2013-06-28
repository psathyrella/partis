package io.github.cmccoy.sam;

import java.util.List;
import java.util.ArrayList;

import net.sf.picard.reference.ReferenceSequence;
import net.sf.picard.reference.ReferenceSequenceFile;
import net.sf.samtools.AlignmentBlock;
import net.sf.samtools.CigarElement;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.util.PeekIterator;

import static com.google.common.base.Preconditions.*;

import io.github.cmccoy.sam.Mutation.MutationType;

public class SAMUtils {

    public static List<Mutation> enumerateMutations(final SAMRecord record,
                                                    final byte[] referenceBases) {
        checkNotNull(record, "null SAM record");
        checkNotNull(referenceBases, "null reference bases");
        // Always be noting 1-based indexing!
        final int start = record.getAlignmentStart() - 1;
        final byte[] qBases = record.getReadBases();
        final List<Mutation> result = new ArrayList<Mutation>();

        int qi = 0, ri = start;

        for(final CigarElement e : record.getCigar().getCigarElements()) {
            final int length = e.getLength();
            switch (e.getOperator()) {
            case H:
            case S:
                break;
            case I:
                result.add(new Mutation(MutationType.INSERTION, ri, "-",
                                        new String(qBases, qi, length)));
                break;
            case D:
                result.add(new Mutation(MutationType.DELETION, ri,
                                        new String(referenceBases, ri, 1),
                                        "-"));
                break;
            case M:
                for(int i = 0; i < length; ++i) {
                    if(qBases[qi+i] != referenceBases[ri+i]) {
                        result.add(new Mutation(MutationType.MUTATION, ri + i,
                                                new String(referenceBases, ri + i, 1),
                                                new String(qBases, qi + i, 1)));
                    }
                }
                break;
            case X:
                // Bases not equal - skip check
                for(int i = 0; i < length; ++i) {
                    result.add(new Mutation(MutationType.MUTATION, ri + i,
                                            new String(referenceBases, ri + i, 1),
                                            new String(qBases, qi + i, 1)));
                }
            case EQ:
                // Bases equal - no action.
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

    /**
     * Calculate the probability that each base matches the reference,
     * given an iterable of equiprobable alignments.
     *
     * \param ref Reference sequences
     * \param records SAM alignments for a given sequence, primary first
     * \returns 100 * p(base matches)
     */
    public static byte[] calculateBaseEqualProbabilities(ReferenceSequenceFile ref,
                                                         Iterable<SAMRecord> records) {
      final PeekIterator<SAMRecord> it = new PeekIterator<SAMRecord>(checkNotNull(records).iterator());
      checkArgument(it.hasNext(), "No records in iterable!");
      final SAMRecord primary = it.peek();
      final byte[] bases = primary.getReadBases();
      checkArgument(!primary.getNotPrimaryAlignmentFlag(), "First record is not primary");

      byte[] matches = new byte[primary.getReadLength()];
      byte[] tot = new byte[primary.getReadLength()];
      for(int i = 0; i < primary.getReadLength(); i++) {
        matches[i] = 0;
        tot[i] = 0;
      }

      while(it.hasNext()) {
        final SAMRecord record = it.next();
        checkState(record.getReadName().equals(primary.getReadName()),
            "Unmatched read names in group (%s, %s)",
            record.getReadName(),
            primary.getReadName());
        final String referenceName = record.getReferenceName();
        final byte[] refBases = ref.getSequence(referenceName).getBases();
        for(final AlignmentBlock b : record.getAlignmentBlocks()) {
          // 0-based index in read, ref
          final int readStart = b.getReadStart() - 1;
          final int refStart = b.getReferenceStart() - 1;
          for(int i = 0; i < b.getLength(); i++) {
            tot[readStart + i]++;
            if(bases[readStart + i] == refBases[refStart + i])
              matches[readStart + i]++;
          }
        }
      }

      byte[] result = new byte[matches.length];
      for(int i = 0; i < tot.length; i++) {
        if(tot[i] > 0) {
          float pMatch = ((float)matches[i]) / tot[i];
          result[i] =  (byte)(100 * pMatch);
        } else {
          result[i] = -1;
        }
      }

      return result;
    }
}
