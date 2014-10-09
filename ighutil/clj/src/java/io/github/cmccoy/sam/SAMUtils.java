package io.github.cmccoy.sam;

import io.github.cmccoy.dna.IUPACUtils;
import io.github.cmccoy.sam.Mutation.MutationType;
import net.sf.samtools.AlignmentBlock;
import net.sf.samtools.CigarElement;
import net.sf.samtools.CigarOperator;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.util.PeekIterator;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import static com.google.common.base.Preconditions.*;

public class SAMUtils {

    public static int countAlignedBases(final SAMRecord record) {
        checkNotNull(record, "null SAM record");
        checkNotNull(record.getCigar(), "Null cigar");
        final List<CigarElement> cigarElements = record.getCigar().getCigarElements();

        int result = 0;
        for (final CigarElement e : cigarElements) {
            final CigarOperator op = e.getOperator();
            if (op == CigarOperator.EQ ||
                op == CigarOperator.M ||
                op == CigarOperator.X)
              result += e.getLength();
        }
        return result;
    }

    public static int countMutations(final SAMRecord record,
                                     final byte[] referenceBases) {
        checkNotNull(record, "null SAM record");
        checkNotNull(referenceBases, "null reference bases");

        int result = 0;
        final byte[] qBases = record.getReadBases();
        for (final AlignmentBlock b : record.getAlignmentBlocks()) {
            final int qStart = b.getReadStart() - 1,
                      rStart = b.getReferenceStart() - 1,
                      length = b.getLength();
            for(int i = 0; i < length; i++) {
                if(referenceBases[rStart + i] != qBases[qStart + i])
                    result++;
            }
        }
        return result;
    }

    /**
     * List all mutations relative to the reference.
     *
     * @param record         A SAM record
     * @param referenceBases Matched reference (all bases)
     * @return
     */
    public static List<Mutation> enumerateMutations(final SAMRecord record,
                                                    final byte[] referenceBases) {
        checkNotNull(record, "null SAM record");
        checkNotNull(referenceBases, "null reference bases");
        // Always be noting 1-based indexing!
        final int start = record.getAlignmentStart() - 1;
        final byte[] qBases = record.getReadBases();
        final List<Mutation> result = new ArrayList<Mutation>();

        int qi = 0, ri = start;

        for (final CigarElement e : record.getCigar().getCigarElements()) {
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
                for (int i = 0; i < length; ++i) {
                    if (qBases[qi + i] != referenceBases[ri + i]) {
                        result.add(new Mutation(MutationType.MUTATION, ri + i,
                                                new String(referenceBases, ri + i, 1),
                                                new String(qBases, qi + i, 1)));
                    }
                }
                break;
            case X:
                // Bases not equal - skip check
                for (int i = 0; i < length; ++i) {
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
            if (e.getOperator().consumesReadBases())
                qi += l;
            if (e.getOperator().consumesReferenceBases())
                ri += l;
        }
        return result;
    }

    /**
     * Calculate the probability that each base matches the reference,
     * given an iterable of equiprobable alignments.
     * <p/>
     * \param ref Reference sequences
     * \param records SAM alignments for a given sequence, primary first
     * \returns 100 * p(base matches)
     */
    public static byte[] calculateBaseEqualProbabilities(final Map<String, byte[]> ref,
                                                         final Iterable<SAMRecord> records) {
        final PeekIterator<SAMRecord> it = new PeekIterator<SAMRecord>(checkNotNull(records).iterator());
        checkArgument(it.hasNext(), "No records in iterable!");
        final SAMRecord primary = it.peek();
        final byte[] bases = primary.getReadBases();
        checkArgument(!primary.getNotPrimaryAlignmentFlag(), "First record is not primary");

        byte[] matches = new byte[primary.getReadLength()];
        byte[] tot = new byte[primary.getReadLength()];

        java.util.Arrays.fill(matches, (byte) 0);
        java.util.Arrays.fill(tot, (byte) 0);

        while (it.hasNext()) {
            final SAMRecord record = it.next();
            checkState(record.getReadName().equals(primary.getReadName()),
                       "Unmatched read names in group (%s, %s)",
                       record.getReadName(),
                       primary.getReadName());
            final String referenceName = record.getReferenceName();
            final byte[] refBases = checkNotNull(ref.get(referenceName),
                                                 "No reference %s", referenceName);
            for (final AlignmentBlock b : record.getAlignmentBlocks()) {
                // 0-based index in read, ref
                final int readStart = b.getReadStart() - 1;
                final int refStart = b.getReferenceStart() - 1;
                for (int i = 0; i < b.getLength(); i++) {
                    tot[readStart + i]++;
                    if (IUPACUtils.isSubset(IUPACUtils.packByte(refBases[refStart + i]),
                                            IUPACUtils.packByte(bases[readStart + i])))
                        matches[readStart + i]++;
                }
            }
        }

        byte[] result = new byte[matches.length];
        for (int i = 0; i < tot.length; i++) {
            if (tot[i] > 0) {
                float pMatch = ((float) matches[i]) / tot[i];
                result[i] = (byte) (100 * pMatch);
            } else {
                result[i] = -1;
            }
        }

        return result;
    }

    public static List<AlignedPair> getAlignedPairs(final SAMRecord record) {
      checkNotNull(record, "null SAM record");
      final byte[] baseEqualProbs = record.getByteArrayAttribute("bq");
      final boolean hasBaseEqualProbs = baseEqualProbs != null;
      checkArgument(!hasBaseEqualProbs || baseEqualProbs.length == record.getReadLength());
      checkNotNull(baseEqualProbs, "Missing bq tag.");
      final List<AlignedPair> result = new ArrayList<AlignedPair>(record.getReadLength());

      // Always be noting 1-based indexing!
      final int start = record.getAlignmentStart() - 1;

      int qi = 0, ri = start;

      for (final CigarElement e : record.getCigar().getCigarElements()) {
        final int length = e.getLength();
        switch (e.getOperator()) {
          case H:
          case S:
            break;
          case I:
            for(int i = 0; i < length; i++)
              result.add(new AlignedPair(-1, qi + i, AlignedPair.MatchesReference.FALSE));
            break;
          case D:
            for(int i = 0; i < length; i++)
              result.add(new AlignedPair(ri + i, -1, AlignedPair.MatchesReference.FALSE));
            break;
          case M:
          case X:
          case EQ:
            // TODO: more sophisticated X / EQ handling? right now: UNKNOWN
            for (int i = 0; i < length; ++i) {
              final AlignedPair.MatchesReference m;
              if (hasBaseEqualProbs)
                  if(baseEqualProbs[qi + i] == 100)
                    m = AlignedPair.MatchesReference.TRUE;
                  else if(baseEqualProbs[qi + i] == 0)
                    m = AlignedPair.MatchesReference.FALSE;
                  else m = AlignedPair.MatchesReference.UNKNOWN;
              else
                m = AlignedPair.MatchesReference.UNKNOWN;
              result.add(new AlignedPair(ri + i, qi + i, m));
            }
            break;
          default:
            throw new IllegalArgumentException("unknown " + e.toString());
        }
        final int l = e.getLength();
        if (e.getOperator().consumesReadBases())
          qi += l;
        if (e.getOperator().consumesReferenceBases())
          ri += l;
      }
      return result;
    }

    /**
     * Calculates the <i>possibility</i> that each base matches the reference.
     * Differs from calculateBaseEqualProbabilities in that if any reference
     * matches the query (including IUPAC), the field is true.
     */
    public static byte[] calculatePossibleMatchToRef(final Map<String, byte[]> ref,
                                                     final Iterable<SAMRecord> records) {
        final byte MISSING = -1;
        final PeekIterator<SAMRecord> it = new PeekIterator<SAMRecord>(checkNotNull(records).iterator());
        checkArgument(it.hasNext(), "No records in iterable!");
        final SAMRecord primary = it.peek();
        final byte[] bases = primary.getReadBases();
        checkArgument(!primary.getNotPrimaryAlignmentFlag(), "First record is not primary");
        byte[] mayMatch = new byte[primary.getReadLength()];
        java.util.Arrays.fill(mayMatch, MISSING);

        final byte[] packedRead = IUPACUtils.packBytes(bases);

        while (it.hasNext()) {
            final SAMRecord record = it.next();
            checkState(record.getReadName().equals(primary.getReadName()),
                       "Unmatched read names in group (%s, %s)",
                       record.getReadName(),
                       primary.getReadName());
            final String referenceName = record.getReferenceName();
            final byte[] refBases = checkNotNull(ref.get(referenceName),
                                                 "No reference %s", referenceName);
            final byte[] packedRef = IUPACUtils.packBytes(refBases);
            for (final AlignmentBlock b : record.getAlignmentBlocks()) {
                // 0-based index in read, ref
                final int readStart = b.getReadStart() - 1;
                final int refStart = b.getReferenceStart() - 1;
                for (int i = 0; i < b.getLength(); i++) {
                    if (mayMatch[readStart + i] == MISSING)
                        mayMatch[readStart + i] = 0;
                    if (IUPACUtils.isSubset(packedRef[refStart + i],
                                            packedRead[readStart + i]))
                        mayMatch[readStart + i] = 1;
                }
            }
        }

        return mayMatch;
    }
}
