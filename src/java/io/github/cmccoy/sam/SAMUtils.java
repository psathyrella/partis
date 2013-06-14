package io.github.cmccoy.sam;

import java.util.List;
import java.util.ArrayList;

import net.sf.samtools.SAMRecord;
import net.sf.samtools.CigarElement;

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
}
