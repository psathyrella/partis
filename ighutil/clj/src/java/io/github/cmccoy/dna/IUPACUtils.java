package io.github.cmccoy.dna;

import java.util.Arrays;
import java.util.Stack;

import static com.google.common.base.Preconditions.*;

/**
 * Utilities for working with IUPAC bases, encoded as bytes.
 * The first 4 bits of each byte are used, in order ACGT.
 * Ambiguous bases are represented as multiple 1 bits.
 */
public class IUPACUtils {
    private final static byte[] LOOKUP =
            {'U', 'T', 'A', 'W', 'C', 'Y', 'M',
                    'H', 'G', 'K', 'R', 'D', 'S', 'B', 'V', 'N'};
    private static int BITS_PER_BASE = 4;
    private static int SINGLE_BASE_MASK = 0xf;

    /**
     * Pack a single byte.
     *
     * @param b An IUPAC character
     * @return byte encoded.
     */
    public static byte packByte(final byte b) {
        switch ((char) b) {
            case 'U':
                return 1;
            case 'T':
                return 1;
            case 'A':
                return 2;
            case 'W':
                return 3;
            case 'C':
                return 4;
            case 'Y':
                return 5;
            case 'M':
                return 6;
            case 'H':
                return 7;
            case 'G':
                return 8;
            case 'K':
                return 9;
            case 'R':
                return 10;
            case 'D':
                return 11;
            case 'S':
                return 12;
            case 'B':
                return 13;
            case 'V':
                return 14;
            case 'N':
                return 15;
            default:
                throw new IllegalArgumentException(String.format("Unknown base: %s", (char)b));
        }
    }

    /**
     * Pack a string into a Byte
     *
     * @param bytes
     * @return the packed equivalent.
     */
    public static byte[] packBytes(final byte[] bytes) {
        checkNotNull(bytes);
        byte[] result = new byte[bytes.length];
        for (int i = 0; i < bytes.length; ++i) {
            result[i] = packByte(bytes[i]);
        }
        return result;
    }

    /**
     * Opposite of packByte.
     *
     * @param s
     * @return
     */
    public static byte unpackByte(final byte s) {
        checkIUPACBase(s);
        return LOOKUP[s];
    }

    public static byte[] unpackBytes(final byte[] packed) {
        checkNotNull(packed);
        byte[] result = new byte[packed.length];
        for (int i = 0; i < packed.length; i++) {
            checkIUPACBase(packed[i]);
            result[i] = unpackByte(packed[i]);
        }
        return result;
    }

    public static boolean areCompatible(final byte s1, final byte s2) {
        checkIUPACBase(s1);
        checkIUPACBase(s2);
        return (s1 & s2) > 0;
    }

    public static boolean areCompatible(final byte[] s1, final byte[] s2) {
        return areCompatible(s1, s2, 0, s1.length);
    }

    public static boolean areCompatible(final byte[] s1, final byte[] s2,
                                        final int start, final int length) {
        checkNotNull(s1);
        checkNotNull(s2);
        checkArgument(length == s2.length,
                "Lengths differ: %d != %d",
                length,
                s2.length);
        checkArgument(start < s1.length,
                "Start after end of array (%d >= %d)",
                start,
                s1.length);
        final int end = start + length;
        checkArgument(end <= s1.length,
                "End after end of array (%d > %d)", end, s1.length);
        for (int i = 0; i < length; ++i) {
            if (!areCompatible(s1[start + i], s2[i]))
                return false;
        }
        return true;
    }

    /**
     * Tests whether s2 is a subset of s1
     */
    public static boolean isSubset(byte s1, byte s2) {
        checkIUPACBase(s1);
        checkIUPACBase(s2);
        return 0 == (s2 & ~s1);
    }

    public static boolean isSubset(final byte[] s1, final byte[] s2,
                                   final int start, final int length) {
        checkNotNull(s1);
        checkNotNull(s2);
        checkArgument(length == s2.length,
                "Lengths differ: %d != %d",
                length,
                s2.length);
        checkArgument(start < s1.length,
                "Start after end of array (%d >= %d)",
                start,
                s1.length);
        final int end = start + length;
        checkArgument(end <= s1.length,
                "End after end of array (%d > %d)", end, s1.length);
        for (int i = 0; i < length; ++i) {
            if (!isSubset(s1[start + i], s2[i]))
                return false;
        }
        return true;
    }

    public static boolean isSubset(final byte[] s1, final byte[] s2) {
        return isSubset(s1, s2, 0, s1.length);
    }

    public static boolean isAmbiguous(final byte b) {
        checkIUPACBase(b);
        return Integer.bitCount((int) b) != 1;
    }

    /**
     * Number of <i>unambiguous</i> bases this character codes
     */
    public static int cardinality(final byte b) {
        checkIUPACBase(b);
        return Integer.bitCount((int) b);
    }

    /**
     * Total number of <i>unambiguous</i> strings this array codes
     */
    public static int cardinality(final byte[] dna) {
        checkNotNull(dna, "Null string passed");
        if (dna.length == 0) return 0;
        int result = 1;
        for (final byte b : dna) {
            result *= cardinality(b);
        }
        return result;
    }

    public static byte[] disambiguate(final byte b) {
        checkIUPACBase(b);
        byte[] result = new byte[cardinality(b)];
        final byte T = 1, A = 2, C = 4, G = 8;
        int i = 0;
        if ((b & A) > 0) result[i++] = 'A';
        if ((b & C) > 0) result[i++] = 'C';
        if ((b & G) > 0) result[i++] = 'G';
        if ((b & T) > 0) result[i++] = 'T';
        return result;
    }

    public static byte[][] disambiguate(final byte[] dna) {
        final int card = cardinality(dna);
        byte[][] result = new byte[card][];
        int k = 0;
        Stack<byte[]> s1 = new Stack<byte[]>(), s2;
        s1.push(new byte[dna.length]);

        for (int i = 0; i < dna.length; i++) {
            s2 = new Stack<byte[]>();
            final byte[] d = disambiguate(dna[i]);
            //	    System.out.format("%s -> %s\n",
            //			      Character.toString((char)unpackByte(dna[i])),
            //			      new String(d));
            while (!s1.empty()) {
                byte[] cur = s1.pop();
                for (int j = 0; j < d.length; j++) {
                    byte[] b = Arrays.copyOf(cur, cur.length);
                    b[i] = d[j];
                    s2.push(b);
                }
            }
            s1 = s2;
        }
        for (int i = 0; i < result.length; i++) {
            checkState(!s1.empty(), "Stack of ambiguous chars empty!");
            result[i] = s1.pop();
        }
        return result;
    }

    public static long packKmer(final byte[] packed, final int start,
                                final int length) {
        checkNotNull(packed);
        checkArgument(start >= 0, "Negative start: %d", start);
        checkArgument(length >= start, "length %d < start %d", length, start);
        long result = 0;
        for (int i = start; i < start + length; i++) {
            result = (result << BITS_PER_BASE) + packed[i];
        }
        return result;
    }

    public static long packKmer(final byte[] packed) {
        return packKmer(packed, 0, packed.length);
    }

    public static byte[] unpackKmer(final long value) {
        final int highest = (int) Long.numberOfTrailingZeros(Long.highestOneBit(value));
        final int n = highest / BITS_PER_BASE + ((highest % BITS_PER_BASE) > 0 ? 1 : 0);
        byte[] result = new byte[n];
        // Fill right-to-left
        for (int i = 0; i < n; i++) {
            // This nastiness just extracts a specific base (4 bits) by
            // masking and shifting.
            final int bitOffset = i * BITS_PER_BASE;
            result[n - 1 - i] = (byte) ((value & (SINGLE_BASE_MASK << bitOffset)) >> bitOffset);
        }
        return result;
    }

    private static void checkIUPACBase(byte b) {
        checkArgument(b < 16, "Invalid encoded IUPAC: %d", b);
        checkArgument(b >= 0, "Invalid encoded IUPAC: %d", b);
    }
}
