package io.github.cmccoy.dna;

import static com.google.common.base.Preconditions.checkArgument;
import static com.google.common.base.Preconditions.checkNotNull;

/**
 * Utilities for working with IUPAC bases, encoded as bytes.
 * The first 4 bits of each byte are used, in order ACGT.
 * Ambiguous bases are represented as multiple "ON" bits.
 */
public class IUPACUtils {
    private final static byte[] LOOKUP =
            {'U', 'T', 'A', 'W', 'C', 'Y', 'M',
                    'H', 'G', 'K', 'R', 'D', 'S', 'B', 'V', 'N'};

    /**
     * Pack a single byte.
     *
     * @param b An IUPAC character
     * @return byte encoded.
     */
    public static byte packByte(final byte b) {
        switch ((char) b) {
            case 'U':
                return 0;
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
                throw new IllegalArgumentException(String.format("Unknown base: %s", b));
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
        return LOOKUP[s];
    }

    public static byte[] unpackBytes(final byte[] packed) {
        checkNotNull(packed);
        byte[] result = new byte[packed.length];
        for (int i = 0; i < packed.length; i++)
            result[i] = unpackByte(packed[i]);
        return result;
    }

    public static boolean areCompatible(final byte s1, final byte s2) {
        return (s1 & s2) > 0;
    }

    public static boolean areCompatible(final byte[] s1, final byte[] s2) {
        return areCompatible(s1, s2, 0, s1.length);
//    checkNotNull(s1);
//    checkNotNull(s2);
//    checkArgument(s1.length == s2.length,
//        "Lengths differ: %d != %d",
//        s1.length,
//        s2.length);
//    for(int i = 0; i < s1.length; i++) {
//      if(!areCompatible(s1[i], s2[i]))
//        return false;
//    }
//    return true;
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


    private static int BITS_PER_BASE = 4;
    private static int SINGLE_BASE_MASK = 0xf;

    public static long packKmer(final byte[] packed, final int start, final int length) {
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
        final int highest = (int)Long.numberOfTrailingZeros(Long.highestOneBit(value));
        final int n = highest / BITS_PER_BASE + ((highest % BITS_PER_BASE) > 0 ? 1 : 0);
        byte[] result = new byte[n];
        // Fill right-to-left
        for(int i = 0; i < n; i++) {
            // This nastiness just extracts a specific base (4 bits) by
            // masking and shifting.
            final int bitOffset = i * BITS_PER_BASE;
            result[n - 1 - i] = (byte) ((value & (SINGLE_BASE_MASK << bitOffset)) >> bitOffset);
        }
        return result;
    }
}
