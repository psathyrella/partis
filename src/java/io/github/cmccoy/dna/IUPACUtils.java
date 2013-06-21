package io.github.cmccoy.dna;

import static com.google.common.base.Preconditions.*;

public class IUPACUtils {
  private final static byte[] LOOKUP =
  {'U', 'T', 'A', 'W', 'C', 'Y', 'M',
    'H', 'G', 'K', 'R', 'D', 'S', 'B', 'V', 'N'};

  public static short packByte(final byte b) {
    switch ((char) b) {
      case 'U': return 0;
      case 'T': return 1;
      case 'A': return 2;
      case 'W': return 3;
      case 'C': return 4;
      case 'Y': return 5;
      case 'M': return 6;
      case 'H': return 7;
      case 'G': return 8;
      case 'K': return 9;
      case 'R': return 10;
      case 'D': return 11;
      case 'S': return 12;
      case 'B': return 13;
      case 'V': return 14;
      case 'N': return 15;
      default:
                throw new IllegalArgumentException(String.format("Unknown base: %s", b));
    }
  }

  public static short[] packBytes(final byte[] bytes) {
    checkNotNull(bytes);
    short[] result = new short[bytes.length];
    for(int i = 0; i < bytes.length; ++i) {
      result[i] = packByte(bytes[i]);
    }
    return result;
  }

  public static byte unpackShort(final short s) {
    return LOOKUP[s];
  }

  public static byte[] unpackShorts(final short[] packed) {
    checkNotNull(packed);
    byte[] result = new byte[packed.length];
    for(int i = 0; i < packed.length; i++)
      result[i] = unpackShort(packed[i]);
    return result;
  }

  public static boolean areCompatible(final short s1, final short s2) {
    return (s1 & s2) > 0;
  }

  public static boolean areCompatible(final short[] s1, final short[] s2) {
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

  public static boolean areCompatible(final short[] s1, final short[] s2,
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
    for(int i = 0; i < length; ++i) {
      if(!areCompatible(s1[start + i], s2[i]))
        return false;
    }
    return true;
  }

  /**
   * Tests whether s2 is a subset of s1
   */
  public static boolean isSubset(short s1, short s2) {
    return 0 == (s2 & ~s1);
  }

  public static boolean isSubset(final short[] s1, final short[] s2,
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
    for(int i = 0; i < length; ++i) {
      if(!isSubset(s1[start + i], s2[i]))
        return false;
    }
    return true;
  }

  public static boolean isSubset(final short[] s1, final short[] s2) {
    return isSubset(s1, s2, 0, s1.length);
  }

}
