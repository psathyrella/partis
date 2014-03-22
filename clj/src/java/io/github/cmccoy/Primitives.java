package io.github.cmccoy;

import static com.google.common.base.Preconditions.*;

/**
 * Functions for working with primitive arrays.
 */
public class Primitives {

    /**
     * Converts a byte array to a double array. Simple.
     */
    public static double[] bytesToDoubles(final byte[] arr) {
        checkNotNull(arr);
        double[] result = new double[arr.length];
        for(int i = 0; i < result.length; i++)
            result[i] = (double)arr[i];
        return result;
    }
}
