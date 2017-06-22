// Copyright (c) 2017 Peter A. Audano III
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published
// by the Free Software Foundation; either version 3 of the License or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Library General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with this program; see the file COPYING.LESSER.  If not, see
// <http://www.gnu.org/licenses/>

package edu.gatech.kestrel.util;

import java.util.Arrays;

/**
 * A collection of numerical utility methods for Kestrel.
 */
public final class NumberUtil {
	
	/** Floating-point values within <code>ZERO_RANGE</code> of 0.0 are considered 0. */
	public static final float ZERO_RANGE = 0.0001F;
	
	/** <code>ZERO_RANGE</code> negated. Saved to speed functions in this class. */
	public static final float N_ZERO_RANGE = -ZERO_RANGE;
	
	/**
	 * Hidden constructor (all methods are static).
	 */
	protected NumberUtil() { 
		return;
	}
	
	/**
	 * Determine if a floating-point value is within <code>ZERO_RANGE</code> of zero,
	 * inclusive (including when <code>f</code> is equal to <code>ZERO_RANGE</code>).
	 * Because of rounding concerns, checking a floating-point number for equality to
	 * zero should be done by a range around zero, which this function does.
	 * 
	 * @param f Floating-point value to check.
	 * 
	 * @return <code>true</code> if <code>f</code> is within <code>ZERO_RANGE</code> of
	 *   zero.
	 */
	public static final boolean isZero(Float f) {
		
		return f <= ZERO_RANGE && f >= N_ZERO_RANGE;
	}
	
	/**
	 * Determine if a floating-point value is within <code>ZERO_RANGE</code> of zero
	 * exclusive (not including when <code>f</code> is equal to <code>ZERO_RANGE</code>).
	 * Because of rounding concerns, checking a floating-point number for equality to
	 * zero should be done by a range around zero, which this function does.
	 * 
	 * @param f Floating-point value to check.
	 * 
	 * @return <code>true</code> if <code>f</code> is within <code>ZERO_RANGE</code> of
	 *   zero.
	 */
	public static final boolean isNotZero(Float f) {
		
		return f > ZERO_RANGE || f < N_ZERO_RANGE;
	}
	
	/**
	 * Get a quantile of k-mer count differences. This is calculated by taking an array
	 * of k-mers and subtracting neighboring k-mers to get a distribution of the
	 * differences between neighbors. The quantile is then calculated on that difference.
	 * 
	 * @param count K-mer count array. If empty or one element, <code>0</code> is returned.
	 *   If it is two elements element, then the difference is returned.
	 * @param quantile Quantile of <code>count</code> to find. Must be between <code>0</code>
	 *   and <code>1</code> (exclusive, cannot be <code>0</code> or <code>1</code>).
	 * 
	 * @return A k-mer count at <code>quantile</code> with the fractional part truncated.
	 * 
	 * @throws NullPointerException If <code>count</code> is <code>null</code>.
	 * @throws IllegalArgumentException If <code>quantile</code> is not strictly between
	 *   <code>0</code> and <code>1</code> (exclusive on both ends).
	 */
	public final static int countDiffQuantile(int[] count, double quantile)
			throws NullPointerException, IllegalArgumentException {
		
		// Declarations
		int lastCount;
		int thisCount;
		int thisCountDiff;
		int loc;
		double offset;
		
		int[] countDiff;
		int nCount;       // Number of elements in count while getting differences, and one less than
		                  // the length of countDiff while computing quantiles
		
		// Check arguments
		if (count == null)
			throw new NullPointerException("Count array is null");
		
		if (quantile <= 0.0F || quantile >= 1.0F)
			throw new IllegalArgumentException("q is out of range (0, 1): " + quantile);
		
		nCount = count.length;
		
		if (nCount < 3) {
			if (nCount < 2)
				return 0;
			
			return Math.abs(count[1] - count[0]);
		}
		
		countDiff = new int[nCount - 1];
		
		// Convert to count differences
		lastCount = count[0];
		
		for (int index = 1; index < nCount; ++index) {
			thisCount = count[index];
			thisCountDiff = lastCount - thisCount;
			
			if (thisCountDiff < 0)
				thisCountDiff = -thisCountDiff;
			
			countDiff[index - 1] = thisCountDiff;
			lastCount = thisCount;
		}
		
		nCount -= 2;  // One less difference than k-mer counts (simplifies quantile calculations)
		
		// Sort
		Arrays.sort(countDiff);
		
		// Calculate quantile
		loc = (int) (nCount * quantile);
		offset = (nCount * quantile) - loc;
		
		// Return quantile
		return (int) (countDiff[loc] * (1 - offset) + countDiff[loc + 1] * offset);
	}
}
