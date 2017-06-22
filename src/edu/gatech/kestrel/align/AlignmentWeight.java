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

package edu.gatech.kestrel.align;

import edu.gatech.kestrel.util.NumberUtil;

/**
 * A structure of sequence alignment weights. The alignment controls assembly over active regions
 * and influences variant calling.
 */
public class AlignmentWeight {
	
	/** Weight for an aligned pair of matching bases (<code>match</code> &gt; 0). */
	public final float match;
	
	/** Weight for an aligned pair of mis-matched bases (<code>mismatch</code> &lt; 0). */
	public final float mismatch;
	
	/** Weight for opening a gap (<code>gapOpen</code> &le; 0). */
	public final float gapOpen;
	
	/** Weight for extending a gap (<code>gapExtend</code> &lt; 0). */
	public final float gapExtend;
	
	/** Initial alignment score. If <code>0</code>, it defaults to <code>match * kSize</code>. */
	public final float initScore;
	
	/** <code>gapOpen + gapExtend</code>. Calculated once to avoid recomputing during alignments. */
	public final float newGap;
	
	/** Default match weight. */
	public static final float DEFAULT_MATCH = 10.0F;
	
	/** Default mismatch weight. */
	public static final float DEFAULT_MISMATCH = -10.0F;
	
	/** Default gap-open weight. */
	public static final float DEFAULT_GAP_OPEN = -40.0F;
	
	/** Default gap-extend weight. */
	public static final float DEFAULT_GAP_EXTEND = -4.0F;
	
	/** Default initial score. */
	public static final float DEFAULT_INIT_SCORE = 0.0F;
	
	
	//
	// Private constructor (use "get()" factory methods
	//
	
	/**
	 * Create a set of alignment weights.
	 * 
	 * @param match Weight for an aligned pair of matching bases.
	 * @param mismatch Weight for an aligned pair of mis-matched bases.
	 * @param gapOpen Weight for opening a gap.
	 * @param gapExtend Weight for extending a gap.
	 * @param initScore Initial alignment score.
	 */
	private AlignmentWeight(float match, float mismatch, float gapOpen, float gapExtend, float initScore) {
		
		// Set fields
		this.match = match;
		this.mismatch = mismatch;
		this.gapOpen = gapOpen;
		this.gapExtend = gapExtend;
		this.initScore = initScore;
		
		newGap = gapOpen + gapExtend;
		
		return;
	}
	
	/**
	 * Get a string representation of this weight structure.
	 * 
	 * @return String representation of this weight structure.
	 */
	@Override
	public String toString() {
		return String.format("(%.4f, %.4f, %.4f, %.4f, %.4f)", match, mismatch, gapOpen, gapExtend, initScore);
	}
	
	/**
	 * Get a string representation of this weight structure.
	 * 
	 * @param precision Floating-point precision of the weights.
	 * 
	 * @return String representation of this weight structure.
	 * 
	 * @throws IllegalArgumentException If <code>precision</code> is out of range [0, 10].
	 */
	public String toString(int precision)
			throws IllegalArgumentException {
		
		if (precision < 0 || precision > 10)
			throw new IllegalArgumentException("Precision out of range [0, 10]: " + precision);
		
		String fmtString = String.format("(%%.%df, %%.%df, %%.%df, %%.%df, %%.%df)", precision, precision, precision, precision, precision);
		
		return String.format(fmtString, match, mismatch, gapOpen, gapExtend, initScore);
	}
	
	
	//
	// AlignmentWeight factory methods
	//
	
	/**
	 * Get a <code>AlignmentWeight</code> object from a string representing the alignment weight as a
	 * comma-separated vector of values. The order of the values is match, mismatch, gap-open, and gap-extend.
	 * If less than four elements are defined or elements in the matrix are blank, they are assigned their default
	 * values (<code>DEFAULT_MATCH</code>, <code>DEFAULT_MISMATCH</code>, <code>DEFAULT_GAP_OPEN</code>,
	 * and <code>DEFAULT_GAP_EXTEND</code>, respectively).
	 * <p/>
	 * The vector may be surrounded by matching parenthesis, angle-braces, square-braces, curly-braces, or
	 * it may not be surrounded at all. Each value may be any valid floating-point number, including ones
	 * in exponential form (e.g. 1.0e2), as well as hexadecimal (0xD2) and octal (0322).
	 * 
	 * @param weightString String of weights. If <code>null</code> or empty, a default set of sequence
	 *   weights is returned.
	 * 
	 * @return A configured <code>AlignmentWeight</code> object.
	 * 
	 * @throws IllegalArgumentException If <code>weightString</code> is not properly formatted or
	 *   contains values that cannot be converted to a floating-point number.
	 */
	public static AlignmentWeight get(String weightString)
			throws IllegalArgumentException {
		
		// Weights
		float match = DEFAULT_MATCH;
		float mismatch = DEFAULT_MISMATCH;
		float gapOpen = DEFAULT_GAP_OPEN;
		float gapExtend = DEFAULT_GAP_EXTEND;
		float initScore = DEFAULT_INIT_SCORE;
		
		// Check arguments
		if (weightString == null)
			weightString = "";
		
		weightString = weightString.trim();
		
		if (weightString.isEmpty())
			return AlignmentWeight.get();
		
		// Remove surrounding braces, parens, etc
		int strLength = weightString.length();
		
		if (strLength >= 2) {
			char[] strArray = weightString.toCharArray();
			
			char firstChar = strArray[0];
			char lastChar = strArray[strArray.length - 1];
			
			switch (firstChar) {
			case '(':
				if (lastChar == ')')
					weightString = weightString.substring(1, strLength - 1);
				else
					throw new IllegalArgumentException("Mis-matched parenthesis in alignment weights: \"" + weightString + "\"");
				
				break;
				
			case '<':
				if (lastChar == '>')
					weightString = weightString.substring(1, strLength - 1);
				else
					throw new IllegalArgumentException("Mis-matched angle-braces in alignment weights: \"" + weightString + "\"");
				
				break;
			
			case '[':
				if (lastChar == ']')
					weightString = weightString.substring(1, strLength - 1);
				else
					throw new IllegalArgumentException("Mis-matched square-brackets in alignment weights: \"" + weightString + "\"");
				
				break;
			
			case '{':
				if (lastChar == '}')
					weightString = weightString.substring(1, strLength - 1);
				else
					throw new IllegalArgumentException("Mis-matched curly-braces in alignment weights: \"" + weightString + "\"");
				
				break;
			
			default:
				if (lastChar == ')' || lastChar == '>' || lastChar == ']' || lastChar == '}')
					throw new IllegalArgumentException("Closing brace or parenthesis in alignment weights with no opening brace or parenthesis: \"" + weightString + "\"");
			}
		}
		
		// Split
		String[] weightTok = weightString.split("\\s*,\\s*");
		
		if (weightTok.length > 5)
			throw new IllegalArgumentException("Weight vector is more than 5 comma-separated values: Values = " + weightTok.length);
		
		// Set weights
		if (weightTok.length > 0 && ! weightTok[0].isEmpty())
			match = normalizeMatch(stringToFloat(weightTok[0]));  // throws IllegalArgumentException
		
		if (weightTok.length > 1 && ! weightTok[1].isEmpty())
			mismatch = normalizeMismatch(stringToFloat(weightTok[1]));  // throws IllegalArgumentException
		
		if (weightTok.length > 2 && ! weightTok[2].isEmpty())
			gapOpen = normalizeGapOpen(stringToFloat(weightTok[2]));  // throws IllegalArgumentException
		
		if (weightTok.length > 3 && ! weightTok[3].isEmpty())
			gapExtend = normalizeGapExtend(stringToFloat(weightTok[3]));  // throws IllegalArgumentException
		
		if (weightTok.length > 4 && ! weightTok[4].isEmpty())
			initScore = normalizeInitScore(stringToFloat(weightTok[4]));  // throws IllegalArgumentException
		
		// Return weights
		return new AlignmentWeight(match, mismatch, gapOpen, gapExtend, initScore);
	}
	
	/**
	 * Create a set of alignment weights.
	 * 
	 * @param match Weight for an aligned pair of matching bases (<code>match</code> &gt; 0).
	 *   If negative, this weight is automatically converted to a positive number.
	 * @param mismatch Weight for an aligned pair of mis-matched bases (<code>mismatch</code> &lt; 0).
	 *   If positive, this weight is automatically converted to a negative number.
	 * @param gapOpen Weight for opening a gap (<code>gapOpen</code> &lt;= 0).
	 *   If positive, this weight is automatically converted to a negative number.
	 * @param gapExtend Weight for extending a gap (<code>gapExtend</code> &lt; 0).
	 *   If positive, this weight is automatically converted to a negative number.
	 * @param initScore Initial alignment score (<code>initScore</code> &gt;= 0).
	 *   If zero, the initial score is <code>match</code> multiplied by the k-mer size
	 *   to obtain the initial score.
	 *   
	 * @return A set of sequence weights.
	 * 
	 * @throws IllegalArgumentException If <code>match</code>, <code>mismatch</code>, or <code>gapExtend</code>
	 *   is within <code>ZERO_RANGE</code> of <code>0</code>. 
	 */
	public static AlignmentWeight get(float match, float mismatch, float gapOpen, float gapExtend, float initScore)
			throws IllegalArgumentException {
		
		return new AlignmentWeight(
				normalizeMatch(match),  // throws IllegalArgumentException
				normalizeMismatch(mismatch),  // throws IllegalArgumentException
				normalizeGapOpen(gapOpen),  // throws IllegalArgumentException
				normalizeGapExtend(gapExtend),  // throws IllegalArgumentException
				normalizeInitScore(initScore)  // throws IllegalArgumentException
		);
	}
	
	/**
	 * Get a set of default sequence weights.
	 * 
	 * @return A set of default sequence weights.
	 */
	public static AlignmentWeight get() {
		return new AlignmentWeight(DEFAULT_MATCH, DEFAULT_MISMATCH, DEFAULT_GAP_OPEN, DEFAULT_GAP_EXTEND, DEFAULT_INIT_SCORE);
	}
	
	
	//
	// Sequence weight factory mutator methods (returns an altered copy)
	//
	
	/**
	 * Get a new sequence weight container with a set match weight. All other weights are
	 * unaltered from this object.
	 * 
	 * @param match Weight for aligned bases that match.
	 * 
	 * @return A reconfigured <code>AlignmentWeight</code> object.
	 * 
	 * @throws IllegalArgumentException If <code>match</code> is within
	 *   <code>ZERO_RANGE</code> of <code>0</code>.
	 */
	public AlignmentWeight getWithMatch(float match)
			throws IllegalArgumentException {
		
		return new AlignmentWeight(normalizeMatch(match), mismatch, gapOpen, gapExtend, initScore);
	}
	
	/**
	 * Get a new sequence weight container with a set mismatch weight. All other weights are
	 * unaltered from this object.
	 * 
	 * @param mismatch Weight for aligned bases that do not match.
	 * 
	 * @return A reconfigured <code>AlignmentWeight</code> object.
	 * 
	 * @throws IllegalArgumentException If <code>mismatch</code> is within
	 *   <code>ZERO_RANGE</code> of <code>0</code>.
	 */
	public AlignmentWeight getWithMismatch(float mismatch)
			throws IllegalArgumentException {
		
		return new AlignmentWeight(match, normalizeMismatch(mismatch), gapOpen, gapExtend, initScore);
	}
	
	/**
	 * Get a new sequence weight container with a set gap-open weight. All other weights are
	 * unaltered from this object.
	 * 
	 * @param gapOpen Weight for opening a gap.
	 * 
	 * @return A reconfigured <code>AlignmentWeight</code> object.
	 */
	public AlignmentWeight getWithGapOpen(float gapOpen) {
		
		return new AlignmentWeight(match, mismatch, normalizeGapOpen(gapOpen), gapExtend, initScore);
	}
	
	/**
	 * Get a new sequence weight container with a set gap-extend weight. All other weights are
	 * unaltered from this object.
	 * 
	 * @param gapExtend Weight for extending a gap.
	 * 
	 * @return A reconfigured <code>AlignmentWeight</code> object.
	 * 
	 * @throws IllegalArgumentException If <code>gapExtend</code> is within
	 *   <code>ZERO_RANGE</code> of <code>0</code>.
	 */
	public AlignmentWeight getWithGapExtend(float gapExtend)
			throws IllegalArgumentException {
		
		return new AlignmentWeight(match, mismatch, gapOpen, normalizeGapExtend(gapExtend), initScore);
	}
	
	/**
	 * Get a new sequence weight container with a set initial score. All other weights are
	 * unaltered from this object.
	 * 
	 * @param initScore Initial score. If negative, the absolute value is used.
	 * 
	 * @return A reconfigured <code>AlignmentWeight</code> object.
	 */
	public AlignmentWeight getWithInitialScore(float initScore) {
		
		return new AlignmentWeight(match, mismatch, gapOpen, gapExtend, normalizeInitScore(initScore));
	}
	
	
	//
	// Score utilities
	//
	
	/**
	 * Get the initial score for alignment starts.
	 * 
	 * @param kSize K-mer size. If the initial score is 0, then <code>match</code> is multiplied
	 *   by this number to obtain the initial score.
	 * 
	 * @return Initial score.
	 * 
	 * @throws IllegalArgumentException If <code>kSize</code> is less than <code>1</code>.
	 */
	public float getInitialScore(int kSize)
			throws IllegalArgumentException {
		
		if (kSize < 1)
			throw new IllegalArgumentException("Cannot get initial score for k-mer size less than 1: " + kSize);
		
		if (NumberUtil.isZero(initScore))
			return match * kSize;
		
		return initScore;
	}
	
	/**
	 * Maximum size of a gap assuming no aligned bases in the active region. This is used as a dynamic
	 * limit by other classes.
	 * 
	 * @param kSize K-mer size for determining the maximum initial score if not explicitly set.
	 * 
	 * @return The maximum gap size assuming no aligned bases in the active region (besides anchor k-mers).
	 */
	public int getMaxExclusiveGapSize(int kSize) {
		
		float initScore = (int) getInitialScore(kSize);
		
		if (initScore > gapOpen)
			return (int) ((initScore + gapOpen) / -gapExtend);
		
		return 0;
	}
	
	//
	// Check and normalize weights
	//
	
	/**
	 * Check weight of aligned matching bases and return a positive value.
	 * 
	 * @param match Weight to normalize.
	 * 
	 * @return A positive value of <code>match</code>.
	 * 
	 * @throws IllegalArgumentException If <code>match</code> is within <code>ZERO_RANGE</code>
	 *   of <code>0</code>.
	 */
	private static float normalizeMatch(float match)
			throws IllegalArgumentException {
		
		// Check match for 0
		if (NumberUtil.isZero(match))
			throw new IllegalArgumentException(String.format(
					"Cannot set weight for matching bases: Weight is zero or too close to zero (%.6f < %.6f < %.6f)",
					NumberUtil.N_ZERO_RANGE,
					match,
					NumberUtil.ZERO_RANGE
					));
		
		// Return positive value
		return Math.abs(match); 
	}
	
	/**
	 * Check weight of aligned mis-matching bases and return a negative value.
	 * 
	 * @param mismatch Weight to normalize.
	 * 
	 * @return A negative value of <code>mismatch</code>.
	 * 
	 * @throws IllegalArgumentException If <code>mismatch</code> is within <code>ZERO_RANGE</code>
	 *   of <code>0</code>.
	 */
	private static float normalizeMismatch(float mismatch)
			throws IllegalArgumentException {
		
		// Check mismatch for 0
		if (NumberUtil.isZero(mismatch))
			throw new IllegalArgumentException(String.format(
					"Cannot set weight for mismatched bases: Weight is zero or too close to zero (%.6f < %.6f < %.6f)",
					NumberUtil.ZERO_RANGE,
					mismatch,
					NumberUtil.N_ZERO_RANGE
					));
		
		// Return positive value
		return -Math.abs(mismatch); 
	}
	
	/**
	 * Check weight for opening a gap and return a negative value.
	 * 
	 * @param gapOpen Weight to normalize.
	 * 
	 * @return A negative value of <code>gapOpen</code>.
	 */
	private static float normalizeGapOpen(float gapOpen) {
		
		// Return positive value
		return -Math.abs(gapOpen); 
	}
	
	/**
	 * Check weight extending a gap and return a negative value.
	 * 
	 * @param gapExtend Weight to normalize.
	 * 
	 * @return A negative value of <code>gapExtend</code>.
	 * 
	 * @throws IllegalArgumentException If <code>gapExtend</code> is within <code>ZERO_RANGE</code>
	 *   of <code>0</code>.
	 */
	private static float normalizeGapExtend(float gapExtend)
			throws IllegalArgumentException {
		
		// Check mismatch for 0
		if (NumberUtil.isZero(gapExtend))
			throw new IllegalArgumentException(String.format(
					"Cannot set weight for mismatched bases: Weight is zero or too close to zero (%.6f < %.6f < %.6f)",
					NumberUtil.ZERO_RANGE,
					gapExtend,
					NumberUtil.N_ZERO_RANGE
					));
		
		// Return positive value
		return -Math.abs(gapExtend); 
	}
	
	/**
	 * Check initial score and return a non-negative number.
	 * 
	 * @param initScore Initial score to normalize.
	 * 
	 * @return A non-negative value of <code>initScore</code>.
	 */
	private static float normalizeInitScore(float initScore) {
		
		// Check match for 0
		if (NumberUtil.isZero(initScore))
			initScore = 0.0F;
		
		// Return positive value
		return Math.abs(initScore);
	}
	
	
	//
	// String processing methods
	//
	
	/**
	 * Convert a weight from a string to a floating-point value. This method supports the typical
	 * string representations of floating point numbers, including exponential form, as well
	 * as octal and hexadecimal integer formats.
	 *  
	 * @param weightString Weight string.
	 * 
	 * @return Floating-point value of the weight.
	 * 
	 * @throws IllegalArgumentException If <code>weightString</code> does not represent a floating
	 *   point number or an integer in hex or octal form.
	 */
	private static float stringToFloat(String weightString)
			throws IllegalArgumentException {
		
		// Check arguments
		assert (weightString != null) :
			"Cannot convert weight: null";
		
		weightString = weightString.trim();
		
		// Convert (float)
		try {
			return Float.parseFloat(weightString);
			
		} catch (NumberFormatException ex) {
			// Ignore, move to integer parsing
		}
		
		// Convert (integer)
		try {
			return Integer.decode(weightString);
			
		} catch (NumberFormatException ex2) {
			// Ignore, throws an exception next
		}
		
		// No valid formats
		throw new IllegalArgumentException("Cannot convert weight to a numeric value: " + weightString);
	}
	
	
	//
	// Equality methods
	//
	
	/**
	 * Determine if the weights in this object are equal to another set of sequence weights. If the
	 * floating point values in the objects (match, mismatch, gapOpen, and gapExtend) are within
	 * <code>ZERO_RANGE</code> of the corresponding value in another set of weights, then they are
	 * considered equal.
	 * 
	 * @param other Object to check. 
	 * 
	 * @return <code>true</code> if the objects represent equal weights, and <code>false</code> if at least
	 *   one weight is not within <code>ZERO_RANGE</code> in the other object. <code>false</code> is
	 *   returned if <code>other</code> is <code>null</code> or not an instance of <code>AlignmentWeight</code>.
	 */
	@Override
	public boolean equals(Object other) {
		
		if (other == null || ! (other instanceof AlignmentWeight))
			return false;
		
		AlignmentWeight otherWeight = (AlignmentWeight) other;
		
		return
				NumberUtil.isZero(match - otherWeight.match) &&
				NumberUtil.isZero(mismatch - otherWeight.mismatch) &&
				NumberUtil.isZero(gapOpen - otherWeight.gapOpen) &&
				NumberUtil.isZero(gapExtend - otherWeight.gapExtend) &&
				NumberUtil.isZero(initScore - otherWeight.initScore)
		;
	}
	
	/**
	 * Determine if this set of weights contains all default values. This object is checked against one
	 * with all default values. For details on how comparisons are made, see <code>equals()</code>.
	 * 
	 * @return <code>true</code> if this set of weights contains default values, and <code>false</code> if
	 *   at least one weight is not default.
	 * 
	 * @see #equals(Object)
	 */
	public boolean isDefault() {
		return this.equals(AlignmentWeight.get());
	}
	
	/**
	 * Get a hash code for this weight object.
	 * 
	 * @return A hash code for this object.
	 */
	@Override
	public int hashCode() {
		
		return Float.valueOf(match).hashCode() ^
				Float.valueOf(mismatch).hashCode() ^
				Float.valueOf(gapOpen).hashCode() ^
				Float.valueOf(gapExtend).hashCode() ^
				Float.valueOf(initScore).hashCode();
	}
}
