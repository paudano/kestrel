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

package edu.gatech.kestrel.activeregion;

import java.util.Arrays;

import edu.gatech.kanalyze.KAnalyzeConstants;
import edu.gatech.kanalyze.util.Base;
import edu.gatech.kanalyze.util.kmer.KmerUtil;
import edu.gatech.kestrel.refreader.ReferenceRegion;
import edu.gatech.kestrel.util.SystemUtil;

/**
 * Defines an active region within a reference sequence. Active regions are regions that
 * may contain variants.
 * <p/>
 * Active regions are defined by the k-mers of a reference sequence. The start and end k-mer index
 * of the active region are k-mers that are present in the reference sequence and the sample
 * (the k-mer counts from the sample are high). The k-mers between these indices (but not
 * including the indices) are suspected to contain at least one variation.
 */
public class ActiveRegion implements Comparable<ActiveRegion>{
	
	/** Reference region where this k-mer occurs. */
	public final ReferenceRegion refRegion;
	
	/**
	 * Location of <code>refSequence.sequence</code> where the active region starts.
	 */
	public final int startIndex;
	
	/**
	 * Location of <code>refSequence.sequence</code> where the active region ends.
	 */
	public final int endIndex;
	
	/**
	 * In an array of k-mers over the reference region, this is the index of the first k-mer in the active
	 * region.
	 */
	public final int startKmerIndex;
	
	/**
	 * In an array of k-mers over the reference region, this is the index of the last k-mer in the active
	 * region.
	 */
	public final int endKmerIndex;
	
	/**
	 * Set to <code>true</code> if this active region extends to the left end of the
	 * reference sequence. This is set only when the active region has no left-end
	 * k-mer anchor. This and <code>rightEnd</code> will never both be <code>true</code>.
	 */
	public final boolean leftEnd;
	
	/**
	 * Set to <code>true</code> if this active region extends to the right end of the
	 * reference sequence. This is set only when the active region has no right-end
	 * k-mer anchor. This and <code>leftEnd</code> will never both be <code>true</code>.
	 */
	public final boolean rightEnd;
	
	/** Summary statistics over this region&apos; wild-type sequence. */
	public final RegionStats stats;
	
	/**
	 * K-mer bordering the left end of this active region or <code>null</code> if <code>leftEnd</code>
	 * is <code>true</code>.
	 */
	private final int[] leftEndKmer;
	
	/**
	 * K-mer bordering the right end of this active region or <code>null</code> if <code>rightEnd</code>
	 * is <code>true</code>.
	 */
	private final int[] rightEndKmer;
	
	/** Translates bytes to bases. */
	private static final Base[] BYTE_TO_BASE = KAnalyzeConstants.getByteToBaseArray();
	
	
	/**
	 * Create a new active region.
	 * 
	 * @param refSequence Reference sequence.
	 * @param startKmerIndex Index of the k-mer array over <code>refSequence.sequence</code> where the
	 *   active region starts, or <code>-1</code> to indicate that the active region reaches the
	 *   left end of the sequence.
	 * @param endKmerIndex Index of <code>refSequence.sequence</code> where the active region ends,
	 *   or <code>-1</code> to indicate that the active region reaches the right end of the sequence.
	 * @param count Array of k-mer counts over the k-mers in <code>refSequence</code>.
	 * @param kUtil K-mer utility.
	 *   
	 * @throws NullPointerException If <code>refSequence</code>, <code>count</code>, or
	 * <code>kUtil</code> is <code>null</code>.
	 * @throws IllegalArgumentException If <code>startKmerIndex</code> and <code>endKmerIndex</code> do not
	 *   define an active region contained within this reference sequence or if both are <code>-1</code>, or
	 *   if the k-mers bordering this active region contain ambiguous bases.
	 */
	public ActiveRegion(ReferenceRegion refSequence, int startKmerIndex, int endKmerIndex, int[] count, KmerUtil kUtil)
			throws NullPointerException, IllegalArgumentException {
		
		// Check arguments
		if (refSequence == null)
			throw new NullPointerException("Cannot create active region on reference sequence: null");
		
		if (kUtil == null)
			throw new NullPointerException("Cannot create active region with k-mer utility: null");
		
		if (count == null)
			throw new NullPointerException("K-mer count array is null");
		
		if (startKmerIndex < -1)
			throw new IllegalArgumentException("Cannot create active region: Start k-mer index is less than -1: " + startKmerIndex);
		
		if (endKmerIndex < -1 || endKmerIndex >= count.length)
			throw new IllegalArgumentException(String.format("Cannot create active region: End k-mer index is not between -1 and the last k-emr in the reference (%d): %d", count.length - 1, endKmerIndex));
		
		if (endKmerIndex > refSequence.size - kUtil.kSize)
			throw new IllegalArgumentException(String.format("Cannot create active region: The end k-mer at index (%d) extends past the right end of the reference sequence (length %d) for the given k-mer size (%d)", endKmerIndex, refSequence.size, kUtil.kSize));
		
		if (startKmerIndex >= 0 && endKmerIndex >= 0) {
			if (startKmerIndex >= endKmerIndex)
				throw new IllegalArgumentException(String.format("Cannot create active region: Start k-mer index (%d) must come before end k-mer index (%d)", startKmerIndex, endKmerIndex));
			
		} else  {
			if (startKmerIndex == -1 && endKmerIndex == -1)
				throw new IllegalArgumentException("Active region may not span the entire reference sequence (start and end k-mer indices are both -1)");
		}
		
		SystemUtil.checkKmerSize(kUtil.kSize, null);  // throws IllegalArgumentException
		
		// Set reference sequence
		this.refRegion = refSequence;
		
		// Translate for reference end cases
		if (startKmerIndex == -1) {
			this.leftEnd = true;
			
			this.startIndex = 0;
			this.startKmerIndex = 0;
			
			leftEndKmer = null;
			
		} else {
			this.leftEnd = false;
			
			this.startIndex = startKmerIndex;
			this.startKmerIndex = startKmerIndex;
			
			leftEndKmer = getKmer(this.startIndex, kUtil);
			
			if (leftEndKmer == null)
				throw new IllegalArgumentException(String.format("Reference sequence (%s) contains an ambiguous base in the left-end k-mer at index %d (kSize = %d)", refSequence.name, this.startIndex, kUtil.kSize));
		}
		
		if (endKmerIndex == -1) {
			this.rightEnd = true;
			
			this.endKmerIndex = count.length - 1;
			this.endIndex = refSequence.size - 1;
			
			rightEndKmer = null;
			
		} else {
			this.rightEnd = false;
			
			this.endKmerIndex = endKmerIndex;
			this.endIndex = endKmerIndex + kUtil.kSize - 1;
			
			rightEndKmer = getKmer(this.endKmerIndex, kUtil);
			
			if (rightEndKmer == null)
				throw new IllegalArgumentException(String.format("Reference sequence (%s) contains an ambiguous base in the right-end k-mer at index %d (kSize = %d)", refSequence.name, this.endKmerIndex, kUtil.kSize));
		}
		
		// Set fields
		this.stats = RegionStats.getStats(count, this.startKmerIndex, this.endKmerIndex);
		
		return;
	}
	
	/**
	 * Get the k-mer bordering the left end of this active region or <code>null</code> if
	 * the active region extends to the left end of the reference with no bordering k-mer.
	 * 
	 * @return K-mer bordering the left end of this active region or <code>null</code> if
	 *   there is no bordering k-mer.
	 */
	public int[] getLeftEndKmer() {
		if (leftEndKmer == null)
			return null;
		
		return Arrays.copyOf(leftEndKmer, leftEndKmer.length);
	}
	
	/**
	 * Get the k-mer bordering the right end of this active region or <code>null</code> if
	 * the active region extends to the right end of the reference with no bordering k-mer.
	 * 
	 * @return K-mer bordering the right end of this active region or <code>null</code> if
	 *   there is no bordering k-mer.
	 */
	public int[] getRightEndKmer() {
		if (rightEndKmer == null)
			return null;
		
		return Arrays.copyOf(rightEndKmer, rightEndKmer.length);
	}
	
	/**
	 * Check for a match between a k-mer and the left-end border of this active region.
	 * 
	 * @param kmer K-mer to check.
	 * @param kUtil K-mer utility for checking.
	 * 
	 * @return <code>true</code> if <code>kmer</code> matches the left-end k-mer of this active
	 *   region or <code>false</code> if they do not match, <code>kmer</code> is <code>null</code>,
	 *   or this region has no left-end bordering k-mer.
	 * 
	 * @throws IllegalArgumentException If the array size of <code>kmer</code> or the left-end
	 *   k-mer is incompatible with <code>kUtil</code>.
	 */
	public boolean matchLeftEnd(int[] kmer, KmerUtil kUtil)
			throws IllegalArgumentException {
		
		// Check arguments
		if (kUtil == null)
			throw new NullPointerException("K-mer utility is null");
		
		if (kmer == null || leftEndKmer == null)
			return false;
		
		// Compare
		try {
			return kUtil.eq(kmer, leftEndKmer);
			
		} catch (IllegalArgumentException ex) {
			
			if (kmer.length < kUtil.wordSize)
				throw new IllegalArgumentException(String.format("K-mer array size mismatch: K-mer size = %d, expected array size >= %d", kmer.length, kUtil.wordSize));
			
			throw new IllegalArgumentException(String.format("K-mer utility mismatch: Left-end k-mer size = %d, expected k-mer array size >= %d (by k-mer utility)", leftEndKmer.length, kUtil.wordSize));
		}
	}
	
	/**
	 * Check for a match between a k-mer and the right-end border of this active region.
	 * 
	 * @param kmer K-mer to check.
	 * @param kUtil K-mer utility for checking.
	 * 
	 * @return <code>true</code> if <code>kmer</code> matches the right-end k-mer of this active
	 *   region or <code>false</code> if they do not match, <code>kmer</code> is <code>null</code>,
	 *   or this region has no right-end bordering k-mer.
	 * 
	 * @throws IllegalArgumentException If the array size of <code>kmer</code> or the right-end
	 *   k-mer is incompatible with <code>kUtil</code>.
	 */
	public boolean matchRightEnd(int[] kmer, KmerUtil kUtil)
			throws IllegalArgumentException {
		
		// Check arguments
		if (kUtil == null)
			throw new NullPointerException("K-mer utility is null");
		
		if (kmer == null || rightEndKmer == null)
			return false;
		
		// Compare
		try {
			return kUtil.eq(kmer, rightEndKmer);
			
		} catch (IllegalArgumentException ex) {
			
			if (kmer.length < kUtil.wordSize)
				throw new IllegalArgumentException(String.format("K-mer array size mismatch: K-mer size = %d, expected array size >= %d", kmer.length, kUtil.wordSize));
			
			throw new IllegalArgumentException(String.format("K-mer utility mismatch: Right-end k-mer size = %d, expected k-mer array size >= %d (by k-mer utility)", rightEndKmer.length, kUtil.wordSize));
		}
	}
	
	/**
	 * Get a k-mer at a given start location of the sequence.
	 * 
	 * @param index Index of the first base of the k-mer.
	 * @param kUtil K-mer utility.
	 * 
	 * @return The k-mer at <code>index</code> or <code>null</code> if the k-mer
	 *   contains ambiguous bases.
	 */
	private int[] getKmer(int index, KmerUtil kUtil) {
		
		// Check arguments
		assert (index >= 0) :
			"index is negative: " + index;
		
		assert (kUtil != null) :
			"kUtil is null";
		
		// Create k-mer and cache sequence
		int[] kmer = new int[kUtil.kmerArraySize];
		byte[] sequence = refRegion.sequence;
		int endIndex = index + kUtil.kSize;
		
		// Check range
		assert (endIndex <= sequence.length) :
			String.format("index (%d) + kSize (%d) > sequence.length (%d)", index, kUtil.kSize, sequence.length);
		
		// Build k-mer
		try {
			for (; index < endIndex; ++index)
				kUtil.append(kmer, BYTE_TO_BASE[sequence[index]]);
			
		} catch (NullPointerException ex) {
			return null;
		}
		
		return kmer;
	}
	
	/**
	 * Get a string representation of this active region.
	 * 
	 * @return A string describing this active region
	 */
	@Override
	public String toString() {
		return "ActiveRegion[name=" + refRegion.name +
				", start=" +
				((startIndex != -1) ? startIndex : "LEFT_END") +
				", end=" +
				((endIndex != -1) ? endIndex : "RIGHT_END") +
				"]";
	}
	
	/**
	 * Compares this active region to another.
	 * 
	 * @param other Other active region this is compared to.
	 * 
	 * @return A negative number, zero, or a positive number if <code>other</code> is less than, equal to,
	 *   or greater than this region (respectively).
	 * 
	 * @throws NullPointerException If <code>other</code> is <code>null</code>.
	 */
	@Override
	public int compareTo(ActiveRegion other)
			throws NullPointerException {
		
		int compareVal;
		
		// Check arguments
		if (other == null)
			throw new NullPointerException("Cannot comapre this region haplotype to null");
		
		// Compare reference region
		compareVal = refRegion.compareTo(other.refRegion);
		
		if (compareVal != 0)
			return compareVal;
		
		// Compare start index
		if (startIndex != other.startIndex)
			return startIndex - other.startIndex;
		
		// Compare end index
		return endIndex - other.endIndex;
	}
}
