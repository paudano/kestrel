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
import edu.gatech.kanalyze.util.StringUtil;
import edu.gatech.kanalyze.util.kmer.KmerUtil;
import edu.gatech.kestrel.counter.CountMap;

/**
 * K-mer count summary statistics over a sequence region. When attached to an active
 * region, it summarizes the wild-type k-mers. When attached to a haplotype, it
 * summarizes the k-mers associated with the k-mers of that haplotype. When applied
 * to a variant, it summarizes the k-mers that touch that variant.
 */
public final class RegionStats {
	
	/** Minimum k-mer count. */
	public final int min;
	
	/** 25th percentile. */
	public final float pct25;
	
	/** 50th percentile. */
	public final float pct50;
	
	/** 75th percentile. */
	public final float pct75;
	
	/** Maximum k-mer count. */
	public final int max;
	
	/** Number of k-mers counted. */
	public final int n;
	
	/** Translates bytes to bases. */
	private static final Base[] BYTE_TO_BASE = KAnalyzeConstants.getByteToBaseArray();
	
	/**
	 * Create a set of summary statistics.
	 * 
	 * @param min Minimum k-mer count.
	 * @param pct25 25th percentile.
	 * @param pct50 50th percentile.
	 * @param pct75 75th percentile.
	 * @param max Maximum k-mer count.
	 * @param n Number of statistics in this set.
	 */
	private RegionStats(int min, float pct25, float pct50, float pct75, int max, int n) {
		
		this.min = min;
		this.pct25 = pct25;
		this.pct50 = pct50;
		this.pct75 = pct75;
		this.max = max;
		this.n = n;
		
		return;
	}
	
	/**
	 * Get summary statistics for a sequence.
	 * 
	 * @param sequence Sequence.
	 * @param start Index of the first base in <code>sequence</code>.
	 * @param end Index of the last base in <code>sequence</code>, exclusive. The range of
	 *   bases characterized is <code>[start, end)</code>.
	 * @param counter K-mer counter.
	 * @param countRev If <code>true</code>, include the count for the reverse-complement
	 *   of each k-mer.
	 * 
	 * @return Summary statistics for the k-mers in <code>sequence</code>.
	 * 
	 * @throws NullPointerException If <code>sequence</code> or <code>counter</code> is <code>null</code>.
	 * @throws IllegalArgumentException If <code>kSize</code> is less than <code>1</code>, <code>start</code> is
	 *   less than <code>0</code>, <code>end</code> is greater than the sequence length, or if <code>end - start</code>
	 *   does not contain at least 2 k-mers.
	 */
	public static final RegionStats getStats(byte[] sequence, int start, int end, CountMap counter, boolean countRev)
			throws NullPointerException, IllegalArgumentException {
		
		KmerUtil kUtil;  // K-mer util
		int[] kmer;      // K-mer
		int[] revKmer;   // Reverse-complement of kmer
		
		Base base;   // A base to be appended to a k-mer
		int nBases;  // Number of bases in the sequence
		
		int seqIndex;    // Current index in sequence
		int countIndex;  // Current index of the count array
		
		int[] count;  // K-mer count array
		int length;   // Length of count
		int nLess1;   // Length - 1 (saved to avoid repetitive calculation)
		
		int loc;       // First of two counts (index of count) for linear interpolation of counts
		float offset;  // Offset for linear interpolation of counts
		
		float pct25;  // 25th percentile
		float pct50;  // 50th percentile
		float pct75;  // 75th percentile
		
		// Check arguments
		if (sequence == null)
			throw new NullPointerException("Sequence is null");
		
		if (counter == null)
			throw new NullPointerException("K-mer counter is null");
		
		kUtil = counter.kUtil;
		
		if (start < 0)
			throw new IllegalArgumentException("start is negative: " + start);
		
		if (end > sequence.length)
			throw new IllegalArgumentException(String.format("end is greater than the sequence length (%d): %d", sequence.length, end));
		
		nBases = end - start;
		length = nBases - kUtil.kSize + 1;
		
		if (length < 2)
			throw new IllegalArgumentException(String.format("sequence must contain 2 kmers (k-mer size = %d): end (%d) - start (%d) = %d", kUtil.kSize, end, start, end - start));
		
		// Set utilities and kmer
		kmer = new int[kUtil.kmerArraySize];
		
		// Create count array
		count = new int[length];
		
		// Fill the first k-mer
		for (seqIndex = 0; seqIndex < kUtil.kSize - 1; ++seqIndex) {
			
			base = BYTE_TO_BASE[sequence[seqIndex]];
			
			if (base == null)
				throw new IllegalArgumentException(String.format("Sequence contains an illegal character at index %d: %s", seqIndex, StringUtil.charDescription(sequence[seqIndex])));
			
			kUtil.append(kmer, base);
		}
		
		if (countRev)
			revKmer = kUtil.revComplement(kmer, null);
		else
			revKmer = null;  // Quell initialization warnings; not accessed if not countRev
		
		countIndex = 0;
		
		// Generate k-mers
		while (seqIndex < nBases) {
			
			base = BYTE_TO_BASE[sequence[seqIndex]];
			
			if (base == null)
				throw new IllegalArgumentException(String.format("Sequence contains an illegal character at index %d: %s", seqIndex, StringUtil.charDescription(sequence[seqIndex])));
			
			kUtil.append(kmer, base);
			
			count[countIndex] = counter.get(kmer);
			
			if (countRev) {
				kUtil.prepend(revKmer, base.getComplement());
				count[countIndex] += counter.get(revKmer);
			}
			
			// Next base
			++seqIndex;
			++countIndex;
		}
		
		Arrays.sort(count);
		
		nLess1 = length - 1;
		
		// 25th percentile
		loc = (int) (nLess1 * .25F);
		offset = (nLess1 * .25F) - loc;
		
		assert (offset <= 1.0F && offset >= 0.0F) :
			String.format("Offset out of bounds (Q(.25), length = %d: %f", length, offset);
		
		pct25 = count[loc] * (1 - offset) + count[loc + 1] * offset;
		
		// 50th percentile
		loc = (int) (nLess1 * .5F);
		offset = (nLess1 * .5F) - loc;
		
		assert (offset <= 1.0F && offset >= 0.0F) :
			String.format("Offset out of bounds (Q(.50), length = %d: %f", length, offset);
		
		pct50 = count[loc] * (1 - offset) + count[loc + 1] * offset;
		
		// 75th percentile
		loc = (int) (nLess1 * .75F);
		offset = (nLess1 * .75F) - loc;
		
		assert (offset <= 1.0F && offset >= 0.0F) :
			String.format("Offset out of bounds (Q(.75), length = %d: %f", length, offset);
		
		pct75 = count[loc] * (1 - offset) + count[loc + 1] * offset;
		
		// Generate summary statistics
		return new RegionStats(
				count[0],
				pct25,
				pct50,
				pct75,
				count[length - 1],
				length
		);
	}
	
	/**
	 * Get summary statistics for a sequence.
	 * 
	 * @param count K-mer counts.
	 * @param start Index of the first element in <code>count</code>.
	 * @param end Index of the last element in <code>count</code>, exclusive. The range of
	 *   bases characterized is <code>[start, end)</code>.
	 * 
	 * @return Summary statistics for the k-mers in <code>sequence</code>.
	 * 
	 * @throws NullPointerException If <code>count</code> is <code>null</code>.
	 * @throws IllegalArgumentException If <code>start</code> is
	 *   less than <code>0</code>, <code>end</code> is greater than the count length, or if <code>end - start</code>
	 *   is less than <code>2</code>.
	 */
	public static final RegionStats getStats(int[] count, int start, int end)
			throws NullPointerException, IllegalArgumentException {
		
		int[] countSlice;    // Array of counts copied from the count array argument
		int countSliceSize;  // Length of countSlice
		
		int loc;       // First of two counts (index of count) for linear interpolation of counts
		float offset;  // Offset for linear interpolation of counts
		
		int nLess1;   // countSliceSize - 1 (saved to avoid repetitive calculation)
		
		float pct25;  // 25th percentile
		float pct50;  // 50th percentile
		float pct75;  // 75th percentile
		
		// Check arguments
		if (count == null)
			throw new NullPointerException("Cannot get stats from count array: null");
		
		if (start < 0)
			throw new IllegalArgumentException("Start index is negative: " + start);
		
		countSliceSize = end - start;
		
		if (end > count.length || countSliceSize < 1)
			throw new IllegalArgumentException(String.format("End index is not less than the array length (%d): the range from start to end is not at least 1 (%d): %d", count.length, countSliceSize, end));
		
		if (countSliceSize == 1)
			return new RegionStats(count[0], count[0], count[0], count[0], count[0], 1);
		
		// Copy and sort counts
		countSlice = Arrays.copyOfRange(count, start, end);
		Arrays.sort(countSlice);
		
		nLess1 = countSliceSize - 1;
		
		// 25th percentile
		loc = (int) (nLess1 * .25F);
		offset = (nLess1 * .25F) - loc;
		
		assert (offset <= 1.0F && offset >= 0.0F) :
			String.format("Offset out of bounds (Q(.25), length = %d: %f", countSliceSize, offset);
		
		pct25 = countSlice[loc] * (1 - offset) + countSlice[loc + 1] * offset;
		
		// 50th percentile
		loc = (int) (nLess1 * .5F);
		offset = (nLess1 * .5F) - loc;
		
		assert (offset <= 1.0F && offset >= 0.0F) :
			String.format("Offset out of bounds (Q(.50), length = %d: %f", countSliceSize, offset);
		
		pct50 = countSlice[loc] * (1 - offset) + countSlice[loc + 1] * offset;
		
		// 75th percentile
		loc = (int) (nLess1 * .75F);
		offset = (nLess1 * .75F) - loc;
		
		assert (offset <= 1.0F && offset >= 0.0F) :
			String.format("Offset out of bounds (Q(.75), length = %d: %f", countSliceSize, offset);
		
		pct75 = countSlice[loc] * (1 - offset) + countSlice[loc + 1] * offset;
		
		// Generate summary statistics
		return new RegionStats(
				countSlice[0],
				pct25,
				pct50,
				pct75,
				countSlice[countSliceSize - 1],
				countSliceSize
		);
	}
	
	/**
	 * Get a string representation of these stats.
	 */
	@Override
	public String toString() {
		return String.format("Stats[min=%d, 25th=%f, 50th=%f, 75th=%f, max=%d, n=%d]", min, pct25, pct50, pct75, max, n);
	}
}
