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

import edu.gatech.kestrel.refreader.ReferenceRegion;

/**
 * Contains results from active region detection.
 */
public class ActiveRegionContainer {
	
	/** Reference region the active region haplotypes were found in. */
	public final ReferenceRegion refRegion;
	
	/** Active regions detected. */
	public final RegionHaplotype[] haplotypes;
	
	/** K-mer frequency stats over the reference region. */
	public final RegionStats stats;
	
	/**
	 * Create a new active region container.
	 * 
	 * @param refRegion Reference region.
	 * @param haplotypes Active regions from the reference. 
	 * @param count An array of k-mer counts over the reference region.
	 * 
	 * @throws NullPointerException If <code>refRegion</code> or <code>count</code> is <code>null</code>. 
	 * @throws IllegalArgumentException If <code>haplotypes</code> is <code>null</code>.
	 */
	public ActiveRegionContainer(ReferenceRegion refRegion, RegionHaplotype[] haplotypes, int[] count)
			throws NullPointerException, IllegalArgumentException {
		
		// Check arguments
		if (refRegion == null)
			throw new NullPointerException("Reference region is null");
		
		if (count == null)
			throw new NullPointerException("Count array is null");
		
		if (haplotypes == null)
			haplotypes = new RegionHaplotype[0];
		
		for (int index = 0; index < haplotypes.length; ++index)
			if (haplotypes[index] == null)
				throw new IllegalArgumentException("Haplotypes contains a null reference at index " + index);
		
		// Get stats
		stats = RegionStats.getStats(count, 0, count.length);
		
		// Set fields
		this.refRegion = refRegion;
		this.haplotypes = Arrays.copyOf(haplotypes, haplotypes.length);
		
		return;
	}
}
