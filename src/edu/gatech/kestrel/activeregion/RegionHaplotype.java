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

/**
 * Associates an active region and haplotypes.
 */
public class RegionHaplotype implements Comparable<RegionHaplotype> {
	
	/** Active region. */
	public final ActiveRegion activeRegion;
	
	/**
	 * Haplotypes covering this active region or an empty array if no haplotypes were
	 * found in this region. Regions with no haplotypes (wild-type only) can be used to
	 * build gVCF files.
	 */
	public final Haplotype[] haplotype;
	
	/** Estimation of the minimum depth. */
	public final int minDepth;
	
	/** <code>true</code> if this region contains only wildtype haplotypes. */
	public final boolean isWildtype;
	
	/**
	 * Create a new active region and associated haplotypes.
	 * 
	 * @param activeRegion Active region.
	 * @param haplotype Array of haplotypes.
	 * 
	 * @throws NullPointerException If <code>activeRegion</code> or <code>haplotype</code>
	 *   is <code>null</code>.
	 * @throws IllegalArgumentException If <code>haplotype</code> contains null references or
	 *   has no elements.
	 */
	public RegionHaplotype(ActiveRegion activeRegion, Haplotype[] haplotype)
			throws NullPointerException, IllegalArgumentException {
		
		int minDepth;
		boolean isWildtype;
		
		// Check arguments
		if (activeRegion == null)
			throw new NullPointerException("Active region is null");
		
		if (haplotype == null)
			throw new NullPointerException("Haplotype array is null");
		
		if (haplotype.length == 0)
			throw new IllegalArgumentException("Haplotype array is empty");
		
		for (int index = 0; index < haplotype.length; ++index)
			if (haplotype[index] == null)
				throw new IllegalArgumentException("Haplotype array contains a null reference at index " + index);
		
		// Assign fields
		this.activeRegion = activeRegion;
		this.haplotype = haplotype;
		
		// Estimate depth and determine wildtype status
		minDepth = activeRegion.stats.min;
		isWildtype = true;
		
		for (Haplotype nextHaplo : haplotype) {
			minDepth += nextHaplo.stats.min;
			
			if (! nextHaplo.isWildtype())
				isWildtype = false;
		}
		
		this.minDepth = minDepth;
		this.isWildtype = isWildtype;
		
		return;
	}
	
	/**
	 * Get a string representing this region-haplotype.
	 * 
	 * @return A string representing this region-haplotype.
	 */
	@Override
	public String toString() {
		return String.format("RegionHaplotype[region=%s, haplotypes=%d, wildtype=%b, mindepth=%d]", activeRegion, haplotype.length, isWildtype, minDepth);
	}
	
	/**
	 * Compare this region to another.
	 * 
	 * @param other Region this region is compared to.
	 * 
	 * @return A negative number, zero, or a positive number if <code>other</code> is less than, equal to,
	 *   or greater than this region (respectively).
	 * 
	 * @throws NullPointerException If <code>other</code> is <code>null</code>.
	 */
	@Override
	public int compareTo(RegionHaplotype other)
			throws NullPointerException {
		
		if (other == null)
			throw new NullPointerException("Cannot comapre this region haplotype to null");
		
		return activeRegion.compareTo(other.activeRegion);
	}
}
