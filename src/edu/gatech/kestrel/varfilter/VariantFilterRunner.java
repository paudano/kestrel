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

package edu.gatech.kestrel.varfilter;

import edu.gatech.kestrel.KestrelConstants;
import edu.gatech.kestrel.variant.VariantCall;

/**
 * Checks variants on a set of filters.
 */
public class VariantFilterRunner {
	
	/** Array of filters. */
	private VariantFilter[] filters;
	
	/** Number of filters. */
	private int nFilter;
	
	/** Capacity of <code>filters</code>. */
	private int filterCapacity;
	
	/** Default capacity. */
	public static final int DEFAULT_FILTER_CAPACITY = 5;
	
	/**
	 * Create a new variant filter.
	 * 
	 * @param filterCapacity Initial filter capacity. If more filters than this are added, the
	 *   capacity is automatically expanded.
	 *   
	 * @throws IllegalArgumentException If <code>filterCapacity</code> is negative.
	 */
	public VariantFilterRunner(int filterCapacity)
			throws IllegalArgumentException {
		
		if (filterCapacity < 0)
			throw new IllegalArgumentException("Cannot set a negative filter capacity: " + filterCapacity);
		
		filters = new VariantFilter[filterCapacity];
		this.filterCapacity = filterCapacity;
		nFilter = 0;
		
		return;
	}
	
	/**
	 * Create a new variant filter.
	 */
	public VariantFilterRunner() {
		this(DEFAULT_FILTER_CAPACITY);
		
		return;
	}
	
	/**
	 * Add a filter to this filter runner.
	 * 
	 * @param variantFilter Variant filter to add. If <code>null</code>, this function
	 *   call has no effect.
	 * 
	 * @throws IllegalStateException If the variant filter capacity is at its maximum value.
	 */
	public void addFilter(VariantFilter variantFilter)
			throws IllegalStateException {
		
		if (variantFilter == null)
			return;
		
		// Expand if needed
		if (nFilter == filterCapacity) {
			
			// Get new capacity
			int newFilterCapacity = (int) (filterCapacity * KestrelConstants.ARRAY_EXPAND_FACTOR);
			
			if (newFilterCapacity < 0) {
				if (filterCapacity == KestrelConstants.MAX_ARRAY_SIZE)
					throw new IllegalStateException("Maximum number of filters reached: " + KestrelConstants.MAX_ARRAY_SIZE);
				
				newFilterCapacity = KestrelConstants.MAX_ARRAY_SIZE;
			}
			
			// Expand and copy
			VariantFilter[] newFilters = new VariantFilter[newFilterCapacity];
			
			for (int index = 0; index < nFilter; ++index)
				newFilters[index] = filters[index];
			
			filters = newFilters;
			filterCapacity = newFilterCapacity;
		}
		
		// Add filter
		filters[nFilter++] = variantFilter;
	}
	
	/**
	 * Add a collection of variant filters.
	 * 
	 * @param variantFilters A collection of variant filters. If <code>null</code>,
	 *   this function call has no effect.
	 * 
	 * @throws IllegalStateException If the variant filter capacity is at its maximum value.
	 */
	public void addFilter(Iterable<VariantFilter> variantFilters)
			throws IllegalStateException {
		
		if (variantFilters == null)
			return;
		
		for (VariantFilter variantFilter : variantFilters)
			addFilter(variantFilter);
		
		return;
	}
	
	/**
	 * Filter against all variant filters.
	 * 
	 * @param variantCall Variant call to filter.
	 * 
	 * @return A variant call or <code>null</code> if the variant was removed by a filter.
	 */
	public VariantCall filter(VariantCall variantCall) {
		
		int filterIndex = 0;
		
		while (filterIndex < nFilter && variantCall != null)
			variantCall = filters[filterIndex++].filter(variantCall);
		
		return variantCall;
	}
}
