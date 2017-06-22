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

package edu.gatech.kestrel.interval;

import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;

import edu.gatech.kestrel.KestrelConstants;

/**
 * Manages a set of region intervals.
 */
public class RegionIntervalContainer {
	
	/**
	 * A map of reference sequence names to a collection of intervals associated with
	 * the reference sequence.
	 */
	private final HashMap<String, ReferenceNode> refTable;
	
	/** The number if intervals in this container. */
	private long size;
	
	/**
	 * Create a new region interval container.
	 */
	public RegionIntervalContainer() {
		
		refTable = new HashMap<>();
		
		return;
	}
	
	/**
	 * Add an interval to this container.
	 * 
	 * @param interval Interval to add.
	 * 
	 * @throws NullPointerException If <code>interval</code> is <code>null</code>.
	 * @throws IllegalArgumentException If <code>interval</code> overlaps with an existing
	 *   interval.
	 * @throws IllegalStateException If the interval array for this reference sequence is
	 *   full and cannot be expanded.
	 */
	public void add(RegionInterval interval)
			throws NullPointerException, IllegalArgumentException, IllegalStateException {
		
		ReferenceNode refNode;
		
		// Check arguments
		if (interval == null)
			throw new NullPointerException("Cannot add interval: null");
		
		// Get reference node
		refNode = refTable.get(interval.sequenceName);
		
		if (refNode == null) {
			refNode = new ReferenceNode();
			refTable.put(interval.sequenceName, refNode);
		}
		
		// Add
		refNode.add(interval);  // throws IllegalArgumentException, IllegalStateException
		
		++size;
		
		return;
	}
	
	/**
	 * Get intervals for a reference by the reference name.
	 * 
	 * @param referenceName Reference name.
	 * 
	 * @return An array of reference intervals or an empty array if there are no intervals
	 *   defined for this reference.
	 *   
	 * @throws NullPointerException If <code>referenceName</code> is <code>null</code>.
	 * @throws IllegalArgumentException If <code>referenceName</code> is empty.
	 */
	public RegionInterval[] getIntervals(String referenceName)
			throws NullPointerException, IllegalArgumentException {
		
		ReferenceNode refNode;
		
		// Check arguments
		if (referenceName == null)
			throw new NullPointerException("Cannot get intervals for reference with name: null");
		
		referenceName = referenceName.trim();
		
		if (referenceName.isEmpty())
			throw new IllegalArgumentException("Cannot get intervals for reference with an empty name");
		
		// Get node
		refNode = refTable.get(referenceName);
		
		if (refNode == null)
			return new RegionInterval[0];
		
		return Arrays.copyOf(refNode.intervals, refNode.intervals.length);
	}
	
	/**
	 * Get an object for mapping a reference sequence name to a sorted array of intervals
	 * associated with that reference. When a reference sequence has no intervals,
	 * <code>null</code> is retured by the map.
	 * 
	 * @return A region map keyed on reference names.
	 */
	public Map<String, RegionInterval[]> getMap() {
		
		HashMap<String, RegionInterval[]> intervalMap;  // Map to be returned
		
		ReferenceNode node;          // Node stored in this container                      
		
		// Init
		intervalMap = new HashMap<>();
		
		// Fill map (copy interval arrays)
		for (String referenceName : refTable.keySet()) {
			node = refTable.get(referenceName);
			intervalMap.put(referenceName, Arrays.copyOf(node.intervals, node.size));
		}
		
		return intervalMap;
	}
	
	/**
	 * Remove all intervals.
	 */
	public void clear() {
		refTable.clear();
		
		return;
	}
	
	/**
	 * Remove all intervals associated with a reference name.
	 * 
	 * @param referenceName Name of the reference sequence. If <code>null</code>, no
	 *   action is taken.
	 */
	public void clear(String referenceName) {
		
		if (referenceName != null)
			refTable.remove(referenceName);
		
		return;
	}
	
	/**
	 * Get the number of intervals in this container.
	 * 
	 * @return The number of intervals in this container.
	 */
	public long getSize() {
		return size;
	}
	
	/**
	 * Determine if this container has no intervals.
	 * 
	 * @return <code>true</code> if this container has no intervals.
	 */
	public boolean isEmpty() {
		return size == 0;
	}
	
	/**
	 * A reference sequence and all intervals on it.
	 */
	private class ReferenceNode {
		
		/** An array of intervals. */
		public RegionInterval[] intervals;
				
		/** Number of objects in <code>intervals</code>. */
		private int size;
		
		/** Capacity of <code>intervals</code>. */
		private int capacity;
		
		/** Interval array size. */
		public static final int DEFAULT_INTERVAL_ARRAY_SIZE = 5;
		
		/**
		 * Create a new reference sequence node.
		 */
		public ReferenceNode() {
			
			capacity = DEFAULT_INTERVAL_ARRAY_SIZE;
			intervals = new RegionInterval[capacity];
			size = 0;
			
			return;
		}
		
		/**
		 * Add an interval to this container.
		 * 
		 * @param interval Interval to add. Must not be <code>null</code>.
		 * 
		 * @throws IllegalArgumentException If <code>interval</code> overlaps with
		 *   another interval in this container.
		 * @throws IllegalStateException If the interval array is full and at its maximum size.
		 */
		public void add(RegionInterval interval)
				throws IllegalArgumentException, IllegalStateException {
			
			assert (interval != null) :
				"interval is null";
			
			// Check capacity and expand array
			if (size == capacity) {
				
				// Get the new capacity
				int newCapacity = (int) (size * KestrelConstants.ARRAY_EXPAND_FACTOR);
				
				if (newCapacity < 0 || newCapacity > KestrelConstants.MAX_ARRAY_SIZE) {
					
					if (capacity == KestrelConstants.MAX_ARRAY_SIZE)
						throw new IllegalStateException("Cannot add interval: Maximum interval array size reached: " + KestrelConstants.MAX_ARRAY_SIZE);
					
					newCapacity = KestrelConstants.MAX_ARRAY_SIZE;
				}
				
				// Create a new array
				RegionInterval[] newIntervals = new RegionInterval[newCapacity];
				
				for (int index = 0; index < size; ++index)
					newIntervals[index] = intervals[index];
				
				capacity = newCapacity;
				intervals = newIntervals;
			}
			
			// Find add location
			int addIndex = size;
			
			while (addIndex > 0 && intervals[addIndex - 1].compareTo(interval) > 0) {
				intervals[addIndex] = intervals[addIndex - 1];  // Shift
				--addIndex;
			}
			
			// Intervals on either side must not overlap on either side
			if (addIndex > 0 && intervals[addIndex - 1].overlaps(interval)) {
				
				IllegalArgumentException ex = new IllegalArgumentException(String.format("New interval overlaps with an existing interval: new=%s, existing=%s", interval.toString(), intervals[addIndex - 1].toString()));
				
				// Unshift
				while (addIndex < size) {
					intervals[addIndex - 1] = intervals[addIndex];
					++addIndex;
				}
				
				throw ex;
			}
			
			if (addIndex < size && intervals[addIndex].overlaps(interval)) {
				IllegalArgumentException ex =  new IllegalArgumentException(String.format("New interval overlaps with an existing interval: new=%s, existing=%s", interval.toString(), intervals[addIndex].toString()));
				
				// Unshift
				while (addIndex < size - 1) {
					intervals[addIndex] = intervals[addIndex + 1];
					++addIndex;
				}
				
				throw ex;
			}
			
			// Add
			intervals[addIndex] = interval;
			++size;
			
			return;
		}
	}
}
