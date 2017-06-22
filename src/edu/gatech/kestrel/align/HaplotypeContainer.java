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

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import edu.gatech.kestrel.activeregion.Haplotype;

/**
 * Manage haplotypes for the aligner.
 */
public class HaplotypeContainer {
	
	/** First haplotype in a haplotype stack. */
	private HaplotypeNode haplotypeHead;
	
	/** Number of haplotypes in this container. */
	private int size;
	
	/** Maximum number of haplotypes this container can accept. */
	public final int limit;
	
	/** Logger. */
	private final Logger logger;
	
	/**
	 * Create a haplotype container.
	 * 
	 * @param limit Haplotype limit.
	 * 
	 * @throws IllegalArgumentException If <code>limit</code> is less than <code>1</code>.
	 */
	public HaplotypeContainer(int limit)
			throws IllegalArgumentException {
		
		logger = LoggerFactory.getLogger(this.getClass());
		
		// Check argument
		if (limit < 1)
			throw new IllegalArgumentException(String.format("Haplotype container limit is less than 1: %d", limit));
		
		this.limit = limit;
		this.size = 0;
		
		return;
	}
	
	/**
	 * Add a haplotype to this container.
	 * 
	 * @param haplotype Haplotype.
	 * 
	 * @throws NullPointerException If <code>haplotype</code> is <code>null</code>.
	 */
	public void add(Haplotype haplotype)
			throws NullPointerException {
		
		Haplotype rmHaplotype;  // Removed haplotype if the container was trimmed to make space
		
		HaplotypeNode haplotypeNode;  // Node to be added.
		
		// Check arguments
		if (haplotype == null)
			throw new NullPointerException("Haplotype container cannot accept add haplotype: null");
		
		// Check capacity
		if (size >= limit) {
			rmHaplotype = removeMinHaplotype(haplotype.stats.min);
			
			if (rmHaplotype == null) {
				logger.trace("Rejecting haplotype: Haplotype container is at capacity: {} [minDepth={}]", haplotype, haplotype.stats.min);
				return;
			}
			
			logger.trace("Rejecting haplotype: Haplotype container is at capacity: {} [minDepth={}]", rmHaplotype, rmHaplotype.stats.min);
		}
		
		// Add haplotype
		haplotypeNode = new HaplotypeNode(haplotype);
		haplotypeNode.next = haplotypeHead;
		haplotypeHead = haplotypeNode;
		
		++size;
		
		return;
	}
	
	/**
	 * Get number of haplotypes in this container.
	 * 
	 * @return Number of haplotypes in this container.
	 */
	public int size() {
		return size;
	}
	
	/**
	 * Get an array of haplotypes.
	 * 
	 * @return Array of haplotypes.
	 */
	public Haplotype[] toArray() {
		
		HaplotypeNode nextNode = haplotypeHead;
		
		Haplotype[] hapArray = new Haplotype[size];
		
		for (int index = 0; index < size; ++index) {
			
			assert(nextNode != null):
				String.format("Haplotype container list ran out while converting to an array at index %d (size=%d)", index, size);
			
			hapArray[index] = nextNode.haplotype;
			nextNode = nextNode.next;
		}
		
		assert(nextNode == null):
			String.format("Haplotype container list contained more haplotypes than it's reported size (size=%d)", size);
		
		return hapArray;
	}
	
	/**
	 * Remove the haplotype with the least minimum k-mer depth if the depth is less than a set
	 * limit. Used to control haplotype container capacity.
	 *  
	 * @param limit Minimum k-mer depth limit.
	 * 
	 * @return The removed haplotype or <code>null<code> if no haplotype was found with a minimum k-mer depth
	 *   less than <code>limit</code>.
	 */
	private Haplotype removeMinHaplotype(int limit) {
		
		assert(limit > 0):
			"Limit must be positive: " + limit;
		
		// Declarations
		HaplotypeNode nextNode = haplotypeHead;  // Current node being analyzed
		
		HaplotypeNode lastNode = null;  // Last node that was nextNode (null when nextNode == haplotypeHead)
		
		HaplotypeNode minUpstream = null; // Node upstream of the minimum, or null if no minimum less than limit or minimum is head
		
		Haplotype rmHaplotype;  // Haplotype removed
		
		int lowestLimit = limit;  // Lowest limit (< limit when a node can be removed)
		
		// Find last min state
		while (nextNode != null) {
			
			if (nextNode.haplotype.stats.min < lowestLimit) {
				minUpstream = lastNode;
				lowestLimit = nextNode.haplotype.stats.min;
			}
			
			lastNode = nextNode;
			nextNode = nextNode.next;
		}
		
		// Do not remove if a lower limit was not found
		if (lowestLimit == limit)
			return null;
		
		// Unlink from head
		if (minUpstream == null) {
			rmHaplotype = haplotypeHead.haplotype;
			
			haplotypeHead = haplotypeHead.next;
			
			--size;
			
			return rmHaplotype;
		}
		
		// Unlink from list
		rmHaplotype = minUpstream.next.haplotype;
		
		minUpstream.next = minUpstream.next.next;
		
		--size;
		
		return rmHaplotype;
	}
	
	/**
	 * One node in a list of haplotype nodes.
	 */
	private class HaplotypeNode {
		
		/** Haplotype. */
		public final Haplotype haplotype;
		
		/** Next node. */
		public HaplotypeNode next;
		
		/**
		 * Create a haplotype node.
		 * 
		 * @param haplotype Haplotype.
		 */
		public HaplotypeNode(Haplotype haplotype) {
			
			assert(haplotype != null):
				"Haplotype is null";
			
			this.haplotype = haplotype;
			this.next = null;
			
			return;
		}
	}
}
