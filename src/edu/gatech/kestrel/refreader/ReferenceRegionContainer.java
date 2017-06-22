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

package edu.gatech.kestrel.refreader;

import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.NoSuchElementException;

import edu.gatech.kestrel.KestrelConstants;
import edu.gatech.kestrel.util.digest.ArrayUtil;

/**
 * Contains reference sequence regions and the reference sequences.
 */
public class ReferenceRegionContainer implements Iterable<ReferenceRegion> {
	
	/** Array of reference sequences. */
	private ReferenceSequence[] referenceSequenceArray;
	
	/** Number of elements in <code>referenceSequenceArray</code>. */
	private int referenceSequenceArraySize;
	
	/**
	 * Maps reference sequence names to a node containing the sequence and reference
	 * regions associated with it.
	 */
	private Map<String, ReferenceSequenceNode> refRegionMap;
	
	/** Number of regions. */
	private int regionCount;
	
	/** Default capacity of sequence and region arrays. */
	private static final int DEFAULT_LIST_CAPACITY = 5;
	
	/**
	 * Create a reference region container.
	 */
	public ReferenceRegionContainer() {
		
		referenceSequenceArray = new ReferenceSequence[DEFAULT_LIST_CAPACITY];
		referenceSequenceArraySize = 0;
		
		refRegionMap = new HashMap<>();
		regionCount = 0;
	}
	
	/**
	 * Add a reference region to this container.
	 * 
	 * @param refRegion Reference region to add.
	 * 
	 * @throws NullPointerException If <code>refRegion</code> is <code>null</code>.
	 * @throws IllegalArgumentException If <code>refRegion</code> is attached to a reference
	 *   sequence with the same name as a reference sequence already added (attached to a
	 *   reference region), but that is different.
	 * @throws IllegalStateException If too many regions attached to the same reference sequence
	 *   are added. The limit is <code>KestrelConstants.MAX_ARRAY_SIZE</code>.
	 * 
	 * @see KestrelConstants#MAX_ARRAY_SIZE
	 */
	public void add(ReferenceRegion refRegion)
			throws NullPointerException, IllegalArgumentException, IllegalStateException {
		
		ReferenceSequenceNode node;
		
		if (refRegion == null)
			throw new NullPointerException("Cannot add reference region: null");
		
		node = refRegionMap.get(refRegion.referenceSequence.name);
		
		if (node == null) {
			node = new ReferenceSequenceNode(refRegion.referenceSequence);
			refRegionMap.put(refRegion.referenceSequence.name, node);
			
			// Add reference sequence
			if (referenceSequenceArraySize == referenceSequenceArray.length)
				referenceSequenceArray = ArrayUtil.expandArray(referenceSequenceArray, ReferenceSequence.class);  // throws IllegalStateException
			
			referenceSequenceArray[referenceSequenceArraySize++] = refRegion.referenceSequence;
			
		} else {
			if (! refRegion.referenceSequence.equals(node.referenceSequence))
				throw new IllegalArgumentException(String.format("Cannot add reference region %s: This region is attached to reference sequence %s, but a reference sequence with the same name has already been added: %s", refRegion.toString(), refRegion.referenceSequence.toString(), node.referenceSequence.toString()));
		}
		
		node.add(refRegion);  // throws IllegalStateException
		++regionCount;
		
		return;
	}
	
	/**
	 * Add an array of reference regions.
	 * 
	 * @param refRegion Reference region array.
	 * 
	 * @param size Number of elements of <code>refRegion</code> to add.
	 * 
	 * @throws NullPointerException If <code>refRegion</code> is <code>null</code>
	 * @throws IllegalArgumentException If an element in <code>refRegion</code> is attached
	 *   to a reference sequence with the same name as a reference sequence already added
	 *   (attached to a reference region), but that is different.
	 * @throws IllegalStateException If too many regions attached to the same reference sequence
	 *   are added. The limit is <code>KestrelConstants.MAX_ARRAY_SIZE</code>.
	 */
	public void addAll(ReferenceRegion[] refRegion, int size)
			throws NullPointerException, IllegalArgumentException, IllegalStateException {
		
		// Check arguments
		if (refRegion == null)
			throw new NullPointerException("Cannot add reference region: null");
		
		if (size > refRegion.length)
			size = refRegion.length;
		
		// Add
		ADD_LOOP:
		for (int index = 0; index < size; ++index) {
			
			if (refRegion[index] == null)
				continue ADD_LOOP;
			
			add(refRegion[index]);  // throws IllegalArgumentException, IllegalStateException
		}
		
		return;
	}
	
	/**
	 * Get reference regions associated with a reference sequence.
	 * 
	 * @param refSequence Reference sequence.
	 * 
	 * @return An array of reference regions.
	 * 
	 * @throws NullPointerException If <code>refSequence</code> is <code>null</code>.
	 * @throws NoSuchElementException If there are no regions associated with <code>refSequence</code>.
	 */
	public ReferenceRegion[] get(ReferenceSequence refSequence)
			throws NullPointerException, NoSuchElementException {
		
		ReferenceSequenceNode node;
		
		if (refSequence == null)
			throw new NullPointerException("Reference sequence is null");
		
		node = refRegionMap.get(refSequence.name);
		
		if (node == null)
			throw new NoSuchElementException("Reference sequence does not have regions in this container: " + refSequence);
		
		return Arrays.copyOf(node.regions, node.size);
	}
	
	/**
	 * Get an iterator for the reference regions in this container.
	 * 
	 * @return A reference region iterator.
	 */
	@Override
	public Iterator<ReferenceRegion> iterator() {
		return new ReferenceRegionIterator();
	}
	
	/**
	 * Get an iterator for the reference sequences in this container.
	 * 
	 * @return A reference sequence iterator.
	 */
	public Iterator<ReferenceSequence> refSequenceIterator() {
		return new ReferenceArrayIterator();
	}
	
	/**
	 * Get an array of reference sequences.
	 * 
	 * @return Array of reference sequences.
	 */
	public ReferenceSequence[] getReferenceSequenceArray() {
		return Arrays.copyOf(referenceSequenceArray, referenceSequenceArraySize);
	}
	
	/**
	 * Sort reference sequences. The iterator will return references in this order. If not
	 * sorted, then the references will be in the order they found as reference regions
	 * were added. Any new reference sequences added after this method call will be appended
	 * to the end in the order they were added.
	 * 
	 * @param comparator Comparator for custom sequence sorting or <code>null</code> to use
	 *   the default sorting by the sequence name.
	 */
	public void sortReferences(Comparator<ReferenceSequence> comparator) {
		
		if (comparator == null)
			Arrays.sort(referenceSequenceArray, 0, referenceSequenceArraySize);
		else
			Arrays.sort(referenceSequenceArray, 0, referenceSequenceArraySize, comparator);
		
		return;
	}
	
	/**
	 * Sort reference sequences. The iterator will return references in this order. If not
	 * sorted, then the references will be in the order they found as reference regions
	 * were added. Any new reference sequences added after this method call will be appended
	 * to the end in the order they were added.
	 */
	public void sortReferences() {
		
		sortReferences(null);
				
		return;
	}
	
	/**
	 * Get the number of reference regions.
	 * 
	 * @return The number of reference regions.
	 */
	public int getSize() {
		return regionCount;
	}
	
	/**
	 * Determine if this container is empty.
	 * 
	 * @return <code>true</code> if the container has no reference regions.
	 */
	public boolean isEmpty() {
		return regionCount == 0;
	}
	
	/**
	 * Get the number of reference sequences.
	 * 
	 * @return The number of reference sequences.
	 */
	public int getReferenceSize() {
		return referenceSequenceArraySize;
	}
	
	/**
	 * Clear all regions from this container.
	 */
	public void clear() {
		
		refRegionMap.clear();
		regionCount = 0;
		
		while (referenceSequenceArraySize > 0)
			referenceSequenceArray[referenceSequenceArraySize--] = null;
		
		referenceSequenceArray[0] = null;
		
		return;
	}
	
	/**
	 * Get a string representation of this container.
	 * 
	 * @return A string representation of this container.
	 */
	@Override
	public String toString() {
		return String.format("ReferenceRegionContainer[references=%d, regions=%d]", referenceSequenceArraySize, regionCount);
	}
		
	/**
	 * Stores a reference sequence and regions associated with it.
	 */
	private class ReferenceSequenceNode {
		
		/** Reference sequence regions are associated with. */
		public final ReferenceSequence referenceSequence;
		
		/** Regions in <code>refSequence</code>. */
		public ReferenceRegion[] regions;
		
		/** Number of elements in <code>regions</code>. */
		public int size;
		
		/**
		 * Create a new reference sequence node to associate a reference sequence with
		 * a set of regions.
		 * 
		 * @param referenceSequence Reference sequence.
		 */
		public ReferenceSequenceNode(ReferenceSequence referenceSequence) {
			
			assert(referenceSequence != null) :
				"Reference sequence is null";
			
			this.referenceSequence = referenceSequence;
			
			regions = new ReferenceRegion[DEFAULT_LIST_CAPACITY];
			size = 0;
			
			return;
		}
		
		/**
		 * Add a region to this container.
		 * 
		 * @param refRegion Region to add. Must not be <code>null</code> and must be attached to
		 *   reference <code>refSequence</code>.
		 * 
		 * @throws IllegalArgumentException If <code>refRegion</code> overlaps with
		 *   another region in this container.
		 * @throws IllegalStateException If the region array is full and at its maximum size.
		 */
		public void add(ReferenceRegion refRegion)
				throws IllegalArgumentException, IllegalStateException {
			
			assert (refRegion != null) :
				"reference region is null";
			
			assert (refRegion.referenceSequence == referenceSequence) :
				String.format("Cannot add region: Mismatching reference sequence: found=%s, expected=%s", refRegion.referenceSequence, referenceSequence);
			
			// Check capacity and expand array
			if (size == regions.length)
				regions = ArrayUtil.expandArray(regions, ReferenceRegion.class);  // throws IllegalStateException
			
			// Find add location
			int addIndex = size;
			
			while (addIndex > 0 && regions[addIndex - 1].compareTo(refRegion) > 0) {
				
				regions[addIndex] = regions[addIndex - 1];  // Shift
				--addIndex;
			}
			
			// Intervals on either side must not overlap
			if (addIndex > 0 && regions[addIndex - 1].interval.overlaps(refRegion.interval)) {
				
				IllegalArgumentException ex = new IllegalArgumentException(String.format("New region overlaps with an existing region: new=%s, existing=%s", refRegion.toString(), regions[addIndex - 1].toString()));
				
				// Unshift
				while (addIndex < size) {
					regions[addIndex - 1] = regions[addIndex];
					++addIndex;
				}
				
				throw ex;
			}
			
			if (addIndex < size && regions[addIndex + 1].interval.overlaps(refRegion.interval)) {
				IllegalArgumentException ex =  new IllegalArgumentException(String.format("New region overlaps with an existing region: new=%s, existing=%s", refRegion.toString(), regions[addIndex].toString()));
				
				// Unshift
				while (addIndex < size) {
					regions[addIndex - 1] = regions[addIndex];
					++addIndex;
				}
				
				throw ex;
			}
			
			// Add
			regions[addIndex] = refRegion;
			++size;
			
			return;
		}
	}
	
	/**
	 * An iterator over the reference regions in this container.
	 */
	private class ReferenceRegionIterator implements Iterator<ReferenceRegion> {
		
		/** Node containing the next reference sequence. */
		private ReferenceSequenceNode node;
		
		/** Next index of <code>refSeqList</code> where the next node is found. */
		private int referenceSequenceListIndex;
		
		/** Index of the next region in <code>node</code>. */
		private int nextRegionIndex;
		
		/**
		 * Get an iterator for reference regions.
		 */
		public ReferenceRegionIterator() {
			
			if (referenceSequenceArraySize == 0) {
				node = null;
				return;
			}
			
			node = refRegionMap.get(referenceSequenceArray[0].name);
			
			assert (node != null) :
				"Node is null for reference sequence: " + referenceSequenceArray[0];
			
			referenceSequenceListIndex = 1;
			nextRegionIndex = 0;
			
			return;
		}
		
		/**
		 * Determine if this region has another element.
		 * 
		 * @return <code>true</code> if this region has another element.
		 */
		@Override
		public boolean hasNext() {
			
			return node != null;
		}
		
		/**
		 * Get the next element from this iterator.
		 * 
		 * @return Next element.
		 * 
		 * @throws NoSuchElementException If there are no more elements.
		 */
		@Override
		public ReferenceRegion next()
				throws NoSuchElementException {
			
			ReferenceRegion refRegion;
			
			if (node == null)
				throw new NoSuchElementException("No more elements");
			
			refRegion = node.regions[nextRegionIndex++];
			
			if (nextRegionIndex >= node.size) {
				
				if (referenceSequenceListIndex < referenceSequenceArraySize) {
					node = refRegionMap.get(referenceSequenceArray[referenceSequenceListIndex++].name);
					
					assert (node != null) :
						"Node is null for reference sequence: " + referenceSequenceArray[referenceSequenceListIndex - 1];
					
				} else {
					node = null;
				}
				
				nextRegionIndex = 0;
			}
			
			return refRegion;
		}
		
		/**
		 * Method disabled.
		 * 
		 * @throws UnsupportedOperationException Always.
		 */
		@Override
		public void remove()
				throws UnsupportedOperationException {
			
			throw new UnsupportedOperationException("Iterator does not support remove()");
			
		}
	}
	
	/**
	 * Iterator for the reference sequences.
	 */
	private class ReferenceArrayIterator implements Iterator<ReferenceSequence> {
		
		/** Next index of the reference sequence array. */
		private int nextIndex;
		
		/**
		 * Create a new iterator.
		 */
		public ReferenceArrayIterator() {
			nextIndex = 0;
			
			return;
		}
		
		/**
		 * Determine if this iterator has another element.
		 * 
		 * @return <code>true</code> if this iterator has another element.
		 */
		@Override
		public boolean hasNext() {
			return nextIndex < referenceSequenceArraySize;
		}
		
		/**
		 * Get the next reference sequence.
		 * 
		 * @return Next element.
		 * 
		 * @throws NoSuchElementException When this iterator is depleted.
		 */
		@Override
		public ReferenceSequence next()
				throws NoSuchElementException {
			
			if (nextIndex >= referenceSequenceArraySize)
				throw new NoSuchElementException("No more reference sequence elements");
			
			return referenceSequenceArray[nextIndex++];
		}
		
		/**
		 * Disabled remove method.
		 * 
		 * @throws UnsupportedOperationException Always.
		 */
		@Override
		public void remove()
				throws UnsupportedOperationException {
			
			throw new UnsupportedOperationException("remove() is not supported by this reference array iterator");
		}
	}
}
