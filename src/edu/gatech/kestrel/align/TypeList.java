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

import edu.gatech.kestrel.KestrelConstants;

/**
 * A list of alignment type bytes used for building a list of alignments.
 */
public final class TypeList implements Cloneable {
	
	/** Type bytes. <code>type[0]</code> is always <code>0</code>. */
	private byte type[];
	
	/**
	 * Number of events of a type at <code>type</code>. <code>n[0]</code> is
	 * always <code>0</code>.
	 */
	private int n[];
	
	/** Number of elements in <code>type</code> and <code>n</code>. */
	private int size;
	
	/** Size of <code>type</code>. */
	private int capacity;
	
	/** Initial size of <code>type</code>. */
	public static final int INITIAL_CAPACITY = 100;
	
	/**
	 * Create a new type list.
	 */
	public TypeList() {
		
		capacity = INITIAL_CAPACITY;
		
		type = new byte[capacity];
		type[0] = 0;
		
		n = new int[capacity];
		n[0] = 0;
		
		size = 1;
		
		return;
	}
	
	/**
	 * Create a new list by deep-copying another.
	 * 
	 * @param o Stack to copy or <code>null</code> to create a default object.
	 */
	public TypeList(TypeList o) {
		
		// Set default if o is null
		if (o == null) {
			
			capacity = INITIAL_CAPACITY;
			
			type = new byte[capacity];
			type[0] = 0;
			
			n = new int[capacity];
			n[0] = 0;
			
			size = 1;
			return;
		}
		
		// Copy
		capacity = o.capacity;
		size = o.size;
		type = new byte[capacity];
		n = new int[capacity];
		
		for (int index = 0; index < size; ++index) {
			type[index] = o.type[index];
			n[index] = o.n[index];
		}
		
		return;
	}
	
	/**
	 * Add a type to the end of this list.
	 * 
	 * @param typeByte Type to add.
	 * 
	 * @throws IllegalStateException If the stack array is at maximum capacity.
	 */
	public final void add(byte typeByte)
			throws IllegalStateException {
		
		// Check arguments
		if (typeByte < AlignNode.MATCH || typeByte > AlignNode.DEL)
			throw new IllegalArgumentException("type is out of range [1, 4]: " + typeByte);
		
		// Add
		if (type[size - 1] == typeByte) {
			// Add to last event
			++n[size - 1];
			
		} else {
			
			// Expand
			if (size == capacity)
				expand();  // throws IllegalStateException
			
			
			// Add
			type[size] = typeByte;
			n[size] = 1;
			
			// Next position
			++size;
		}
		
		return;
	}
	
	/**
	 * Get an alignment from this list.
	 * 
	 * @param reverse If <code>true</code>, return the list in the reverse from the order
	 *   types were added.
	 *   
	 * @return A list of alignments.
	 */
	public final AlignNode toAlignment(boolean reverse) {
		
		// Declarations
		AlignNode alignHead;  // Head of the list
		
		// Init
		alignHead = null;
		
		// Build list (skipping index=0)
		if (reverse) {
			for (int index = 1; index < size; ++index)
				alignHead = new AlignNode(type[index], n[index], alignHead);
			
		} else {
			for (int index = size - 1; index > 0; --index)
				alignHead = new AlignNode(type[index], n[index], alignHead);
		}
		
		return alignHead;
	}
	
	/**
	 * Clear elements in this list.
	 */
	public final void clear() {
		size = 1;
		
		return;
	}
	
	/**
	 * Copy elements from another list into this list and remove the contents of this list.
	 * 
	 * @param o Type to copy from. If <code>null</code>, this list is cleared.
	 */
	public final void cloneFrom(TypeList o) {
		
		// Check arguments
		if (o == null) {
			this.size = 1;
			return;
		}
		
		// Check capacity
		if (capacity < o.size) {
			type = new byte[o.capacity];
			n = new int[o.capacity];
			
			this.capacity = o.capacity;
		}
		
		// Copy
		size = o.size;
		
		for (int index = 0; index < size; ++index) {
			type[index] = o.type[index];
			n[index] = o.n[index];
		}
		
		return;
	}
	
	/**
	 * Clone this type stack.
	 * 
	 * @return A deep-copy clone of this stack.
	 */
	@Override
	public final TypeList clone() {
		return new TypeList(this);
	}
	
	/**
	 * Expand this stack array.
	 * 
	 * @throws IllegalStateException If the stack array is at maximum capacity.
	 */
	private void expand()
			throws IllegalStateException {
		
		int newCapacity;  // New capacity
		byte[] newType;   // New array of types
		int[] newN;       // New array of n
		
		// Get new capacity
		newCapacity = (int) (capacity * KestrelConstants.ARRAY_EXPAND_FACTOR);
		
		if (newCapacity < 0) {
			
			if (capacity == KestrelConstants.MAX_ARRAY_SIZE)
				throw new IllegalStateException("Cannot expand type list: Reached maximum capacity: " + KestrelConstants.MAX_ARRAY_SIZE);
			
			capacity = KestrelConstants.MAX_ARRAY_SIZE;
		}
		
		newType = new byte[newCapacity];
		newN = new int[newCapacity];
		
		// Copy elements
		for (int index = 0; index < size; ++index) {
			newType[index] = type[index];
			newN[index] = n[index];
		}
		
		// Assign fields
		type = newType;
		n = newN;
		capacity = newCapacity;
		
		return;
	}
}
