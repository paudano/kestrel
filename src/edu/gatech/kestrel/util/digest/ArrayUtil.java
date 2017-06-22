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

package edu.gatech.kestrel.util.digest;

import java.lang.reflect.Array;

import edu.gatech.kestrel.KestrelConstants;

/**
 * A collection of utility methods for operating on arrays.
 */
public class ArrayUtil {
	
	/**
	 * Hide the default constructor. All methods are static.
	 */
	private ArrayUtil() {
		return;
	}
	
	/**
	 * Expand an array of elements and copy existing elements to the
	 * expanded array.
	 * 
	 * @param <T> Type of elements in the array to be expanded.
	 * @param arr Array to be expanded.
	 * @param arrClass A class object for <code>T</code>. This is required to
	 *   create a new array of type <code>T</code>.
	 * 
	 * @return An expanded array with elements from <code>arr</code> copied
	 *   into it.
	 * 
	 * @throws NullPointerException If <code>arr</code> or <code>arrClass</code>
	 *   is <code>null</code>.
	 * @throws IllegalStateException <code>arr</code> is already at the maximum array
	 *   size, <code>KestrelConstants.MAX_ARRAY_SIZE</code>.
	 */
	@SuppressWarnings("unchecked")
	public static final <T> T[] expandArray(T[] arr, Class<T> arrClass)
			throws NullPointerException, IllegalStateException {
		
		T[] newArr;  // New array
		int newCapacity;  // Capacity of the new array
		int capacity;     // Capacity of the existing array 
		
		// Check arguments
		if (arr == null)
			throw new NullPointerException("Cannot expand array: null");
		
		if (arrClass == null)
			throw new NullPointerException("Cannot expand array with class: null");
		
		// Set capacity
		capacity = arr.length;
		newCapacity = (int) (capacity * KestrelConstants.ARRAY_EXPAND_FACTOR);
		
		if (newCapacity < 0 || newCapacity > KestrelConstants.MAX_ARRAY_SIZE) {
			
			if (capacity == KestrelConstants.ARRAY_EXPAND_FACTOR)
				throw new IllegalStateException("Cannot expand array: Array is at its maximum size: " + KestrelConstants.MAX_ARRAY_SIZE);
			
			newCapacity = KestrelConstants.MAX_ARRAY_SIZE;
		}
		
		// Create array and copy
		newArr = (T[]) Array.newInstance(arrClass, newCapacity);
		
		--capacity;  // Start at capacity - 1
		
		while (capacity >= 0) {
			newArr[capacity] = arr[capacity];
			--capacity;
		}
		
		return newArr;
	}
}
