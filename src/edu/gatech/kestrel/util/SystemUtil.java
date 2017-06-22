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

package edu.gatech.kestrel.util;

import edu.gatech.kestrel.KestrelConstants;

/**
 * A collection of system utilities.
 */
public final class SystemUtil {
	
	/**
	 * Hidden default constructor (all methods are static).
	 */
	private SystemUtil() {
		return;
	}
	
	/**
	 * Emulate <code>Object.toString</code> on an object even if <code>toString()</code> was
	 * overridden.
	 * 
	 * @param o Object.
	 * 
	 * @return A string similar to <code>o.toString()</code> with the <code>Object.toString()</code>
	 *   version of the method, or &quot;null&quot; if <code>o</code> is <code>null</code>.
	 */
	public static String objectToString(Object o) {
		
		if (o == null)
			return "null";
		
		return o.getClass().getName() + "@" + Integer.toHexString(System.identityHashCode(o));
	}
	
	/**
	 * Verify a k-mer size is valid for Kestrel and throw and exception if it
	 * is not.
	 * 
	 * @param kSize K-mer size.
	 * @param msg Additional information prepended to the error message or
	 *   <code>null</code> if there is no additional information to prepend.
	 *   
	 * @throws IllegalArgumentException If <code>kSize</code> is less than
	 *   <code>KestrelConstants.MIN_KMER_SIZE</code>.
	 */
	public static void checkKmerSize(int kSize, String msg)
			throws IllegalArgumentException {
		
		if (kSize < KestrelConstants.MIN_KMER_SIZE) {
			
			if (msg != null)
				throw new IllegalArgumentException(String.format("%s: K-mer size is less than the minimum (%d): %d", msg, KestrelConstants.MIN_KMER_SIZE, kSize)); 
			
			throw new IllegalArgumentException(String.format("K-mer size is less than the minimum (%d): %d", KestrelConstants.MIN_KMER_SIZE, kSize));
		}
		
		return;
	}
	
	/**
	 * Check if a k-mer size is valid for Kestrel.
	 * 
	 * @param kSize K-mer size.
	 * 
	 * @return <code>true</code> if <code>kSize</code> is a valid k-mer size for
	 *   Kestrel.
	 */
	public boolean isValidKmerSize(int kSize) {
		return kSize >= KestrelConstants.MIN_KMER_SIZE;
	}
}
