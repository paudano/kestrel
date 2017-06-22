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

import java.util.Arrays;

/**
 * Represents a cryptographic digest.
 */
public class Digest implements Comparable<Digest> {
	
	/** Digest algorithm. */
	public final String algorithm;
	
	/** Digest bytes. */
	private final byte[] digestBytes;
	
	/** Translates one half of a byte to a hex character. */
	private static final char[] NYBLE_TO_HEX_CHAR = new char[] {
			'0', '1', '2', '3',
			'4', '5', '6', '7',
			'8', '9', 'a', 'b',
			'c', 'd', 'e', 'f'
	};
	
	/**
	 * Create a new digest.
	 * 
	 * @param digestBytes Digest bytes.
	 * @param algorithm Algorithm.
	 * 
	 * @throws NullPointerException If <code>digestBytes</code> or <code>algorithm</code>
	 *   is <code>null</code>.
	 * @throws IllegalArgumentException If <code>digestBytes</code> or <code>algorithm</code>
	 *   is empty.
	 */
	public Digest(byte[] digestBytes, String algorithm)
			throws NullPointerException, IllegalArgumentException {
		
		// Check arguments
		if (digestBytes == null)
			throw new NullPointerException("Digest bytes is a null array");
		
		if (digestBytes.length == 0)
			throw new IllegalArgumentException("Digest bytes is an empty array");
		
		if (algorithm == null)
			throw new NullPointerException("Digest algorithm is null");
		
		algorithm = algorithm.trim();
		
		if (algorithm.isEmpty())
			throw new IllegalArgumentException("Digest algorithm is empty");
		
		// Set fields
		this.digestBytes = Arrays.copyOf(digestBytes, digestBytes.length);
		this.algorithm = algorithm;
		
		return;
	}
	
	/**
	 * Get a string of the digest bytes.
	 * 
	 * @return A string of the digest bytes.
	 */
	@Override
	public String toString() {
		StringBuilder builder = new StringBuilder();
		
		for (byte digestByte : digestBytes) {
			builder.append(NYBLE_TO_HEX_CHAR[(digestByte >>> 4) & 0xf]);
			builder.append(NYBLE_TO_HEX_CHAR[digestByte & 0xf]);
		}
		
		return builder.toString();
	}
	
	/**
	 * Deterimine if an object is equal to this digest.
	 * 
	 * @param o Other object.
	 * 
	 * @return <code>true</code> if <code>o</code> is a digest object and contains
	 *   equivalent information. Returns <code>false</code> if <code>o</code> is
	 *   <code>null</code>. 
	 */
	@Override
	public boolean equals(Object o) {
		
		Digest oDigest;
		
		if (o == this)
			return true;
		
		if (o == null)
			return false;
		
		if (! (o instanceof Digest))
			return false;
		
		oDigest = (Digest) o;
		
		if (! algorithm.equals(oDigest.algorithm))
			return false;
		
		if (digestBytes.length != oDigest.digestBytes.length)
			return false;
		
		for (int index = 0; index < digestBytes.length; ++index)
			if (digestBytes[index] != oDigest.digestBytes[index])
				return false;
		
		return true;
	}
	
	/**
	 * Get a hash code for this digest.
	 * 
	 * @return hash code for this digest.
	 */
	@Override
	public int hashCode() {
		int code = algorithm.hashCode();
		
		for (int index = 0; index < digestBytes.length; ++index)
			code ^= digestBytes[index];
		
		return code;
	}
	
	/**
	 * Compare digest.
	 * 
	 * @param o Other digest.
	 * 
	 * @return A negative number, zero, or a positive number if this digest is
	 *   less than, equal to, or greater than <code>o</code>.
	 * 
	 * @throws NullPointerException If <code>o</code> is <code>null</code>.
	 */
	@Override
	public int compareTo(Digest o)
			throws NullPointerException {
		
		int cmpVal;
		int lastIndex;
		
		if (o == null)
			throw new NullPointerException("Cannot compare to digest: null");
		
		// Algorithm
		cmpVal = algorithm.compareTo(o.algorithm);
		
		if (cmpVal != 0)
			return cmpVal;
		
		// Digest
		lastIndex = digestBytes.length;
		
		if (lastIndex > o.digestBytes.length)
			lastIndex = o.digestBytes.length;
		
		for (int index = 0; index < lastIndex; ++index) {
			if (digestBytes[index] != o.digestBytes[index])
				return digestBytes[index] - o.digestBytes[index];
		}
		
		// Compare length
		if (digestBytes.length == o.digestBytes.length)
			return 0;
		else if (digestBytes.length < o.digestBytes.length)
			return -1;
		else
			return 1;
	}
}
