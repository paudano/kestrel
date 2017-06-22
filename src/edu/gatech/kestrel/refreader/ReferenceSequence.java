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

import edu.gatech.kestrel.interval.RegionInterval;
import edu.gatech.kestrel.util.digest.Digest;

/**
 * Represents one reference sequence. Variants are called on the whole reference or
 * regions within the reference.
 */
public class ReferenceSequence implements Comparable<ReferenceSequence> {
	
	/** Sequence name. */
	public final String name;
	
	/** Size of this sequence. */
	public final int size;
	
	/** Cryptographic hash of the reference sequence. This is usually an MD5 hash. */
	public final Digest digest;
	
	/**
	 * A string describing where this sequence came from, such as the file name. This may be useful for error
	 * reporting after the sequences are read.
	 */
	public final String sourceName;
	
	/**
	 * Create a new reference sequence.
	 * 
	 * @param name Name of this sequence.
	 * @param size Size of this sequence.
	 * @param digest Digest of the original bases of this sequence or <code>null</code>
	 *   if no digest was performed. When present, this is normally an MD5 crytographic hash.
	 * @param sourceName Name describing where the sequence source came from, such as the
	 *   name of the file it was read from. This is used for error reporting. If <code>null</code>
	 *   or empty, &quot;<UNKNOWN_SOURCE>&quot; is assigned.
	 * 
	 * @throws NullPointerException If <code>name</code> is <code>null</code>.
	 * @throws IllegalArgumentException If <code>name</code> is empty or
	 *   <code>size</code> is less than <code>1</code>.
	 */
	public ReferenceSequence(String name, int size, Digest digest, String sourceName)
			throws NullPointerException, IllegalArgumentException {
		
		// Check arguments
		if (name == null)
			throw new NullPointerException("Sequence name is null");
		
		name = name.trim();
		
		if (name.isEmpty())
			throw new IllegalArgumentException("Sequence name is empty");
		
		if (size < 1)
			throw new IllegalArgumentException("Sequence size is less than 1: " + size);
		
		if (sourceName == null)
			sourceName = "";
		
		sourceName = sourceName.trim();
		
		if (sourceName.isEmpty())
			sourceName = "<UNKNOWN_SOURCE>";
		
		// Assign fields
		this.name = name;
		this.size = size;
		this.digest = digest;
		this.sourceName = sourceName;
	}
	
	/**
	 * Get an interval that spans this reference sequence.
	 * 
	 * @return Region interval that spans this reference sequence.
	 */
	public RegionInterval getReferenceInterval() {
		return new RegionInterval(null, name, 1, size, true);
	}
	
	/**
	 * Get a string representation of this sequence.
	 * 
	 * @return A string representation of this sequence.
	 */
	@Override
	public String toString() {
		return String.format("ReferenceSequence[name=%s, size=%d%s]",
				name,
				size,
				(digest != null) ? ", digest=" + digest.toString() : ""
		);
	}
	
	/**
	 * Determine if this object is equal to another. This test skips the
	 * <code>sourceName</code> field.
	 * 
	 * @param o Object to test.
	 * 
	 * @return <code>true</code> if <code>o</code> is an instance of
	 *   <code>ReferenceSequence</code> and contains equivalent information.
	 */
	@Override
	public boolean equals(Object o) {
		
		ReferenceSequence oRef;
		
		if (o == this)
			return true;
		
		if (o == null)
			return false;
		
		if (! (o instanceof ReferenceSequence))
			return false;
		
		oRef = (ReferenceSequence) o;
		
		return
				name.equals(oRef.name) &&
				size == oRef.size &&
				digest.equals(digest);
	}
	
	/**
	 * Get a hash code for this sequence.
	 * 
	 * @return Hash code for this sequence.
	 */
	@Override
	public int hashCode() {
		return name.hashCode() ^ size ^ digest.hashCode();
	}
	
	/**
	 * Compare this reference sequence to another.
	 * 
	 * @return A negative number, zero, or a positive number if this reference is
	 *   less than, equal to, or greater than <code>o</code>.
	 * 
	 * @throws NullPointerException If <code>o</code> is <code>null</code>.
	 */
	@Override
	public int compareTo(ReferenceSequence o)
			throws NullPointerException {
		
		int cmpVal;
		
		if (o == null)
			throw new NullPointerException("Cannot compare to reference sequence: null");
		
		// Name
		cmpVal = name.compareTo(o.name);
		
		if (cmpVal != 0)
			return cmpVal;
		
		// Size
		cmpVal = size - o.size;
		
		if (cmpVal != 0)
			return cmpVal;
		
		// Digest
		return digest.compareTo(o.digest);
	}
}
