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

/**
 * Describes a region of a reference sequence.
 */
public class RegionInterval implements Comparable<RegionInterval> {
	
	/**
	 * Name of this region. This name will never be <code>null</code>, empty, or contain
	 * tab characters.
	 */
	public final String name;
	
	/**
	 * Start position of this interval (inclusive) where the first base of the sequence
	 * is <code>1</code>.
	 */
	public final int start;
	
	/**
	 * End position of this interval (inclusive) where the first base of the sequence
	 * is <code>1</code>.
	 */
	public final int end;
	
	/**
	 * Name of the sequence this region is in. This name will never be <code>null</code>,
	 * empty, or contain tab characters.
	 */
	public final String sequenceName;
	
	/**
	 * <code>true</code> if the region interval is in the reference orientation.
	 */
	public final boolean isFwd;
	
	/**
	 * Create an interval of a reference region.
	 * 
	 * @param name Name of this region. If <code>null</code> or empty, the sequence name is used
	 *   by replacing whitespace with underscores and appending the start and end interval with
	 *   an underscore and separating start and end with a dash.
	 * @param sequenceName Reference sequence name.
	 * @param start Start location relative to the reference sequence where the first base of the
	 *   reference is <code>1</code>. This location is inclusive (the base is in this region).
	 * @param end End location relative to the reference sequence where the first base of the
	 *   reference is <code>1</code>. This location is inclusive (the base is in this region).
	 * @param isFwd <code>true</code> if this interval is in the same orientation as the reference
	 *   sequence.
	 * 
	 * @throws NullPointerException If <code>sequenceName</code> is <code>null</code>.
	 * @throws IllegalArgumentException If <code>sequenceName</code> is empty, if
	 *   <code>end</code> is not greater than <code>start</code>, or if <code>start</code> or
	 *   <code>end</code> is less than <code>1</code>.
	 */
	public RegionInterval(String name, String sequenceName, int start, int end, boolean isFwd)
			throws NullPointerException, IllegalArgumentException {
		
		// Check arguments
		if (sequenceName == null)
			throw new NullPointerException("Sequence name is null");
		
		sequenceName = sequenceName.trim();
		
		if (sequenceName.isEmpty())
			throw new IllegalArgumentException("Sequence name is empty");
		
		else if (sequenceName.contains("\t"))
			throw new IllegalArgumentException("Sequence name contains tab characters: " + sequenceName);
		
		if (start < 1)
			throw new IllegalArgumentException("Start position is not a positive number: " + start);
		
		if (end < 1)
			throw new IllegalArgumentException("End position is not a positive number: " + end);
		
		if (start > end)
			throw new IllegalArgumentException(String.format("Start position (%d) is greater than the end position (%d)", start, end));
		
		if (name == null)
			name = "";
		
		name = name.trim();
		
		if (name.isEmpty())
			name = String.format("%s_%d-%d", sequenceName.replaceAll("\\s+", "_"), start, end);
		else
			if (name.contains("\t"))
				throw new IllegalArgumentException("Interval name contains tab characters: " + name);
		
		// Set fields
		this.name = name;
		this.sequenceName = sequenceName;
		this.start = start;
		this.end = end;
		this.isFwd = isFwd;
		
		return;
	}
	
	/**
	 * Create an interval of a reference region. This constructor assumes that the reference is in
	 * the same orientation as the reference.
	 * 
	 * @param name Name of this region. If <code>null</code> or empty, the sequence name is used
	 *   by replacing whitespace with underscores and appending the start and end interval with
	 *   an underscore and separating start and end with a dash.
	 * @param sequenceName Reference sequence name.
	 * @param start Start location.
	 * @param end End location.
	 * 
	 * @throws NullPointerException If <code>sequenceName</code> is <code>null</code>.
	 * @throws IllegalArgumentException If <code>sequenceName</code> is empty, if
	 *   <code>end</code> is not greater than <code>start</code>, or if <code>start</code> or
	 *   <code>end</code> is less than <code>1</code>.
	 */
	public RegionInterval(String name, String sequenceName, int start, int end)
			throws NullPointerException, IllegalArgumentException {
		
		this(name, sequenceName, start, end, true);
		
		return;
	}
	
	/**
	 * Create an interval of a reference region and use the values of <code>start</code> and
	 * <code>end</code> to determine if this region is in the same orientation as the reference
	 * (<code>start &lt; end</code>) or a reverse-complement of the reference region
	 * (<code>start &gt; end</code>).
	 * 
	 * @param name Name of this region. If <code>null</code> or empty, the sequence name is used
	 *   by replacing whitespace with underscores and appending the start and end interval with
	 *   an underscore and separating start and end with a dash.
	 * @param sequenceName Reference sequence name.
	 * @param start Start location.
	 * @param end End location.
	 * 
	 * @return A region interval object.
	 * 
	 * @throws NullPointerException If <code>sequenceName</code> is <code>null</code>.
	 * @throws IllegalArgumentException If <code>sequenceName</code> is empty, if
	 *   <code>end</code> is equal to <code>start</code>, or if <code>start</code> or
	 *   <code>end</code> is less than <code>1</code>.
	 */
	public static RegionInterval getAutoInterval(String name, String sequenceName, int start, int end)
			throws NullPointerException, IllegalArgumentException {
		
		if (start > end)
			return new RegionInterval(name, sequenceName, end, start, false);
		
		return new RegionInterval(name, sequenceName, start, end, true);
	}
	
	/**
	 * Get a string representation of this interval.
	 * 
	 * @return A string representation of this interval.
	 */
	@Override
	public String toString() {
		return String.format("RegionInterval[name=%s, start=%d, end=%d, fwd=%b, seq=%s]", name, start, end, isFwd, sequenceName);
	}
	
	/**
	 * Compares the sequence name, start, and end positions (in that order).
	 * 
	 * @return A negative value, <code>0</code>, or a positive value if this interval is
	 *   less than, equal to, or greater than <code>interval</code> (respectively). 
	 */
	@Override
	public int compareTo(RegionInterval interval)
			throws NullPointerException {
		
		// Declarations
		int cmpVal;  // Value from comparing two fields
		
		// Check arguments
		if (interval == null)
			throw new NullPointerException("Cannot compare to interval: null");
		
		// Compare name
		cmpVal = sequenceName.compareTo(interval.sequenceName);
		
		if (cmpVal != 0)
			return cmpVal;
		
		// Compare start
		cmpVal = start - interval.start;
		
		if (cmpVal != 0)
			return cmpVal;
		
		// Compare end
		return end - interval.end;
	}
	
	/**
	 * Determine if this interval overlaps with another.
	 * 
	 * @param interval Interval to test.
	 * 
	 * @return <code>true</code> if both intervals are on the same reference sequence and
	 *   overlap.
	 * 
	 * @throws NullPointerException If <code>interval</code> is <code>null</code>.
	 */
	public boolean overlaps(RegionInterval interval)
			throws NullPointerException {
		
		
		// Check arguments
		if (interval == null)
			throw new NullPointerException("Cannot check for overlaps with interval: null");
		
		return
				sequenceName.equals(interval.sequenceName) &&
				start < interval.end &&
				interval.start < end;
	}
}
