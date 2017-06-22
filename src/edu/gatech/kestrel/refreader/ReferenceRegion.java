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

import edu.gatech.kanalyze.util.RBTree;
import edu.gatech.kanalyze.util.StringUtil;
import edu.gatech.kestrel.interval.RegionInterval;

/**
 * Represents one reference sequence region where variants are called. Variants are called against reference
 * sequences, and this object contains the sequence, sequence name, and other useful metadata.
 */
public class ReferenceRegion implements Comparable<ReferenceRegion> {
	
	/** Name of this reference region. This is equivalent to <code>interval.name</code>. */
	public final String name;
	
	/** Reference sequence this region belongs to. */
	public final ReferenceSequence referenceSequence;
	
	/** Reference sequence interval this region spans. */
	public final RegionInterval interval;
	
	/** The reference sequence bytes. */
	public final byte[] sequence;
	
	/** Size of this sequence (same as <code>sequence.length</code>). */
	public final int size;
	
	/**
	 * The length of the left flank. If no left flank was added to the interval, this value is <code>0</code>.
	 * Flanks can improve variant calling near the end of intervals when the interval does not reach the end
	 * of the reference sequence. However, the variants that extend into these flanks must be discarded and
	 * the location of the variants must be shifted so they are relative to the start of the reference or the
	 * start of the interval (otherwise, they are relative to <code>sequence</code>, which has these flanks).
	 */
	public final int leftFlankLength;
	
	/**
	 * The index of <code>sequence</code> where the right flank begins. If no right flank was added to
	 * the interval, then this value is <code>size</code>. Flanks can improve variant calling near the end
	 * of intervals when the interval does not reach the end of the reference sequence. However, the variants
	 * that extend into these flanks must be discarded and the location of the variants must be shifted so they
	 * are relative to the start of the reference or the start of the interval (otherwise, they are relative to
	 * <code>sequence</code>, which has these flanks).
	 */
	public final int rightFlankIndex;
	
	/**
	 * Add this number to a one-based coordinate into the reference sequence to obtain the position in
	 * <code>sequence</code> for that base.
	 */
	public final int sequenceOffset;
	
	/** Contains a record of regions with ambiguous bases in the reference region. */
	private RBTree<AmbiguousRegion, Boolean> ambiRegions;  // Tree of regions with ambiguous bases
	
	/** Normalizes a base to a capital IUPAC code (omitting gaps, which are not allowed). */
	private static final byte[] NORM_BASE = new byte[256];
	
	/** Complement a normalized base including ambiguous bases. */
	private static final byte[] COMPL_BASE = new byte[256];
	
	/** Checks a byte value to see if it is ambiguous. */
	private static final boolean[] IS_AMBIGUOUS = new boolean[256];
	
	static {
		
		// NORM_BASE: Single base
		NORM_BASE['A'] = 'A'; NORM_BASE['a'] = 'A';
		NORM_BASE['C'] = 'C'; NORM_BASE['c'] = 'C';
		NORM_BASE['G'] = 'G'; NORM_BASE['g'] = 'G';
		NORM_BASE['T'] = 'T'; NORM_BASE['t'] = 'T';
		NORM_BASE['U'] = 'U'; NORM_BASE['u'] = 'u';
		
		// NORM_BASE: Two bases
		NORM_BASE['R'] = 'R'; NORM_BASE['r'] = 'R';
		NORM_BASE['Y'] = 'Y'; NORM_BASE['y'] = 'Y';
		NORM_BASE['S'] = 'S'; NORM_BASE['s'] = 'S';
		NORM_BASE['W'] = 'W'; NORM_BASE['w'] = 'W';
		NORM_BASE['K'] = 'K'; NORM_BASE['k'] = 'K';
		NORM_BASE['M'] = 'M'; NORM_BASE['m'] = 'M';
		
		// NORM_BASE: Three bases
		NORM_BASE['B'] = 'B'; NORM_BASE['b'] = 'B';
		NORM_BASE['D'] = 'D'; NORM_BASE['d'] = 'D';
		NORM_BASE['H'] = 'H'; NORM_BASE['h'] = 'H';
		NORM_BASE['V'] = 'V'; NORM_BASE['v'] = 'V';
		
		// NORM_BASE: Any base
		NORM_BASE['N'] = 'N'; NORM_BASE['n'] = 'N';
		
		// IS_NOT_AMBIGUOUS
		for (int index = 0; index < 256; ++index)
			IS_AMBIGUOUS[index] = true;
		
		IS_AMBIGUOUS['A'] = false;
		IS_AMBIGUOUS['C'] = false;
		IS_AMBIGUOUS['G'] = false;
		IS_AMBIGUOUS['T'] = false;
		IS_AMBIGUOUS['U'] = false;
	}
	
	static {
		
		// Initialize all bases to N
		for (int index = 0; index < COMPL_BASE.length; ++index)
			COMPL_BASE[index] = 'N';
		
		// COMPL_BASE: Upper-case
		COMPL_BASE['A'] = 'T';
		COMPL_BASE['C'] = 'G';
		COMPL_BASE['G'] = 'C';
		COMPL_BASE['T'] = 'A';
		COMPL_BASE['U'] = 'A';
		COMPL_BASE['R'] = 'Y';
		COMPL_BASE['Y'] = 'R';
		COMPL_BASE['S'] = 'S';
		COMPL_BASE['W'] = 'W';
		COMPL_BASE['K'] = 'M';
		COMPL_BASE['M'] = 'K';
		COMPL_BASE['B'] = 'V';
		COMPL_BASE['D'] = 'H';
		COMPL_BASE['H'] = 'D';
		COMPL_BASE['V'] = 'B';
		COMPL_BASE['N'] = 'N';
		
		// COMPL_BASE: Lower-case
		COMPL_BASE['a'] = 'T';
		COMPL_BASE['c'] = 'G';
		COMPL_BASE['g'] = 'C';
		COMPL_BASE['t'] = 'A';
		COMPL_BASE['u'] = 'A';
		COMPL_BASE['r'] = 'Y';
		COMPL_BASE['y'] = 'R';
		COMPL_BASE['s'] = 'S';
		COMPL_BASE['w'] = 'W';
		COMPL_BASE['k'] = 'M';
		COMPL_BASE['m'] = 'K';
		COMPL_BASE['b'] = 'V';
		COMPL_BASE['d'] = 'H';
		COMPL_BASE['h'] = 'D';
		COMPL_BASE['v'] = 'B';
		COMPL_BASE['n'] = 'N';
	}
	
	/**
	 * Create a new reference region defined as an interval over a reference sequence.
	 * 
	 * @param referenceSequence Reference sequence this region belongs to.
	 * @param interval Interval of the reference sequence is region spans or <code>null</code>
	 *   if this region spans the entire reference sequence.
	 * @param buffer Buffer containing the sequence bytes.
	 * @param bufferOffset The offset from the start of <code>buffer</code> to start copying
	 *   bases.
	 * @param leftFlankLength Number of bases before <code>interval</code> copied to the sequence.
	 *   All bytes of both flanks and the region interval must be in <code>buffer</code> starting
	 *   at <code>bufferOffset</code>.
	 * @param rightFlankLength Number of bases after <code>interval</code> copied to the sequence.
	 *   All bytes of both flanks and the region interval must be in <code>buffer</code> starting
	 *   at <code>bufferOffset</code>.
	 *   
	 * @throws NullPointerException If <code>referenceSequence</code> or <code>buffer</code>
	 *   is <code>null</code>.
	 * @throws IllegalArgumentException If <code>bufferOffset</code> is negative, either flank length is
	 *   negative, or the buffer is too small to contain this reference region and both flanks starting
	 *   from <code>bufferOffset</code>.
	 */
	public ReferenceRegion(ReferenceSequence referenceSequence, RegionInterval interval, byte[] buffer, int bufferOffset, int leftFlankLength, int rightFlankLength)
			throws NullPointerException, IllegalArgumentException {
		
		long size;  // Size of this region
		
		// Check arguments
		if (referenceSequence == null)
			throw new NullPointerException("Reference sequence is null");
		
		if (interval == null)
			interval = referenceSequence.getReferenceInterval();
		
		if (buffer == null)
			throw new NullPointerException("Cannot copy reference sequence from buffer: null");
		
		if (bufferOffset < 0)
			throw new IllegalArgumentException("Buffer offset is negative: " + bufferOffset);
		
		if (leftFlankLength < 0)
			throw new IllegalArgumentException("Left flank length is negative: " + leftFlankLength);
		
		if (rightFlankLength < 0)
			throw new IllegalArgumentException("Right flank length is negative: " + rightFlankLength);
		
		size = (long) interval.end - interval.start + 1 + leftFlankLength + rightFlankLength;
		
		if (size > Integer.MAX_VALUE)
			throw new IllegalArgumentException(String.format("Reference region is too large (max = %d): %d", Integer.MAX_VALUE, size));
		
		if (buffer.length < bufferOffset + size)
			throw new IllegalArgumentException(String.format("Buffer is not large enough to contain this reference region: Buffer size = %d, offset = %d, region size = %d, flank size = %d", buffer.length, bufferOffset, size, (leftFlankLength + rightFlankLength)));
		
		// Set fields
		this.name = interval.name;
		this.referenceSequence = referenceSequence;
		this.interval = interval;
		this.size = (int) size;
		this.leftFlankLength = leftFlankLength;
		this.rightFlankIndex = (int) size - rightFlankLength;
		
		// Copy buffer
		this.sequence = copySequenceFromBuffer(interval, buffer, bufferOffset, (int) size);  // throws IllegalArgumentException
		
		// Set ambiguous regions (must be called after sequence and size fields are set)
		setAmbiguousRegions();
		
		// Set offset
		sequenceOffset = getSequenceOffset();
		
		return;
	}
	
	/**
	 * Create a new reference sequence that spans an entire reference sequence.
	 * 
	 * @param referenceSequence Reference sequence this region belongs to.
	 * @param buffer Buffer containing the sequence bytes.
	 * @param bufferOffset The offset from the start of <code>buffer</code> to start copying
	 *   bases.
	 *   
	 * @throws NullPointerException If <code>referenceSequence</code> or <code>buffer</code>
	 *   is <code>null</code>.
	 * @throws IllegalArgumentException If <code>bufferOffset</code> is negative, either flank length is
	 *   negative, or the buffer is too small to contain this reference region and both flanks starting
	 *   from <code>bufferOffset</code>.
	 */
	public ReferenceRegion(ReferenceSequence referenceSequence, byte[] buffer, int bufferOffset)
			throws NullPointerException, IllegalArgumentException {
		
		long size;  // Size of this region
		
		// Check arguments
		if (referenceSequence == null)
			throw new NullPointerException("Reference sequence is null");
		
		if (buffer == null)
			throw new NullPointerException("Cannot copy reference sequence from buffer: null");
		
		if (bufferOffset < 0)
			throw new IllegalArgumentException("Buffer offset is negative: " + bufferOffset);
		
		size = (long) referenceSequence.size + bufferOffset;
		
		if (buffer.length < bufferOffset + size)
			throw new IllegalArgumentException(String.format("Buffer is not large enough to contain this reference sequence: Buffer size = %d, offset = %d: %s", buffer.length, bufferOffset, referenceSequence.toString()));
		
		// Set fields
		this.name = referenceSequence.name;
		this.referenceSequence = referenceSequence;
		this.interval = referenceSequence.getReferenceInterval();
		this.size = (int) size;
		this.leftFlankLength = 0;
		this.rightFlankIndex = (int) size;
		
		// Copy buffer
		this.sequence = copySequenceFromBuffer(interval, buffer, bufferOffset, (int) size);  // throws IllegalArgumentException
		
		// Set ambiguous regions (must be called after sequence and size fields are set)
		setAmbiguousRegions();
		
		// Set offset
		sequenceOffset = getSequenceOffset();
		
		return;
	}
	
	/**
	 * Create a new reference sequence from a reference interval.
	 * 
	 * @param incompleteRegion Incomplete region object.
	 * @param referenceSequence Reference sequence.
	 * @param sequence Sequence bytes.
	 */
	private ReferenceRegion(ReferenceRegion.IncompleteRegion incompleteRegion, ReferenceSequence referenceSequence, byte[] sequence) {
		
		// Check arguments
		assert (incompleteRegion != null) :
			"Incomplete region object is null";
		
		assert (referenceSequence != null) :
			"Reference sequence is null";
		
		assert (sequence != null) :
			"Sequence is null";
		
		// Set fields
		this.name = incompleteRegion.interval.name;
		this.referenceSequence = referenceSequence;
		this.interval = incompleteRegion.interval;
		this.sequence = sequence;
		this.leftFlankLength = incompleteRegion.leftFlankLength;
		this.rightFlankIndex = incompleteRegion.rightFlankIndex;
		this.size = incompleteRegion.size;
		
		// Set ambiguous regions (must be called after sequence and size fields are set)
		setAmbiguousRegions();
		
		// Set offset
		sequenceOffset = getSequenceOffset();
		
		return;
	}
	
	/**
	 * Reverse complement the sequence in this reference.
	 */
	public void reverseComplement() {
		
		int index = 0;
		int size = this.size;
		int lastIndex = size / 2;
		byte temp;
		
		// Reverse complement
		while (index < lastIndex) {
			temp = COMPL_BASE[sequence[index]];
			sequence[index] = COMPL_BASE[sequence[size - index - 1]];
			sequence[size - index - 1] = temp;
			
			++index;
		}
		
		// Complement the middle base (did not change position)
		if (size % 2 != 0)
			sequence[lastIndex] = COMPL_BASE[sequence[lastIndex]];
		
		return;
	}
	
	/**
	 * Get a string description of this reference sequence.
	 * 
	 * @return A string description of this reference sequence.
	 */
	@Override
	public String toString() {
		return String.format("ReferenceRegion[name=%s, size=%d]", name, size);
	}
	
	/**
	 * Determine if a range within this region contains ambiguous bases using indices
	 * in the sequence field where the first base is this region <code>0</code>.
	 * 
	 * @param startIndex Start index in <code>sequence</code>.
	 * @param endIndex End index in <code>sequence</code>.
	 * 
	 * @return <code>true</code> if the region contains at least one ambiguous base.
	 * 
	 * @throws IllegalArgumentException If <code>startIndex</code> is less than <code>0</code>
	 *   or <code>endIndex</code> is less than <code>startIndex</code>.
	 */
	public boolean containsAmbiguousByIndex(int startIndex, int endIndex)
			throws IllegalArgumentException {
		
		if (startIndex < 0 || endIndex < startIndex)
			throw new IllegalArgumentException(String.format("startTndex must be greater than 0, and endIndex must not be lest than startIndex: startIndex=%d, endIndex=%d", startIndex, endIndex));
		
		return ambiRegions.containsKey(new AmbiguousRegion(startIndex, endIndex));
	}
	
	/**
	 * Determine if a range within this region contains ambiguous bases using coordinates
	 * of the bases in the sequence where the first base is <code>1</code>.
	 * 
	 * @param startCoordinate Start coordinate relative to this region.
	 * @param endCoordinate End coordinate relative to this region.
	 * 
	 * @return <code>true</code> if the region contains at least one ambiguous base.
	 * 
	 * @throws IllegalArgumentException If <code>startCoordinate</code> is less than <code>1</code>
	 *   or <code>endCoordinate</code> is less than <code>startCoordinate</code>.
	 */
	public boolean containsAmbiguousByBaseCoordinate(int startCoordinate, int endCoordinate)
			throws IllegalArgumentException {
		
		if (startCoordinate < 0 || endCoordinate < startCoordinate)
			throw new IllegalArgumentException(String.format("startCoordinate must be greater than 0, and endCoordinate must not be lest than startCoordinate: startCoordinate=%d, endCoordinate=%d", startCoordinate, endCoordinate));
		
		return ambiRegions.containsKey(new AmbiguousRegion(startCoordinate - 1, endCoordinate - 1));
	}
	
	/**
	 * Determine if a variant beginning on index <code>start</code> and ending on index
	 * <code>end</code> spans into the left or right flank of this reference region.
	 * 
	 * @param start Start index of the variant where the first base is <code>0</code> (inclusive).
	 * @param end End index of the variant where the first base is <code>0</code> (inclusive).
	 * 
	 * @return <code>true</code> if this variant spans into a flank of this reference region.
	 * 
	 * @throws IllegalArgumentException If <code>start</code> is less than <code>0</code> or
	 *   <code>end</code> is less than <code>start</code>.
	 */
	public boolean isFlankByIndex(int start, int end)
			throws IllegalArgumentException {
		
		if (start < 0 || end < start)
			throw new IllegalArgumentException(String.format("Start and end locations must be 0 or greater and end must not be less than start: start=%d, end=%d", start, end));
		
		return (start < leftFlankLength || end >= rightFlankIndex);
	}
	
	/**
	 * Determine if a variant beginning on base <code>start</code> and ending on base
	 * <code>end</code> spans into the left or right flank of this reference region.
	 * 
	 * @param start Start coordinate of the variant where the first base is <code>1</code> (inclusive).
	 * @param end End coordinate of the variant where the first base is <code>1</code> (inclusive).
	 * 
	 * @return <code>true</code> if this variant spans into a flank of this reference region.
	 * 
	 * @throws IllegalArgumentException If <code>start</code> is less than <code>1</code> or
	 *   <code>end</code> is less than <code>start</code>.
	 */
	public boolean isFlankByCoordinate(int start, int end)
			throws IllegalArgumentException {
		
		if (start < 1 || end < start)
			throw new IllegalArgumentException(String.format("Start and end locations must be 1 or greater and end must not be less than start: start=%d, end=%d", start, end));
		
		return (start <= leftFlankLength || (end - 1) >= rightFlankIndex);
	}
	
	/**
	 * Compares this reference sequence on sequence name.
	 */
	@Override
	public int compareTo(ReferenceRegion refRegion) {
		
		if (refRegion == null)
			throw new NullPointerException("Cannot compare reference region: null");
		
		return interval.compareTo(refRegion.interval);
	}
	
	/**
	 * Get a base from this reference region.
	 * 
	 * @param location Location of the base in the reference sequence coordinate system (1-based).
	 * 
	 * @return Base at <code>location</code> in the reference sequence.
	 * 
	 * @throws IllegalArgumentException If <code>location</code> is not contained in this reference region. A location
	 *   in the flanking regions is allowed since that base is in <code>sequence</code>.
	 */
	public byte getBase(int location)
			throws IllegalArgumentException {
		
		// Get sequence index
		int seqIndex = location + sequenceOffset;
		
		// Check arguments
		if (seqIndex < 0 || seqIndex >= size) {
			throw new IllegalArgumentException(String.format(
					"Base location %d is out of bounds for interval %s with flanks (l=%d, r=%d)",
					location, interval, leftFlankLength, size - rightFlankIndex
			));
		}
		
		return sequence[seqIndex];
	}
	
	/**
	 * Get the offset of sequences in the reference coordinate (1-based) to the index of
	 * <code>sequence</code> of that base (0-based).
	 *  
	 * @return Offset to translate a reference coordinate to a sequence index.
	 */
	private int getSequenceOffset() {
		
		return -(interval.start - leftFlankLength);
	}
	
	/**
	 * Find and mark ambiguous regions in the reference sequence. This must be called
	 * after the <code>sequence</code> and <code>size</code> fields are set.
	 */
	private void setAmbiguousRegions() {
		
		// Declarations
		int seqIndex;
		int ambiEndIndex;
		int seqSize;
		byte[] sequence;
		
		// Check state
		assert (this.sequence != null) :
			"Reference sequence was not set before setAmbiguousRegions() was callled";
		
		// Create ambiguous region structure
		ambiRegions = new RBTree<>();
		sequence = this.sequence;
		seqSize = this.size;
		
		// Scan
		seqIndex = 0;
		
		while (seqIndex < seqSize) {
			
			if (IS_AMBIGUOUS[sequence[seqIndex]]) {
				
				ambiEndIndex = seqIndex;
				
				while (ambiEndIndex < seqSize && IS_AMBIGUOUS[sequence[ambiEndIndex]])
					++ambiEndIndex;
				
				if (ambiEndIndex == seqSize)
					ambiEndIndex -= 1;
				
				ambiRegions.put(new AmbiguousRegion(seqIndex, ambiEndIndex), Boolean.TRUE);
				
				seqIndex = ambiEndIndex;
			}
			
			++seqIndex;
		}
		
		return;
	}
	
	/**
	 * Copy a sequence from a buffer.
	 * 
	 * @param interval Interval defining this region.
	 * @param buffer Buffer.
	 * @param offset Start offset.
	 * @param size Number of bytes to copy.
	 * 
	 * @return An array of sequence bytes containing normalized IUPAC bases (all upper-case).
	 * 
	 * @throws IllegalArgumentException If the sequence bytes contains gaps or non-IUPAC
	 *   bases.
	 */
	private static byte[] copySequenceFromBuffer(RegionInterval interval, byte[] buffer, int offset, int size)
			throws IllegalArgumentException {
		
		byte[] sequence;
		byte baseByte;
		
		// Check arguments
		assert (interval != null) :
			"interval is null";
		
		assert (buffer != null) :
			"buffer is null";
		
		assert (offset >= 0) :
			"Offset is negative: " + offset;
		
		assert (size > 0) :
			"size is negative: " + size;
		
		// Init
		sequence = new byte[size];
		
		// Copy buffer
		for (int index = 0; index < size; ++index) {
			baseByte = NORM_BASE[buffer[index + offset]];
			
			// Check base (must be IUPAC and not a gap)
			if (baseByte == 0) {
				
				if (buffer[index] == '.' || buffer[index] == '-')
					throw new IllegalArgumentException(String.format("Found a gap in reference region %s (reference sequence %s) at location %d (reference offset %d): %c", interval.name, interval.sequenceName, (index + 1), (index + interval.start), buffer[index]));
				
				throw new IllegalArgumentException(String.format("Found non-IUPAC character in reference region %s (reference sequence %s) at location %d (reference offset %d): %s", interval.name, interval.sequenceName, (index + 1), (index + interval.start), StringUtil.charDescription(buffer[index])));
			}
			
			sequence[index] = baseByte;
		}
		
		return sequence;
	}
	
	/**
	 * Represents an ambiguous region designed to be stored in a tree structure. When building
	 * the tree, it is assumed the reference regions will not overlap. For querying reference
	 * regions, this object will claim it is equal to any other where the regions overlap.
	 */
	private static class AmbiguousRegion implements Comparable<AmbiguousRegion> {
		
		/** Start position of this region (inclusive) where the first base is <code>0</code>. */
		public final int start;
		
		/** End position of this region (inclusive) where the first base is <code>0</code>. */
		public final int end;
		
		/**
		 * Create a reference region.
		 * 
		 * @param start Start position of this region (inclusive) where the first base
		 *   is <code>0</code>.
		 * @param end End position of this region (inclusive) where the first base is
		 *   <code>0</code>.
		 */
		public AmbiguousRegion(int start, int end) {
			
			assert (start >= 0 && end >= start) :
				String.format("Bad start and end positions in ReferenceRegion(): start=%d, end=%d", start, end);
			
			this.start = start;
			this.end = end;
			
			return;
		}
		
		/**
		 * Compare this region to another. If the regions overlap, this method claims they are the same
		 * even if the start and end positions are different.
		 * 
		 * @param o Other region.
		 * 
		 * @return <code>-1</code> if this region comes before <code>o</code>, <code>1</code> if it
		 *   comes after, and <code>0</code> if the regions overlap.
		 */
		@Override
		public int compareTo(AmbiguousRegion o) {
			
			assert (o != null) :
				"Cannot compare to object: null";
			
			if (o.end < start)
				return 1;
			
			if (o.start > end)
				return -1;
			
			return 0;
		}
	}
	
	/**
	 * Stores data for a reference region and saves it until the reference region can be created
	 * from it. Only one reference region may be created from one of these objects since the
	 * internal buffers are transferred to the region object.
	 * <p/>
	 * Using these incomplete regions, the sequence from the reference regions may be stored
	 * while the reference sequence is being read. If a digest is used on the reference sequence,
	 * the reference sequence object cannot be created until the entire sequence is read. This
	 * object caches the sequence data until the reference sequence object is available.  
	 */
	public static class IncompleteRegion {
		
		/** Interval defining this region. */
		public final RegionInterval interval;
		
		/** Size of the sequence. */
		public final int size;
		
		/**
		 * The length of the right flank. If no right flank was added to the interval, this value is <code>0</code>.
		 * Flanks can improve variant calling near the end of intervals when the interval does not reach the end
		 * of the reference sequence. However, the variants that extend into these flanks must be discarded and
		 * the location of the variants must be shifted so they are relative to the start of the reference or the
		 * start of the interval (otherwise, they are relative to <code>sequence</code>, which has these flanks).
		 */
		public final int leftFlankLength;
		
		/**
		 * The index of <code>sequence</code> where the right flank begins. If no right flank was added to
		 * the interval, then this value is <code>size</code>. Flanks can improve variant calling near the end
		 * of intervals when the interval does not reach the end of the reference sequence. However, the variants
		 * that extend into these flanks must be discarded and the location of the variants must be shifted so they
		 * are relative to the start of the reference or the start of the interval (otherwise, they are relative to
		 * <code>sequence</code>, which has these flanks).
		 */
		public final int rightFlankIndex;
		
		/** The reference sequence bytes. */
		private byte[] sequence;
		
		/**
		 * Create an interval for a reference region to store interval data before the reference sequence
		 * object is created. This object may be used to create one (but only one) reference region since
		 * the internal buffers are transfered to the region.
		 * 
		 * @param interval Interval that defines this region.
		 * @param buffer A buffer to read the sequence bytes from. The entire interval must be in this
		 *   region.
		 * @param bufferOffset Offset in <code>buffer</code> to start reading from.
		 * @param leftFlankLength Number of bases before <code>interval</code> copied to the sequence.
		 *   All bytes of both flanks and the region interval must be in <code>buffer</code> starting
		 *   at <code>bufferOffset</code>.
		 * @param rightFlankLength Number of bases after <code>interval</code> copied to the sequence.
		 *   All bytes of both flanks and the region interval must be in <code>buffer</code> starting
		 *   at <code>bufferOffset</code>.
		 * 
		 * @throws NullPointerException If <code>interval</code> or <code>buffer</code> is <code>null</code>.
		 * @throws IllegalArgumentException If <code>bufferOffset</code> is negative or the buffer is not large
		 *   enough to read the interval starting at <code>bufferOffset</code>.
		 */
		public IncompleteRegion(RegionInterval interval, byte[] buffer, int bufferOffset, int leftFlankLength, int rightFlankLength)
				throws NullPointerException, IllegalArgumentException {
			
			long size;         // Length of this interval
			byte[] sequence;  // Sequence bytes
			
			// Check arguments
			if (interval == null)
				throw new NullPointerException("Interval is null");
			
			if (buffer == null)
				throw new NullPointerException("Sequence buffer is null");
			
			if (bufferOffset < 0)
				throw new IllegalArgumentException("Buffer offset is negative: " + bufferOffset);
			
			if (leftFlankLength < 0)
				throw new IllegalArgumentException("Left flank length is negative: " + leftFlankLength);
			
			if (rightFlankLength < 0)
				throw new IllegalArgumentException("Right flank length is negative: " + rightFlankLength);
			
			size = (long) interval.end - interval.start + 1 + leftFlankLength + rightFlankLength;
			
			if (size > Integer.MAX_VALUE)
				throw new IllegalArgumentException(String.format("Reference region is too large (max = %d): %d", Integer.MAX_VALUE, size));
			
			if (buffer.length < bufferOffset + size)
				throw new IllegalArgumentException(String.format("Buffer is not large enough to contain this reference region: Buffer size = %d, offset = %d, region size = %d, flank size = %d", buffer.length, bufferOffset, size, (leftFlankLength + rightFlankLength)));
			
			// Copy sequence
			sequence = copySequenceFromBuffer(interval, buffer, bufferOffset, (int) size);
			
			// Set fields
			this.interval = interval;
			this.sequence = sequence;
			this.leftFlankLength = leftFlankLength;
			this.rightFlankIndex = (int) size - rightFlankLength;
			this.size = (int) size;
			
			return;
		}
		
		/**
		 * Get a reference region from this interval object. This method can be called only once
		 * per object.
		 * 
		 * @param referenceSequence Reference sequence.
		 * 
		 * @return A reference interval object.
		 * 
		 * @throws NullPointerException If <code>referenceSequence</code> is <code>null</code>.
		 * @throws IllegalStateException If this method has already been called on this object.
		 */
		public ReferenceRegion getRegion(ReferenceSequence referenceSequence)
				throws NullPointerException, IllegalStateException {
			
			ReferenceRegion region;
			
			// Check arguments and state
			if (sequence == null)
				throw new IllegalStateException("This interval has been used to create a reference sequence and cannot be used to create another");
			
			if (referenceSequence == null)
				throw new NullPointerException("Reference sequence is null");
			
			// Create region
			region = new ReferenceRegion(this, referenceSequence, sequence);
			
			// Free sequence
			sequence = null;
			
			// Return this region
			return region;
		}
	}
}
