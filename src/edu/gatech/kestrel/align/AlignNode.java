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

/**
 * Represents one element of an alignment where multiple bases with the same event
 * type are merged into one <code>AlignNode</code> object. For example, 10 aligned
 * bases followed by a mismatch and 5 more aligned bases (16 total) are represented
 * as 3 events (10 matched, 1 mismatch, 5 matched).
 * <p/>
 * The natural representation of this structure is a CIGAR string (see the SAM file
 * specification). In the example above, the cooresponding CIGAR string would be "10=1X5=".
 * Because Kestrel relies on this structure to call variants, the newer standard where
 * "=" and "X" distinguish between matched and mismatched bases is required (the older standard
 * allowed "M" for any base that is not aligned with a gap).
 */
public final class AlignNode implements Comparable<AlignNode> {
	
	/** Type of this alignment event. */
	public final byte type;
	
	/** CIGAR string character for this event. */
	public final char cigarChar;
	
	/** Number of bases in this event. */
	public final int n;
	
	/** Next node in the alignment or <code>null</code> if this is the last element. */
	public final AlignNode next;
	
	/** Assign to <code>type</code> to indicate the two current bases are aligned and they match. */
	public static final byte MATCH = TraceNode.TYPE_MATCH;
	
	/** Assign to <code>type</code> to indicate the two current bases are aligned and they do not match. */
	public static final byte MISMATCH = TraceNode.TYPE_MISMATCH;
	
	/** Assign to <code>type</code> to indicate a gap in the reference sequence at this position. */
	public static final byte INS = TraceNode.TYPE_GAP_REF;
	
	/** Assign to <code>type</code> to indicate a gap in the consensus sequence at this position. */
	public static final byte DEL = TraceNode.TYPE_GAP_CON;
	
	/** Converts type constants to a string. */
	public static final char[] CIGAR_CHARS = TraceNode.getCigarArray();
	
	/**
	 * Create a new alignment node.
	 * 
	 * @param type Type of this alignment event.
	 * @param n The number of bases in this event.
	 * @param next Next node in the alignment.
	 * 
	 * @throws IllegalArgumentException If <code>type</code> is note one of the TYPE
	 *   constants defined in this class.
	 */
	public AlignNode(byte type, int n, AlignNode next)
			throws IllegalArgumentException {
		
		if (type < MATCH || type > DEL)
			throw new IllegalArgumentException("type is out of range [1, 4]: " + type);
		
		this.type = type;
		this.n = n;
		this.cigarChar = CIGAR_CHARS[type];
		this.next = next;
		
		return;
	}
	
	/**
	 * Get a CIGAR string representing an alignment starting from this node.
	 * 
	 * @return A string representing this alignment.
	 */
	public String getCigarString() {
		
		StringBuilder builder = new StringBuilder();
		AlignNode node = this;
		
		while (node != null) {
			builder.append(node.n);
			builder.append(node.cigarChar);
			
			node = node.next;
		}
		
		return builder.toString();
	}
	
	/**
	 * Compare an alignment starting at this node to another alignment starting at another
	 * node, <code>oNode</code>. Returns a value less than <code>0</code> if the alignment
	 * on this node sorts before the alignment on the other node, a value greater than
	 * <code>0</code> if the alignment on the other node sorts before the alignment on this
	 * node, and 0 if the alignments are equal. Note that only the alignment events are
	 * compared, and so they must come from the same haplotype (same reference and consensus
	 * sequences).
	 * <p/>
	 * This comparison represents the canonical order of alternate alignments. It is possible
	 * that a reference and consensus sequence pair has several maximum-score alignments that
	 * explain their relationship. For example, a single base deletion in a homopolymer repeat
	 * (e.g. &quot;TTTTT&quot;) could be explained by a deletion of any one of the bases (any
	 * of the 5 T&apos;s in this example).
	 * <p/>
	 * The first disagreement of the alignments determines which comes first. The order is
	 * mismatch, reference gap, consensus gap, and match. For example, if one sequence contains
	 * a match where another contains a reference gap, the node with the reference gap comes first.
	 * In the homopolymer repeat example, a deletion of the first T comes before the other alignments.
	 * Since Kestrel uses the first alignment as the canonical alignment, it will always report the deletion
	 * of the first T.
	 */
	@Override
	public final int compareTo(AlignNode oNode) {
		
		AlignNode tNode = this;      // This node
		AlignNode oNodeOrg = oNode;  // Save for error reporting
		
		byte tType;  // A type from this alignment
		byte oType;  // A type from the other alignment
		int nDiff;   // Length difference between two nodes
		
		while (oNode != null && tNode != null) {
			
			tType = tNode.type;
			oType = oNode.type;
			nDiff = tNode.n - oNode.n;
			
			// Mismatch type or length, compare first non-matching base
			if (tType != oType || nDiff != 0) {
				
				// Adjust types if one alignment is longer than the other
				if (tType == oType) {
					
					if (nDiff > 0)
						oType = (oNode.next != null) ? oNode.next.type : 0;
					else
						tType = (tNode.next != null) ? tNode.next.type : 0;
					
					assert (oType != tType) :
						"Two sequential alignment nodes have the same type: " + ((nDiff > 0) ? oNodeOrg.getCigarString() : this.getCigarString());
				}

				return ((tType != MATCH) ? tType : 5) -
					   ((oType != MATCH) ? oType : 5);
			}
			
			// Advance
			oNode = oNode.next;
			tNode = tNode.next;
		}
		
		// The shorter of otherwise equal lists sorts first
		if (tNode == null)
			return -1;
		
		if (oNode == null)
			return 1;
		
		// Same alignment
		return 0;
	}
}
