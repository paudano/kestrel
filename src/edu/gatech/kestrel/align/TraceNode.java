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

import java.util.Arrays;

/**
 * The alignment trace-back matrix is represented as a linked list wap of
 * <code>TraceNode</code> objects. This map is build from the root out to the branches,
 * and it is traversed in reverse order from the branches back to the root. A matrix
 * is normally used to implement the Smith-Waterman and Needleman-Wunsch alignment algorithms.
 * By storing it as a linked map, empty cells with no transitions do not take memory.
 * <p/>
 * For alignments that allow an affine-gap (different gap open and gap extend penalties), three
 * score matrices are used to build the alignment, and the traceback matrix must represent transitions
 * among these matrices. For example, transitioning from the alignment matrix to the reference-gap
 * matrix represents the boundary of a deletion in the reference sequence.
 * <p/>
 * This node may be associated with any of the three score matrices (align, reference-gap, and
 * consensus-gap), and it may transition to one or more of these three matrices.  
 */
public final class TraceNode {
	
	/** Score of the alignment at this point. */
	public final float score;
	
	/** Type of this node (align, gap reference, or gap consensus). */
	public final byte type;
	
	/** Next node. */
	public final TraceNode nextNode;
	
	/**
	 * If the map branches at this position in the trace-back matrix on its traversal
	 * back to the root, then this points to the node it branches back to. The node this
	 * link points to is in the same position in the matrix as this (current object) node. 
	 */
	public final TraceNode branchNode;
	
	/**
	 * Assign to <code>type</code> to indicate that this node goes nowhere. This should only be
	 * used by <code>ZERO_NODE</code>.
	 */
	public static final byte TYPE_NONE = 0;
	
	/** Assign to <code>type</code> to indicate the two current bases are aligned and they match. */
	public static final byte TYPE_MATCH = 1;
	
	/** Assign to <code>type</code> to indicate the two current bases are aligned and they do not match. */
	public static final byte TYPE_MISMATCH = 2;
	
	/** Assign to <code>type</code> to indicate a gap in the reference sequence at this position. */
	public static final byte TYPE_GAP_REF = 3;
	
	/** Assign to <code>type</code> to indicate a gap in the consensus sequence at this position. */
	public static final byte TYPE_GAP_CON = 4;
	
	/** Maps type constants to strings for <code>toString()</code>. */
	private static final String[] TYPE_STRING = new String[] {"NONE", "ALIGN_MATCH", "ALIGN_MISMATCH", "GAP_REFERENCE", "GAP_CONSENSUS"};
	
	/**
	 * Type constant to CIGAR. Note that TYPE_CIGAR[0] (*) should never show up in a CIGAR string, but it is
	 * in this array as a placeholder since TYPE_NONE = 0.
	 */
	private static final char[] CIGAR_CHARS = new char[] {'*', '=', 'X', 'I', 'D'};
	
	
	
	/**
	 * Indicates a dead-end in the map. This node is used by dynamic alignments in matrix positions where
	 * there is no valid node, and it avoids checking for <code>null</code> references.
	 */
	public static final TraceNode ZERO_NODE = new TraceNode(0.0F, TYPE_NONE, null, null);
	
	/**
	 * Create a new trace node.
	 * 
	 * @param score Score of the alignment at this point. This is used while building the alignment.
	 * @param type Type of node (align or gap).
	 * @param nextNode Next node in the map after applying <code>type</code>.
	 * @param branchNode Branch to another node at this position (do not apply <code>type</code>).
	 * 
	 * @throws IllegalArgumentException If <code>score</code> is negative or <code>type</code> is
	 *   out of range.
	 */
	public TraceNode(float score, byte type, TraceNode nextNode, TraceNode branchNode)
			throws IllegalArgumentException {
		
		if (score < 0.0F)
			throw new IllegalArgumentException("score is negative");
		
		if (type < TYPE_NONE || type > TYPE_GAP_CON)
			throw new IllegalArgumentException(String.format("type is out of range [%d, %d]: %d", TYPE_NONE, TYPE_GAP_CON, type));
		
		this.score = score;
		this.type = type;
		this.nextNode = nextNode;
		this.branchNode = branchNode;
		
		return;
	}
	
	/**
	 * Create a new trace node.
	 * 
	 * @param score Score of the alignment at this point. This is used while building the alignment.
	 * @param type Type of node (align or gap).
	 * @param nextNode Next node in the map after applying <code>type</code>.
	 */
	public TraceNode(float score, byte type, TraceNode nextNode) {
		
		this.score = score;
		this.type = type;
		this.nextNode = nextNode;
		this.branchNode = null;
		
		return;
	}
	
	/**
	 * Get an array of CIGAR characters indexed by the <code>TYPE</code> constants defined
	 * on this class.
	 * 
	 * @return Array of CIGAR characters.
	 */
	public static final char[] getCigarArray() {
		return Arrays.copyOf(CIGAR_CHARS, CIGAR_CHARS.length);
	}
	
	/**
	 * Get a string description of this trace node.
	 * 
	 * @return A string description of this trace node.
	 */
	@Override
	public String toString() {
		return String.format("TraceNode[score=%f, type=%s]", score, TYPE_STRING[type]);
	}
}
