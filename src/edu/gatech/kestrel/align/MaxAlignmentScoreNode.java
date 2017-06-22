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
 * Tracks one or more maximum alignment scores and the information needed to trace
 * back from the location where the maximum score occurred.
 */
public class MaxAlignmentScoreNode {
	
	/** Node where trace-back begins. */
	public final TraceNode traceNode;
	
	/** Number of consensus bases. */
	public final int nConsensusBases;
	
	/** If the maximum score occurred more than once, then this links to the next one. */
	public MaxAlignmentScoreNode next;
	
	/**
	 * Set to <code>true</code> when a haplotype is built from this node. This flag avoids
	 * duplicate haplotypes arising from attempts to save state when a haplotype appears to split
	 * and the state is restored.
	 */
	public boolean haplotypeBuilt;
	
	/**
	 * Create a max alignment score node.
	 * 
	 * @param traceNode Node where trace-back begins.
	 * @param nConsensusBases Number of consensus bases.
	 * @param next If the maximum score occurred more than once, then this links to the next one.
	 * 
	 * @throws NullPointerException If <code>traceNode</code> is <code>null</code>.
	 * @throws IllegalArgumentException if <code>nConsensusBases</code> is less than
	 *   <code>1</code>.
	 */
	public MaxAlignmentScoreNode(TraceNode traceNode, int nConsensusBases, MaxAlignmentScoreNode next) {
		
		// Check arguments
		if (traceNode == null)
			throw new NullPointerException("traceNode is null");
		
		if (nConsensusBases < 1)
			throw new IllegalArgumentException("nConsensusBases is less than 1: " + nConsensusBases);
		
		// Assign fields
		this.traceNode = traceNode;
		this.nConsensusBases = nConsensusBases;
		this.next = next;
		
		haplotypeBuilt = false;
		
		return;
	}
	
	/**
	 * Get a string representation of this max score.
	 * 
	 * @return A string representation of this max score.
	 */
	@Override
	public String toString() {
		return String.format("MaxAlignment[len=%d, last=%s, trace=%s]", nConsensusBases, traceNode.type, traceNode.toString());
	}
}
