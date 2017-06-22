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

package edu.gatech.kestrel.align.state;

import edu.gatech.kestrel.align.TraceNode;

/**
 * A node in a linked list of cached trace nodes. The alignment, reference gap,
 * and consensus gap matrices are each cached with one of these lists. Links to
 * the default trace node (zero score) are not stored.
 */
public final class TraceNodeContainer {
	
	/** Index in the matrix column where this node belongs. */
	public final int index;
	
	/** Trace node. */
	public final TraceNode node;
	
	/** Next container in this list. */
	public final TraceNodeContainer next;
	
	/**
	 * Create a new node container for one <code>TraceNode</code>.
	 * 
	 * @param index Index in the matrix column where this node is restored to.
	 * @param node Node.
	 * @param next Next node in this list.
	 */
	public TraceNodeContainer(int index, TraceNode node, TraceNodeContainer next) {
		
		// Check arguments
		assert (index >= 0) :
			"TraceNodeContainer(): Negative index: " + index;
		
		assert (node != null) :
			"TraceNodeContainer(): node is null";
		
		// Assign fields
		this.index = index;
		this.node = node;
		this.next = next;
		
		return;
	}
}
