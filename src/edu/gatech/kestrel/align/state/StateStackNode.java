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

import edu.gatech.kanalyze.util.Base;
import edu.gatech.kestrel.align.MaxAlignmentScoreNode;

/**
 * When an alignment splits, this node stores information needed to restore the
 * state of this alignment for backtracking. This occurs when more than one k-mer
 * may be added to the alignment at once (for example, evidence was found for two
 * variants at the same position, or one variant and the wild-type are both present).
 * <p/>
 * This is stored as a stack because when multiple states are stored, back-tracking
 * most go back to the last node found. Doing this out of order will result in incorrect
 * data (for example, part the consensus sequence may have been overwritten).  
 */
public final class StateStackNode {
	
	/**
	 * K-mer that the aligner should start with for the next base, and the next base
	 * to be added must be the last base of this k-mer.
	 */
	public final int[] kmer;
	
	/** Base that should be added when this state is restored. */
	public final Base nextBase;
	
	/** Number of bases in the consensus at this point of the alignment. */
	public final int consensusSize;
	
	/** A list of nodes in the alignment column. */
	public final TraceNodeContainer alignContainer;
	
	/** A list of nodes in the reference-gap column. */
	public final TraceNodeContainer gapRefContainer;
	
	/** A list of nodes in the consensus-gap column. */
	public final TraceNodeContainer gapConContainer;
	
	/** The maximum score of this alignment. */
	public final float maxAlignmentScore;
	
	/** Links to the start of maximum-score alignments. */ 
	public final MaxAlignmentScoreNode maxAlignmentScoreNode;
	
	/** Minimum depth for all k-mers in the alignment up to this state. */
	public final int minDepth;
	
	/**
	 * Links to the next node down the stack. This is the node that was saved before this one,
	 * and it will be the next node to be restored after this one.
	 */
	public StateStackNode nextNodeDown;
	
	/** Links to the next node up the stack. This is the node that was saved after this one. */
	public StateStackNode nextNodeUp;
	
	/**
	 * @param kmer K-mer that the aligner should start with for the next base, and the next base
	 *   to be added must be the last base of this k-mer. This constructor only saves a reference
	 *   to the k-mer, so it <strong>MUST</code> be copied (see
	 *   <code>edu.gatech.kanalyze.util.KmerUtil.copy()</code>).
	 * @param nextBase Base that should be added when this state is restored.
	 * @param consensusSize Number of bases in the consensus at this point of the alignment.
	 * @param alignContainer A list of nodes in the alignment column.
	 * @param gapRefContainer A list of nodes in the reference-gap column.
	 * @param gapConContainer A list of nodes in the consensus-gap column.
	 * @param maxAlignmentScore The maximum score of this alignment.
	 * @param maxAlignmentScoreNode Links to the start of maximum-score alignments.
	 * @param minDepth Minimum depth for all k-emrs in the alignment up to this state.
	 * @param nextNodeDown If more than one stack node is present, this links to the next one down
	 *   the stack.
	 */
	public StateStackNode(
			int[] kmer,
			Base nextBase,
			int consensusSize,
			TraceNodeContainer alignContainer,
			TraceNodeContainer gapRefContainer,
			TraceNodeContainer gapConContainer,
			float maxAlignmentScore,
			MaxAlignmentScoreNode maxAlignmentScoreNode,
			int minDepth,
			StateStackNode nextNodeDown) {
		
		// Check arguments
		assert (kmer != null) :
			"StateStackNode(): kmer is null";
		
		assert (nextBase != null) :
			"StateStackNode(): nextbase is null";
		
		assert (consensusSize > 0) :
			"StateStackNode(): consensusSize is less than 1: " + consensusSize;
		
		assert (maxAlignmentScore >= 0.0F) :
			"StateStackNode(): maxAlignmentScore is negative: " + maxAlignmentScore;
		
		assert (minDepth > 0) :
			"StateStackNode(): minDepth is less than 1: " + minDepth;
		
		// Assign fields
		this.kmer = kmer;
		this.nextBase = nextBase;
		this.consensusSize = consensusSize;
		this.alignContainer = alignContainer;
		this.gapRefContainer = gapRefContainer;
		this.gapConContainer = gapConContainer;
		this.maxAlignmentScore = maxAlignmentScore;
		this.maxAlignmentScoreNode = maxAlignmentScoreNode;
		this.minDepth = minDepth;
		this.nextNodeDown = nextNodeDown;
		this.nextNodeUp = null;
		
		return;
	}
	
	/**
	 * Get a restored state object from this object. This allows vital information to
	 * be returned from this node without exposing references to other nodes (not yet restored)
	 * or holding onto references for data that can be garbage-collected.
	 * 
	 * @return A restored state object.
	 */
	public RestoredState getRestoredState() {
		return new RestoredState(kmer, consensusSize, minDepth);
	}
}
