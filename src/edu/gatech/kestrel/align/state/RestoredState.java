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

import edu.gatech.kanalyze.util.KmerHashSet;

/**
 * Contains references to select fields of the state stack. This can be sent outside of the
 * aligner after a state is restored without exposing references to other states.
 */
public class RestoredState {
	
	/**
	 * K-mer that the aligner should start with for the next base, and the next base
	 * to be added must be the last base of this k-mer.
	 */
	public final int[] kmer;
	
	/** Number of bases in the consensus at this point of the alignment. */
	public final int consensusSize;
	
	/** Minimum depth for all k-mers in the alignment up to this state. */
	public final int minDepth;
	
	/** K-mer hash for cycle detection. */
	public final KmerHashSet kmerHash;
	
	/** Maximum number of repeated k-mers. */
	public final int repeatCount;
	
	/**
	 * Create a restored state object.
	 * 
	 * @param kmer Kmer the aligner should start with for the next base.
	 * @param consensusSize Size of the consensus sequence.
	 * @param minDepth Minimum depth.
	 */
	public RestoredState(int[] kmer, int consensusSize, int minDepth, KmerHashSet kmerHash, int repeatCount) {
		
		this.kmer = kmer;
		this.consensusSize = consensusSize;
		this.minDepth = minDepth;
		this.kmerHash = kmerHash;
		this.repeatCount = repeatCount;
		
		return;
	}
}
