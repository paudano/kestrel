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

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import edu.gatech.kanalyze.util.Base;
import edu.gatech.kanalyze.util.KmerHashSet;
import edu.gatech.kanalyze.util.kmer.KmerUtil;
import edu.gatech.kestrel.KestrelConstants;
import edu.gatech.kestrel.activeregion.ActiveRegion;
import edu.gatech.kestrel.activeregion.Haplotype;
import edu.gatech.kestrel.activeregion.RegionStats;
import edu.gatech.kestrel.align.state.RestoredState;
import edu.gatech.kestrel.align.state.StateStackNode;
import edu.gatech.kestrel.align.state.TraceNodeContainer;
import edu.gatech.kestrel.counter.CountMap;

/**
 * Represents one consensus sequence and the alignment of that sequence with the reference
 * sequence over an active region. This object is built one base at a time by adding
 * consensus bases to the alignment.
 */
public final class KmerAligner {
	
	/** K-mer utility. */
	public final KmerUtil kUtil;
	
	/** K-mer size. */
	private final int kSize;
	
	/** Logger. */
	private final Logger logger;
	
	/** Alignment weights. */
	public final AlignmentWeight alnWeight;
	
	/** If <code>true</code>, build the trace matrix. */
	public final boolean trace;
	
	/** <code>true</code> until <code>init()</code> is called at least once. */
	private boolean notInit;
	
	/** Active region the progressive alignment extends. */
	private ActiveRegion activeRegion;
	
	/** Index of <code>activeRegion.refSequence.sequence</code> where the alignment begins. */
	private int refBaseStart;
	
	/** Index of <code>activeRegion.refSequence.sequence</code> where the alignment ends. */
	private int refBaseEnd;
	
	/** Length of the reference sequence (<code>refBaseEnd - refBaseStart + 1</code>. */
	private int refLength;
	
	/**
	 * If <code>true</code>, this alignment is build from right to left. This is set if the active
	 * region reaches the left end of the sequence and it the alignment must anchor on the right
	 * end of the active region and traverse from right to left.
	 */
	private boolean reverse;
	
	/**
	 * If <code>true</code>, this alignment allows deletions on the ends of the alignment. This is
	 * set if the active region reaches either the left or right end of the reference sequence. 
	 */
	private boolean allowEndDeletion;
	
	/**
	 * The last column of the alignment score matrix for aligned bases. Always length
	 * <code>refLength</code>.
	 */
	private TraceNode[] matrixColAlign;
	
	/**
	 * The last column of the alignment score matrix for gaps in the reference sequence (insertions).
	 * Always length <code>refLength</code>.
	 */
	private TraceNode[] matrixColGapRef;
	
	/**
	 * The last column of the alignment score matrix for gaps in the consensus sequence (deletions).
	 * Always length <code>refLength</code>.
	 */
	private TraceNode[] matrixColGapCon;
	
	/**
	 * As a base is added to the alignment, this is the new <code>matrixColAlign</code> as it
	 * is being built.
	 */
	private TraceNode[] matrixColAlignNext;
	
	/**
	 * As a base is added to the alignment, this is the new <code>matrixColGapRef</code> as it
	 * is being built.
	 */
	private TraceNode[] matrixColGapRefNext;
	
	/**
	 * As a base is added to the alignment, this is the new <code>matrixColGapCon</code> as it
	 * is being built.
	 */
	private TraceNode[] matrixColGapConNext;
	
	/** Size of column matrices. May be larger than the reference length (aligner is reusable). */
	private int referenceCapacity;
	
	/** The consensus sequence. */
	private byte[] consensus;
	
	/** Number of positions in <code>consensus</code> with data. */
	private int consensusSize;
	
	/** <code>consensus.length</code> stored for quicker access. */
	private int consensusCapacity;
	
	/** Maximum alignment score. */
	private float maxAlignmentScore;
	
	/** Maximum alignment score node where traceback begins. */
	private MaxAlignmentScoreNode maxAlignmentScoreNode;
	
	/** The initial consensus sequence array is the active region size multiplied by this factor. */
	private float CONS_SIZE_MULTIPLIER = 1.2F;
	
	/** Trace matrix if <code>trace</code> is <code>true</code>, else, <code>null</code>. */
	private TraceMatrix traceMatrix;
	
	/**
	 * If multiple states are cached (several consensus sequences for one active region),
	 * then they are stored here.
	 */
	private StateStackNode stateStack;
	
	/** The number of states in <code>stateStack</code>. */
	private int nState;
	
	/**
	 * The maximum number of saved alignment states allowed. For some genomic regions where k-mers are
	 * shared with other regions, saving and restoring states can result in a loop that takes a very long
	 * time to resolve. This parameter restricts the number of stored states, and when it is reached, the
	 * least likely stored states are discarded.
	 */
	private int maxState;
	
	/**
	 * A cache of type lists used for building alignments.
	 */
	private TypeList[] typeListCache;
	
	/** Number of elements in <code>typeListCache</code>. */
	private int typeListCacheSize;
	
	/**
	 * Number of elements in <code>typeListCache</code>.
	 */
	private static final int TYPE_LIST_CACHE_CAPACITY = 0;  // DBGTMP: Set to 0 from 10 for troubleshooting
	
	/** Default maximum number of stored states. */
	public static final int DEFAULT_MAX_STATE = 10;
	
	
	/**
	 * Create a new alignment object.
	 * 
	 * @param kUtil K-mer utility.
	 * @param alnWeight Alignment weights.
	 * @param trace If <code>true</code>, record the trace matrix in matrix form. This
	 *   option will have a significant performance impact, and it should only be used
	 *   for troubleshooting and debugging.
	 * 
	 * @throws NullPointerException If <code>kUtil</code> or <code>alnWeight</code> is
	 *   <code>null</code>.
	 * @throws IllegalArgumentException If <code>kUtil.kSize</code> is less than
	 *   <code>KestrelConstants.MIN_KMER_SIZE</code>.
	 */
	public KmerAligner(KmerUtil kUtil, AlignmentWeight alnWeight, boolean trace)
			throws NullPointerException, IllegalArgumentException {
		
		logger = LoggerFactory.getLogger(KmerAligner.class);
		
		// Check arguments
		if (kUtil == null)
			throw new NullPointerException("K-mer utility is null");
		
		if (kUtil.kSize < KestrelConstants.MIN_KMER_SIZE)
			throw new IllegalArgumentException(String.format("K-mer size is less than %d: %d", KestrelConstants.MIN_KMER_SIZE, kUtil.kSize));
		
		if (alnWeight == null)
			throw new NullPointerException("Cannot create alignment with alignment weights: null");
		
		this.kUtil = kUtil;
		this.alnWeight = alnWeight;
		this.trace = trace;
		
		notInit = true;
		
		kSize = kUtil.kSize;
		
		// Create alignment matrix columns
		matrixColAlign = new TraceNode[0];
		matrixColGapRef = new TraceNode[0];
		matrixColGapCon = new TraceNode[0];
		
		matrixColAlignNext = new TraceNode[0];
		matrixColGapRefNext = new TraceNode[0];
		matrixColGapConNext = new TraceNode[0];
		
		referenceCapacity = 0;
		
		// Initialize the type list cache
		typeListCache = new TypeList[TYPE_LIST_CACHE_CAPACITY];
		typeListCacheSize = 0;
		
		// Initialize other fields
		reverse = false;
		allowEndDeletion = false;
		
		// Initialize saved states
		stateStack = null;
		nState = 0;
		
		maxState = DEFAULT_MAX_STATE;
		
		return;
	}
	
	/**
	 * Initialize the alignment.
	 * 
	 * @param activeRegion Active region to align over.
	 * 
	 * @throws NullPointerException If <code>activeRegion</code> is <code>null</code>.
	 */
	public void init(ActiveRegion activeRegion)
			throws NullPointerException {
		
		int newConsensusCapacity;  // New consensus sequence capacity
		
		// Check and set active region
		if (activeRegion == null)
			throw new NullPointerException("Cannot initialize alignment with active region: null");
		
		this.activeRegion = activeRegion;
		
		// Reset fields
		consensusSize = 0;
		stateStack = null;
		
		maxAlignmentScore = 0.0F;
		maxAlignmentScoreNode = null;
		
		// Get sequence start and end positions
		if (activeRegion.leftEnd) {
			reverse = true;
			allowEndDeletion = true;
			
		} else {
			reverse = false;
			allowEndDeletion = activeRegion.rightEnd;
		}
		
		// Assign sequence start and end fields
		refBaseStart = activeRegion.startIndex;
		refBaseEnd = activeRegion.endIndex;
		refLength = refBaseEnd - refBaseStart + 1;
		
		// Expand reference matrix columns
		if (refLength > referenceCapacity) {
			referenceCapacity = refLength + 100;  // Add 100 bases to avoid re-expanding
			
			if (referenceCapacity < 0 || referenceCapacity > KestrelConstants.MAX_ARRAY_SIZE)
				referenceCapacity = KestrelConstants.MAX_ARRAY_SIZE;
			
			matrixColAlign = new TraceNode[referenceCapacity];
			matrixColGapRef = new TraceNode[referenceCapacity];
			matrixColGapCon = new TraceNode[referenceCapacity];
			
			matrixColAlignNext = new TraceNode[referenceCapacity];
			matrixColGapRefNext = new TraceNode[referenceCapacity];
			matrixColGapConNext = new TraceNode[referenceCapacity];
		}
		
		// Create consensus sequence
		newConsensusCapacity = (int) (refLength * CONS_SIZE_MULTIPLIER);
		
		if (newConsensusCapacity < 0 || newConsensusCapacity > KestrelConstants.MAX_ARRAY_SIZE)
			newConsensusCapacity = KestrelConstants.MAX_ARRAY_SIZE;
		
		if (newConsensusCapacity > consensusCapacity) {
			
			consensusCapacity = newConsensusCapacity;
			
			if (consensusCapacity == KestrelConstants.MAX_ARRAY_SIZE)
				logger.trace("Aligner creating consensus sequence at the maximum array size: " + consensusCapacity);
			
			consensus = new byte[consensusCapacity];
			consensusCapacity = newConsensusCapacity;
		}
		
		consensusSize = 0;
		
		// Initialize max alignment
		maxAlignmentScore = 0.0F;
		maxAlignmentScoreNode = null;
		
		// Set other fields
		stateStack = null;
		nState = 0;
		
		// Initialize alignment
		initAlignment();
		
		notInit = false;
		
		return;
	}
	
	/**
	 * Initialize this alignment and prepare it for adding bases.
	 */
	private void initAlignment() {
		
		TraceNode lastNode;  // Last trace-node object
		float initScore;     // Initial alignment score
		
		byte[] sequence;   // Reference sequence
		byte[] consensus;  // Consensus sequence
		
		// Init
		lastNode = null;
		initScore = alnWeight.getInitialScore(kSize);
		
		sequence = activeRegion.refRegion.sequence;
		consensus = this.consensus;
		
		// Assign default values to matrix columns
		for (int index = 0; index < refLength; ++index) {
			matrixColAlign[index] = TraceNode.ZERO_NODE;
			matrixColGapRef[index] = TraceNode.ZERO_NODE;
			matrixColGapCon[index] = TraceNode.ZERO_NODE;
		}
		
		matrixColAlignNext[0] = TraceNode.ZERO_NODE;
		matrixColGapRefNext[0] = TraceNode.ZERO_NODE;
		matrixColGapConNext[0] = TraceNode.ZERO_NODE;
		
		// Begin match by aligning the k-mer that seeded the alignment
		if (trace) {
			traceMatrix = new TraceMatrix(refLength);
			
			for (int count = 0; count < kSize; ++count) {
				lastNode = new TraceNode(initScore, TraceNode.TYPE_MATCH, lastNode);
				traceMatrix.nextCol();
				traceMatrix.set(count, TraceNode.TYPE_MATCH, TraceNode.TYPE_MATCH);
			}
			
		} else {
			for (int count = 0; count < kSize; ++count)
				lastNode = new TraceNode(initScore, TraceNode.TYPE_MATCH, lastNode);
		}
		
		// Assign match matrix column
		matrixColAlign[kSize - 1] = lastNode;
		
		// Allow a gap to open in the haplotype at this position
		if (initScore + alnWeight.newGap > 0) {
			matrixColGapCon[kSize] = new TraceNode(initScore + alnWeight.newGap, TraceNode.TYPE_GAP_CON, lastNode);
			
			lastNode = matrixColGapCon[kSize];
			
			for (int index = kSize + 1; index < refLength && lastNode.score + alnWeight.gapExtend > 0; ++index) {
				matrixColGapCon[index] = new TraceNode(lastNode.score + alnWeight.gapExtend, TraceNode.TYPE_GAP_CON, lastNode);
				lastNode = matrixColGapCon[index];
			}
		}
		
		// Store bases in the consensus sequence
		if (reverse) {
			for (int index = 0; index < kSize; ++index)
				consensus[index] = sequence[refBaseEnd - index];
			
		} else {
			for (int index = 0; index < kSize; ++index)
				consensus[index] = sequence[refBaseStart + index];
		}
		
		consensusSize = kSize;
		
		// Clear the type list
		clearTypeList();
		
		return;
	}
	
	/**
	 * Add a base to this alignment.
	 * 
	 * @param base Base to add.
	 * 
	 * @return <code>true</code> if the alignment score might be improved by adding
	 *   another base, or <code>false</code> if the maximum alignment score has been
	 *   reached.
	 * 
	 * @throws NullPointerException If <code>base</code> is <code>null</code>.
	 * @throws IllegalStateException If the consensus sequence has reached its
	 *   maximum length, <code>KestrelConstants.MAX_ARRAY_SIZE</code>, or if
	 *   <code>init()</code> was not called at least once.
	 * 
	 * @see KestrelConstants#MAX_ARRAY_SIZE
	 */
	public final boolean addBase(Base base)
			throws NullPointerException, IllegalStateException {
		
		float alignScore;   // Score if bases are aligned
		float refGapScore;  // Score if a gap is inserted in the reference sequence
		float conGapScore;  // Score if a gap is inserted in the consensus sequence
		float maxScore;     // Maximum of the align and gap scores
		
		int refLength;     // Length of the reference sequence
		int refIndex;      // Index of reference sequence at the current base
		
		float addAlignScore;  // alnWeight.match if bases match, and alnWeight.mismatch if they do not
		byte alignType;       // ALN_MATCH if bases match, or ALN_MISMATCH if they do not
		
		TraceNode lastNode;
		
		float maxPotScore;     // Maximum potential score
		float newMaxPotScore;  // New maximum potential score (calculated and compared to maxPotScore
		
		TraceNode[] tnSwap;  // Matrix column swap
		
		// Check arguments and state
		if (notInit)
			throw new IllegalStateException("Aligner was not initialized (must call init())");
		
		if (base == null)
			throw new NullPointerException("Cannot add base to alignment: null");
		
		// Check capacity and expand
		if (consensusSize == consensusCapacity)
			expandConsensus();  // throws IllegalStateException
		
		// Init
		refLength = this.refLength;
		
		maxPotScore = 0.0F;
		
		// Add base to the consensus sequence
		consensus[consensusSize++] = base.baseCharByte;
		
		//
		// Align table
		//
		for (int index = 1; index < refLength; ++index) { // Start at 1, can never align from the top of the matrix
			
			// Get index of the reference base
			if (reverse)
				refIndex = refBaseEnd - index;
			else
				refIndex = index + refBaseStart;
			
			// Set score (match or mismatch)
			if (activeRegion.refRegion.sequence[refIndex] == base.baseCharByte) {
				addAlignScore = alnWeight.match;
				alignType = TraceNode.TYPE_MATCH;
				
			} else {
				addAlignScore = alnWeight.mismatch;
				alignType = TraceNode.TYPE_MISMATCH;
			}
			
			alignScore = 0.0F;
			refGapScore = 0.0F;
			conGapScore = 0.0F;
			maxScore = 0.0F;
			
			lastNode = null;
			
			// Align from align table
			if (matrixColAlign[index - 1] != TraceNode.ZERO_NODE) {
				
				alignScore = matrixColAlign[index - 1].score + addAlignScore;
				
				if (alignScore > maxScore)
					maxScore = alignScore;
			}
			
			// Align from ref-gap table
			if (matrixColGapRef[index - 1] != TraceNode.ZERO_NODE) {
				
				refGapScore = matrixColGapRef[index  - 1].score + addAlignScore;
				
				if (refGapScore > maxScore)
					maxScore = refGapScore;
			}
			
			// Align from con-gap table
			if (matrixColGapCon[index - 1] != TraceNode.ZERO_NODE) {
				
				conGapScore = matrixColGapCon[index - 1].score + addAlignScore;
				
				if (conGapScore > maxScore)
					maxScore = conGapScore;
			}
			
			// Add node(s)
			if (maxScore > 0.0F) {
				
				if (alignScore == maxScore) {
					lastNode = new TraceNode(maxScore, alignType, matrixColAlign[index - 1], lastNode);
					matrixColAlignNext[index] = lastNode;
				}
				
				if (refGapScore == maxScore) {
					lastNode = new TraceNode(maxScore, alignType, matrixColGapRef[index - 1], lastNode);
					matrixColAlignNext[index] = lastNode;
				}
				
				if (conGapScore == maxScore) {
					lastNode = new TraceNode(maxScore, alignType, matrixColGapCon[index - 1], lastNode);
					matrixColAlignNext[index] = lastNode;
				}
				
				// Calculate maximum potential score from this location
				newMaxPotScore = maxScore + (refLength - index - 1) * alnWeight.match;
				
				if (newMaxPotScore > maxPotScore)
					maxPotScore = newMaxPotScore;
				
			} else {
				matrixColAlignNext[index] = TraceNode.ZERO_NODE;
			}
		}
		
		// If score of the bottom element is the maximum, record this as a start position
		maxScore = matrixColAlignNext[refLength - 1].score;
		
		if (maxScore >= maxAlignmentScore && maxScore > 0) {
			
			// Replace if maxScore > maxAlignmentScore, and append if maxScore == maxAlignmentScore
			maxAlignmentScoreNode = new MaxAlignmentScoreNode(
					matrixColAlignNext[refLength - 1], consensusSize,
					maxScore > maxAlignmentScore ? null : maxAlignmentScoreNode  // Replace list if greater, prepend if equal
			);
			
			maxAlignmentScore = maxScore;
		}
		
		
		//
		// Insertion table (reference gap)
		//
		
		for (int index = 0; index < refLength; ++index) {
			
			alignScore = 0.0F;
			refGapScore = 0.0F;
			conGapScore = 0.0F;
			maxScore = 0.0F;
			
			lastNode = null;
			
			// Gap from align table
			if (matrixColAlign[index] != TraceNode.ZERO_NODE) {
				alignScore = matrixColAlign[index].score + alnWeight.newGap;
				
				if (alignScore > maxScore)
					maxScore = alignScore;
			}
			
			// Gap from ref-gap table
			if (matrixColGapRef[index] != TraceNode.ZERO_NODE) {
				refGapScore = matrixColGapRef[index].score + alnWeight.gapExtend;
				
				if (refGapScore > maxScore)
					maxScore = refGapScore;
			}
			
			// Gap from con-gap table
			if (matrixColGapCon[index] != TraceNode.ZERO_NODE) {
				conGapScore = matrixColGapCon[index].score + alnWeight.newGap;
				
				if (conGapScore > maxScore)
					maxScore = conGapScore;
			}
			
			// Add node(s)
			if (maxScore > 0.0F) {
				
				if (alignScore == maxScore) {
					lastNode = new TraceNode(maxScore, TraceNode.TYPE_GAP_REF, matrixColAlign[index], lastNode);
					matrixColGapRefNext[index] = lastNode;
				}
				
				if (refGapScore == maxScore) {
					lastNode = new TraceNode(maxScore, TraceNode.TYPE_GAP_REF, matrixColGapRef[index], lastNode);
					matrixColGapRefNext[index] = lastNode;
				}
				
				if (conGapScore == maxScore) {
					lastNode = new TraceNode(maxScore, TraceNode.TYPE_GAP_REF, matrixColGapCon[index], lastNode);
					matrixColGapRefNext[index] = lastNode;
				}
				
			} else {
				matrixColGapRefNext[index] = TraceNode.ZERO_NODE;
			}
		}
		
		//
		// Deletion table (haplotype gap)
		//
		
		for (int index = 1; index < refLength; ++index) { // Start at 1, can never align from the top of the matrix
			
			alignScore = 0.0F;
			refGapScore = 0.0F;
			conGapScore = 0.0F;
			maxScore = 0.0F;
			
			lastNode = null;
			
			// Align from align table
			if (matrixColAlignNext[index - 1] != TraceNode.ZERO_NODE) {
				
				alignScore = matrixColAlignNext[index - 1].score + alnWeight.newGap;
				
				if (alignScore > maxScore)
					maxScore = alignScore;
			}
			
			// Align from ref-gap table
			if (matrixColGapRefNext[index - 1] != TraceNode.ZERO_NODE) {
				
				refGapScore = matrixColGapRefNext[index  - 1].score + alnWeight.newGap;
				
				if (refGapScore > maxScore)
					maxScore = refGapScore;
			}
			
			// Align from con-gap table
			if (matrixColGapConNext[index - 1] != TraceNode.ZERO_NODE) {
				
				conGapScore = matrixColGapConNext[index - 1].score + alnWeight.gapExtend;
				
				if (conGapScore > maxScore)
					maxScore = conGapScore;
			}
			
			// Add node(s)
			if (maxScore > 0.0F) {
				
				if (alignScore == maxScore) {
					lastNode = new TraceNode(maxScore, TraceNode.TYPE_GAP_CON, matrixColAlignNext[index - 1], lastNode);
					matrixColGapConNext[index] = lastNode;
				}
				
				if (refGapScore == maxScore) {
					lastNode = new TraceNode(maxScore, TraceNode.TYPE_GAP_CON, matrixColGapRefNext[index - 1], lastNode);
					matrixColGapConNext[index] = lastNode;
				}
				
				if (conGapScore == maxScore) {
					lastNode = new TraceNode(maxScore, TraceNode.TYPE_GAP_CON, matrixColGapConNext[index - 1], lastNode);
					matrixColGapConNext[index] = lastNode;
				}
				
				// Calculate maximum potential score from this location
				if (allowEndDeletion) {
					newMaxPotScore = maxScore + (refLength - index - 1) * alnWeight.match;
					
					if (newMaxPotScore > maxPotScore)
						maxPotScore = newMaxPotScore;
				}
				
			} else {
				matrixColGapConNext[index] = TraceNode.ZERO_NODE;
			}
		}
		
		if (allowEndDeletion) {
			
			// If score of the bottom element is the maximum, record this as a start position
			maxScore = matrixColGapConNext[refLength - 1].score;
			
			if (maxScore >= maxAlignmentScore && maxScore > 0) {
				
				// Replace if maxScore > maxAlignmentScore, and append if maxScore == maxAlignmentScore
				maxAlignmentScoreNode = new MaxAlignmentScoreNode(
						matrixColGapConNext[refLength - 1], consensusSize,
						maxScore > maxAlignmentScore ? null : maxAlignmentScoreNode  // Replace list if greater, append if equal
				);
				
				maxAlignmentScore = maxScore;
			}
		}
		
		// Set trace matrix if enabled
		if (trace) {
			
			TraceNode node;
			
			traceMatrix.nextCol();
			
			for (int index = 0; index < refLength; ++index) {
				
				// Align table
				node = matrixColAlignNext[index];
				
				if (node != TraceNode.ZERO_NODE) {
					
					while (node != null) {
						
						traceMatrix.set(index, node.nextNode.type, node.type);
						
						node = node.branchNode;
					}
				}
				
				// Insertion (reference gap)
				node = matrixColGapRefNext[index];
				
				if (node != TraceNode.ZERO_NODE) {
					
					while (node != null) {
						
						traceMatrix.set(index, node.nextNode.type, node.type);
						
						node = node.branchNode;
					}
				}
				
				// Deletion (consensus gap)
				node = matrixColGapConNext[index];
				
				if (node != TraceNode.ZERO_NODE) {
					
					while (node != null) {
						
						traceMatrix.set(index, node.nextNode.type, node.type);
						
						node = node.branchNode;
					}
				}
			}
		}
		
		// Advance matrix columns
		tnSwap = matrixColAlign;
		matrixColAlign = matrixColAlignNext;
		matrixColAlignNext = tnSwap;
		
		tnSwap = matrixColGapRef;
		matrixColGapRef = matrixColGapRefNext;
		matrixColGapRefNext = tnSwap;
		
		tnSwap = matrixColGapCon;
		matrixColGapCon = matrixColGapConNext;
		matrixColGapConNext = tnSwap;
		
		return maxPotScore >= maxAlignmentScore && maxPotScore > 0.0F;
	}
			
	/**
	 * Save an alignment state.
	 * 
	 * @param kmer The next k-mer in this alignment. The last base of the k-mer is the first
	 *   to be added when the alignment restores this state.
	 * @param nextBase Base that should be added when this state is restored.
	 * @param minDepth Minimum depth for all k-mers in the alignment up to this point.
	 * 
	 * @throws NullPointerException If <code>kmer</code> or <code>nextBase</code> is
	 *   <code>null</code>.
	 * @throws IllegalArgumentException If <code>minDepth</code> is less than <code>1</code>.
	 */
	public final void saveState(int[] kmer, Base nextBase, int minDepth, KmerHashSet kmerHash, int repeatCount)
			throws NullPointerException, IllegalArgumentException {
		
		TraceNodeContainer alignContainer;
		TraceNodeContainer gapRefContainer;
		TraceNodeContainer gapConContainer;
		
		int index;
		
		// Check arguments
		if (kmer == null)
			throw new NullPointerException("Cannot save state for k-mer: null");
		
		if (nextBase == null)
			throw new NullPointerException("Cannot save state with next base: null");
		
		if (minDepth < 1)
			throw new IllegalArgumentException("Minimum depth may not be less than 1: " + minDepth);
		
		if (kmerHash == null)
			throw new NullPointerException("Connot save state with k-mer hash: null");
		
		if (repeatCount < 0)
			throw new IllegalArgumentException("Repeat count must be non-negative: " + repeatCount);
		
		// Check max depth
		if (nState == maxState) {
			
			if (! removeLastMinState(minDepth)) {
				logger.trace("Rejecting state save: State stack is at capacity: {} [minDepth={}]", kUtil.toBaseString(kmer), minDepth);
				return;
			}
		}
		
		// Store alignment matrix column
		alignContainer = null;
		index = 0;
		
		while (index < refLength) {
			
			if (matrixColAlign[index] != TraceNode.ZERO_NODE)
				alignContainer = new TraceNodeContainer(index, matrixColAlign[index], alignContainer);
			
			++index;
		}
		
		// Store reference-gap matrix column
		gapRefContainer = null;
		index = 0;
		
		while (index < refLength) {
			if (matrixColGapRef[index] != TraceNode.ZERO_NODE)
				gapRefContainer = new TraceNodeContainer(index, matrixColGapRef[index], gapRefContainer);
			
			++index;
		}
		
		// Store consensus-gap matrix column
		gapConContainer = null;
		index = 0;
		
		while (index < refLength) {
			if (matrixColGapCon[index] != TraceNode.ZERO_NODE)
				gapConContainer = new TraceNodeContainer(index, matrixColGapCon[index], gapConContainer);
			
			++index;
		}
		
		// Save state
		stateStack = new StateStackNode(
				kmer,
				nextBase,
				consensusSize,
				alignContainer,
				gapRefContainer,
				gapConContainer,
				maxAlignmentScore,
				maxAlignmentScoreNode,
				minDepth,
				stateStack,
				new KmerHashSet(kmerHash),
				repeatCount
		);
		
		if (stateStack.nextNodeDown != null)  // Set upward link on previous state
			stateStack.nextNodeDown.nextNodeUp = stateStack;
		
		++nState;
		
		return;
	}
	
	/**
	 * Restore the last state cached.
	 * 
	 * @return The state being resumed or <code>null</code> if there are no
	 *   states to restore.
	 */
	public final RestoredState restoreState() {
		
		TraceNodeContainer container;  
		StateStackNode thisState;      // State to be restored
		
		int index;  // Index for traversing tables
		
		// Check if there is a state to restore
		if (stateStack == null)
			return null;
		
		thisState = stateStack;
		
		// Restore fields
		consensusSize = thisState.consensusSize;
		maxAlignmentScore = thisState.maxAlignmentScore;
		maxAlignmentScoreNode = thisState.maxAlignmentScoreNode;
		
		// Restore alignment table
		container = stateStack.alignContainer;
		index = refLength - 1;
		
		while (container != null) {
			
			while (index > container.index)
				matrixColAlign[index--] = TraceNode.ZERO_NODE;
			
			matrixColAlign[index--] = container.node;
			container = container.next;
		}
		
		while (index >= 0)
			matrixColAlign[index--] = TraceNode.ZERO_NODE;
		
		// Restore reference gap table
		container = thisState.gapRefContainer;
		index = refLength - 1;
		
		while (container != null) {
			
			while (index > container.index)
				matrixColGapRef[index--] = TraceNode.ZERO_NODE;
			
			matrixColGapRef[index--] = container.node;
			container = container.next;
		}
		
		while (index >= 0)
			matrixColGapRef[index--] = TraceNode.ZERO_NODE;
		
		// Restore consensus gap table
		container = thisState.gapConContainer;
		index = refLength - 1;
		
		while (container != null) {
			
			while (index > container.index)
				matrixColGapCon[index--] = TraceNode.ZERO_NODE;
			
			matrixColGapCon[index--] = container.node;
			container = container.next;
		}
		
		while (index >= 0)
			matrixColGapCon[index--] = TraceNode.ZERO_NODE;
		
		// Add base
		addBase(thisState.nextBase);
		
		// Remove from state stack
		stateStack = stateStack.nextNodeDown;
		
		// Return k-mer
		return thisState.getRestoredState();
	}
	
	/**
	 * Determine if this alignment has saved states to be restored.
	 * 
	 * @return <code>true</code> if this alignment has saved states.
	 */
	public final boolean hasCachedStates() {
		return stateStack != null;
	}
	
	/**
	 * Get an array of haplotypes from this aligner. If multiple maximum-score consensus sequences
	 * were found, then this list may contain more than one element.
	 * 
	 * @param counter K-mer counter for generating summary statistics over haplotypes.
	 * @param countReverseKmers <code>true</code> if k-mers and their reverse complement are
	 *   counted.
	 * 
	 * @return An array of haplotypes. If no alignments were found, this array will be empty.
	 * 
	 * @throws NullPointerException If <code>counter</code> is <code>null</code>.
	 */
	public Haplotype[] getHaplotypes(CountMap counter, boolean countReverseKmers) {
		
		// Declarations
		MaxAlignmentScoreNode scoreNode;  // Current node
		
		Haplotype[] haplotypeList;  // List of haplotypes
		int index;                  // Index of haplotypeList
		
		byte[] consensus = this.consensus;  // Saved for performance (many reads to copy)
		
		byte[] haploConsensus;   // Consensus sequence for a haplotype
		int haploConsensusSize;  // Size of haploConsenus
		
		RegionStats stats;  // Summary statistics for the haplotype k-mer counts
		
		// Check arguments
		if (counter == null)
			throw new NullPointerException("K-mer counter is null");
		
		// Trim haplotypes that do not match this active region
		trimHaplotypes();
		
		// Get list length
		scoreNode = maxAlignmentScoreNode;
		index = 0;
		
		while (scoreNode != null) {
			
			if (! scoreNode.haplotypeBuilt)
				++index;
			
			scoreNode = scoreNode.next;
		}
		
		// Create haplotype list
		haplotypeList = new Haplotype[index];
		scoreNode = maxAlignmentScoreNode;
		index = 0;
		
		// Fill haplotypes
		while (scoreNode != null) {
			
			// Skip haplotypes that have already been constructed
			if (scoreNode.haplotypeBuilt) {
				scoreNode = scoreNode.next;
				
				continue;
			}
			
			// Create consensus for haplotype
			haploConsensusSize = scoreNode.nConsensusBases;
			haploConsensus = new byte[haploConsensusSize];
			
			if (reverse) {
				for (int copyIndex = 0; copyIndex < haploConsensusSize; ++copyIndex)
					haploConsensus[copyIndex] = consensus[haploConsensusSize - copyIndex - 1];
				
			} else {
				for (int copyIndex = 0; copyIndex < haploConsensusSize; ++copyIndex)
					haploConsensus[copyIndex] = consensus[copyIndex];
			}
			
			// Create haplotype
			stats = RegionStats.getStats(haploConsensus, 0, haploConsensusSize, counter, countReverseKmers);
			haplotypeList[index] = new Haplotype(haploConsensus, activeRegion, getAlignment(scoreNode.traceNode), scoreNode.traceNode.score, traceMatrix, stats);
			
			// Flag
			scoreNode.haplotypeBuilt = true;
			
			// Next max score node
			scoreNode = scoreNode.next;
			++index;
		}
		
		// Return haplotypes
		return haplotypeList;
	}
	
	/**
	 * Get reverse state. If <code>reverse</code>, then the alignment is build in reverse. This happens
	 * when the active region reaches the left end of the reference and it must build from a k-mer in
	 * the reference and extend toward the left end.
	 * 
	 * @return <code>true</code> if this alignment is built in reverse.
	 */
	public boolean isReverse() {
		return reverse;
	}
	
	/**
	 * If <code>true</code>, then this alignment allows a deletion event on the ends of the alignment. This
	 * happens when the active region reaches either the left or right end of the reference and there is no
	 * k-mer to anchor the k-mer alignment to.
	 * 
	 * @return <code>true</code> if the alignment may end in a deletion.
	 */
	public boolean isAllowEndDeletion() {
		return allowEndDeletion;
	}
	
	/**
	 * Expand the consensus sequence array.
	 * 
	 * @throws IllegalStateException If the capacity is already at its maximum value,
	 *   <code>KestrelConstants.MAX_ARRAY_SIZE</code>.
	 * 
	 * @see KestrelConstants#MAX_ARRAY_SIZE
	 */
	private final void expandConsensus()
			throws IllegalStateException {
		
		int newCapacity;
		byte[] newCons;
		
		int conSize = consensusSize;
		
		// Check arguments
		if (consensusCapacity == KestrelConstants.MAX_ARRAY_SIZE)
			throw new IllegalStateException("Cannod expand consensus sequence: Already at maximum size: " + KestrelConstants.MAX_ARRAY_SIZE);
		
		// Get new capacity
		newCapacity = (int) (consensusCapacity * KestrelConstants.ARRAY_EXPAND_FACTOR);
		
		if (newCapacity < 0)
			newCapacity = KestrelConstants.MAX_ARRAY_SIZE;
		
		if (newCapacity == consensusCapacity)
			++newCapacity;
		
		// Create consensus array and copy
		newCons = new byte[newCapacity];
		
		for (int index = 0; index < conSize; ++index)
			newCons[index] = consensus[index];
		
		consensus = newCons;
		consensusCapacity = newCapacity;
		
		return;
	}
	
	/**
	 * Set the maximum number of states to be saved. The least-likely states are trimmed when
	 * this threshold is reached.
	 * 
	 * @param maxState Maximum number of states.
	 * 
	 * @throws IllegalArgumentException If <code>maxState</code> is less than <code>1</code>.
	 * 
	 * @see #DEFAULT_MAX_STATE
	 */
	public void setMaxState(int maxState)
		throws IllegalArgumentException {
		
		if (maxState < 1)
			throw new IllegalArgumentException("Maximum number of states must not be less than 1: " + maxState);
		
		this.maxState = maxState;
		
		return;
	}
	
	/**
	 * Get the maximum number of states that may be saved.
	 * 
	 * @return The maximum number of states that may be saved.
	 * 
	 * @see #setMaxState(int)
	 */
	public int getMaxState() {
		return maxState;
	}
	
	/**
	 * Remove haplotypes that do not end in a deletion (if <code>allowEndDeletion</code>
	 * is <code>true</code> and that do not end in a k-mer that matches the corresponding
	 * end of the active region.
	 */
	private void trimHaplotypes() {
		
		// Declarations
		MaxAlignmentScoreNode lastNode;   // Nodes are added here after they are checked
		MaxAlignmentScoreNode nextNode;   // Node currently being analyzed
		
		byte[] refSeq;  // Reference sequence
		byte[] conSeq;  // Consensus sequence
		
		int refIndex;     // Location in the reference sequence
		int conIndex;     // Location in the consensus sequence
		
		int kSize;  // K-mer size
		
		boolean remove;  // Set to true when a node is to be removed
		
		// Do not trim regions that reach ends (may not end in a full k-mer)
		if (activeRegion.leftEnd || activeRegion.rightEnd)
			return;
		
		// Init
		kSize = this.kSize;
		
		refSeq = activeRegion.refRegion.sequence;
		conSeq = consensus;
		
		// Check nodes
		nextNode = maxAlignmentScoreNode;
		lastNode = null;
		
		while (nextNode != null) {
			remove = false;

			// Start at the index of the last k-mer int the consensus sequence
			conIndex = nextNode.nConsensusBases - kSize;
			
			// Search for k matches
			if (reverse) {
				refIndex = activeRegion.startIndex + kSize - 1;
				
				while (conIndex < nextNode.nConsensusBases) {
					if (refSeq[refIndex--] != conSeq[conIndex++]) {
						remove = true;
						break;
					}
				}
				
			} else {
				refIndex = activeRegion.endIndex - kSize + 1;
				
				while (conIndex < nextNode.nConsensusBases) {
					if (refSeq[refIndex++] != conSeq[conIndex++]) {
						remove = true;
						break;
					}
				}
			}
			
			if (remove) {
				logger.trace("Trimming alignment: Mismatch end k-mer for active region: {}: {}", activeRegion.toString(), nextNode.toString());
				
				// Unlink
				if (lastNode == null) {
					maxAlignmentScoreNode = maxAlignmentScoreNode.next;
					
				} else {
					lastNode.next = nextNode.next;
				}
				
			} else {
				// Record last accepted node
				lastNode = nextNode;
			}
			
			nextNode = nextNode.next;
		}
		
		return;
	}
	
	/**
	 * Get alignments from this aligner in random order (<code>Haplotype</code> will sort them).
	 * 
	 * @param traceNode Maximum-score trace node.
	 * 
	 * @return An array of alignments characterizing the relationship between the reference and
	 *   consensus sequences.
	 */
	private AlignNode[] getAlignment(TraceNode traceNode) {
		
		int nAlign;  // Number of alignments
		
		AlignStart alignmentStart;  // A linked list to the start of alignments
		AlignStart thisStart;       // Current element of alignment start being extended
		AlignStart lastStart;       // Last element where duplicated starts are added
		
		TypeList typeList;  // TypeList object in the current alignment
		
		// Init
		nAlign = 1;
		alignmentStart = new AlignStart(traceNode, getTypeList(null), null);
		thisStart = alignmentStart;
		lastStart = alignmentStart;
		
		// Trace
		while (thisStart != null) {
			
			// Get trace node
			traceNode = thisStart.traceNode;
			typeList = thisStart.typeList;
			
			while (traceNode != null) {
				
				// Handle branches
				if (traceNode.branchNode != null) {
					lastStart.next = new AlignStart(traceNode.branchNode, getTypeList(typeList), null);
					lastStart = lastStart.next;
					++nAlign;
				}
				
				typeList.add(traceNode.type);
				traceNode = traceNode.nextNode;
			}
			
			// Move to next start and extend
			thisStart = thisStart.next;
		}
		
		// Create a list of alignments
		AlignNode[] alignList = new AlignNode[nAlign];
		
		for (int index = 0; index < nAlign; ++index) {
			alignList[index] = alignmentStart.typeList.toAlignment(! reverse);
			returnTypeList(alignmentStart.typeList);
			alignmentStart = alignmentStart.next;
		}
		
		// Return alignments
		return alignList;
	}
	
	/**
	 * Get a type list from the cache of unused list objects or create a new list
	 * if the cache is empty.
	 * 
	 * @param oList Copy contents of this list into the type list returned by this
	 *   method. If <code>null</code>, a default empty list is returned.
	 * 
	 * @return A type list object.
	 */
	private TypeList getTypeList(TypeList oList) {
		
		TypeList tList;
		
		// Get existing list and initialize
		if (typeListCacheSize > 0) {
			tList = typeListCache[--typeListCacheSize];
			tList.cloneFrom(oList);
			
			typeListCache[typeListCacheSize] = null;  // Free for GC
			
			return tList;
		}
		
		// Return a new list copying oList
		return new TypeList(oList);
	}
	
	/**
	 * Return an unused type list to the cache. If the cache is full, the type
	 * list object is discarded.
	 * 
	 * @param list List to return. If <code>null<code>, nothing is done.
	 */
	private void returnTypeList(TypeList list) {
		
		// Do not add if null or the cache is full
		if (list == null || typeListCacheSize == TYPE_LIST_CACHE_CAPACITY)
			return;
		
		// Cache list
		typeListCache[typeListCacheSize++] = list;
		
		return;
	}
	
	/**
	 * Clear the type-list cache used for processing alignments.
	 */
	public void clearTypeList() {
		
		while (typeListCacheSize > 0)
			typeListCache[--typeListCacheSize] = null;
		
		return;
	}
	
	/**
	 * Remove the last minimum state less than <code>minDepthLimit</code>.
	 * 
	 * @param minDepthLimit The removed state must have a limit less than this value.
	 * 
	 * @return <code>true</code> if a state with a minimum depth less than
	 *   <code>minDepthLimit</code> was found and removed and <code>false</code> if
	 *   a state with a minimum depth of less than this value was found.
	 */
	private boolean removeLastMinState(int minDepthLimit) {
		
		// Declarations
		StateStackNode stateNode = stateStack;
		StateStackNode minState;
		
		// Check arguments
		assert (nState == maxState) :
			String.format("removeLastMinState(): Called when nState (%d) != maxState (%d)", nState, maxState);
		
		assert (minDepthLimit > 0) :
			"removeLastMinState(): Called with minDepthLimit < 0: " + minDepthLimit;
		
		// Search for a minimum state
		minState = null;
		
		while (stateNode != null) {
			
			if (stateNode.minDepth < minDepthLimit)
				if (minState == null || stateNode.minDepth < minState.minDepth)
					minState = stateNode;
			
			stateNode = stateNode.nextNodeDown;
		}
		
		// Return if a minimum state was not found
		if (minState == null)
			return false;
		
		logger.trace("Removing saved state: State stack is at capacity: {} [minDepth={}]", kUtil.toBaseString(minState.kmer), minState.minDepth);
		
		// Unlink state
		if (minState.nextNodeUp == null) {
			stateStack = minState.nextNodeDown;
			
		} else {
			minState.nextNodeUp.nextNodeDown = minState.nextNodeDown;
		}
		
		if (minState.nextNodeDown != null)
			minState.nextNodeDown.nextNodeUp = minState.nextNodeUp;
		
		// State was removed
		--nState;
		
		return true;
	}
	
	/**
	 * Tracks alignment start positions for building a list of possible alignments from the
	 * trace matrix.
	 */
	private class AlignStart {
		
		/** List of types. */
		public TypeList typeList;
		
		/** Node at this start. */
		public TraceNode traceNode;
		
		/**
		 * The start of the next alignment or <code>null</code> if this is the last
		 * alignment start position.
		 */
		public AlignStart next;
		
		/**
		 * Create a new alignment start.
		 * 
		 * @param traceNode Trace node at this alignment start.
		 * @param typeList A list of type events in this alignment.
		 * @param next Next alignment start.
		 */
		public AlignStart(TraceNode traceNode, TypeList typeList, AlignStart next) {
			
			assert (traceNode != null) :
				"Trace node is null";
			
			this.traceNode = traceNode;
			this.next = next;
			
			if (typeList == null)
				this.typeList = new TypeList();
			else
				this.typeList = typeList;
			
			return;
		}
	}
}
