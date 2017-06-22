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

package edu.gatech.kestrel.activeregion;

import java.util.Arrays;

import edu.gatech.kestrel.align.AlignNode;
import edu.gatech.kestrel.align.TraceMatrix;

/**
 * Represents one haplotype over an active region.
 */
public class Haplotype {
	
	/** Active region this haplotype is defined over. */
	public final ActiveRegion activeRegion;
	
	/**
	 * The consensus sequence of this haplotype. This is the consensus sequence from
	 * a k-mer guided progressive alignment over the active region.
	 */
	public final byte[] sequence;
	
	/** Length of the consensus sequence. */
	public final int length;
	
	/** Summary statistics for the bases in this haplotype. */
	public final RegionStats stats;
	
	/** Canonical alignment. */
	public final AlignNode alignment;
	
	/**
	 * Alignment of this haplotype with the active region. With the initial bases at
	 * the first base of the reference and consensus sequences, this series of alignment
	 * events describes how the sequences are aligned (e.g. align bases at the current position
	 * and advance both, insert a gap in the reference sequence and advance only the consensus,
	 * or insert a gap in the consensus and advance only the reference). After this process, all
	 * bases of the reference and consensus sequences must be used, and it cannot extend past
	 * the last base of either sequence.
	 * <p/>
	 * This list of alignments is sorted in its canonical order. For two alignments, the first
	 * disagreement determines their order, which is prioritized by mismatch, reference gap,
	 * consensus gap, and match. For example, if the first disagreement between two alignments
	 * yields a reference gap in one and a match in the other, the one with a reference gap is
	 * sorted first. The alignment used to call variants is the one that is first in this list,
	 * and it is also copied to the <code>alignment</code> field.
	 */
	private final AlignNode[] alignmentList;
	
	/** Number of alignments between the reference and consensus sequences. */
	public final int nAlign;
	
	/** Alignment score for the alignments in this haplotype. */
	public final float alignmentScore;
	
	/**
	 * If the aligner was set to preserve the traceback matrix, then it is stored here. Otherwise,
	 * this field is <code>null</code>.
	 */
	public final TraceMatrix traceMatrix;
	
	/**
	 * Define a new haplotype.
	 * 
	 * @param sequence Sequence of consensus bases.
	 * @param activeRegion Active region this haplotype covers.
	 * @param alignmentList Alignment of the reference and consensus sequences.
	 * @param alignmentScore Alignment score.
	 * @param traceMatrix Alignment traceback matrix or <code>null</code> if the aligner was not set to
	 *   record it.
	 * @param stats Summary statistics for the k-mer counts in <code>sequence</code>.
	 * 
	 * @throws NullPointerException If <code>sequence</code> or <code>activeRegion</code>
	 *   is <code>null</code>.
	 */
	public Haplotype(byte[] sequence, ActiveRegion activeRegion, AlignNode[] alignmentList, float alignmentScore, TraceMatrix traceMatrix, RegionStats stats)
			throws NullPointerException {
		
		// Check arguments
		if (sequence == null)
			throw new NullPointerException("Cannot create haplotype with sequence: null");
		
		if (activeRegion == null)
			throw new NullPointerException("Cannot create haplotype with active region: null");
		
		if (alignmentList == null)
			throw new NullPointerException("Cannot create haplotype with alignmentList: null");
		
		if (alignmentList.length == 0)
			throw new IllegalArgumentException("Cannot create haplotype: alignmentList is an empty array");
		
		if (stats == null)
			throw new NullPointerException("Cannot create haplotype with summary statistics: null");
		
		// Assign fields
		this.sequence = sequence;
		this.activeRegion = activeRegion;
		this.alignmentScore = alignmentScore;
		this.traceMatrix = traceMatrix;
		this.stats = stats;
		
		length = sequence.length;
		
		// Copy, check, and sort alignments
		this.alignmentList = Arrays.copyOf(alignmentList, alignmentList.length);
		
		for (int index = 0; index < this.alignmentList.length; ++index)
			if (this.alignmentList[index] == null)
				throw new IllegalArgumentException("Cannot create haplotype: alignmentList contains a null reference at index " + index);
		
		Arrays.sort(this.alignmentList);
		
		// Save canonical alignment
		this.alignment = alignmentList[0];
		this.nAlign = this.alignmentList.length;
		
		return;
	}
	
	/**
	 * Get the alignment strings of one alignment where <code>alignmentIndex</code> is the
	 * alignment in a sorted list (starting at 0 for the canonical alignment).
	 * 
	 * @param alignmentIndex Index of the alignment in the list of alternate alignments where
	 *   <code>0</code> is the canonical alignment used for variant calls and
	 *   <code>this.nAlign - 1</code> is the last alternate alignment.
	 * 
	 * @return An array of three strings depicting the alignment of the reference and consensus
	 *   sequence with gaps inserted appropriately. The first string is the reference sequence,
	 *   the second is vertical bars (where bases match) or spaces (mismatches or gaps), and
	 *   the third is the consensus sequence.
	 * 
	 * @throws IllegalArgumentException If <code>alignmentIndex</code> is negative or greater than
	 *   the last alternate alignment index.
	 */
	public String[] getAlignmentString(int alignmentIndex)
			throws IllegalArgumentException {
		
		// Declarations
		StringBuilder refBuilder;  // Reference string
		StringBuilder barBuilder;  // Bar (matched bases) or spaces between the reference and consensus strings
		StringBuilder conBuilder;  // Consensus string
		
		AlignNode node;  // Current alignment node
		
		int refIndex;  // Reference sequence index
		int conIndex;  // Consensus sequence index
		
		byte[] reference;  // Reference sequence
		byte[] consensus;  // Consensus sequence
		
		// Check arguments
		if (alignmentIndex >= nAlign || alignmentIndex < 0)
			throw new IllegalArgumentException(String.format("altAlign is less than 0 or greater than the last alternate alignment index (%d): %d", nAlign - 1, alignmentIndex));
		
		// Initialize
		refBuilder = new StringBuilder();
		barBuilder = new StringBuilder();
		conBuilder = new StringBuilder();
		
		node = alignmentList[alignmentIndex];
		
		refIndex = activeRegion.startIndex;
		conIndex = 0;
		
		reference = activeRegion.refRegion.sequence;
		consensus = this.sequence;
		
		// Traverse
		while (node != null) {
			
			if (node.type == AlignNode.MATCH) {
				for (int index = 0; index < node.n; ++index) {
					refBuilder.append((char) reference[refIndex++]);
					conBuilder.append((char) consensus[conIndex++]);
					barBuilder.append('|');
				}
				
			} else if (node.type == AlignNode.MISMATCH) {
				for (int index = 0; index < node.n; ++index) {
					refBuilder.append((char) reference[refIndex++]);
					conBuilder.append((char) consensus[conIndex++]);
					barBuilder.append(' ');
				}
				
			} else if (node.type == AlignNode.INS) {
				for (int index = 0; index < node.n; ++index) {
					refBuilder.append('-');
					conBuilder.append((char) consensus[conIndex++]);
					barBuilder.append(' ');
				}
				
			} else if (node.type == AlignNode.DEL) {
				for (int index = 0; index < node.n; ++index) {
					refBuilder.append((char) reference[refIndex++]);
					conBuilder.append('-');
					barBuilder.append(' ');
				}
				
			} else {
				assert (false) :
					"Unknown node type: " + node.type;
			}
			
			node = node.next;
		}
		
		return new String[] {refBuilder.toString(), barBuilder.toString(), conBuilder.toString()};
	}
	
	/**
	 * Determine if this alignment represents the wildtype sequence.
	 * 
	 * @return <code>true</code> if this alignment represents the wildtype sequence.
	 */
	public boolean isWildtype() {
		return alignmentList[0].type == AlignNode.MATCH && alignmentList[0].next == null;
	}
	
	/**
	 * Get a string version of this haplotype.
	 * 
	 * @return String version of this haplotype.
	 */
	@Override
	public String toString() {
		return String.format("Haplotype[length=%d, min_depth=%d, region=%s]", length, stats.min, activeRegion);
	}
}
