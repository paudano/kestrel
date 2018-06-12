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
import edu.gatech.kestrel.activeregion.ActiveRegion;
import edu.gatech.kestrel.activeregion.Haplotype;
import edu.gatech.kestrel.align.state.RestoredState;
import edu.gatech.kestrel.counter.CountMap;
import edu.gatech.kanalyze.util.kmer.KmerUtil;
import edu.gatech.kanalyze.util.KmerHashSet;

/**
 * Builds consensus sequences for active regions using k-mer evidence to choose bases and
 * an aligner to guide the building process.
 */
public class KmerAlignmentBuilder {
	
	/** K-mer size. */
	public final KmerUtil kUtil;
	
	/** Logger. */
	private final Logger logger;
	
	/** K-mer aligner. */
	private final KmerAligner aligner;
	
	/** Alignment weights. */
	public final AlignmentWeight alnWeight;
	
	/** Count reverse complement k-mers in region statistics if <code>true</code>. */
	public final boolean countReverseKmers;
	
	/** K-mer counter. */
	private final CountMap counter;
	
	/** Maximum number of haplotypes that may be accepted. */
	private int maxHaplotypes;
	
	/** Maximum repeat count. */
	public final int maxRepeatCount;
	
	/** Default maximum number of haplotypes. */
	public static final int DEFAULT_MAX_HAPLOTYPES = 15;
	
	/**
	 * Create a new aligner.
	 * 
	 * @param kUtil K-mer utility.
	 * @param counter K-mer counter.
	 * @param alnWeight Alignment weights.
	 * @param countReverseKmers Count reverse complement k-mers in region statistics if
	 *   <code>true</code>.
	 * @param maxRepeatCount Maximum number of times k-mers may be repeated within an assembly
	 *    before aborting the assembly.
	 * @param trace If <code>true</code>, record the trace matrix in matrix form. This
	 *   option will have a significant performance impact, and it should only be used
	 *   for troubleshooting and debugging.
	 * 
	 * @throws NullPointerException If <code>kUtil</code>, <code>counter</code>, or
	 *   <code>alnWeight</code> is <code>null</code>.
	 * @throws IllegalArgumentException If <code>kUtil.kSize</code> is less than
	 *   <code>KestrelConstants.MIN_KMER_SIZE</code> or if <code>maxRepeatCount</code> is negative.
	 */
	public KmerAlignmentBuilder(KmerUtil kUtil, CountMap counter, AlignmentWeight alnWeight, boolean countReverseKmers, int maxRepeatCount, boolean trace)
			throws NullPointerException, IllegalArgumentException {
		
		logger = LoggerFactory.getLogger(KmerAlignmentBuilder.class);
		
		// Check arguments
		if (kUtil == null)
			throw new NullPointerException("K-mer utility is null");
		
		if (counter == null)
			throw new NullPointerException("K-mer counter is null");
		
		if (alnWeight == null)
			alnWeight = AlignmentWeight.get(null);  // Get default alignment weight
		
		if (maxRepeatCount < 0)
			throw new IllegalArgumentException("Maximum repeat count must be non-negative: " + maxRepeatCount);
		
		// Assign fields
		this.kUtil = kUtil;
		this.counter = counter;
		this.alnWeight = alnWeight;
		this.countReverseKmers = countReverseKmers;
		this.maxRepeatCount = maxRepeatCount;
		
		maxHaplotypes = DEFAULT_MAX_HAPLOTYPES;
		
		// Create aligner
		aligner = new KmerAligner(kUtil, alnWeight, trace);  // throws IllegalArgumentException
		
		return;
	}
	
	/**
	 * Find haplotypes in an active region.
	 * 
	 * @param activeRegion Active region.
	 * 
	 * @return Haplotypes in the region defined by <code>activeRegion</code> or an
	 *   empty array if no haplotypes were found.
	 * 
	 * @throws NullPointerException If <code>activeRegion</code> is <code>null</code>.
	 */
	public Haplotype[] getHaplotypes(ActiveRegion activeRegion)
			throws NullPointerException {
		
		// Check arguments
		if (activeRegion == null)
			throw new NullPointerException("Cannot find haplotypes in active region: null");
		
		logger.trace("Building haplotypes: {}", activeRegion.toString());
		
		// Create alignment
		aligner.init(activeRegion);
		
		// Get first kmer (alignment is already seeded with this k-mer)
		if (aligner.isReverse())
			return buildRev(activeRegion);
		
		return buildFwd(activeRegion);
	}
	
	/**
	 * Build haplotypes from left to right.
	 * 
	 * @param activeRegion Active region.
	 * 
	 * @return An array of haplotypes or an empty array if none are found.
	 */
	private Haplotype[] buildFwd(ActiveRegion activeRegion) {
		
		KmerAligner aligner = this.aligner;
		
		int[] kmer;       // Current k-mer
		int[] revKmer;    // Reverse k-mer for counting both strands
		int[] kmerCache;  // A temporary k-mer cache
		
		int minDepth;                 // Minimum depth of the haplotype currently being assembled
		RestoredState restoredState;  // Saved state that was restored
		
		Base base;      // Base to be added
		int revBase;    // A shifted complement-base of base to be restored to revKmer after queries
		int baseCount;  // Count with a query base
		int maxCount;   // Maximum count of all possible k-mer queries
		
		int lastKmerWord = kUtil.wordSize - 1;  // Last word of the k-mer
		int lastWordMask = (~ 0) << 2;          // Mask for removing the last base of the last k-mer word 
		
		int mswMask = kUtil.mswMask >>> 2;  // Mask for removing the first base of the first k-mer word
		
		int shiftC = Base.C.intVal << kUtil.mswMaskFirstShift;  // C shifted into the first position of the first word
		int shiftG = Base.G.intVal << kUtil.mswMaskFirstShift;  // G shifted into the first position of the first word
		int shiftT = Base.T.intVal << kUtil.mswMaskFirstShift;  // T shifted into the first position of the first word
		
		HaplotypeContainer haplotypeList;
		
		KmerHashSet kmerHash = new KmerHashSet(kUtil.kSize);  // Hash set for cycle detection
		int repeatCount = 0;  // Number of repeated k-mers
		
		// Init
		haplotypeList = new HaplotypeContainer(maxHaplotypes);
		
		kmer = activeRegion.getLeftEndKmer();
		revKmer = new int[kUtil.wordSize];
		
		// Set initial minimum depth
		minDepth = counter.get(kmer);
		
		if (countReverseKmers)
			minDepth += counter.get(revKmer);
				
		// Iterate until all possible paths from the initial k-mer are explored
		ITER_LOOP:
		while (true) {
			
			kUtil.revComplement(kmer, revKmer);
			
			BASE_LOOP:
			do {
				// A
				kUtil.append(kmer, Base.A);
				
				maxCount = counter.get(kmer);
				base = Base.A;
				revBase = shiftT;
				
				if (countReverseKmers) {
					kUtil.prepend(revKmer, Base.T);
					maxCount += counter.get(revKmer);
				}
				
				// C
				kmer[lastKmerWord] = kmer[lastKmerWord] & lastWordMask | Base.C.intVal;
				baseCount = counter.get(kmer);
				
				if (countReverseKmers) {
					revKmer[0] = (revKmer[0] & mswMask) | shiftG;
					baseCount += counter.get(revKmer);
				}
				
				if (baseCount > 0) {
					
					// Save the least count state and continue with the highest count state
					if (maxCount > baseCount) {
						
						if (logger.isTraceEnabled())
							logger.trace("Align split: Saving state {} (count={}, added=C, saving=this)", kUtil.toBaseString(kmer), baseCount);
						
						aligner.saveState(kUtil.copy(kmer, null), Base.C, (baseCount > minDepth && minDepth > 0) ? minDepth : baseCount, kmerHash, repeatCount);
					
					} else {
						
						if (maxCount > 0) {
							kmerCache = kUtil.copy(kmer, null);
							kmerCache[lastKmerWord] = kmerCache[lastKmerWord] & lastWordMask | base.intVal;
							
							if (logger.isTraceEnabled())
								logger.trace("Align split: Saving state {} (count={}, added=C, saving=cache)", kUtil.toBaseString(kmerCache), maxCount);
							
							aligner.saveState(kmerCache, base, (maxCount > minDepth && minDepth > 0) ? minDepth : maxCount, kmerHash, repeatCount);
						}
						
						maxCount = baseCount;
						base = Base.C;
						revBase = shiftG;
					}
				}
				
				// G
				kmer[lastKmerWord] = kmer[lastKmerWord] & lastWordMask | Base.G.intVal;
				baseCount = counter.get(kmer);
				
				if (countReverseKmers) {
					revKmer[0] = (revKmer[0] & mswMask) | shiftC;
					baseCount += counter.get(revKmer);
				}
				
				if (baseCount > 0) {
					
					// Save the least count state and continue with the highest count state
					if (maxCount > baseCount) {
						
						if (logger.isTraceEnabled())
							logger.trace("Align split: Saving state {} (count={}, added=G, saving=this)", kUtil.toBaseString(kmer), baseCount);
						
						aligner.saveState(kUtil.copy(kmer, null), Base.G, (baseCount > minDepth && minDepth > 0) ? minDepth : baseCount, kmerHash, repeatCount);
						
					} else {
						
						if (maxCount > 0) {
							kmerCache = kUtil.copy(kmer, null);
							kmerCache[lastKmerWord] = kmerCache[lastKmerWord] & lastWordMask | base.intVal;
							
							if (logger.isTraceEnabled())
								logger.trace("Align split: Saving state {} (count={}, added=G, saving=cache)", kUtil.toBaseString(kmerCache), maxCount);
							
							aligner.saveState(kmerCache, base, (maxCount > minDepth && minDepth > 0) ? minDepth : maxCount, kmerHash, repeatCount);
						}
						
						maxCount = baseCount;
						base = Base.G;
						revBase = shiftC;
					}
				}
				
				// T
				kmer[lastKmerWord] = kmer[lastKmerWord] & lastWordMask | Base.T.intVal;
				baseCount = counter.get(kmer);
				
				if (countReverseKmers) {
					revKmer[0] = (revKmer[0] & mswMask);
					baseCount += counter.get(revKmer);
				}
				
				if (baseCount > 0) {
					
					// Save the least count state and continue with the highest count state
					if (maxCount > baseCount) {
						
						if (logger.isTraceEnabled())
							logger.trace("Align split: Saving state {} (count={}, added=T, saving=this)", kUtil.toBaseString(kmer), baseCount);
						
						aligner.saveState(kUtil.copy(kmer, null), Base.T, (baseCount > minDepth && minDepth > 0) ? minDepth : baseCount, kmerHash, repeatCount);
					
					} else {
						
						if (maxCount > 0) {
							kmerCache = kUtil.copy(kmer, null);
							kmerCache[lastKmerWord] = kmerCache[lastKmerWord] & lastWordMask | base.intVal;
							
							if (logger.isTraceEnabled())
								logger.trace("Align split: Saving state {} (count={}, added=T, saving=cache)", kUtil.toBaseString(kmerCache), maxCount);
							
							aligner.saveState(kmerCache, base, (maxCount > minDepth && minDepth > 0) ? minDepth : maxCount, kmerHash, repeatCount);
						}
						
						maxCount = baseCount;
						base = Base.T;
						revBase = 0x0;
					}
				}
				
				// Stop if max count is 0
				if (maxCount == 0)
					break;
				
				// Replace max base
				kmer[lastKmerWord] = kmer[lastKmerWord] & lastWordMask | base.intVal;
				revKmer[0] = revKmer[0] & mswMask | revBase;
				
				// Cycle detection
				if (! kmerHash.add(kmer)) {
					++repeatCount;
					
					if (repeatCount > maxRepeatCount) {
						logger.trace("Cycle detected: {} (repeated k-mers = {}): Breaking forward assembly", kUtil.toBaseString(kmer), repeatCount);
						
						break BASE_LOOP;
					}
				}
				
				// Update depth
				if (maxCount < minDepth)
					minDepth = maxCount;
				
			} while (aligner.addBase(base));
			
			// Add haplotypes
			if (minDepth > 0) {
				for (Haplotype haplotype : aligner.getHaplotypes(counter, countReverseKmers)) {
					logger.trace("Adding haplotype: {}", haplotype);
					haplotypeList.add(haplotype);
				}
			}
			
			// Restore last state if the trace split
			restoredState = aligner.restoreState();
			
			if (restoredState == null)
				break ITER_LOOP;
			
			kmer = restoredState.kmer;
			minDepth = restoredState.minDepth;
			
			kmerHash = restoredState.kmerHash;
			repeatCount = restoredState.repeatCount;
		}
		
		logger.trace("Built {} haplotypes (fwd): {}", haplotypeList.size(), activeRegion.toString());
		
		// Return haplotypes
		return haplotypeList.toArray();
	}
	
	/**
	 * Set the maximum number of states to be saved. The least-likely states are trimmed when
	 * this threshold is reached.
	 * 
	 * @param maxState Maximum number of states.
	 * 
	 * @throws IllegalArgumentException If <code>maxState</code> is less than <code>1</code>.
	 * 
	 * @see KmerAligner#DEFAULT_MAX_STATE
	 */
	public void setMaxState(int maxState)
		throws IllegalArgumentException {
		
		if (maxState < 1)
			throw new IllegalArgumentException("Maximum number of states must not be less than 1: " + maxState);
		
		aligner.setMaxState(maxState);
		
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
		return aligner.getMaxState();
	}
	
	/**
	 * Set the maximum number of haplotypes to be saved. The least-likely haplotypes are trimmed when
	 * this threshold is reached.
	 * 
	 * @param maxHaplotypes Maximum number of haplotypes.
	 * 
	 * @throws IllegalArgumentException If <code>maxHaplotypes</code> is less than <code>1</code>.
	 * 
	 * @see #DEFAULT_MAX_HAPLOTYPES
	 */
	public void setMaxHaplotypes(int maxHaplotypes)
		throws IllegalArgumentException {
		
		if (maxHaplotypes < 1)
			throw new IllegalArgumentException("Maximum number of haplotypes must not be less than 1: " + maxHaplotypes);
		
		this.maxHaplotypes = maxHaplotypes;
		
		return;
	}
	
	/**
	 * Get the maximum number of haplotypes that may be saved.
	 * 
	 * @return The maximum number of haplotypes that may be saved.
	 * 
	 * @see #setMaxHaplotypes(int)
	 */
	public int getMaxHaplotypes() {
		return maxHaplotypes;
	}
	
	/**
	 * Build haplotypes from left to right.
	 * 
	 * @param activeRegion Active region.
	 * 
	 * @return An array of haplotypes or an empty array if none are found.
	 */
	private Haplotype[] buildRev(ActiveRegion activeRegion) {
		
		KmerAligner aligner = this.aligner;
		
		int[] kmer;       // Current k-mer
		int[] revKmer;    // Reverse k-mer for counting both strands
		int[] kmerCache;  // A temporary k-mer cache
		
		int minDepth;                 // Minimum depth of the haplotype currently being assembled
		RestoredState restoredState;  // Saved state that was restored
		
		Base base;      // Base to be added
		int revBase;    // A shifted complement-base of base to be restored to revKmer after queries
		int baseCount;  // Count with a query base
		int maxCount;   // Maximum count of all possible k-mer queries
		
		int lastKmerWord = kUtil.wordSize - 1;  // Last word of the k-mer
		int lastWordMask = (~ 0) << 2;          // Mask for removing the last base of the last k-mer word 
		
		int mswMask = kUtil.mswMask >>> 2;  // Mask for removing the first base of the first k-mer word
		
		int shiftC = Base.C.intVal << kUtil.mswMaskFirstShift;  // C shifted into the first position of the first word
		int shiftG = Base.G.intVal << kUtil.mswMaskFirstShift;  // G shifted into the first position of the first word
		int shiftT = Base.T.intVal << kUtil.mswMaskFirstShift;  // T shifted into the first position of the first word
		
		HaplotypeContainer haplotypeList;
		
		KmerHashSet kmerHash = new KmerHashSet(kUtil.kSize);  // Hash set for cycle detection
		int repeatCount = 0;  // Number of repeated k-mers
		
		// Init
		haplotypeList = new HaplotypeContainer(maxHaplotypes);
		
		kmer = activeRegion.getRightEndKmer();
		revKmer = new int[kUtil.wordSize];
		
		// Set initial minimum depth
		minDepth = counter.get(kmer);
		
		if (countReverseKmers)
			minDepth += counter.get(revKmer);
		
		// Iterate until all possible paths from the initial k-mer are explored
		ITER_LOOP:
		while (kmer != null) {
			
			kUtil.revComplement(kmer, revKmer);
			
			BASE_LOOP:
			do {
				// A
				kUtil.prepend(kmer, Base.A);
				
				maxCount = counter.get(kmer);
				base = Base.A;
				revBase = Base.T.intVal;
				
				if (countReverseKmers) {
					kUtil.append(revKmer, Base.T);
					maxCount += counter.get(revKmer);
				}
				
				// C
				kmer[0] = kmer[0] & mswMask | shiftC;
				baseCount = counter.get(kmer);
				
				if (countReverseKmers) {
					revKmer[lastKmerWord] = (revKmer[lastKmerWord] & lastWordMask) | Base.G.intVal;
					baseCount += counter.get(revKmer);
				}
				
				if (baseCount > 0) {
					
					// Save the least count state and continue with the highest count state
					if (maxCount > baseCount) {
						
						if (logger.isTraceEnabled())
							logger.trace("Align split: Saving state {} (count={}, added=C, saving=this)", kUtil.toBaseString(kmer), baseCount);
						
						aligner.saveState(kUtil.copy(kmer, null), Base.C, (baseCount > minDepth && minDepth > 0) ? minDepth : baseCount, kmerHash, repeatCount);
					
					} else {
						
						if (maxCount > 0) {
							kmerCache = kUtil.copy(kmer, null);
							kmerCache[0] = kmerCache[0] & mswMask | (base.intVal << kUtil.mswMaskFirstShift);
							
							if (logger.isTraceEnabled())
								logger.trace("Align split: Saving state {} (count={}, added=C, saving=cache)", kUtil.toBaseString(kmerCache), maxCount);
							
							aligner.saveState(kmerCache, base, (maxCount > minDepth && minDepth > 0) ? minDepth : maxCount, kmerHash, repeatCount);
						}
						
						maxCount = baseCount;
						base = Base.C;
						revBase = Base.G.intVal;
					}
				}
				
				// G
				kmer[0] = kmer[0] & mswMask | shiftG;
				baseCount = counter.get(kmer);
				
				if (countReverseKmers) {
					revKmer[lastKmerWord] = (revKmer[lastKmerWord] & lastWordMask) | Base.C.intVal;
					baseCount += counter.get(revKmer);
				}
				
				if (baseCount > 0) {
					
					// Save the least count state and continue with the highest count state
					if (maxCount > baseCount) {
						
						if (logger.isTraceEnabled())
							logger.trace("Align split: Saving state {} (count={}, added=G, saving=this)", kUtil.toBaseString(kmer), baseCount);
						
						aligner.saveState(kUtil.copy(kmer, null), Base.G, (baseCount > minDepth && minDepth > 0) ? minDepth : baseCount, kmerHash, repeatCount);
						
					} else {
						
						if (maxCount > 0) {
							kmerCache = kUtil.copy(kmer, null);
							kmerCache[0] = kmerCache[0] & mswMask | (base.intVal << kUtil.mswMaskFirstShift);
							
							if (logger.isTraceEnabled())
								logger.trace("Align split: Saving state {} (count={}, added=G, saving=cache)", kUtil.toBaseString(kmerCache), maxCount);
							
							aligner.saveState(kmerCache, base, (maxCount > minDepth && minDepth > 0) ? minDepth : maxCount, kmerHash, repeatCount);
						}
						
						maxCount = baseCount;
						base = Base.G;
						revBase = Base.C.intVal;
					}
				}
				
				// T
				kmer[0] = kmer[0] & mswMask | shiftT;
				baseCount = counter.get(kmer);
				
				if (countReverseKmers) {
					revKmer[lastKmerWord] = (revKmer[lastKmerWord] & lastWordMask);
					baseCount += counter.get(revKmer);
				}
				
				if (baseCount > 0) {
					
					// Save the least count state and continue with the highest count state
					if (maxCount > baseCount) {
						
						if (logger.isTraceEnabled())
							logger.trace("Align split: Saving state {} (count={}, added=T, saving=this)", kUtil.toBaseString(kmer), baseCount);
						
						aligner.saveState(kUtil.copy(kmer, null), Base.T, (baseCount > minDepth && minDepth > 0) ? minDepth : baseCount, kmerHash, repeatCount);
					
					} else {
						
						if (maxCount > 0) {
							kmerCache = kUtil.copy(kmer, null);
							kmerCache[0] = kmerCache[0] & mswMask | (base.intVal << kUtil.mswMaskFirstShift);
							
							if (logger.isTraceEnabled())
								logger.trace("Align split: Saving state {} (count={}, added=T, saving=cache)", kUtil.toBaseString(kmerCache), maxCount);
							
							aligner.saveState(kmerCache, base, (maxCount > minDepth && minDepth > 0) ? minDepth : maxCount, kmerHash, repeatCount);
						}
						
						maxCount = baseCount;
						base = Base.T;
						revBase = 0x0;
					}
				}
				
				// Stop if max count is 0
				if (maxCount == 0)
					break;
				
				// Replace max base
				kmer[0] = kmer[0] & mswMask | (base.intVal << kUtil.mswMaskFirstShift);
				revKmer[lastKmerWord] = revKmer[lastKmerWord] & lastWordMask | revBase;
				
				// Cycle detection
				if (! kmerHash.add(kmer)) {
					++repeatCount;
					
					if (repeatCount > maxRepeatCount) {
						logger.trace("Cycle detected: {}: Breaking reverse assembly");
						
						break BASE_LOOP;
					}
				}
				
			} while (aligner.addBase(base));
			
			for (Haplotype haplotype : aligner.getHaplotypes(counter, countReverseKmers)) {
				logger.trace("Adding haplotype: {}", haplotype);
				haplotypeList.add(haplotype);
			}
			
			// Restore last state if the trace split
			restoredState = aligner.restoreState();
			
			if (restoredState == null)
				break ITER_LOOP;
			
			kmer = restoredState.kmer;
			minDepth = restoredState.minDepth;
		}
		
		logger.trace("Built {} haplotypes (rev): {}", haplotypeList.size(), activeRegion.toString());
		
		// Return haplotypes
		return haplotypeList.toArray();
	}
}
