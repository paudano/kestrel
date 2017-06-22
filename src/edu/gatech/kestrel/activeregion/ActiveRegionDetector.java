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

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import edu.gatech.kanalyze.KAnalyzeConstants;
import edu.gatech.kanalyze.util.Base;
import edu.gatech.kanalyze.util.kmer.KmerUtil;
import edu.gatech.kestrel.align.AlignmentWeight;
import edu.gatech.kestrel.align.KmerAligner;
import edu.gatech.kestrel.align.KmerAlignmentBuilder;
import edu.gatech.kestrel.counter.CountMap;
import edu.gatech.kestrel.refreader.ReferenceRegion;

/**
 * Detects active regions in reference sequence by searching for areas where the k-mer count
 * changes abruptly. These active regions become the targets of variant detection.
 * <p/>
 * Active regions are identified by searching for a significant loss in counts over the expected
 * k-mers in a reference sequence. When the count difference between neighboring k-mers is
 * positive (a k-mer has a lower count than its left-neighbor), then a right-scan is initiated.
 * This scan is anchored on the left neighbor, which is the last k-mer with a high count before the
 * loss. The scan terminates when the count for some downstream k-mer reaches a threshold close
 * to the count of the left-anchor k-mer, and this k-mer becomes the right anchor. The active
 * region is defined as all bases from the start of the left anchor to the end of the right
 * anchor.
 * <p/>
 * Different datasets have different characteristics. The read depth may vary significantly
 * among samples and among reference within the same sample. In order to adapt to the data,
 * a difference threshold is found dynamically using the data and parameters
 * to guide threshold selection. The threshold is determined by taking all k-mer counts in a
 * reference, computing the absolute value of the k-mer difference
 * between each neighboring k-mer, and choosing a quantile from that set of differences. For
 * example, if the quantile is <code>0.90</code>, then at most 10% of the neighboring k-mers
 * will have a difference large enough to trigger an active region scan. The higher the quantile,
 * the more stringent the threshold becomes, and a larger difference will be required to trigger
 * an active region scan. If this threshold is less than the minimum difference, then the minimum
 * difference is used instead.
 * <p/>
 * Some variants may occur close enough to an end of the reference such that there is no anchor
 * k-mer on the right or left side. A right scan, as described above, may not find a k-mer to
 * anchor the right side into. When this occurs, an active region that extends to the end of the
 * reference may be defined, and the aligner may treat these differently. Since active regions with
 * anchors on both sides provide the strongest evidence for variants called within them, the scanner
 * may be configured to ignore regions without an anchor on one side. However, this will result in a
 * failure to identify true variants near the ends of the reference. If a variant is <code>k</code>
 * bases or less from the end, all variants near it are missed until a variant is more than
 * <code>k</code> bases away from another. In this case, un-anchored ends or a smaller k-mer size
 * is necessary to recover them.
 * <p/>
 * If a variant occurs near the left end of the reference, this detector will see a large negative
 * count difference between neighbors for the first k-mer to leave the variant region. In this case,
 * the right k-mer will have a higher count than the left. As with a right scan that reaches the end
 * with no anchor, an active region may also extend to the left end.
 * <p/>
 * If a k-mer in the reference is shared with another part of the genome, k-mer counts from both
 * regions are mixed together in k-mer space (the collection of k-mers generated from the sample).
 * This will cause a peak where the k-mer count increases dramatically and decreases again. The
 * edges of this peak will look like a large negative change in counts between two neighbors followed
 * by a large positive change. This detector employs peak detection to step around them. The number
 * of k-mers scanned while searching for a peak is tunable.
 * <p/>
 * Active regions missing an anchor on one end may be limited in length to prevent a set of high
 * k-mer counts from arbitrarily running the active region to the end. This could occur in very noisy
 * data, a sharp decline in read depth, or a large area of homology. The parameter
 * <code>endScanLimit</code> is a multiple of the maximum gap length allowable by the alignment
 * weight vector (assuming no other variants). This gap length is calculated as
 * <code>endScanLimit * ((init_score - gap_open) / gap_extend + 1)</code>. If <code>endScanLimit</code>
 * is <code>0.0</code>, then end calling is effectively disabled. It may be set to a high number to
 * allow end scans that are arbitrarily long. Overflow detection is employed, so this value may be set
 * arbitrarily high.
 * <p/>
 * A k-mer count is a surrogate for read depth, and read depth is not typically uniform over the
 * reference region. When the read depth declines over an active region, the counts may not recover
 * and the right end of an active region may be missed. This can result in active regions that are
 * far too long, and it will cause the anchored k-mer aligner to slow down or consume more memory
 * than necessary. To account for this, an exponential decay function is applied to the recovery
 * threshold. As the scan extends to the right, the recovery threshold is reduced at each step.
 * <p/>
 * The exponential decay function can be tuned by two parameters. The first, <code>expDecayMin</code>,
 * defines a lower asymptotic bound of the function. <code>expDecayMin</code> is a proportion (a value
 * between <code>0.0</code> and <code>1.0</code> that is multiplied to the original k-mer count. The
 * decay value will start at the initial recovery threshold and decrease approaching the minimum count,
 * but never passing it. If this value is <code>0.0</code>, k-mer count recovery threshold may decline
 * to <code>1</code>. If this value is <code>1.0</code>, the decay function is not used and the detector
 * falls back to finding a k-mer with a count within the difference threshold of the anchor k-mer count.
 * <p/>
 * The second parameter, <code>expDecayAlpha</code> determines how quickly the recovery value decays. It is
 * defined as the proportion of decay from the left anchor k-mer count to the minimum decay value
 * at <code>k</code> k-mers from the left anchor. Since this is a proportion, its value is always
 * between <code>0.0</code> and <code>1.0</code> (exclusive for both value; may not be <code>0.0</code>
 * or <code>1.0</code>). In other words, at <code>k</code> k-mers from the left anchor, the
 * recovery value will be <code>expDecayAlpha</code> multiplied by the left anchor k-mer count. At
 * <code>2 * k</code> k-mers, it will be <code>expDecayAlpha^2</code>. This is true for all values
 * <code>n * k</code>. A higher value causes the decay to occur more slowly, and a lower value will
 * cause the decay to occur rapidly.
 * <p/>
 * Putting this together, assume the left anchor k-mer has a count of 200 and the minimum is set to
 * <code>0.50</code>. The exponential decay will start at <code>200</code> and decline to no less
 * than <code>100</code>. The range the function decays over is <code>100</code>
 * (<code>200 - 100</code>). If <code>expDecayAlpha</code> is <code>0.80</code>, and the k-mer size is
 * <code>48</code>, then at <code>48</code> bases from the anchor, the recovery value will be
 * <code>180</code> (<code>0.80 * 100 + 100</code>). At <code>96</code> k-mers, it will be <code>164</code>
 * (<code>100 * 0.80^2 + 100</code>), and so on. 
 * <p/>
 * In an exponential decay function, lambda controls the rate of decay. It is the rate of change of the
 * function as a proportion of the amount (a larger amount will decline more than a smaller amount in a
 * given time). Alpha is a similar measure, however, it is related to the k-mer size for convenience.
 * It is much easier to intuitively choose the alpha value than it is a lambda value. Alpha and the
 * k-mer size are used to compute lambda (see below), and lambda is never set directly.
 * <p/>
 * The exponential decay function is defined as:<br/>
 * <code>f(x) = range * exp(-x * lambda) + min</code><br/>
 * <ol>
 *   <li><code>f(x)</code> := The recovery threshold <code>x</code> k-mers from the anchor</li>
 *   <li><code>range</code> := The length from the initial recovery value to the minimum value</li>
 *   <li><code>exp</code> := The exponential function (<em>e</em> raised to some power)</li>
 *   <li><code>min</code> := <code>expDecayMin</code> * left-anchor count</li>
 *   <li><code>lambda</code> := <code>-log(expDecayAlpha) / k</code></li>
 *   <li><code>expDecayAlpha</code> := The proportion of decay at f(k) k-mers from the start decay (<code>f(k) = expDecayAlpha * f(0)</code>)</li>
 *   <li><code>k</code> := The k-mer size</li>
 * </ol>
 */
public class ActiveRegionDetector {
	
	/** Logger object. */
	private static final Logger logger = LoggerFactory.getLogger(ActiveRegionDetector.class);
	
	/** Size of k-mers in this detector. */
	public final int kSize;
	
	/** K-mer counts used to identify active regions. */
	private final CountMap counter;
	
	/** K-mer utility without minimizers. */
	private final KmerUtil kUtil;
	
	/**
	 * If <code>true</code>, active regions must be bordered on both sides by unaltered k-mers. The
	 * variants called in these active regions are supported by better evidence than those called
	 * in active regions that reach an end, however, variants close to the end of a reference may
	 * be missed.
	 */
	private boolean anchorBothEnds;
	
	/** If <code>true</code>, count reverse k-mers in region statistics. */
	private boolean countReverseKmers;
	
	/** Minimum k-mer difference to flag a potential active region. */
	private int minimumDifference;
	
	/**
	 * K-mer count differences must be at least this quantile of total difference to trigger an
	 * active region scan over the reference sequence. <code>0.0</code> disables picking the
	 * threshold by quantile and <code>minimumDifference</code> is used as the threshold. If
	 * the quantile is less than <code>minimumDifference</code>, then the minimum is used. The
	 * difference of neighboring k-mers is defined as their distance apart (the absolute value
	 * of one k-mer count subtracted by the count of its neighbor).
	 */
	private double differenceQuantile;
	
	/**
	 * For peak detection, scan forward over a suspected peak this number of k-mers. If
	 * <code>0</code>, disable peak detection.
	 */
	private int peakScanLength;
	
	/**
	 * The maximum length of a scan over the reference
	 * is determined by taking the maximum length of a gap (assuming no aligned bases)
	 * with one added (to support end calling when gaps are disabled) as defined by the alignment
	 * score vector and multiplying it by this value.
	 */
	private double scanLimitFactor;
	
	/**
	 * The end scan limit length after computing it from the alignment weight vector and
	 * <code>scanLimitFactor</code>.
	 * 
	 * @see #scanLimitFactor
	 */
	private int scanLimit;
	
	/**
	 * Lower asymptotic bound for the exponential decay function applied the recovery threshold
	 * in an active region scan. The minimum decay value will not be less than this proportion
	 * of the anchor k-mer count. This value is between <code>0.0</code> and <code>1.0</code>.
	 * If <code>0.0</code>, the count may decline to <code>1</code>, and if <code>1.0</code>,
	 * the exponential decay heuristic is effectively disabled.
	 */
	private double expDecayMin;
		
	/**
	 * The decay function is applied to a range of k-mer counts from the initial recovery threshold
	 * to the minimum decay (defined by <code>expDecayMin</code>). At <code>k</code> k-mers, the
	 * function will have declined proportion of the decay range. This controls the rate of the decay,
	 * and relates it to the k-mer size. 
	 */
	private double expDecayAlpha;
	
	/** Exponential decay lambda. */
	private double expDecayLambda;
	
	/** Alignment weights. */
	private AlignmentWeight alignmentWeight;
	
	/**
	 * If <code>true</code>, then allow haplotypes and variants to occur over ambiguous regions.
	 */
	private boolean callAmbiguousRegions;
	
	/** If the k-mer count does not recover, search for a region where the k-mer count increases abruptly. */ 
	private boolean recoverRightAnchor;
	
	/**
	 * When set, emit active regions with no variant haplotypes. These will be inserted in the gaps
	 * between active regions with variants.
	 */
	private boolean emitWildtypeActiveRegions;
	
	/**
	 * Generate a full alignment matrix for each haplotype and attach it to the haplotype object.
	 * This can be useful for debugging, but it has significant performance and memory consumption
	 * implications. 
	 */
	private boolean traceHaplotypeAlignment;
	
	/** Builds haplotypes from active regions. */
	private KmerAlignmentBuilder alignmentBuilder;
	
	/** Get the complement of a base. */
	private static final Base[] BASE_COMPLEMENT = new Base[] {Base.T, Base.G, Base.C, Base.A};
	
	/** Translate a nucleotide from a byte to a base. */
	private static final Base[] BYTE_TO_BASE = KAnalyzeConstants.getByteToBaseArray();
	
	
	/** Default minimum k-mer difference. */
	public static final int DEFAULT_MINIMUM_DIFFERENCE = 5;
	
	/** Default difference quantile. */
	public static final double DEFAULT_DIFFERENCE_QUANTILE = 0.90;
	
	/** Default anchor property. */
	public static final boolean DEFAULT_ANCHOR_BOTH_ENDS = true;
	
	/** Default reverse-complement counting property. */
	public static final boolean DEFAULT_COUNT_REVERSE_KMERS = true;
	
	/** Default property to call variants in regions with ambiguous bases. */
	public static final boolean DEFAULT_CALL_AMBIGUOUS_REGIONS = true;
	
	/** Default number of k-mers to scan for peak detection. */
	public static final int DEFAULT_PEAK_SCAN_LENGTH = 7;
	
	/** Default end-scan limit factor. */
	public static final double DEFAULT_SCAN_LIMIT_FACTOR = 5.0;
	
	/** Default alignment trace attribute. */
	public static final boolean DEFAULT_TRACE_HAPLOTYPE_ALIGNMENT = false;
	
	/** Default exponential decay minimum. */
	public static final double DEFAULT_EXP_MIN = 0.55;
	
	/** Default exponential decay alpha. */
	public static final double DEFAULT_EXP_ALPHA = 0.80;
	
	/** Default value for right anchor recovery. */
	public static final boolean DEFAULT_RECOVER_RIGHT_ANCHOR = true;
	
	/** Default value for emitting wildtype active regions. */
	public static final boolean DEFAULT_EMIT_WILDTYPE_ACTIVE_REGIONS = false;
	
	/**
	 * Create a new active region detector.
	 * 
	 * @param counter Counter with k-mers and associated k-mer counts in a data set.
	 * 
	 * @throws NullPointerException If <code>counter</code> is <code>null</code>.
	 */
	public ActiveRegionDetector(CountMap counter)
			throws NullPointerException {
		
		// Check arguments
		if (counter == null)
			throw new NullPointerException("K-mer counter is null");
		
		this.kSize = counter.kUtil.kSize;
		this.counter = counter;
		
		kUtil = KmerUtil.get(kSize);
		
		setAlignmentWeight(AlignmentWeight.get(null));
		
		setCountReverseKmers(DEFAULT_COUNT_REVERSE_KMERS);
		setMinimumDifference(DEFAULT_MINIMUM_DIFFERENCE);
		setDifferenceQuantile(DEFAULT_DIFFERENCE_QUANTILE);
		setAnchorBothEnds(DEFAULT_ANCHOR_BOTH_ENDS);
		setCallAmbiguousRegions(DEFAULT_CALL_AMBIGUOUS_REGIONS);
		setRecoverRightAnchor(DEFAULT_RECOVER_RIGHT_ANCHOR);
		
		setPeakScanLength(DEFAULT_PEAK_SCAN_LENGTH);
		setDecayMinimum(DEFAULT_EXP_MIN);
		setDecayAlpha(DEFAULT_EXP_ALPHA);
		setScanLimitFactor(DEFAULT_SCAN_LIMIT_FACTOR);
		setEmitWildtypeActiveRegions(DEFAULT_EMIT_WILDTYPE_ACTIVE_REGIONS);
		
		setTraceHaplotypeAlignment(DEFAULT_TRACE_HAPLOTYPE_ALIGNMENT);
		
		initAlignmentBuilder();  // Call after all other fields are initialized
		
		return;
	}
	
	/**
	 * Get active regions.
	 * 
	 * @param refRegion Reference region.
	 * 
	 * @return An array of active regions.
	 * 
	 * @throws NullPointerException If <code>refRegion</code> is <code>null</code>.
	 */
	public ActiveRegionContainer getActiveRegions(ReferenceRegion refRegion)
			throws NullPointerException {
		
		// Declarations
		List<RegionHaplotype> haplotypeList;  // List of haplotypes and their active regions
		
		int diffThresholdR;  // Difference threshold for a forward scan
		int diffThresholdL;  // Difference threshold for a reverse scan
		
		int[] refCount;     // Array of counts for each k-mer
		int refCountIndex;  // Current element of refCount
		int refCountSize;   // Length of refCount
		
		int recoveryValue;     // Value computed to end a specific active region
		int scanEndIndex;      // refCountIndex where an active region scan ends
		int lastScanIndex;     // Scanning will not pass this value
		
		int peakScanIndex;    // Index of refCount for a peak scan
		int peakScanLimit;    // Max value of peakScanIndex
		int nPeak;            // The number of peaks in an active region
		int lastValleyIndex;  // The end location of a valley (at least k k-mers between peaks)
		
		int countL;     // Left count (refCount[refCountIndex - 1])
		int countR;     // Right count (refCount[refCountIndex])
		int countDiff;  // countL - countR
		
		boolean noAmbiguousBases;  // Set to true if no ambiguous bases are allowed in active regions
		int lastRegionEnd;         // Index of the last region found; Avoids region collisions
		
		double expDecayRange;     // Range of exponential decay from anchor k-mer count to the lower bound
		double expDecayMinValue;  // Exponential decay asymptotic lower bound
		
		KmerAlignmentBuilder alignmentBuilder;  // Builds haplotypes from active regions
		
		ActiveRegion activeRegion;  // An active region (variant calling will be attempted over this region)
		Haplotype[] haplotypes;     // Haplotypes called from an active region
		
		// Check arguments
		if (refRegion == null)
			throw new NullPointerException("Reference sequence is null");
		
		// Init
		alignmentBuilder = this.alignmentBuilder;
		
		haplotypeList = new ArrayList<>();
		
		noAmbiguousBases = ! callAmbiguousRegions;
		lastRegionEnd = 0;
		
		// Get counts
		refCount = getCounts(refRegion);
		refCountSize = refCount.length;
		
		if (refCountSize < 2)
			return new ActiveRegionContainer(refRegion, null, refCount);
		
		// Set the scan and recovery thresholds
		diffThresholdR = getDifferenceThreshold(refCount, 0, refCountSize) - 1;  // Subtract 1: Code uses < and >, not <= and >=
		diffThresholdL = -diffThresholdR;
		
		logger.trace("Difference threshold: {}", diffThresholdR + 1);
		
		// Scan counts for variants
		refCountIndex = 1;
		countL = refCount[0];
				
		// Search for active regions over the reference
		REF_SEARCH:
		while (refCountIndex < refCountSize) {
			
			countR = refCount[refCountIndex];
			
			countDiff = countL - countR;
			
			if (countDiff > diffThresholdR) {
				// Scan right
				
				scanEndIndex = refCountIndex + 1;
				
				// Reset peak stats
				nPeak = 0;
				peakScanIndex = 0;
				lastValleyIndex = 0;
				
				lastScanIndex = refCountIndex + scanLimit;
				
				if (lastScanIndex < 0)
					lastScanIndex = Integer.MAX_VALUE;
				
				// Do scan
				SCAN_LOOP:
				while (scanEndIndex <= lastScanIndex) {
					
					if (expDecayMin == 1.0) {
						// Constant threshold scan
						
						recoveryValue = countL - diffThresholdR;
						
						if (recoveryValue < 1)
							recoveryValue = 1;
						
						logger.trace("Right scan: Constant search: index={}, threshold={}, value={}", refCountIndex, diffThresholdR, recoveryValue);
						
						while (scanEndIndex < refCountSize && refCount[scanEndIndex] < recoveryValue)
							++scanEndIndex;
						
					} else {
						// Exponential-decay threshold scan
						
						// Exponential decay recovery value
						expDecayMinValue = countL * expDecayMin;
						
						if (expDecayMinValue < 1)
							expDecayMinValue = 1;
						
						expDecayRange = countL - expDecayMinValue;
						
						logger.trace("Right scan: Exponential search: index={}, range={}, min={}, lambda={}", refCountIndex, expDecayRange, expDecayMinValue, expDecayLambda);
						
						while (
								scanEndIndex < refCountSize &&
								refCount[scanEndIndex] <
									expDecayRange * Math.exp((refCountIndex - scanEndIndex) * expDecayLambda) + expDecayMinValue  // range * e^{-(scanEndIndex - refCountIndex) * lambda} + min
							) {
							
							++scanEndIndex;
						}
						
						// Set recovery value for peak detection
						recoveryValue = (int) (expDecayRange * Math.exp((refCountIndex - scanEndIndex) * expDecayLambda) + expDecayMinValue);
					}
					
					// Setup peak detection (ensures a stable recovery)
					if (peakScanLength == 0)
						break SCAN_LOOP;
					
					
					// Detect valley (at least k k-mers between peaks), peakScanIndex points to the first
					// low-count k-mer. A valley is at least k k-mers without a peak
					if (peakScanIndex > 0 && scanEndIndex - peakScanIndex >= kSize)  // Not first peak
						lastValleyIndex = scanEndIndex; 
					
					else if (peakScanIndex == 0 && scanEndIndex - refCountIndex >= kSize)  // First peak
						lastValleyIndex = scanEndIndex;
					
					peakScanIndex = scanEndIndex;
					peakScanLimit = scanEndIndex + peakScanLength;
					
					if (peakScanLimit > refCountSize)
						peakScanLimit = refCountSize;
					
					// Perform peak detection
					while (peakScanIndex < peakScanLimit) {
						
						if (refCount[peakScanIndex] < recoveryValue) {
							
							++nPeak;
							
							scanEndIndex = peakScanIndex;
							
							// Check for an excessive number of peaks
							if (nPeak > 3) {
								
								if ((scanEndIndex - refCountIndex) / nPeak < kSize) {
									
									if (lastValleyIndex > 0) {  // Go back to the end of the last region where at least kSize k-mers were low between peaks
										logger.trace("Right-scan (from={}, to={}, left={}, right={}, diff={}): Peak scan falling-back to last valley location due to frequent peaks (n={}, valley={})", refCountIndex - 1, scanEndIndex, countL, countR, countDiff, nPeak, lastValleyIndex);
										scanEndIndex = lastValleyIndex;
										break SCAN_LOOP;
									}
									
									// Excessive peaks
									
									logger.trace("Right-scan (from={}, to={}, left={}, right={}, diff={}): Peaks too frequent (n={}, interavg={})", refCountIndex - 1, scanEndIndex, countL, countR, countDiff, nPeak, (scanEndIndex - refCountIndex) / nPeak);
									
									countL = countR;
									++refCountIndex;
									continue REF_SEARCH;
								}	
							}
							
							continue SCAN_LOOP;
						}
						
						++peakScanIndex;
					}
					
					// Check for peaks at the end of the reference
					if (peakScanIndex == refCountSize && lastValleyIndex > 0) {  // Go back to the end of the last region where at least kSize k-mers were low between peaks
						scanEndIndex = lastValleyIndex;
						break SCAN_LOOP;
					}
					
					// No peak found, end scan
					break SCAN_LOOP;
				}
				
				// Check for scan limit
				if (scanEndIndex > lastScanIndex) {
					logger.trace("Right-scan (from={}, to={}, left={}, right={}, diff={}): Discarded: Region length reached the scan limit length ({}) before a recovery k-mer was found", refCountIndex - 1, scanEndIndex, countL, countR, countDiff, scanLimit);
					
					countL = countR;
					++refCountIndex;
					continue REF_SEARCH;
				}
				
				// Check anchors
				if (scanEndIndex < refCountSize) {
					// Active region is right-anchored
					
					// Ignore if region is too short.
					if (scanEndIndex - refCountIndex < kSize - 1) {
						
						logger.trace("Right-scan (from={}, to={}, left={}, right={}, diff={}): Discarded: Region length is less than k - 1 ({})", refCountIndex - 1, scanEndIndex, countL, countR, countDiff, kUtil.kSize - 1);
						
						countL = countR;
						++refCountIndex;
						continue REF_SEARCH;
					}
					
					// Ignore if regions contains ambiguous bases and it is not allowed
					if (noAmbiguousBases && refRegion.containsAmbiguousByIndex(refCountIndex, scanEndIndex + kSize - 1)) {
						
						logger.trace("Right-scan (from={}, to={}, left={}, right={}, diff={}): Discarded: Ambiguous regions are disabled", refCountIndex - 1, scanEndIndex, countL, countR, countDiff);
						
						countL = countR;
						++refCountIndex;
						continue REF_SEARCH;
					}
					
					activeRegion = new ActiveRegion(refRegion, refCountIndex - 1, scanEndIndex, refCount, kUtil);
					
					// Get haplotypes
					haplotypes = alignmentBuilder.getHaplotypes(activeRegion);
					
					if (haplotypes.length == 0 || (haplotypes.length == 1 && haplotypes[0].isWildtype())) {
						
						if (haplotypes.length == 0)
							logger.trace("Right-scan (from={}, to={}, left={}, right={}, diff={}): Discarded: No haplotypes", refCountIndex - 1, scanEndIndex, countL, countR, countDiff);
						else
							logger.trace("Right-scan (from={}, to={}, left={}, right={}, diff={}): Discarded: Found only the reference haplotype", refCountIndex - 1, scanEndIndex, countL, countR, countDiff);
						
						countL = countR;
						++refCountIndex;
						continue REF_SEARCH;
					}
					
					logger.trace("Right-scan (from={}, to={}, left={}, right={}, diff={}): Found {} haplotype(s)", refCountIndex - 1, scanEndIndex, countL, countR, countDiff, haplotypes.length);
					
					haplotypeList.add(new RegionHaplotype(activeRegion, haplotypes));
					
					countL = refCount[scanEndIndex];
					refCountIndex = scanEndIndex + 1;
					lastRegionEnd = scanEndIndex;
					
					
				} else {
					// Active region is not right-anchored
					
					if (recoverRightAnchor) {
						// Attempt to find the most likely right anchor even though the k-mer count did not recover
						
						scanEndIndex = refCountIndex + kSize;
						
						RECOVERY_SCAN:
						while (scanEndIndex < refCountSize) {
							
							if (refCount[scanEndIndex] - refCount[scanEndIndex - 1] > diffThresholdR) {
								logger.trace("Right-scan (from={}, to=RIGHT_END, left={}, right={}, diff={}): Recovery scan found k-mer increase at {}", refCountIndex - 1, countL, countR, countDiff, scanEndIndex);
								break RECOVERY_SCAN;
							}
							
							++scanEndIndex;
						}
					}
					
					if (scanEndIndex >= refCountSize) {
						
						scanEndIndex = -1;
						
						// Ignore if all active regions must be anchored or if the end-scan limit is exceeded
						if (anchorBothEnds) {
							logger.trace("Right-scan (from={}, to=RIGHT_END, left={}, right={}, diff={}): Discarded: Anchor-both-ends is set", refCountIndex - 1, countL, countR, countDiff);
							
							countL = countR;
							++refCountIndex;
							
							continue REF_SEARCH;
						}
						
						// Check end-scan limit
						if ((refCountSize - refCountIndex) > scanLimit) {
							logger.trace("Right-scan (from={}, to=RIGHT_END, left={}, right={}, diff={}): Discarded: Scan limit reached ({})", refCountIndex - 1, countL, countR, countDiff, scanLimit);
							
							countL = countR;
							++refCountIndex;
							continue REF_SEARCH;
						}
						
						// Ignore if regions contains ambiguous bases and it is not allowed
						if (noAmbiguousBases && refRegion.containsAmbiguousByIndex(refCountIndex, refCountSize - 1)) {
							
							logger.trace("Right-scan (from={}, to=RIGHT_END, left={}, right={}, diff={}): Discarded: Ambiguous regions are disabled", refCountIndex - 1, countL, countR, countDiff);
							
							countL = countR;
							++refCountIndex;
							continue REF_SEARCH;
						}
						
					}
					
					activeRegion = new ActiveRegion(refRegion, refCountIndex - 1, scanEndIndex, refCount, kUtil);
					
					// Get haplotypes
					haplotypes = alignmentBuilder.getHaplotypes(activeRegion);
					
					if (haplotypes.length == 0 || (haplotypes.length == 1 && haplotypes[0].isWildtype())) {
						
						if (scanEndIndex == -1) {
							
							if (haplotypes.length == 0)
								logger.trace("Right-scan (from={}, to=RIGHT_END, left={}, right={}, diff={}): Discarded: No haplotypes", refCountIndex - 1, countL, countR, countDiff);
							else
								logger.trace("Right-scan (from={}, to=RIGHT_END, left={}, right={}, diff={}): Discarded: Found only the reference haplotype", refCountIndex - 1, countL, countR, countDiff);
							
						} else {
							if (haplotypes.length == 0)
								logger.trace("Right-scan (from={}, to={}, left={}, right={}, diff={}): Discarded: No haplotypes", refCountIndex - 1, scanEndIndex, countL, countR, countDiff);
							else
								logger.trace("Right-scan (from={}, to={}, left={}, right={}, diff={}): Discarded: Found only the reference haplotype", refCountIndex - 1, scanEndIndex, countL, countR, countDiff);
						}
						
						countL = countR;
						++refCountIndex;
						continue REF_SEARCH;
					}
					
					if (scanEndIndex == -1) {
						logger.trace("Right-scan (from={}, to=RIGHT_END, left={}, right={}, diff={}): Found {} haplotype(s)", refCountIndex - 1, countL, countR, countDiff, haplotypes.length);
						
						refCountIndex = refCountSize;
						
					} else {
						logger.trace("Right-scan (from={}, to={}, left={}, right={}, diff={}): Found {} haplotype(s)", refCountIndex - 1, scanEndIndex, countL, countR, countDiff, haplotypes.length);
						
						refCountIndex = scanEndIndex + 1;
						countL = refCount[scanEndIndex];
					}
					
					
					haplotypeList.add(new RegionHaplotype(activeRegion, haplotypes));
					lastRegionEnd = scanEndIndex;
				}

				
			} else if (countDiff < diffThresholdL) {
				
				///////////////////////////////////////////////////
				///////////////////////////////////////////////////
				////                 Scan left                 ////
				///////////////////////////////////////////////////
				///////////////////////////////////////////////////
				
				logger.trace("Left-scan threshold at (from={}, to=LEFT_END, left={}, right={}, diff={})", refCountIndex, countL, countR, countDiff);
				
				// Check for a peak before scanning
				if (peakScanLength > 0) {
					recoveryValue = countL + diffThresholdR;
					scanEndIndex = refCountIndex + 1;
					lastScanIndex = refCountIndex + peakScanLength;
					
					if (lastScanIndex >= refCountSize)
						lastScanIndex = refCountSize;
					
					while (scanEndIndex < lastScanIndex) {
						
						if (refCount[scanEndIndex] <= recoveryValue &&  // Peak count lowers again
							refCount[refCountIndex] - refCount[scanEndIndex] < diffThresholdR)  // Right side of peak must not drop into a valley
						{ 
							
							logger.trace("Skipping peak from {} to {}", refCountIndex, scanEndIndex);
							
							// Peak detected, skip it
							countL = refCount[scanEndIndex];
							refCountIndex = scanEndIndex + 1;
							continue REF_SEARCH;
						}
						
						++scanEndIndex;
					}
				}
				
				// Stop if the length is longer than the end scan limit
				if (refCountIndex > scanLimit) {
					logger.trace("Left-scan threshold at (from={}, to=LEFT_END, left={}, right={}, diff={}): Discarded: Scan limit reached ({})", refCountIndex, countL, countR, countDiff, scanLimit);
					
					countL = countR;
					++refCountIndex;
					continue REF_SEARCH;
				}
				
				// Scan left
				scanEndIndex = refCountIndex - 1;
				
				// Reset peak stats
				nPeak = 0;
				peakScanIndex = 0;
				lastValleyIndex = 0;
				
				// Do scan
				SCAN_LOOP:
				while (true) {
					
					if (expDecayMin == 1.0) {
						// Constant threshold scan
						
						// Set recovery value
						recoveryValue = countR - diffThresholdR;
						
						if (recoveryValue < 1)
							recoveryValue = 1;
						
						// Set last scan index
						lastScanIndex = lastRegionEnd;
						
						logger.trace("Left scan: Constant search: index={}, threshold={}, value={}", refCountIndex, diffThresholdR, recoveryValue);
						
						while (scanEndIndex >= lastScanIndex && refCount[scanEndIndex] < recoveryValue)
							--scanEndIndex;
						
					} else {
						// Exponential-decay threshold scan
						
						// Exponential decay recovery value
						expDecayMinValue = countR * expDecayMin;
						
						if (expDecayMinValue < 1)
							expDecayMinValue = 1;
						
						expDecayRange = countR - expDecayMinValue;
						
						// Set last scan index
						lastScanIndex = lastRegionEnd;
						
						logger.trace("Left scan: Exponential search: index={}, range={}, min={}, lambda={}", refCountIndex, expDecayRange, expDecayMinValue, expDecayLambda);
						
						while (
								scanEndIndex >= lastScanIndex &&
								refCount[scanEndIndex] <
									expDecayRange * Math.exp((refCountIndex - scanEndIndex) * expDecayLambda) + expDecayMinValue  // range * e^{-(scanEndIndex - refCountIndex) * lambda} + min
							) {
							
							--scanEndIndex;
						}
						
						// Set recovery value for peak detection
						recoveryValue = (int) (expDecayRange * Math.exp((refCountIndex - scanEndIndex) * expDecayLambda) + expDecayMinValue);
					}
					
					// Setup peak detection (ensures a stable recovery)
					if (peakScanLength == 0)
						break SCAN_LOOP;
					
					peakScanIndex = scanEndIndex;
					peakScanLimit = scanEndIndex - peakScanLength;
					
					if (peakScanLimit < -1)
						peakScanLimit = -1;
					
					// Perform peak detection
					while (peakScanIndex > peakScanLimit) {
						
						if (refCount[peakScanIndex] > recoveryValue) {
							
							++nPeak;
							
							scanEndIndex = peakScanIndex;
							
							// Check for an excessive number of peaks
							if (nPeak > 3) {
								
								if ((refCountIndex - scanEndIndex) / nPeak < kSize) {
									
									// Excessive peaks
									
									logger.trace("Left-scan (from={}, to={}, left={}, right={}, diff={}): Peaks too frequent (n={}, interavg={})", refCountIndex - 1, scanEndIndex, countL, countR, countDiff, nPeak, (refCountIndex - scanEndIndex) / nPeak);
									
									countL = countR;
									++refCountIndex;
									continue REF_SEARCH;
								}
							}
							
							continue SCAN_LOOP;
						}
						
						--peakScanIndex;
					}
					
					// No peak found, end scan
					break SCAN_LOOP;
				}
				
				// Scan must reach the left end of the reference
				if (scanEndIndex > 0) {
					logger.trace("Left-scan threshold at (from={}, to={}, left={}, right={}, diff={}): Discarded: Scan does not reach left end", refCountIndex, scanEndIndex, countL, countR, countDiff);
					
					countL = countR;
					++refCountIndex;
					continue REF_SEARCH;
				}
				
				// Recover by searching for a k-mer count change
				if (recoverRightAnchor && refCountIndex - scanLimit < 0) {  // TODO: Change recoverRightAnchor to a var that controls right and left
					// Attempt to find the most likely left anchor even though the k-mer count did not recover
					
					scanEndIndex = refCountIndex - kSize;
					
					if (scanEndIndex < 0)
						scanEndIndex = 0;
					
					RECOVERY_SCAN:
					while (scanEndIndex > 0) {
						
						if (refCount[scanEndIndex - 1] - refCount[scanEndIndex] > diffThresholdR) {
							logger.trace("Left-scan (from={}, to=LEFT_END, left={}, right={}, diff={}): Recovery scan found k-mer increase at {}", refCountIndex - 1, countL, countR, countDiff, scanEndIndex);
							break RECOVERY_SCAN;
						}
						
						--scanEndIndex;
					}
					
					if (scanEndIndex == 0)
						scanEndIndex = -1;
					
				} else {
					scanEndIndex = -1;
				}
				
				// Stop if variants are already called to the left
				if (scanEndIndex < lastRegionEnd && lastRegionEnd > 0) {
					
					if (scanEndIndex == -1)
						logger.trace("Left-scan threshold at (from={}, to=LEFT_END, left={}, right={}, diff={}): Discarded: Cannot overrun active regions to the left (last region ends at {})", refCountIndex, countL, countR, countDiff, lastRegionEnd);
					else
						logger.trace("Left-scan threshold at (from={}, to={}, left={}, right={}, diff={}): Discarded: Cannot overrun active regions to the left (last region ends at {})", refCountIndex, scanEndIndex, countL, countR, countDiff, lastRegionEnd);
					
					countL = countR;
					++refCountIndex;
					continue REF_SEARCH;
				}
				
				// Stop if scan goes to the left end and end calling is disabled
				if (scanEndIndex < 0 && anchorBothEnds) {
					logger.trace("Left-scan threshold at (from={}, to=LEFT_END, left={}, right={}, diff={}): Discarded: Anchor-both-ends is set", refCountIndex, countL, countR, countDiff);
					
					countL = countR;
					++refCountIndex;
					continue REF_SEARCH;
				}
				
				// Stop if the region contains ambiguous bases and calling ambiguous regions is not allowed
				if (noAmbiguousBases && refRegion.containsAmbiguousByIndex(scanEndIndex, refCountIndex)) {
					
					if (scanEndIndex == 0)
						logger.trace("Left-scan threshold at (from={}, to=LEFT_END, left={}, right={}, diff={}): Discarded: Ambiguous regions are disabled", refCountIndex, countL, countR, countDiff);
					else
						logger.trace("Left-scan threshold at (from={}, to={}, left={}, right={}, diff={}): Discarded: Ambiguous regions are disabled", refCountIndex, scanEndIndex, countL, countR, countDiff);
					
					countL = countR;
					++refCountIndex;
					continue REF_SEARCH;
				}
				
				// Create active region
				activeRegion = new ActiveRegion(refRegion, scanEndIndex, refCountIndex, refCount, kUtil);
				
				// Get haplotypes
				haplotypes = alignmentBuilder.getHaplotypes(activeRegion);
				
				if (haplotypes.length == 0 || (haplotypes.length == 1 && haplotypes[0].isWildtype())) {
					
					// Log
					if (haplotypes.length == 0)
						if (scanEndIndex >= 0)
							logger.trace("Left-scan threshold at (from={}, to=LEFT_END, left={}, right={}, diff={}): Discarded: No haplotypes", refCountIndex, countL, countR, countDiff);
						else
							logger.trace("Left-scan threshold at (from={}, to={}, left={}, right={}, diff={}): Discarded: No haplotypes", refCountIndex, scanEndIndex, countL, countR, countDiff);
					else
						if (scanEndIndex >= 0)
							logger.trace("Left-scan threshold at (from={}, to=LEFT_END, left={}, right={}, diff={}): Discarded: Found only the reference haplotype", refCountIndex, countL, countR, countDiff);
						else
							logger.trace("Left-scan threshold at (from={}, to={}, left={}, right={}, diff={}): Discarded: Found only the reference haplotype", refCountIndex, scanEndIndex, countL, countR, countDiff);
										
					countL = countR;
					++refCountIndex;
					continue REF_SEARCH;
				}
				
				if (scanEndIndex >= 0)
					logger.trace("Left-scan threshold at (from={}, to={}, left={}, right={}, diff={}): Found {} haplotype(s)", refCountIndex, scanEndIndex, countL, countR, countDiff, haplotypes.length);
				else
					logger.trace("Left-scan threshold at (from={}, to=LEFT_END, left={}, right={}, diff={}): Found {} haplotype(s)", refCountIndex, countL, countR, countDiff, haplotypes.length);
				
				haplotypeList.add(new RegionHaplotype(activeRegion, haplotypes));
				
				countL = countR;
				++refCountIndex;
				lastRegionEnd = refCountIndex;
				
			} else {
				// No scan, move to the next base
				countL = countR;
				++refCountIndex;
			}
		}
		
		// Add wildtype active regions
		if (emitWildtypeActiveRegions) {
			
			RegionHaplotype[] haploArray;  // Array of haplotype regions
			int haploArrayIndex;           // Current element in halpoArray
			int nextRegionIndex;                // Index of the next region
			
			// Get a sorted array of haplotype regions
			haploArray = haplotypeList.toArray(new RegionHaplotype[0]);
			
			Arrays.sort(haploArray);
			
			haplotypeList.clear();
			
			// Initialize traverse
			if (haploArray.length > 0)
				nextRegionIndex = haploArray[0].activeRegion.startIndex;
			else
				nextRegionIndex = refCountSize;
			
			refCountIndex = 0;
			haploArrayIndex = 0;
			
			// Traverse to end of k-mer array
			while (nextRegionIndex < refCountSize) {
				
				if (refCountIndex == nextRegionIndex) {
					// Active region with variants
					
					// Add
					haplotypeList.add(haploArray[haploArrayIndex]);
					
					// Increment indices
					if (haploArray[haploArrayIndex].activeRegion.rightEnd)
						refCountIndex = nextRegionIndex = refCountSize;  // Stop traversal
					else
						refCountIndex = haploArray[haploArrayIndex].activeRegion.endKmerIndex;
					
					++haploArrayIndex;
					
					// Check for the end of the haplotype array
					if (haploArrayIndex < haploArray.length)
						nextRegionIndex = haploArray[haploArrayIndex].activeRegion.startKmerIndex;
					else
						nextRegionIndex = refCountSize;
					
				} else {
					// Wildtype region
					
					// Get index for the end of this active region
					if (nextRegionIndex < refCountSize)
						scanEndIndex = nextRegionIndex;
					else
						scanEndIndex = refCountSize - 1;
					
					// Create active region
					activeRegion = new ActiveRegion(refRegion, refCountIndex, scanEndIndex, refCount, kUtil);
					
					// Get wildtype haplotype
					haplotypes = alignmentBuilder.getHaplotypes(activeRegion);
					
					// Add active region
					haplotypeList.add(new RegionHaplotype(activeRegion, haplotypes));
					
					// Increment index
					refCountIndex = nextRegionIndex;
				}
			}
		}
		
		// Return active regions
		return new ActiveRegionContainer(refRegion, haplotypeList.toArray(new RegionHaplotype[0]), refCount);
	}
	
	/**
	 * When set to <code>true</code>, active regions must be anchored on both ends
	 * by a k-mer that is presumably not altered. This prevents variant calls on the
	 * ends of the reference regions, but ensures that variants are called with the
	 * best possible evidence to support them. If set to <code>false</code>, variant
	 * regions may extend to one end of the reference, but must be anchored on at
	 * least one side by a presumably unaltered k-mer.
	 * 
	 * @param anchorBothEnds Set to <code>true</code> to anchor all active regions
	 *   on both sides by an unaltered k-mer or <code>false</code> to allow an active
	 *   region to extend to one end of the reference.
	 * 
	 * @see #DEFAULT_ANCHOR_BOTH_ENDS
	 */
	public void setAnchorBothEnds(boolean anchorBothEnds) {
		this.anchorBothEnds = anchorBothEnds;
		
		return;
	}
	
	/**
	 * Get the both-end-anchor property.
	 * 
	 * @return <code>true</code> if active regions are anchored on both ends by an
	 *   unaltered k-mer.
	 * 
	 * @see #setAnchorBothEnds(boolean)
	 */
	public boolean getAnchorBothEnds() {
		return anchorBothEnds;
	}
	
	/**
	 * While scanning for variants, use each reference k-mer and its reverse complement
	 * if <code>true</code>. If the sequence data is in one orientation, then this
	 * should be set to false.
	 * 
	 * @param countReverseKmers Count the reverse complement k-mers in the active region
	 *   with the forward k-mers if <code>true</code>. 
	 * 
	 * @see #DEFAULT_COUNT_REVERSE_KMERS
	 */
	public void setCountReverseKmers(boolean countReverseKmers) {
		this.countReverseKmers = countReverseKmers;
		
		initAlignmentBuilder();
		
		return;
	}
	
	/**
	 * Determine if the reverse complement k-mers are counted with the forward k-mers
	 * in the active region.
	 * 
	 * @return <code>true</code> if reverse complement k-mers are counted.
	 * 
	 * @see #setCountReverseKmers(boolean)
	 */
	public boolean getCountReverseKmers() {
		return countReverseKmers;
	}
	
	/**
	 * Set the minimum difference between neighboring k-mers to trigger a scan for an
	 * active region.
	 *  
	 * @param minimumDifference Minimum difference between two k-mers to trigger a scan
	 *   for an active region. 
	 * 
	 * @throws IllegalArgumentException If <code>minimumDifference</code> is less than
	 *   <code>1</code>.
	 * 
	 * @see #DEFAULT_MINIMUM_DIFFERENCE
	 */
	public void setMinimumDifference(int minimumDifference)
			throws IllegalArgumentException {
		
		if (minimumDifference < 1)
			throw new IllegalArgumentException("Cannotset minimum difference (to trigger a correction attempt) to a value less than 1: " + minimumDifference);
		
		this.minimumDifference = minimumDifference;
		
		return;
	}
	
	/**
	 * Get the minimum difference between neighboring k-mers to trigger a scan for an
	 * active region.
	 * 
	 * @return Minimum difference between neighboring k-mers to trigger a scan for an
	 *   active region.
	 * 
	 * @see #setMinimumDifference(int)
	 */
	public int getMinimumDifference() {
		return minimumDifference;
	}
	
	/**
	 * If this is set to a value greater than <code>0.0</code>, then the difference between
	 * neighboring k-mers that triggers a scan for an active region is determined by the
	 * distribution of k-mer differences over the reference for each sample. This is
	 * determined by taking the absolute value of the difference between each neighboring
	 * k-mer count and computing the quantile, <code>differenceQuantile</code>, over that
	 * distribution. If the computed quantile is less than the minimum difference, then the
	 * minimum difference is used. If reverse complement k-mers are counted, then this quantile
	 * is computed with the reverse complement k-mer count.
	 * <p/>
	 * Read depth varies significantly in sequence data, and setting a hard limit can work well
	 * for one reference and fail for another. This method of quantiles makes choosing the limit
	 * dynamic based on the expected distribution of the data.
	 * 
	 * @param differenceQuantile Quantile computed over the distribution of count differences
	 *   between neighbors.
	 * 
	 * @throws IllegalArgumentException If <code>differenceQuantile</code> is negative or greater
	 *   than <code>1.0</code>.
	 * 
	 * @see #DEFAULT_DIFFERENCE_QUANTILE
	 * @see #setMinimumDifference(int)
	 */
	public void setDifferenceQuantile(double differenceQuantile)
			throws IllegalArgumentException {
		
		if (differenceQuantile < 0.0 || differenceQuantile >= 1.0)
			throw new IllegalArgumentException("Cannot set minimum difference quantile (dynamic difference per reference to trigger a correction attempt): Quantile must be a positive number less than 1.0 (or 0.0 to disable): " + differenceQuantile);
		
		this.differenceQuantile = differenceQuantile;
		
		return;
	}
	
	/**
	 * Get the difference quantile or <code>0.0</code> if computing based on quantile is disabled.
	 * 
	 * @return Difference quantile.
	 * 
	 * @see #setDifferenceQuantile(double)
	 */
	public double getDifferenceQuantile() {
		return differenceQuantile;
	}
	
	/**
	 * Set the exponential decay minimum. This is the minimum value (asymptotic lower bound)
	 * of the exponential decay function as a proportion of the anchor k-mer count. If this
	 * value is <code>0.0</code>, k-mer count recovery threshold may decline to <code>1</code>.
	 * If this value is <code>1.0</code>, the decay function is not used and the detector falls
	 * back to finding a k-mer with a count within the difference threshold of the anchor k-mer
	 * count.
	 * 
	 * @param expDecayMin Exponential decay minimum. Must be between <code>0.0</code>
	 *   and <code>1.0</code> (inclusive).
	 * 
	 * @throws IllegalArgumentException If <code>expDecayMin</code> is less than <code>0.0</code>
	 *   or greater than <code>1.0</code>.
	 * 
	 * @see #DEFAULT_EXP_MIN
	 */
	public void setDecayMinimum(double expDecayMin)
			throws IllegalArgumentException {
		
		if (expDecayMin < 0.0 || expDecayMin > 1.0)
			throw new IllegalArgumentException("Exponential decay minimum must not be negative or greater than 1.0: " + expDecayMin);
		
		this.expDecayMin = expDecayMin;
		
		return;
	}
	
	/**
	 * Get the exponential decay minimum as a proportion of the anchor k-mer count.
	 * 
	 * @return Exponential decay minimum.
	 * 
	 * @see #setDecayMinimum(double)
	 */
	public double getDecayMinimum() {
		return expDecayMin;
	}
	
	/**
	 * Set the exponential decay alpha, which controls how quickly the recovery threshold
	 * declines. Alpha is defined as the proportion of decay from the lower bound to the
	 * minimum count at <code>k</code> k-mers from the anchor. A lower value causes a steeper
	 * decline toward the minimum threshold value.
	 * 
	 * @param expDecayAlpha Exponential decay alpha. Must be between <code>0.0</code>
	 *   and <code>1.0</code> (exclusive).
	 * 
	 * @throws IllegalArgumentException If <code>expDecayMin</code> is not between <code>0.0</code>
	 *   and <code>1.0</code> (exclusive).
	 * 
	 * @see #DEFAULT_EXP_ALPHA
	 */
	public void setDecayAlpha(double expDecayAlpha)
			throws IllegalArgumentException {
		
		if (expDecayAlpha <= 0.0 || expDecayAlpha >= 1.0)
			throw new IllegalArgumentException("Exponential decay alpha must be between 0.0 and 1.0 (exclusive): " + expDecayAlpha);
		
		this.expDecayAlpha = expDecayAlpha;
		
		expDecayLambda = -Math.log(expDecayAlpha) / kUtil.kSize;
		
		return;
	}
	
	/**
	 * Get the exponential decay alpha.
	 * 
	 * @return Exponential decay alpha.
	 * 
	 * @see #setDecayAlpha(double)
	 */
	public double getDecayAlpha() {
		return expDecayAlpha;
	}
	
	/**
	 * Get the exponential decay lambda calculated from the exponential decay alpha
	 * and the k-mer size. This is the lambda of the exponential decay function. 
	 * <p/>
	 * <code>lambda</code> := <code>-log(alpha) / k</code>
	 * 
	 * @return The lambda of the exponential decay function.
	 */
	public double getDecayLambda() {
		return expDecayLambda;
	}
	
	/**
	 * If set to <code>true</code>, allow active regions to include ambiguous bases.
	 * 
	 * @param callAmbiguousRegions Allow ambiguous bases in active regions if
	 *   <code>true</code>.
	 * 
	 * @see #DEFAULT_CALL_AMBIGUOUS_REGIONS
	 */
	public void setCallAmbiguousRegions(boolean callAmbiguousRegions) {
		this.callAmbiguousRegions = callAmbiguousRegions;
		
		return;
	}
	
	/**
	 * Get the property allow ambiguous bases in active regions.
	 *  
	 * @return <code>true</code> if active regions may contain ambiguous bases.
	 * 
	 * @see #setCallAmbiguousRegions(boolean)
	 */
	public boolean getCallAmbiguousRegions() {
		return callAmbiguousRegions;
	}
	
	/**
	 * Set the weight vector used for alignments. If a vector is never set, the
	 * default alignment weight vector provided by <code>AlignmentWeight.get()</code>
	 * is used.
	 * 
	 * @param alignmentWeight Alignment weight vector.
	 * 
	 * @throws NullPointerException If <code>alignmentWeight</code> is
	 *   <code>null</code>.
	 */
	public void setAlignmentWeight(AlignmentWeight alignmentWeight)
			throws NullPointerException {
		
		if (alignmentWeight == null)
			throw new NullPointerException("Alignment weight may not be null");
		
		this.alignmentWeight = alignmentWeight;
		
		// Update state and values computed off of this alignment weight
		setScanLimitFactor(scanLimitFactor);
		initAlignmentBuilder();
		
		return;
	}
	
	/**
	 * Get the alignment weight vector used for alignments.
	 *  
	 * @return The alignment weight vector used for alignments.
	 * 
	 * @see #setAlignmentWeight(AlignmentWeight)
	 */
	public AlignmentWeight getAlignmentWeight() {
		return alignmentWeight;
	}
	
	/**
	 * Set the number of k-mers to scan when detecting and stepping over a suspected peak.
	 * A peak can occur when a k-mer matches another region of the genome and the k-mer counts
	 * from that region are added. If set to <code>0</code>, peak detection is disabled.
	 * 
	 * @param peakScanLength Number of k-mers to scan when detecting a peak or <code>0</code>
	 *   to disable peak detection.
	 * 
	 * @throws IllegalArgumentException If <code>peakScanLength</code> is negative.
	 * 
	 * @see #DEFAULT_PEAK_SCAN_LENGTH
	 */
	public void setPeakScanLength(int peakScanLength)
			throws IllegalArgumentException {
		
		if (peakScanLength < 0)
			throw new IllegalArgumentException("Peak scan length is negative: " + peakScanLength);
		
		this.peakScanLength = peakScanLength;
		
		return;
	}
	
	/**
	 * Get the number of k-mers to scan during peak detection.
	 * 
	 * @return The number of k-mers to scan during peak detection.
	 * 
	 * @see #setPeakScanLength(int)
	 */
	public int getPeakScanLength() {
		return peakScanLength;
	}
	
	/**
	 * Set the limit on how many k-mers a scan to over the reference is allowed to
	 * traverse. The maximum number of bases is determined by multiplying this factor by
	 * the k-mer size and adding the maximum length of a gap. The computed limit will
	 * adjusted so that it will never be less than the k-mer size or greater than
	 * <code>Integer.MAX_VALUE</code>.
	 * 
	 * @param scanLimitFactor The scan limit factor.
	 * 
	 * @throws IllegalArgumentException If <code>scanLimitFactor</code> is negative.
	 */
	public void setScanLimitFactor(double scanLimitFactor)
			throws IllegalArgumentException {
		
		int maxGapSize = alignmentWeight.getMaxExclusiveGapSize(kSize);
		double scanLimitDb;  // Scan limit as a double
		
		// Check arguments
		if (scanLimitFactor < 0.0)
			throw new IllegalArgumentException("Scan limit factor must not be negative: " + scanLimitFactor);
		
		// Set limit
		this.scanLimitFactor = scanLimitFactor;
		
		// Safe set (avoid integer overflow)
		scanLimitDb = maxGapSize + scanLimitFactor * kSize;
		
		if (scanLimitDb > Integer.MAX_VALUE)
			scanLimit = Integer.MAX_VALUE;
		else
			scanLimit = (int) scanLimitDb;
		
		// Cannot be less than the k-mer size
		if (scanLimit < kSize)
			scanLimit = kSize;
		
		return;
	}
	
	/**
	 * Get the scan limit.
	 * 
	 * @return Scan limit.
	 * 
	 * @see #setScanLimitFactor(double)
	 */
	public double getScanLimitFactor() {
		return scanLimitFactor;
	}
	
	/**
	 * Get the end-scan limit length computed by multiplying the end-scan limit by
	 * the maximum gap length (assuming no aligned bases).
	 * 
	 * @return End-scan limit length.
	 * 
	 * @see #setScanLimitFactor(double)
	 * @see #getScanLimitFactor()
	 */
	public int getScanLimitLength() {
		return scanLimit;
	}
	
	/**
	 * Attache a trace matrix to each haplotype if set to <code>true</code>. Enabling this option
	 * will have significant performance and memory consumption implications, and it should only
	 * be used for targeted debugging efforts.
	 * 
	 * @param traceHaplotypeAlignment Haplotype trace option.
	 */
	public void setTraceHaplotypeAlignment(boolean traceHaplotypeAlignment) {
		this.traceHaplotypeAlignment = traceHaplotypeAlignment;
		
		initAlignmentBuilder();
		
		return;
	}
	
	/**
	 * Get the haplotype trace option.
	 * 
	 * @return Haplotype trace option.
	 * 
	 * @see #setTraceHaplotypeAlignment(boolean)
	 */
	public boolean getTraceHaplotypeAlignment() {
		return traceHaplotypeAlignment;
	}
	
	/**
	 * Set the maximum number of haplotypes the aligner may save. The least-likely states
	 * are trimmed when this threshold is reached.
	 * 
	 * @param maxAlignerState Maximum number of states.
	 * 
	 * @throws IllegalArgumentException If <code>maxAlignerState</code> is less than <code>1</code>.
	 * 
	 * @see KmerAligner#DEFAULT_MAX_STATE
	 */
	public void setMaxAlignerState(int maxAlignerState)
		throws IllegalArgumentException {
		
		alignmentBuilder.setMaxState(maxAlignerState); // Throws IllegalArgumentException
		
		return;
	}
	
	/**
	 * Get the maximum number of states that may be saved.
	 * 
	 * @return The maximum number of states that may be saved.
	 * 
	 * @see #setMaxAlignerState(int)
	 */
	public int getMaxAlignerState() {
		return alignmentBuilder.getMaxState();
	}
	
	/**
	 * Set the maximum number of haplotypes to be saved. The least-likely haplotypes are trimmed when
	 * this threshold is reached.
	 * 
	 * @param maxHaplotypes Maximum number of haplotypes.
	 * 
	 * @throws IllegalArgumentException If <code>maxHaplotypes</code> is less than <code>1</code>.
	 * 
	 * @see KmerAligner#DEFAULT_MAX_HAPLOTYPES
	 */
	public void setMaxHaplotypes(int maxHaplotypes)
		throws IllegalArgumentException {
		
		alignmentBuilder.setMaxHaplotypes(maxHaplotypes);
		
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
		return alignmentBuilder.getMaxHaplotypes();
	}
	
	/**
	 * Get k-mer counts from the reference sequence and populate the tree of ambiguous regions.
	 * 
	 * @param refSequence Sequence to get k-mer counts from.
	 * 
	 * @return An array of k-mer counts of an empty array if the reference region is not long
	 *   enough to contain a k-mer.
	 */
	private int[] getCounts(ReferenceRegion refSequence) {
		
		// Declarations
		int[] count;      // K-mer counts from this reference region
		int countIndex;   // Location in counts where the next k-mer count is written
		int countLength;  // Length of counts
		
		CountMap counter;  // K-mer counter
		
		byte[] sequence;  // Sequence (refSequence.sequence)
		int seqIndex;     // Current location in sequence
		int seqLength;    // Size of sequence (refSequence.size)
		
		KmerUtil kUtil;  // Kmer utility with no minimizer characteristics
		int kSize;       // K-mer size
		int[] fwdKmer;   // Forward k-mer
		int[] revKmer;   // Reverse-complement of fwdKmer
		
		boolean countReverseKmers;  // Count reverse-complement of k-mers if true
		int lastBaseLoad;           // seqIndex when fwdKmer (and revKmer) is fully loaded 
		Base base;                  // Current nucleotide base
		
		// Check arguments
		assert (refSequence != null) :
			"refSequence is null";
		
		// Get k-mer size and utility object
		kUtil =this.kUtil;
		kSize = kUtil.kSize;
		
		// Check length
		countLength = refSequence.size - kSize + 1;
		
		if (countLength <= 0)
			return new int[0];
		
		// Init
		fwdKmer = new int[kUtil.kmerArraySize];
		revKmer = new int[kUtil.kmerArraySize];
		
		count = new int[countLength];
		countIndex = 0;
		
		counter = this.counter;
		
		sequence = refSequence.sequence;
		seqIndex = 0;
		seqLength = refSequence.size;
		
		countReverseKmers = this.countReverseKmers;
		
		// Read bases and get k-mer counts
		REF_LOOP:
		while (seqIndex < seqLength) {
			
			// Get load limit
			lastBaseLoad = seqIndex + kSize - 1;
			
			// Stop if no more k-mers can be built from this sequence
			if (lastBaseLoad > seqLength) {
				
				// Write 0 to end of count array
				while (countIndex < countLength)
					count[countIndex++] = 0;
				
				break REF_LOOP;
			}
			
			// Load k-mers and get counts
			while (seqIndex < lastBaseLoad) {
				base = BYTE_TO_BASE[sequence[seqIndex]];
				
				// Handle ambiguous bases
				if (base == null) {  // Ambiguous base
					
					while (countIndex <= seqIndex)
						count[countIndex++] = 0;
					
					++seqIndex;
					
					continue REF_LOOP;
				}
				
				// Add to k-mers
				kUtil.append(fwdKmer, base);
				kUtil.prepend(revKmer, BASE_COMPLEMENT[base.intVal]);
				
				// Advance load indices
				++seqIndex;
			}
			
			// Get k-mer counts
			while (seqIndex < seqLength) {
				
				base = BYTE_TO_BASE[sequence[seqIndex]];
				
				if (base == null)
					continue REF_LOOP;  // Let the k-mer load loop deal with ambiguous bases
				
				kUtil.append(fwdKmer, base);
				kUtil.prepend(revKmer, BASE_COMPLEMENT[base.intVal]);
				
				count[countIndex++] = counter.get(fwdKmer) + (countReverseKmers ? counter.get(revKmer) : 0);
				
				++seqIndex;
			}
		}
		
		// Verify all counts were loaded
		assert (countIndex == count.length) :
			String.format("countIndex (%d) != count.length (%d)", countIndex, count.length);
		
		// Return count array
		return count;
	}
	
	/**
	 * Get the difference threshold by computing a quantile over the set of
	 * differences between each neighboring k-mer.
	 * 
	 * @param count Count array to compute quantiles from.
	 * @param start Start index of the count range to analyze (inclusive).
	 * @param end End index of the count range to analyze (exclusive).
	 * 
	 * @return An array with the difference and recovery thresholds.
	 * 
	 * @see #setDifferenceQuantile(double)
	 */
	private int getDifferenceThreshold(int[] count, int start, int end) {
		
		// Declarations
		int lastCount;
		int thisCount;
		int thisCountDiff;
		int loc;
		double offset;
		
		int[] countDiff;
		int nCount;       // Number of elements in count while getting differences, and one less than
		                  // the length of countDiff while computing quantiles
		
		int diffQuantileValue;  // Difference threshold quantile value
		
		// Check arguments
		assert (count != null) :
			"Count array is null";
		
		nCount = end - start;
		
		assert (nCount > 2) :
			"Count array length is less than 3: " + nCount;
		
		nCount -= 1;  // Convert to the number of count differences (1 less than the number of counts) 
		countDiff = new int[nCount];
		
		// Convert to count differences
		lastCount = count[start];
		
		for (int index = 0; index < nCount; ++index) {
			thisCount = count[start + index];
			thisCountDiff = lastCount - thisCount;
			
			countDiff[index] = (thisCountDiff < 0) ? -thisCountDiff : thisCountDiff;  // Absolute value of diff
			
			lastCount = thisCount;
		}
		
		nCount -= 1;  // One less difference than k-mer count differences (simplifies quantile calculations)
		
		// Sort
		Arrays.sort(countDiff);
		
		// Calculate difference quantile
		if (differenceQuantile > 0.0) {
			loc = (int) (nCount * differenceQuantile);
			offset = (nCount * differenceQuantile) - loc;
			
			diffQuantileValue = (int) (countDiff[loc] * (1 - offset) + countDiff[loc + 1] * offset);
			
			if (diffQuantileValue < minimumDifference)
				diffQuantileValue = minimumDifference;
			
		} else {
			// Use minimum difference if disabled
			diffQuantileValue = minimumDifference;
		}
		
		// Return quantile
		return diffQuantileValue;
	}
	
	/**
	 * Set the right anchor recovery property. If the k-mer count does not recover to a value close
	 * to the pre-variant value, then search for a location where the count increases abruptly and
	 * attempt to create a haplotype.
	 * 
	 * @param recoverRightAnchor Right anchor recovery property.
	 * 
	 * @see #DEFAULT_RECOVER_RIGHT_ANCHOR
	 */
	public void setRecoverRightAnchor(boolean recoverRightAnchor) {
		this.recoverRightAnchor = recoverRightAnchor;
	}
	
	/**
	 * Get the right anchor recovery property.
	 * 
	 * @return Right anchor recovery property.
	 * 
	 * @see #setRecoverRightAnchor(boolean)
	 */
	public boolean getRecoverRightAnchor() {
		return recoverRightAnchor;
	}
	
	/**
	 * Set the property to emit active region objects for wildtype regions.
	 * 
	 * @param emitWildtypeActiveRegions <code>true</code> if wildtype regions should be
	 *   emitted.
	 */
	public void setEmitWildtypeActiveRegions(boolean emitWildtypeActiveRegions) {
		this.emitWildtypeActiveRegions = emitWildtypeActiveRegions;
		
		return;
	}
	
	/**
	 * Get the property to emit active region objects for wildtype regions.
	 * 
	 * @return Emit active region property.
	 */
	public boolean getEmitWildtypeActiveRegions() {
		return emitWildtypeActiveRegions;
	}
	
	/**
	 * Initialize the alignment builder. This must be called whenever <code>alignmentWeight</code>,
	 * <code>countReverseKmers</code>, or <code>traceHaplotypeAlignment</code> in changed.
	 */
	private void initAlignmentBuilder() {
		alignmentBuilder = new KmerAlignmentBuilder(kUtil, counter, alignmentWeight, countReverseKmers, traceHaplotypeAlignment);
		
		return;
	}
}
