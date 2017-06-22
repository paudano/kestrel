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

package edu.gatech.kestrel.refreader;

import java.io.IOException;
import java.security.MessageDigest;
import java.security.NoSuchAlgorithmException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import edu.gatech.kanalyze.batch.BatchCache;
import edu.gatech.kanalyze.batch.SequenceBatch;
import edu.gatech.kanalyze.util.BoundedQueue;
import edu.gatech.kanalyze.util.SequenceNameTable;
import edu.gatech.kestrel.interval.RegionInterval;
import edu.gatech.kestrel.runner.KestrelRunnerBase;
import edu.gatech.kestrel.util.digest.Digest;
import edu.gatech.kestrel.util.digest.NullMessageDigest;

/**
 * Collects the output from a KAnalyze sequence reader and transforms them into a set
 * of <code>ReferenceSequence</code> objects. 
 */
public class ReadCollectorRunner implements Runnable {
	
	/** Logger. */
	private final Logger logger = LoggerFactory.getLogger(ReadCollectorRunner.class);
	
	/**
	 * Contains a set of intervals to be extracted from the reference sequence, or is <code>null</code>
	 * if variants are called over the reference sequences.
	 */
	private Map<String, RegionInterval[]> intervalMap;
	
	/** Container of reference regions produced by this reader. */
	private final ReferenceRegionContainer refRegionContainer;
	
	/** Queue of sequence reads. */
	private final BoundedQueue<SequenceBatch> sequenceQueue;
	
	/** Cache for recycling sequence batches. */
	private final BatchCache<SequenceBatch> sequenceBatchCache;
	
	/** Table of sequence names written by the sequence reader. */
	private final SequenceNameTable nameTable;
	
	/** The length of flanks to add to intervals when extracting regions of a reference sequence. */
	private final int flankLength;
	
	/** <code>true</code> to extract regions on intervals. */
	private boolean doInterval;
	
	/** <code>true</code> unless this collector is stopped. */
	private boolean isActive;
	
	/** Engine for computing digests. */
	private final MessageDigest digestEngine;
	
	/** Algorithm used in <code>digestEngine</code>. */
	private final String digestAlgorithm;
	
	
	
	// For processing intervals
	
	/**
	 * Name of the reference sequence currently loaded into <code>refIntervals</code> or
	 * <code>null</code> if no sequence has been loaded.
	 */
	private String refSequenceName;
	
	/** Sequence source name. */
	private String sequenceSourceName;
	
	/** Intervals for the current reference sequence. */
	private RegionInterval[] refInterval;
	
	/** Interval currently being read. */
	private RegionInterval currentInterval;
	
	/** Interval currently being read. */
	private int refIntervalIndex;
	
	/**
	 * When sequence length reaches this value, the current region is done. This value has the left flank length
	 * added to it.
	 */
	private int intervalEndSeqLength;
	
	
	/** A list of incomplete reference regions. */
	private List<ReferenceRegion.IncompleteRegion> incompleteRegionList;
	
	/** Length of the left flank. */
	private int leftFlankLength;
	
	/** Size of the reference sequence. */
	private int refSequenceSize;
	
	/** <code>true</code> if the sequence description should be removed from a sequence name. */
	private boolean removeSequenceDescription;
	
	/** <code>true</code> if regions on the negative strand should be reverse complemented before variant calling. */
	private boolean reverseComplementNegativeStrand;
	
	
	// Constants
	
	/** Name of the cryptographic hash algorithm applied to sequences. */
	public static final String DIGEST_ALGORITHM = "MD5";
	
	/** Default remove sequence description property. */
	public static final boolean DEFAULT_REMOVE_SEQUENCE_DESCRIPTION = true;
	
	/** Default reverse-complement negative strand regions. */
	public static final boolean DEFAULT_REVERSE_COMPLEMENT_NEGATIVE_STRAND = false;
	
	/** An empty set of intervals used by the sequence reader when a sequence has no intervals. */
	private static final RegionInterval[] EMPTY_INTERVAL = new RegionInterval[0];
	
	
	
	/**
	 * Create a read collector.
	 * 
	 * @param refRegionContainer Reference region container.
	 * @param nameTable Table of sequence and sequence source names.
	 * @param sequenceQueue Sequence queue.
	 * @param sequenceBatchCache Sequence batch cache.
	 * @param intervalMap Map of intervals for each reference sequence. If <code>null</code>, all
	 *   reference sequences are extracted and become the reference regions variants are called on
	 *   (i.e. variants are called on entire reference sequences rather than intervals within those
	 *   reference sequences).
	 * @param flankLength Length of flanks.
	 * 
	 * @throws NullPointerException If any argument other than <code>intervalMap</code> is <code>null</code>.
	 * @throws IllegalArgumentException If <code>flankLength</code> is negative.
	 */
	public ReadCollectorRunner(ReferenceRegionContainer refRegionContainer, SequenceNameTable nameTable, BoundedQueue<SequenceBatch> sequenceQueue, BatchCache<SequenceBatch> sequenceBatchCache, Map<String, RegionInterval[]> intervalMap, int flankLength)
			throws NullPointerException, IllegalArgumentException {
		
		MessageDigest digestEngine;
		String digestAlgorithm;
		
		// Check arguments
		if (refRegionContainer == null)
			throw new NullPointerException("refRegionContainer is null");
		
		if (nameTable == null)
			throw new NullPointerException("Name table is null");
		
		if (sequenceQueue == null)
			throw new NullPointerException("Sequence queue is null");
		
		if (sequenceBatchCache == null)
			throw new NullPointerException("Sequence batch cache is null");
		
		if (flankLength < 0)
			throw new IllegalArgumentException("Flank length is null");
		
		if (intervalMap == null) {
			intervalMap = new HashMap<>();
			doInterval = false;
			
		} else {
			doInterval = true;
		}
		
		// Set fields
		this.refRegionContainer = refRegionContainer;
		this.nameTable = nameTable;
		this.sequenceQueue = sequenceQueue;
		this.sequenceBatchCache = sequenceBatchCache;
		
		this.intervalMap = intervalMap;
		this.flankLength = flankLength;
		
		// Set sequence processing state
		refInterval = EMPTY_INTERVAL;
		refIntervalIndex = 0;
		
		incompleteRegionList = new ArrayList<>();
		
		leftFlankLength = 0;
		refSequenceSize = 0;
		
		refSequenceName = null;
		sequenceSourceName = null;
		refSequenceSize = 0;
		
		isActive = true;
		
		removeSequenceDescription = DEFAULT_REMOVE_SEQUENCE_DESCRIPTION;
		reverseComplementNegativeStrand = DEFAULT_REVERSE_COMPLEMENT_NEGATIVE_STRAND;
		
		// Set digest engine
		try {
			digestEngine = MessageDigest.getInstance(DIGEST_ALGORITHM);
			digestAlgorithm = DIGEST_ALGORITHM;
			
		} catch (NoSuchAlgorithmException ex) {
			logger.warn("Digest implementation could not be found for message digest algorithm: {}", DIGEST_ALGORITHM);
			
			digestEngine = new NullMessageDigest();
			digestAlgorithm = NullMessageDigest.ALGORITHM;
		}
		
		this.digestEngine = digestEngine;
		this.digestAlgorithm = digestAlgorithm;
		
		return;
	}
	
	/**
	 * Set the property to remove the sequence description from each sequence name. The
	 * sequence name is defined as everything up to the first whitespace character, and the
	 * description follows. When this option is enabled, the whitespace and description are
	 * removed.
	 * 
	 * @param removeSequenceDescription Remove sequence description parameter.
	 * 
	 * @see #DEFAULT_REMOVE_SEQUENCE_DESCRIPTION
	 */
	public void setRemoveDescription(boolean removeSequenceDescription) {
		this.removeSequenceDescription = removeSequenceDescription;
		
		return;
	}
	
	/**
	 * Get the property to remove sequence descriptions.
	 * 
	 * @return The remove sequence description property.
	 * 
	 * @see #setRemoveDescription(boolean)
	 */
	public boolean getRemoveDescription() {
		return removeSequenceDescription;
	}
	
	/**
	 * When set to <code>true</code>, the &quot;reverse complement negative strand&quot; property reverse
	 * complements any reference regions defined on a negative strand interval. The interval
	 * object has an attribute delineating strandedness, and this is usually set by reading
	 * a file with this information (e.g. BED file).
	 * 
	 * @param reverseComplementNegativeStrand Reverse complement negative strand regions if <code>true</code>.
	 *   Otherwise, regions are left in the orientation found in the reference sequence.
	 * 
	 * @see #DEFAULT_REVERSE_COMPLEMENT_NEGATIVE_STRAND
	 */
	public void setRevComplementNegStrand(boolean reverseComplementNegativeStrand) {
		this.reverseComplementNegativeStrand = reverseComplementNegativeStrand;
	}
	
	/**
	 * Get the property to reverse complement negative strand regions.
	 * 
	 * @return The &quot;Reverse complement negative strand&quot; property.
	 * 
	 * @see #setRevComplementNegStrand(boolean)
	 */
	public boolean getRevComplementNegStrand() {
		return reverseComplementNegativeStrand;
	}
	
	/**
	 * Run.
	 * 
	 * @throws ReadCollectorRunnerException If any error occurs.
	 */
	@Override
	public void run()
			throws ReadCollectorRunnerException {
		
		// Buffer
		byte[] seqBuffer;  // Sequence buffer
		int bufferIndex;   // Current position in seqBuffer
		
		// Batch
		SequenceBatch batch;  // Sequence batch bases are retrieved from
		int batchIndex;       // Current position in batch.sequence
		int seqStartIndex;    // Current position of the sequence start array in the current batch
		int seqStartPos;      // Position of the start of the next sequence
		int lastBatchIndex;   // Last index of the batch sequence in a copy operation
		
		// Sequence ID
		int sourceId;    // ID of the sequence source
		int sequenceId;  // ID of the sequence currently being read
		
		// Init
		seqBuffer = new byte[KestrelRunnerBase.DEFAULT_READER_SEQUENCE_BUFFER_SIZE];
		
		sourceId = 0;  // Flag no sequences read when 0
		sequenceId = 1;
		
		batchIndex = 0;
		bufferIndex = 0;
		
		incompleteRegionList = new ArrayList<>();
		
		// Process batches
		BATCH_LOOP:
		while (isActive) {
			
			// Get next batch
			try {
				batch = sequenceQueue.take();
				
			} catch (InterruptedException ex) {
				
				continue BATCH_LOOP;
			}
			
			if (batch == null)
				break BATCH_LOOP;
			
			// Reset batch index and sequence start position
			batchIndex = 0;
			
			seqStartIndex = 0;
			seqStartPos = batch.getSequenceStartLocation(seqStartIndex);
			
			// Check for sequence source ID change (or first batch)
			if (sourceId != batch.sourceId) {
				
				// Add last sequence if this is not the first batch
				if (sourceId != 0) {
					
					try {
						bufferIndex = flushInterval(seqBuffer, bufferIndex);
						
						flushIncompleteReferenceRegions();
						refSequenceSize = 0;
						
						assert (bufferIndex == 0) :
							"flushInterval() did not return 0 when switching sources";
						
					} catch (IllegalStateException ex) {
						err(ex.getMessage(), ex);
						return;
					}
				}
				
				bufferIndex = 0;
				refSequenceSize = 0;
				sourceId = batch.sourceId;
				sequenceId = 1;
				
				// Set sequence and source names
				refSequenceName = nameTable.getSequenceNameWithDefault(sourceId, sequenceId);
				sequenceSourceName = nameTable.getSourceNameWithDefault(sourceId);
				
				// Get the set of intervals for this reference
				setIntervalsForReference();
				
				logger.trace("Processing reference: sequence = {}, source = {}", refSequenceName, sequenceSourceName);
			}
			
			// Check: Batch sequence ID must not be less than sequenceId
			if (batch.sequenceId < sequenceId) {
				err(String.format("Sequence ID in batch is less than the last sequence ID: Batch (src=%d, seq=%d), last sequence ID = %d", batch.sourceId, batch.sequenceId, sequenceId), null);
				return;
			}
			
			// Add sequences until sequenceId matches batch.sequenceId
			while (sequenceId < batch.sequenceId) {
				
				try {
					bufferIndex = flushInterval(seqBuffer, bufferIndex);
					
					flushIncompleteReferenceRegions();
					
					assert (bufferIndex == 0) :
						"flushInterval() did not return 0 when switching sequences within a source";
					
				} catch (IllegalStateException ex) {
					err(ex.getMessage(), ex);
					return;
				}
				
				// Set next sequence
				++sequenceId;
				refSequenceName = nameTable.getSequenceNameWithDefault(sourceId, sequenceId);
				refSequenceSize = 0;
				bufferIndex = 0;
				
				// Get the set of intervals for this reference
				setIntervalsForReference();
				
				logger.trace("Processing reference: sequence = {}, source = {}", refSequenceName, sequenceSourceName);
			}
			
			// Read sequence from batch
			SEQUENCE_LOOP:
			while (batchIndex < batch.size) {
				
				// At sequence start within the batch (several sequences may be packed in the buffer)
				if (batchIndex == seqStartPos) {
					
					// Flush reference sequence sequence intervals
					try {
						bufferIndex = flushInterval(seqBuffer, bufferIndex);
						
						flushIncompleteReferenceRegions();
						
					} catch (IllegalStateException ex) {
						err(ex.getMessage(), ex);
						return;
					}
					
					// Set next sequence
					refSequenceSize = 0;
					++sequenceId;
					refSequenceName = nameTable.getSequenceNameWithDefault(sourceId, sequenceId);
					
					// Get the set of intervals for this reference
					setIntervalsForReference();
					
					// Set next start position
					++seqStartIndex;
					seqStartPos = batch.getSequenceStartLocation(seqStartIndex);
					
					continue SEQUENCE_LOOP;
				}
				
				// Check for the end of an interval or inter-interval space
				if (refSequenceSize == intervalEndSeqLength) {
					bufferIndex = flushInterval(seqBuffer, bufferIndex);
					
					continue SEQUENCE_LOOP;
				}
				
				// Find copy limit
				lastBatchIndex = batch.size;
				
				if (lastBatchIndex > seqStartPos)
					lastBatchIndex = seqStartPos;
				
				if (refSequenceSize + lastBatchIndex - batchIndex > intervalEndSeqLength)
					lastBatchIndex = batchIndex + intervalEndSeqLength - refSequenceSize;
				
				refSequenceSize += lastBatchIndex - batchIndex;
				digestEngine.update(batch.sequence, batchIndex, lastBatchIndex - batchIndex);
				
				// Copy if an interval is active
				if (currentInterval != null || ! doInterval) {
					
					// Expand seqBuffer if the remaining capacity is too small to hold this write operation
					if (lastBatchIndex - batchIndex > seqBuffer.length - bufferIndex) {
						
						try {
							seqBuffer = expandSeqBuffer(seqBuffer, lastBatchIndex - batchIndex);
							
						} catch (IllegalArgumentException ex) {
							err(ex.getMessage(), null);  // Just use the exception message and not the exception itself
							return;
						}
					}
					
					// Copy
					while (batchIndex < lastBatchIndex)
						seqBuffer[bufferIndex++] = batch.sequence[batchIndex++];
					
				} else {
					batchIndex = lastBatchIndex;
				}
				
				// Save reference intervals
				while (refSequenceSize > intervalEndSeqLength) {
					
					try {
						bufferIndex = flushInterval(seqBuffer, bufferIndex);
						
					} catch (IllegalStateException ex) {
						err(ex.getMessage(), ex);
						return;
					}
				}
			}
			
			// Account for sequence starts at the end of this batch
			while (batchIndex == seqStartPos) {
				
				// Flush reference sequence sequence intervals
				try {
					bufferIndex = flushInterval(seqBuffer, bufferIndex);
					
					flushIncompleteReferenceRegions();
					
				} catch (IllegalStateException ex) {
					err(ex.getMessage(), ex);
					return;
				}
				
				// Set next sequence
				refSequenceSize = 0;
				++sequenceId;
				refSequenceName = nameTable.getSequenceNameWithDefault(sourceId, sequenceId);
				
				// Get the set of intervals for this reference
				setIntervalsForReference();
				
				// Set next start position
				++seqStartIndex;
				seqStartPos = batch.getSequenceStartLocation(seqStartIndex);
			}
			
			// Return batch to cache
			sequenceBatchCache.add(batch);
		}
		
		// Add last sequence
		if (sourceId > 0) {
			// Flush reference sequence sequence intervals
			try {
				bufferIndex = flushInterval(seqBuffer, bufferIndex);
				
				flushIncompleteReferenceRegions();
				
			} catch (IllegalStateException ex) {
				err(ex.getMessage(), ex);
				return;
			}
		}
		
		return;
	}
	
	/**
	 * Set the next interval. This is called after flushing an interval.
	 */
	private void setNextInterval() {
		
		int intervalStart;
		
		// Check state
		assert (currentInterval == null) :
			"Cannot set next interval: Current interval has not been flushed: " + currentInterval;
		
		if (refIntervalIndex >= refInterval.length || ! doInterval) {
			intervalEndSeqLength = Integer.MAX_VALUE;
			
			return;
		}
		
		// Check: The next interval is out of range (first base not yet read)
		intervalStart = refInterval[refIntervalIndex].start - flankLength - 1;  // -1: Intervals start at 1, arrays at 0
		
		if (intervalStart < 0)
			intervalStart = 0;
		
		if (intervalStart > refSequenceSize) {
			intervalEndSeqLength = intervalStart;
			
			return;
		}
		
		// Next interval is in range
		currentInterval = refInterval[refIntervalIndex++];
		
		intervalEndSeqLength = currentInterval.end + flankLength;  // No -1 (end is exclusive)
		
		if (intervalEndSeqLength < 0)
			intervalEndSeqLength = Integer.MAX_VALUE;
		
		leftFlankLength = currentInterval.start - (intervalStart + 1);
		
		return;
	}
	
	/**
	 * Flush the current interval.
	 * 
	 * @param buffer Buffer.
	 * @param size buffer size;
	 * 
	 * @return The new buffer size.
	 * 
	 * @throws IllegalStateException If the sequence ends before the end of the interval.
	 */
	private int flushInterval(byte[] buffer, int size)
			throws IllegalStateException {
		
		// If not handling intervals, the whole reference sequence is the interval
		if (! doInterval) {
			
			assert (size == refSequenceSize) :
				String.format("Number of bytes in buffer (%d) must be the reference sequence size (%d) when extracting references without interval regions", size, refSequenceSize);
			
			currentInterval = new RegionInterval(refSequenceName, refSequenceName, 1, refSequenceSize, true);
		}
		
		// Do not flush if a current interval is not present
		if (currentInterval != null) {
			
			int rightFlankLength = refSequenceSize - currentInterval.end;
			
			if (rightFlankLength < 0)
				throw new IllegalStateException(String.format("Reference sequence ends with length %d before reference region was read: %s", refSequenceSize, currentInterval));
			
			if (rightFlankLength > flankLength)
				rightFlankLength = flankLength;
			
			// Add
			incompleteRegionList.add(new ReferenceRegion.IncompleteRegion(
					currentInterval, buffer, 0, leftFlankLength, rightFlankLength
			));
			
			currentInterval = null;
		}
		
		// Set the next interval
		setNextInterval();
		
		if (currentInterval == null)
			return 0;
		
		// Find the next start location
		int intervalStart = currentInterval.start - leftFlankLength - 1;
		
		// Shift buffer
		if (intervalStart < refSequenceSize) {
			int nBufferCopy = refSequenceSize - intervalStart;
			int copyStart = size - nBufferCopy;
			
			for (int copyIndex = 0; copyIndex < nBufferCopy; ++copyIndex)
				buffer[copyIndex] = buffer[copyIndex + copyStart];
			
			return nBufferCopy;
		}
		
		return 0;
	}
	
	/**
	 * Flush the list of incomplete reference regions to the set of fully initialized reference
	 * regions. This method must be called after the entire reference sequence is read. This
	 * flush also computes the sequence digest and resets the digest engine for the next
	 * reference sequence.
	 * 
	 * @throws IllegalStateException If the reference sequence is too short and at least one reference
	 *   region is incomplete.
	 */
	private void flushIncompleteReferenceRegions()
			throws IllegalStateException {
		
		ReferenceSequence referenceSequence;
		ReferenceRegion referenceRegion;
		
		// Check for missing and incomplete intervals
		if (currentInterval != null)
			throw new IllegalStateException(String.format("Missing intervals for reference sequence %s: Reference was too short to reach the end of interval: %s (reference size = %d)", refSequenceName, currentInterval, refSequenceSize));
		
		if (refIntervalIndex < refInterval.length)
			throw new IllegalStateException(String.format("Missing intervals for reference sequence %s: Reference was too short to reach the start of interval: %s (reference size = %d)", refSequenceName, refInterval[refIntervalIndex], refSequenceSize));
		
		// Complete reference regions
		if (! incompleteRegionList.isEmpty()) {
			referenceSequence = new ReferenceSequence(refSequenceName, refSequenceSize, new Digest(digestEngine.digest(), digestAlgorithm), sequenceSourceName);
			
			logger.info("Reference sequence {} (length={}, {}={})", referenceSequence.name, referenceSequence.size, referenceSequence.digest.algorithm, referenceSequence.digest.toString());
			
			for (ReferenceRegion.IncompleteRegion region : incompleteRegionList) {
				referenceRegion = region.getRegion(referenceSequence);
				
				if (reverseComplementNegativeStrand && ! referenceRegion.interval.isFwd)
					referenceRegion.reverseComplement();
				
				refRegionContainer.add(referenceRegion);
			}
		}
		
		// Reset
		incompleteRegionList.clear();
		
		return;
	}
	
	/**
	 * Set intervals for a new reference sequence.
	 */
	private void setIntervalsForReference() {
		
		// Check state
		assert (currentInterval == null) :
			"Current interval is not null";
		
		if (removeSequenceDescription) {  // Sequence name is everything before the first whitespace character
			refSequenceName = refSequenceName.split("\\s+", 2)[0];
		}
		
		// Get the set of intervals for this reference
		refInterval = intervalMap.get(refSequenceName);
		
		if (refInterval == null)
			refInterval = EMPTY_INTERVAL;
		
		refIntervalIndex = 0;
		
		// Set interval state
		setNextInterval();
	}
	
	/**
	 * Expand the sequence buffer size.
	 * 
	 * @param seqBuffer Sequence buffer.
	 * @param minExpandSize The sequence buffer must expand by at least this many bases.
	 * 
	 * @return An expanded sequence buffer with elements copied from the old buffer.
	 * 
	 * @throws IllegalArgumentException If <code>seqBuffer.length</code> is
	 *   <code>Integer.MAX_VALUE</code>.
	 */
	private byte[] expandSeqBuffer(byte[] seqBuffer, int minExpandSize)
			throws IllegalArgumentException {
		
		assert (seqBuffer != null) :
			"Cannot expand sequence buffer: null";
		
		assert (minExpandSize >= 0) :
			"minExpandSize is negative: " + minExpandSize;
		
		int size = seqBuffer.length;
		int newSize = seqBuffer.length * 2;
		byte[] newBuffer;
		
		// Add the minimum expand size to the current buffer length to get the final buffer length
		minExpandSize += seqBuffer.length;
		
		if (minExpandSize < 0)
			throw new IllegalArgumentException("Sequence is longer than the maximum size: " + Integer.MAX_VALUE);
		
		// Set the new size
		newSize = seqBuffer.length;
		
		do {
			newSize *= 2;
			
			if (newSize < 0)
				newSize = Integer.MAX_VALUE;
			
		} while (newSize < minExpandSize);
		
		// Create a larger buffer
		newBuffer = new byte[newSize];
		
		// Copy
		for (int index = 0; index < size; ++index)
			newBuffer[index] = seqBuffer[index];
		
		return newBuffer;
	}
	
	/**
	 * Signal this reader to stop.
	 */
	public void stop() {
		isActive = false;
		
		return;
	}
	
	/**
	 * Record an error to the <code>errThrowable</code> field to signal to the caller that this
	 * reader collector failed.
	 * 
	 * @param errMsg Error message.
	 * @param cause Cause of this error, or <code>null</code> if there is no throwable cause.
	 * 
	 * @throws ReadCollectorRunnerException Always.
	 */
	private void err(String errMsg, Throwable cause)
			throws ReadCollectorRunnerException {
		
		assert (errMsg != null && ! errMsg.isEmpty()) :
			"Error message is null or empty";
		
		cause = new IOException("Error parsing reference sequence: " + errMsg, cause);
		
		throw new ReadCollectorRunnerException(cause);
	}
}
