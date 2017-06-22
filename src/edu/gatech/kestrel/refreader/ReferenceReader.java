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
import java.util.ConcurrentModificationException;
import java.util.Map;
import java.util.concurrent.locks.Lock;
import java.util.concurrent.locks.ReentrantLock;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import edu.gatech.kanalyze.KAnalyzeConstants;
import edu.gatech.kanalyze.batch.BatchCache;
import edu.gatech.kanalyze.batch.SequenceBatch;
import edu.gatech.kanalyze.batch.SequenceBatchCache;
import edu.gatech.kanalyze.comp.reader.SequenceSource;
import edu.gatech.kanalyze.util.BoundedQueue;
import edu.gatech.kanalyze.util.SequenceNameTable;
import edu.gatech.kanalyze.util.kmer.KmerUtil;
import edu.gatech.kestrel.interval.RegionInterval;


/**
 * Base class for all reference sequence readers.
 */
public class ReferenceReader {
	
	/** K-mer size. */
	public final KmerUtil kUtil;
	
	/** Logger. */
	private final Logger logger = LoggerFactory.getLogger(ReferenceReader.class);
	
	/** The length of flanks to add to intervals when extracting regions of a reference sequence. */
	private int flankLength;
	
	/** Loader used for finding reader classes. */
	private ClassLoader loader;
	
	/** Run lock. Prevents race conditions if more than one thread attempts to use this reader. */
	private Lock runLock;
	
	/** Thread that has acquired <code>runLock</code>. */
	private Thread currentThread;
	
	/** Executes threads. */
	private ReaderThreadRunner threadRunner;
	
	/** <code>true</code> if the sequence description should be removed from a sequence name. */
	private boolean removeSequenceDescription;
	
	/** <code>true</code> if regions on the negative strand should be reverse complemented before variant calling. */
	private boolean reverseComplementNegativeStrand;
	
	/** Default remove sequence description property. */
	public static final boolean DEFAULT_REMOVE_SEQUENCE_DESCRIPTION = ReadCollectorRunner.DEFAULT_REMOVE_SEQUENCE_DESCRIPTION;
	
	/** Default reverse-complement negative strand regions. */
	public static final boolean DEFAULT_REVERSE_COMPLEMENT_NEGATIVE_STRAND = false;
	
	
	/**
	 * Create a new reference reader.
	 * 
	 * @param kUtil K-mer utility.
	 * @param sequenceBatchSize Size of sequence batches returned by the sequence reader.
	 * @param sequenceCacheSize Number of sequence batches to cache as they are returned to the
	 *   sequence reader.
	 * @param sequenceQueueSize Number of sequence batches to queue that holds sequence batches from the
	 *   sequence reader. If less than 1, the KAnalyze default sequence queue size is used.
	 * @param loader Class loader for finding sequence readers. If <code>null</code>, the default class
	 *   loader is used.
	 * 
	 * @throws NullPointerException If <code>kUtil</code> is <code>null</code>.
	 * @throws IllegalArgumentException If <code>sequenceBatchSize</code> or <code>sequenceCacheSize</code>
	 *   is too small for the sequence batch cache component.
	 */
	public ReferenceReader (KmerUtil kUtil, int sequenceBatchSize, int sequenceCacheSize, int sequenceQueueSize, ClassLoader loader)
			throws NullPointerException, IllegalArgumentException {
		
		if (kUtil == null)
			throw new NullPointerException("K-mer utility is null");
		
		if (loader == null)
			loader = ReferenceReader.class.getClassLoader();
		
		this.kUtil = kUtil;
		this.loader = loader;
		
		flankLength = 0;
		
		removeSequenceDescription = DEFAULT_REMOVE_SEQUENCE_DESCRIPTION;
		reverseComplementNegativeStrand = DEFAULT_REVERSE_COMPLEMENT_NEGATIVE_STRAND;
		
		runLock = new ReentrantLock();
		
		return;
	}
	
	/**
	 * Create a new reference reader.
	 * 
	 * @param kUtil K-mer utility.
	 * @param loader Class loader for finding sequence readers. If <code>null</code>, the default system class
	 *   loader is used.
	 * 
	 * @throws NullPointerException If <code>kUtil</code> is <code>null</code>.
	 */
	public ReferenceReader (KmerUtil kUtil, ClassLoader loader)
			throws NullPointerException {
		
		this(kUtil, KAnalyzeConstants.DEFAULT_SEQUENCE_BATCH_SIZE, KAnalyzeConstants.DEFAULT_CACHE_SIZE, KAnalyzeConstants.DEFAULT_SEQUENCE_QUEUE_SIZE, loader);
	}
	
	/**
	 * Create a new reference reader.
	 * 
	 * @param kUtil K-mer utility.
	 * 
	 * @throws NullPointerException If <code>kUtil</code> is <code>null</code>.
	 */
	public ReferenceReader(KmerUtil kUtil)
			throws NullPointerException {
		
		this(kUtil, KAnalyzeConstants.DEFAULT_SEQUENCE_BATCH_SIZE, KAnalyzeConstants.DEFAULT_CACHE_SIZE, KAnalyzeConstants.DEFAULT_SEQUENCE_QUEUE_SIZE, null);
		
		return;
	}
	
	/**
	 * Set the length of flanks for interval regions. Adding flanks helps active region
	 * detection, but it does not otherwise affect the variant calls. If there are no
	 * intervals and the reader extracts full reference sequences, the length of flanks
	 * has no effect.
	 * 
	 * @param flankLength Length of flanks or <code>0</code> to disable flanks.
	 * 
	 * @throws IllegalArgumentException If <code>flankLength</code> is negative.
	 */
	public void setFlankLength(int flankLength)
			throws IllegalArgumentException {
		
		if (flankLength < 0)
			throw new IllegalArgumentException("Flank length is negative: " + flankLength);
		
		this.flankLength = flankLength;
	}
	
	/**
	 * Get the length of flanks to add to intervals when extracting regions from reference
	 * sequences. Adding flanks helps active region detection, but it does not otherwise affect
	 * the variant calls.
	 * 
	 * @return Length of flanks or <code>0</code> if flanks are disabled.
	 */
	public int getFlankLength() {
		return flankLength;
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
	 * Read a sequence source and return a list of reference sequences.
	 * 
	 * @param sources Sequence sources. If <code>null</code> or empty, an empty reference
	 *   container is returned.
	 * @param intervalMap A map of intervals where regions should be extracted for variant
	 *   calling, or <code>null</code> to extract the whole of each reference sequence. 
	 * 
	 * @return A container with reference regions where variants should be called.
	 * 
	 * @throws IOException If there is any error reading the file including problems with the file format.
	 * @throws ConcurrentModificationException If more than one thread attempts to execute this runner.
	 */
	public ReferenceRegionContainer read(SequenceSource[] sources, Map<String, RegionInterval[]> intervalMap)
		throws IOException, ConcurrentModificationException {
		
		if (! runLock.tryLock()) {
			Thread currentThreadSafe = currentThread;  // Avoid race conditions
			
			throw new ConcurrentModificationException("Another thread is executing this reference runner: " + ((currentThreadSafe != null) ? currentThreadSafe.getName() : "Unknown Thread"));
		}
		
		try {
			currentThread = Thread.currentThread();
			
			// Region container
			ReferenceRegionContainer refRegionContainer;  // Reference regions read
			
			// Reader pipeline
			ReaderRunner readerRunner;
			ReadCollectorRunner collectorRunner;
			
			SequenceNameTable nameTable;
			BoundedQueue<SequenceBatch> sequenceQueue;
			BatchCache<SequenceBatch> sequenceBatchCache;
			
			// Runner
			Throwable errThrowable;
			
			// Check sources
			if (sources == null || sources.length == 0) {
				logger.trace("Reading 0 sequence sources (source array is {})", ((sources == null) ? "null" : "empty"));
				return new ReferenceRegionContainer();
			}
			
			// Init data structures
			refRegionContainer = new ReferenceRegionContainer();
			nameTable = new SequenceNameTable();
			
			// Setup reader components
			sequenceQueue = new BoundedQueue<>();
			sequenceBatchCache = new SequenceBatchCache(KmerUtil.get(1), -1, KAnalyzeConstants.DEFAULT_SEQUENCE_BATCH_SIZE, 0);
			
			readerRunner = new ReaderRunner(nameTable, sequenceQueue, sequenceBatchCache, loader);
			collectorRunner = new ReadCollectorRunner(refRegionContainer, nameTable, sequenceQueue, sequenceBatchCache, intervalMap, flankLength);
			
			collectorRunner.setRemoveDescription(removeSequenceDescription);
			collectorRunner.setRevComplementNegStrand(reverseComplementNegativeStrand);
			
			readerRunner.setSource(sources);
			sequenceQueue.noShutdownIgnoreInterrupt();
			
			logger.trace("Reading {} sequence sources with {} intervals (flank length = {})", sources.length, (intervalMap != null) ? intervalMap.size() : 0, flankLength);
			
			// Create callback object
			threadRunner = new ReaderThreadRunner(readerRunner, collectorRunner, sequenceQueue);
			
			// Prepare to run
			logger.trace("Reading {} sequence source(s)", sources.length);
			
			threadRunner.run();
			
			errThrowable = threadRunner.getThrowable();
			
			// Check for errors
			if (errThrowable != null) {
				
				if (errThrowable instanceof ReadCollectorRunnerException)
					errThrowable = ((ReadCollectorRunnerException) errThrowable).cause;
				
				if (errThrowable instanceof IOException)
					throw (IOException) errThrowable;
				
				throw new IOException(errThrowable.getMessage(), errThrowable);
			}
			
			logger.trace("Done reading sequence sources");
			
			// Return reference sequences
			return refRegionContainer;
			
		} finally {
			runLock.unlock();
		}
	}
	
	/**
	 * Signal this reader to stop.
	 */
	public void stop() {
		
		ReaderThreadRunner threadRunner = this.threadRunner;  // Avoid race conditions on the object field
		
		if (threadRunner != null)
			threadRunner.stop();
		
		return;
	}
}
