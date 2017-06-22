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
import java.util.Hashtable;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import edu.gatech.kanalyze.batch.BatchCache;
import edu.gatech.kanalyze.batch.SequenceBatch;
import edu.gatech.kanalyze.comp.reader.ReaderInitException;
import edu.gatech.kanalyze.comp.reader.SequenceReader;
import edu.gatech.kanalyze.comp.reader.SequenceSource;
import edu.gatech.kanalyze.util.BoundedQueue;
import edu.gatech.kanalyze.util.SequenceNameTable;
import edu.gatech.kanalyze.util.kmer.KmerUtil;

/**
 * Runs KAnalyze sequence readers.
 */
public class ReaderRunner implements Runnable {
	
	/** Logger. */
	private final Logger logger = LoggerFactory.getLogger(ReaderRunner.class);
	
	/** Queue of sequence reads. */
	private final BoundedQueue<SequenceBatch> sequenceQueue;
	
	/** Cache for recycling sequence batches. */
	private final BatchCache<SequenceBatch> sequenceBatchCache;
	
	/** Table of sequence names written by the sequence reader. */
	private final SequenceNameTable nameTable;
	
	/** <code>true</code> unless this collector is stopped. */
	private boolean isActive;
	
	/** List of sources to read. */
	private SequenceSource[] sources;
	
	/** Set of sequence readers allocated. */
	private Hashtable<String, SequenceReader> readerTable;
	
	/** Loader used for finding reader classes. */
	private ClassLoader loader;
	
	/** Reader. */
	private SequenceReader reader;
	
	
	/**
	 * Create a sequence read runner.
	 * 
	 * @param nameTable Table of sequence and sequence source names.
	 * @param sequenceQueue Sequence queue.
	 * @param sequenceBatchCache Sequence batch cache.
	 * @param loader Class loader or <code>null</code> to use the default loader.
	 */
	public ReaderRunner(SequenceNameTable nameTable, BoundedQueue<SequenceBatch> sequenceQueue, BatchCache<SequenceBatch> sequenceBatchCache, ClassLoader loader) {
		
		// Check arguments
		if (nameTable == null)
			throw new NullPointerException("Name table is null");
		
		if (sequenceQueue == null)
			throw new NullPointerException("Sequence queue is null");
		
		if (sequenceBatchCache == null)
			throw new NullPointerException("Sequence batch cache is null");
		
		if (loader == null)
			loader = ReaderRunner.class.getClassLoader();
		
		// Set fields
		this.nameTable = nameTable;
		this.sequenceQueue = sequenceQueue;
		this.sequenceBatchCache = sequenceBatchCache;
		this.loader = loader;
		
		readerTable = new Hashtable<String, SequenceReader>();
		
		sources = null;
		
		isActive = true;
		
		return;
	}
	
	/**
	 * Run this reader.
	 * 
	 * @throws ReadCollectorRunnerException If any error occurs.
	 */
	@Override
	public void run()
			throws ReadCollectorRunnerException {
		
		if (sources == null)
			return;
		
		// Read each source
		for (SequenceSource source : sources) {
			
			logger.trace("Reading source: {}", source);
			
			if (! isActive)
				break;
			
			if (source == null)
				continue;
			
			// Get reader from cache
			reader = readerTable.get(source.formatType);
			
			// Create reader if not in cache
			if (reader == null) {
				
				try {
					reader = SequenceReader.getReader(source.formatType, loader);
					
				} catch (ReaderInitException ex) {
					err(ex);
					return;
				}
				
				if (reader == null) {
					err("No sequence reader for format: " + source.formatType, null);
					return;
				}
				
				logger.trace("Loaded reader: {} ({})", reader.name, reader.getClass().getName());
			}
			
			// Initialize reader
			try {
				reader.init(KmerUtil.get(1), sequenceQueue, sequenceBatchCache, nameTable, null, loader);
				
			} catch (IllegalArgumentException ex) {
				err(ex);
				return;
				
			} catch (IllegalStateException ex) {
				err(ex);
				return;
			}
			
			// Read sequences
			try {
				reader.read(source);
				
			} catch (IllegalArgumentException ex) {
				err(ex);
				return;
				
			} catch (IOException ex) {
				err(ex);
				return;
				
			} catch (SecurityException ex) {
				err(ex);
				return;
				
			} catch (IllegalStateException ex) {
				err(ex);
				return;
			}
		}
		
		return;
	}
	
	/**
	 * Set sequence sources for this reader. This should not be called while the runner
	 * is executing. Any sources added that have not been processed will be discarded.
	 * 
	 * @param sources A list of sequence sources.
	 */
	public void setSource(SequenceSource[] sources) {
		
		if (sources == null)
			sources = new SequenceSource[0];
		
		this.sources = sources;
	}
	
	/**
	 * Signal this reader to stop.
	 */
	public void stop() {
		
		SequenceReader reader = this.reader;
		
		isActive = false;
		
		if (reader != null)
			reader.stop();
		
		return;
	}
	
	/**
	 * Generate an error in <code>run</code>.
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
		
		cause = new IOException("Error reading reference sequence: " + errMsg, cause);
		
		throw new ReadCollectorRunnerException(cause);
	}
	
	/**
	 * Generate an error in <code>run</code>.
	 * 
	 * @param cause Cause of this error, or <code>null</code> if there is no throwable cause.
	 * 
	 * @throws ReadCollectorRunnerException Always.
	 */
	private void err(Throwable cause)
			throws ReadCollectorRunnerException {
		
		assert (cause != null) :
			"Error cause is null or empty";
		
		if (cause == null)
			cause = new IOException("Unknown error");
		
		cause = new IOException("Error reading reference sequence: " + cause.getMessage(), cause);
		
		throw new ReadCollectorRunnerException(cause);
	}
}