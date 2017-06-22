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

package edu.gatech.kestrel.counter;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ConcurrentModificationException;
import java.util.concurrent.locks.ReentrantLock;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import edu.gatech.kanalyze.comp.reader.SequenceSource;
import edu.gatech.kanalyze.module.count.CountModule;
import edu.gatech.kanalyze.util.kmer.KmerUtil;
import edu.gatech.kestrel.io.InputSample;

/**
 * An interface for classes that manage k-mer counts for Kestrel.
 */
public abstract class CountMap {
	
	/** K-mer utility for k-mers managed by this map. */
	public final KmerUtil kUtil;
	
	/** Logger object. */
	private final Logger logger;
	
	/** K-mer count module. */
	protected final CountModule countModule;
	
	/**
	 * Sample this map is currently processing or <code>null</code> before
	 * <code>set()</code> or if <code>set()</code> does not complete normally.
	 */
	protected InputSample sample;
	
	/** Lock for map operations. */
	private final ReentrantLock mapLock;
	
	/** Thread that has mapLock. */
	private Thread mapLockThread;
	
	/** Set to <code>true</code> if the map was aborted. */
	private boolean isAborted;
	
	/**
	 * Create a new count map.
	 * 
	 * @param kUtil K-mer utility.
	 * @param countModule A configured count module.
	 * 
	 * @throws NullPointerException If <code>kUtil</code> or <code>countModule</code>
	 *   is <code>null</code>.
	 */
	protected CountMap(KmerUtil kUtil, CountModule countModule)
			throws NullPointerException {
		
		// Check arguments
		if (kUtil == null)
			throw new NullPointerException("K-mer utility is null");
		
		if (countModule == null)
			throw new NullPointerException("Count module is null");
		
		// Set fields
		logger = LoggerFactory.getLogger(CountMap.class);
		
		mapLock = new ReentrantLock();
		mapLockThread = null;
		isAborted = false;
		
		this.kUtil = kUtil;
		this.countModule = countModule;
		
		sample = null;
		
		return;
	}
	
	
	//
	// Implementation defined methods
	//
	
	/**
	 * Get a k-mer from this map. This method must be called after <code>set()</code>
	 * or the results are undefined.
	 * 
	 * @param kmer K-mer to get.
	 * 
	 * @return The k-mer count or <code>0</code> if the k-mer is not in this
	 *   map.
	 * 
	 * @throws NullPointerException If <code>kmer</code> is <code>null</code..
	 * @throws IllegalArgumentException If <code>kmer</code> length is less than
	 *   <code>kUtil.wordSize</code>.
	 * @throws IllegalStateException The implementation may throw this exception if
	 *   <code>set()</code> was never called.
	 */
	public abstract int get(int[] kmer)
			throws NullPointerException, IllegalArgumentException, IllegalStateException;
	
	/**
	 * Called after the sample is set on the count module, but before it is run.
	 * 
	 * @return <code>true</code> if the count module must be run. Some maps can read directly
	 *   from a file and do not need the count module. Returning <code>false</code> does not
	 *   indicate an error.
	 *   
	 * @throws Exception May throw any exception.
	 */
	protected boolean preModuleRun()
			throws Exception {
		
		return true;
	}
	
	/**
	 * Called after the count module is run, the count module fails, or after
	 * <code>preModuleRun()</code> fails.
	 * 
	 * @param onError <code>true</code> if the module or <code>preModuleRun()</code>
	 *   threw an exception.
	 * @param aborted <code>true</code> if the pipeline is stopped because <code>abort()</code>
	 *   was called.
	 * 
	 * @throws Exception May throw any exception.
	 */
	protected void postModuleRun(boolean onError, boolean aborted)
			throws Exception {
		
		return;
	}
	
	
	//
	// Public interface
	//
	
	/**
	 * Set a new sample for this count map.
	 * 
	 * @param sample Sample to set.
	 * 
	 * @throws NullPointerException If <code>sample</code> is <code>null</code>.
	 * @throws FileNotFoundException If <code>sample</code> contains a file that cannot be found.
	 * @throws IOException If any error occurs setting the sample.
	 * @throws ConcurrentModificationException If another thread is attempting to run this method.
	 */
	public final void set(InputSample sample)
			throws NullPointerException, FileNotFoundException, IOException, ConcurrentModificationException {
		
		boolean runCountPipeline;
		
		boolean noChain;  // Used by exception handlers to avoid chaining the same exception more than once
		
		// Check arguments
		if (sample == null)
			throw new NullPointerException("Sample is null");
		
		// Check lock
		if (! mapLock.tryLock())
			throw new ConcurrentModificationException("Cannot set sample: Another thread is attempting to set the sample at the same time: " + mapLockThread.getName());
		
		mapLockThread = Thread.currentThread();
		isAborted = false;
		
		noChain = false;
		
		try {
			
			// Set sources
			this.sample = sample;
			
			countModule.clearSources();
			
			for (SequenceSource seqSource : sample.sources)
				countModule.addSource(seqSource);
			
			if (isAborted)
				return;
			
			// Implementation-defined pre-count run
			try {
				runCountPipeline = preModuleRun();
				
			} catch (Exception ex) {
				logger.error("Error setting k-mer counts for sample {}: {} ({})", sample.name, ex.getMessage(), ex.getClass().getSimpleName());
				
				try {
					postModuleRun(true, false);
					
				} catch (Exception ex2) {
					logger.warn("Error in post-module run: {} ({})", ex2.getMessage(), ex2.getClass().getSimpleName());
					return;
				}
	
				this.sample = null;
				
				throw new IOException(String.format("Error setting k-mer counts (in map pre-run) for sample %s: %s (%s)", sample.name, ex.getMessage(), ex.getClass().getSimpleName()));
			}
			
			if (isAborted)
				return;
			
			try {
				
				if (runCountPipeline) {
					// Run count module
					logger.trace("Getting k-mer counts: {}", sample.name);
					
					countModule.run();
				}
				
				// Implementation-defined post-count run
				try {
					postModuleRun(false, false);
					
				} catch (Exception ex2) {
					logger.warn("Error in post-module run: {} ({})", ex2.getMessage(), ex2.getClass().getSimpleName());
					
					noChain = true;
					
					throw new IOException(String.format("Error setting k-mer counts (in map post-run) for sample %s: %s (%s)", sample.name, ex2.getMessage(), ex2.getClass().getSimpleName()));
				}
				
			} catch (Throwable ex) {
				logger.error("Unexpected error setting sample {}: {} ({})", sample.name, ex.getMessage(), ex.getClass().getSimpleName());
				
				try {
					postModuleRun(true, false);
					
				} catch (Exception ex2) {
					logger.warn("Error in post-module run: {} ({})", ex2.getMessage(), ex2.getClass().getSimpleName());
				}
				
				this.sample = null;
				
				if (noChain)
					throw ex;
				
				throw new IOException(String.format("Error setting k-mer counts (in count module) for sample %s: %s (%s)", sample.name, ex.getMessage(), ex.getClass().getSimpleName()));
			}
			
		} finally {
			mapLockThread = null;
			mapLock.unlock();
		}
		
		return;
	}
	
	/**
	 * Abort the count pipeline if it is running.
	 */
	public final void abort() {
		
		Thread mapLockThread;
		
		logger.info("Aborting count pipeline");
		
		// Check lock
		if (! mapLock.isLocked())
			return;
		
		// Prepare
		mapLockThread = this.mapLockThread;  // Save to avoid a race condition
		isAborted = true;
		
		// Abort pipeline
		countModule.abort();
		
		// Interrupt thread
		if (mapLockThread != null)
			mapLockThread.interrupt();
		
		// Implementation-defined post-count run
		try {
			postModuleRun(false, true);
			
		} catch (Exception ex2) {
			logger.warn("Error in post-module run (run aborted): {} ({})", ex2.getMessage(), ex2.getClass().getSimpleName());
		}
		
		return;
	}
}
