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

import edu.gatech.kanalyze.util.BoundedQueue;

/**
 * Executes readers.
 */
public class ReaderThreadRunner
		implements Thread.UncaughtExceptionHandler, Runnable {
	
	/** Error if one was thrown. */
	private Throwable errThrowable;
	
	/** Reader. */
	private final ReaderRunner readerRunner;
	
	/** Read collector. */
	private final ReadCollectorRunner collectorRunner;
	
	/** Queue to be stopped. */
	private final BoundedQueue<?> sequenceQueue;
	
	/** Reader thread. */
	private final Thread readerThread;
	
	/** Read collector thread. */
	private final Thread collectorThread;
	
	/** <code>true</code> until this runner is stopped. */
	private boolean isActive;
	
	/**
	 * Create a new runner object.
	 * 
	 * @param readerRunner Reader.
	 * @param collectorRunner Read collector.
	 * @param sequenceQueue Sequence queue between the runners. This is stopped by
	 *   <code>run()</code> when the reader runner stops.
	 * 
	 * @throws NullPointerException If any argument is <code>null</code>.
	 */
	public ReaderThreadRunner(ReaderRunner readerRunner, ReadCollectorRunner collectorRunner, BoundedQueue<?> sequenceQueue)
			throws NullPointerException {
		
		// Check arguments
		if (readerRunner == null)
			throw new NullPointerException("Cannot create reader thread callback with reader: null");
		
		if (collectorRunner == null)
			throw new NullPointerException("Cannot create reader thread callback with read collector: null");
		
		// Set runners
		this.readerRunner = readerRunner;
		this.collectorRunner = collectorRunner;
		
		// Create threads
		readerThread = new Thread(readerRunner);
		collectorThread = new Thread(collectorRunner);
		
		readerThread.setName("ReferenceReader");
		readerThread.setUncaughtExceptionHandler(this);
		readerThread.setDaemon(true);
		
		collectorThread.setName("ReferenceReadCollector");
		collectorThread.setUncaughtExceptionHandler(this);
		collectorThread.setDaemon(true);
		
		this.sequenceQueue = sequenceQueue;
		
		// Set state
		errThrowable = null;
		isActive = true;
		
		return;
	}
	
	/**
	 * Throw an error and stop this reader.
	 * 
	 * @param errMsg Error message.
	 * @param cause Cause of this error, or <code>null</code> if there is no throwable cause.
	 */
	public void error(String errMsg, Throwable cause) {
		
		// Check arguments
		if (errMsg == null)
			errMsg = "";
		
		if (errMsg.isEmpty())
			errMsg = "Unknown error (software bug)";
		
		if (cause == null)
			cause = new IOException(errMsg);
		
		// Set error
		if (errThrowable == null)
			errThrowable = cause;
		
		stop();
		
		return;
	}
	
	/**
	 * Get a throwable object if either component failed.
	 * 
	 * @return Throwable object or <code>null</code> if no errors were thrown.
	 */
	public Throwable getThrowable() {
		return errThrowable;
	}
	
	/**
	 * Catch exception and save in <code>threadThrowable[threadIndex]</code> if another
	 * throwable is not already saved there.
	 * 
	 * @param t Thread throwable came from.
	 * @param e Throwable.
	 */
	@Override
	public void uncaughtException(Thread t, Throwable e) {
		
		// Assign if a throwable was not already assigned
		if (errThrowable == null)
			errThrowable = e;
		
		stop();
		
		return;
	}
	
	/**
	 * Run the reader threads.
	 */
	@Override
	public void run() {
		
		// Start threads
		collectorThread.start();
		readerThread.start();
		
		// Join reader thread
		while (isActive) {
			
			try {
				readerThread.join();
				break;
				
			} catch (InterruptedException ex) {
				// Ignore
			}
		}
		
		// Stop the sequence queue
		sequenceQueue.shutdownIgnoreInterrupt();
		
		// Join collector thread
		while (isActive) {
			
			try {
				collectorThread.join();
				break;
				
			} catch (InterruptedException ex) {
				// Ignore
			}
		}
		
		return;
	}
	
	/**
	 * Stop the running threads.
	 */
	public void stop() {
		
		isActive = false;
		
		sequenceQueue.shutdownIgnoreInterrupt();
		
		readerRunner.stop();
		collectorRunner.stop();
		
		if (readerThread.isAlive())
			readerThread.interrupt();
		
		if (collectorThread.isAlive())
			collectorThread.interrupt();
	}
}
