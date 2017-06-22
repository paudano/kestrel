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

package edu.gatech.kestrel.runner;

import java.io.PrintStream;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import edu.gatech.kanalyze.KAnalyzeConstants;
import edu.gatech.kanalyze.KAnalyzeRunnable;
import edu.gatech.kanalyze.condition.ConditionEvent;
import edu.gatech.kanalyze.condition.ConditionListener;

/**
* Listen for all conditions and output a formatted message to a specified print stream.
* If an error occurs, terminate the JVM this listener is running in.
*/
public class RunnerConditionListener implements ConditionListener {
	
	/** Logger. */
	private final Logger logger = LoggerFactory.getLogger(RunnerConditionListener.class);
	
	/** Name of the application reporting the error. */
	private String appName;
	
	/** Output stream. */
	private PrintStream stream;
	
	/**
	 * Runnable to abort if and error condition occurs and <code>abortOnError</code>
	 * is <code>true</code>.
	 */
	private KAnalyzeRunnable runnable;
	
	/** Thread to interrupt if <code>runnable</code> is stopped. */
	private Thread runnableThread;
	
	/** Write a stack trace with messages when <code>true</code>. */
	private boolean stackTrace;
	
	/**
	 * Ignore the <code>trace</code> flag of conditions. If <code>false</code>, a stack trace is
	 * reported for error conditions with a non-<code>null</code> cause.
	 */
	private boolean ignoreTraceFlag;
	
	/** Default application name. */
	public static final String DEFAULT_APP_NAME = "UnknownRunnable";
	
	/**
	 * Create a new stream condition listener.
	 * 
	 * @param appName Name of the application reporting errors. This name is
	 *   used for error messages. If <code>null</code> or empty,
	 *   <code>DEFAULT_APP_NAME</code> is used.
	 * @param stream Stream output is set to. If <code>null</code>, output is
	 *   written to <code>System.err</code>.
	 * @param runnable Runnable to stop if an error condition is received. If
	 *   <code>null</code>, no module will be aborted.
	 * @param runnableThread Thread to be interrputed if <code>runnable</code> is
	 *   stopped or <code>null</code> if no thread should be interrupted.
	 * @param stackTrace If <code>true</code>, output the stack trace.
	 *   
	 * @see #DEFAULT_APP_NAME
	 */
	public RunnerConditionListener(String appName, PrintStream stream, KAnalyzeRunnable runnable, Thread runnableThread, boolean stackTrace) {
		
		if (appName == null)
			appName = "";
		
		appName = appName.trim();
		
		if (appName.isEmpty())
			appName = DEFAULT_APP_NAME;
		
		if (stream == null)
			stream = System.err;
		
		this.appName = appName;
		this.stream = stream;
		this.runnable = runnable;
		this.runnableThread = runnableThread;
		this.stackTrace = stackTrace;
		
		ignoreTraceFlag = false;
		
		return;
	}
	
	/**
	 * Create a new stream condition listener that will write messages to
	 * <code>System.err</code>.
	 * 
	 * @param runnable Runnable to stop if an error condition is received. If
	 *   <code>null</code>, no module will be aborted.
	 * @param runnableThread Thread to be interrputed if <code>runnable</code> is
	 *   stopped or <code>null</code> if no thread should be interrupted.
	 * @param stackTrace If <code>true</code>, output the stack trace.
	 */
	public RunnerConditionListener(KAnalyzeRunnable runnable, Thread runnableThread, boolean stackTrace) {
		
		this (null, null, runnable, runnableThread, stackTrace);
		
		return;
	}
	
	/**
	 * Create a new stream condition listener that will write messages to
	 * <code>System.err</code>.
	 * 
	 * @param runnable Runnable to stop if an error condition is received. If
	 *   <code>null</code>, no module will be aborted.
	 * @param runnableThread Thread to be interrputed if <code>runnable</code> is
	 *   stopped or <code>null</code> if no thread should be interrupted.
	 */
	public RunnerConditionListener(KAnalyzeRunnable runnable, Thread runnableThread) {
		
		this (null, null, runnable, runnableThread, false);
		
		return;
	}
	
	
	/**
	 * Determine if this listener will output a stack trace with all condition
	 * messages. This is a debugging tool and should not be generally be enabled
	 * 
	 * @return <code>true</code> if a module would be output a stack trace
	 *   with condition messages.
	 */
	public boolean isStackTrace() {
		return stackTrace;
	}
	
	/**
	 * Set to <code>true</code> to print a stack trace with all conditions.
	 * 
	 * @param stackTrace <code>true</code> to print a stack trace with
	 *   conditions.
	 */
	public void setStackTrace(boolean stackTrace) {
		this.stackTrace = stackTrace;
	}
	
	/**
	 * Set to <code>true</code> to ignore the <code>trace</code> flag in conditions.
	 * If <code>false</code> (default), then a stack trace is reported when the
	 * condition is an error, the cause is not <code>null</code>, and the condition&apos;s
	 * <code>trace</code> flag is <code>true</code>.
	 * 
	 * @return Ignore trace flag.
	 */
	public boolean isIgnoreTraceFlag() {
		return ignoreTraceFlag;
	}
	
	/**
	 * Set to <code>true</code> to ignore the <code>trace</code> flag in conditions.
	 * If <code>false</code> (default), then a stack trace is reported when the
	 * condition is an error, the cause is not <code>null</code>, and the condition&apos;s
	 * <code>trace</code> flag is <code>true</code>.
	 * 
	 * @param ignoreTraceFlag Flag to be set.
	 */
	public void setIgnoreTraceFlag(boolean ignoreTraceFlag) {
		this.ignoreTraceFlag = ignoreTraceFlag;
		
		return;
	}

	/**
	 * Handle condition.
	 * 
	 * @param condition Condition to be handled.
	 */
	@Override
	public void conditionOccurred(ConditionEvent condition) {
		
		if (condition == null)
			return;
		
		boolean isError = condition.isError();
		
		logger.trace("Caught condition: {} (error={})", condition, isError);
		
		// Output error
		stream.println(
			appName + ": " +
			condition.toString()
		);
		
		if (condition.returnCode == KAnalyzeConstants.ERR_USAGE)
			stream.println(appName + ": Try the -h option for help");
		
		// Write stack trace
		if (condition.cause != null) {
			if (stackTrace || (condition.isError() && condition.trace))
				condition.cause.printStackTrace(stream);
		}
		
		// Abort module
		if (isError && runnable != null) {
			try {
				runnable.stop();
				
			} catch (UnsupportedOperationException ex) {
				// Ignore (module does not support abort())
			}
			
			if (runnableThread != null)
				runnableThread.interrupt();
		}
		
		System.exit(condition.returnCode);
		
		return;
	}
}
