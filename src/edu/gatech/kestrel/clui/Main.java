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

package edu.gatech.kestrel.clui;

import java.io.FileNotFoundException;
import java.io.OutputStream;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import ch.qos.logback.classic.LoggerContext;
import ch.qos.logback.classic.encoder.PatternLayoutEncoder;
import ch.qos.logback.classic.spi.ILoggingEvent;
import ch.qos.logback.core.OutputStreamAppender;
import edu.gatech.kanalyze.condition.ConditionEvent;
import edu.gatech.kanalyze.util.argparse.HelpRunException;
import edu.gatech.kestrel.KestrelConstants;
import edu.gatech.kestrel.LogLevel;
import edu.gatech.kestrel.io.StreamableOutput;
import edu.gatech.kestrel.runner.ConfigurationErrorException;
import edu.gatech.kestrel.runner.KestrelRunner;
import edu.gatech.kestrel.runner.RunnerConditionListener;

/**
 * Runs the command-line user interface (CLUI). Processes command-line arguments and drives the
 * Kestrel application programming interface (API) to execute the command.
 */
public class Main {
	
	/** Thread name for the Kestrel runner. */
	public static final String RUNNER_THREAD_NAME = "KestrelRunner";
	
	/**
	 * Run Kestrel from the command-line.
	 * 
	 * @param args Command-line arguments. If <code>null</code>, an empty set of
	 *   arugments is used.
	 */
	public static void main(String[] args) {
		
		// Declarations
		Logger logger;  // Logger object
		KestrelRunner runner;  // Kestrel runner
		Thread kestrelThread;  // Thread for runner
		Throwable runnerThrowable;  // An throwable generated in the runner thread
		
		// Get logger
		logger = LoggerFactory.getLogger(Main.class);
		
		// Create the runner and argument parser
		runner = new KestrelRunner();
		
		kestrelThread = new Thread(runner);
		kestrelThread.setName(RUNNER_THREAD_NAME);
		
		runner.addListener(new RunnerConditionListener("Kestrel", System.err, runner, kestrelThread, true));
		
		// Set logging implementation for Apache Commons Logging. Uses the standard Java logging
		// first implemented in Java 1.4 (Jdk14Logger).
		System.setProperty("org.apache.commons.logging.Log", "org.apache.commons.logging.impl.Jdk14Logger");
		
		try {
			runner.configure(args);
			
		} catch (HelpRunException ex) {
			// Help was run, exit normally
			System.exit(KestrelConstants.ERR_NONE);
			
		} catch (ConfigurationErrorException ex) {
			err(ex.getMessage(), ex.errCode);  // calls System.exit()
		}
		
		// Report warnings
		for (ConditionEvent condition : runner.getWarningConditions())
			warn(condition.message);
		
		// Setup logging
		configureLogger(runner);  // Errors call System.exit()
		
		logger.info("Log level: {}", runner.getLogLevel());
		
		// Run Kestrel
		logger.trace("Starting thread: {}", RUNNER_THREAD_NAME);
		
		kestrelThread.start();
		
		while (true) {
			
			try {
				kestrelThread.join();
				break;
				
			} catch (InterruptedException ex) {
				logger.trace("Caught InterruptedException waiting for thread {}: Ignoring", RUNNER_THREAD_NAME);
			}
		}
		
		// Check state
		runnerThrowable = runner.getThrowable();
		
		if (runnerThrowable != null) {
			
			StreamableOutput logOutput = runner.getLogFile();
			
			if (! logOutput.isScreen()) {
				System.err.println(runnerThrowable.getMessage());
				runnerThrowable.printStackTrace(System.err);
			}
			
			logger.trace("main(): Complete with errors");
			
		} else {
			logger.trace("main(): Complete");
		}
		
		
		
		return;
	}
	
	/**
	 * Configure logging according to the runner's configuration. If there is an error
	 * configuring the logger, this method calls <code>System.exit()</code>.
	 * 
	 * @param runner Configure logging according to the parameters set on this runner.
	 */
	private static void configureLogger(KestrelRunner runner) {
		
		// Check arguments
		assert (runner != null) :
			"Cannot configure logger for runner: null";
		
		// Get log level
		LogLevel level = runner.getLogLevel();
		
		if (level == LogLevel.OFF)
			return;
		
		// Get logger
		org.slf4j.Logger slfRootLogger = LoggerFactory.getLogger(org.slf4j.Logger.ROOT_LOGGER_NAME);
		
		// Logger must be an instance of ch.qos.logback.classic.Logger
		// If it is not, then something strange happened (perhaps another logger implementation on the classpath?)
		if (! (slfRootLogger instanceof ch.qos.logback.classic.Logger))
			err(
					"Logging system expected a logback logger, but received: " + slfRootLogger.getClass().getName() +
					": Another SLF4J binding may be in the classpath: " + System.getProperty("java.class.path"),
					KestrelConstants.ERR_SYSTEM
			);
		
		ch.qos.logback.classic.Logger rootLogger = (ch.qos.logback.classic.Logger) slfRootLogger;
		
		// Get file
		StreamableOutput logFile = runner.getLogFile();
		OutputStream loggerOutputStream = null;
		
		try {
			loggerOutputStream = logFile.getStream();
			
		} catch (FileNotFoundException ex) {
			err("File not found while opening log file: " + logFile.name + ": " + ex.getMessage(), KestrelConstants.ERR_FILENOTFOUND);
		}
		
		// Get context
		org.slf4j.ILoggerFactory iLoggerFactory = LoggerFactory.getILoggerFactory();
		
		if (! (iLoggerFactory instanceof LoggerContext))
			err(
					"Logging system expected a logback context, but received: " + iLoggerFactory.getClass().getName() +
					"\n" +
					"Another SLF4J binding may be in the classpath: " + System.getProperty("java.class.path"),
					KestrelConstants.ERR_SYSTEM
			);
		
		LoggerContext context = (LoggerContext) iLoggerFactory;
		
		// Set encoder
		PatternLayoutEncoder encoder = new PatternLayoutEncoder();
		encoder.setContext(context);
		encoder.setPattern("%d{HH:mm:ss} [%thread] %-5level %logger{36} - %msg%n");
		encoder.start();
		
		// Create appender
		OutputStreamAppender<ILoggingEvent> appender = new OutputStreamAppender<ILoggingEvent>();
		appender.setName("KestrelLogger");
		appender.setContext(context);
		appender.setEncoder(encoder);
		appender.setOutputStream(loggerOutputStream);
		appender.start();
		
		// Add to root logger
		rootLogger.detachAndStopAllAppenders();
		rootLogger.addAppender(appender);
		
		// Set log level
		rootLogger.setLevel(level.level);
		
		return;
	}
	
	/**
	 * Report one or more error messages separated by a newline.
	 * 
	 * @param msg Message or messages to report. Must not be <code>null</code> or empty.
	 * @param errCode Error code. This must not be <code>Constants.ERR_NONE</code>.
	 */
	private static void err(String msg, int errCode) {
		
		// Check arguments
		assert (msg != null) :
			"Attempted error message with message string: null";
		
		msg = msg.trim();
		
		assert (! msg.isEmpty()) :
			"Attempted error message with an empty string";
		
		// Split into multiple messages.
		String[] msgTok = msg.split("\\s*,\\s*");
		
		// Output messages
		for (String message : msgTok) {
			
			if (message.isEmpty())
				continue;
			
			System.err.println(KestrelConstants.PROG_NAME + ": ERROR: " + message);
		}
		
		// Print -h option
		if (errCode == KestrelConstants.ERR_USAGE)
			System.err.println(KestrelConstants.PROG_NAME + ": Try the -h option for help");
		
		// Terminate
		if (errCode == KestrelConstants.ERR_NONE) {
			System.err.println(KestrelConstants.PROG_NAME + ": Program attempted to return non-error code" + KestrelConstants.ERR_NONE + " for an error condition (please report this as a program bug and include the command that was run)");
			errCode = KestrelConstants.ERR_SYSTEM;
		}
		
		System.exit(errCode);
	}
	
	/**
	 * Report one or more warning messages separated by a newline.
	 * 
	 * @param msg Message or messages to report. Must not be <code>null</code> or empty.
	 */
	private static void warn(String msg) {
		
		// Check arguments
		assert (msg != null) :
			"Attempted warning message with message string: null";
		
		msg = msg.trim();
		
		assert (! msg.isEmpty()) :
			"Attempted warning message with an empty string";
		
		// Split into multiple messages.
		String[] msgTok = msg.split("\\s*,\\s*");
		
		// Output messages
		for (String message : msgTok) {
			
			if (message.isEmpty())
				continue;
			
			System.err.println(KestrelConstants.PROG_NAME + ": WARNING: " + message);
		}
		
		return;
	}
}
