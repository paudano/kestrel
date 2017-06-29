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

package edu.gatech.kestrel;

import ch.qos.logback.classic.Level;

/**
 * Organizes logging levels for this program and facilitates translating names to log levels. 
 */
public enum LogLevel {
	
	/** Log all messages. */
	ALL (Level.ALL),
	
	/** Log all up to program trace messages. */
	TRACE (Level.TRACE),
	
	/** Log all up to debug messages. */
	DEBUG (Level.DEBUG),
	
	/** Log all up to information messages. */
	INFO (Level.INFO),
	
	/** Log all up to warning messages. */
	WARN (Level.WARN),
	
	/** Log only error messages. */
	ERROR (Level.ERROR),
	
	/** Disable logging. */
	OFF (Level.OFF);
	
	/** Logback level associated with this log level element. */
	public final Level level;
	
	/**
	 * Create a new log level.
	 * 
	 * @param level Logback level associated with this level.
	 */
	private LogLevel(Level level) {
		
		this.level = level;
		
		return;
	}
	
	/**
	 * Get a log level associated with a level name. <code>Level.toLevel(String)</code> does
	 * not report errors (just returns <code>Level.DEBUG</code>, so this method checks the
	 * string for a known level before returning one.
	 * 
	 * @param levelName Level name (case insensitive).
	 * 
	 * @return The log level identified by <code>levelName</code>.
	 * 
	 * @throws NullPointerException If <code>levelName</code> is <code>null</code>.
	 * @throws IllegalArgumentException If <code>levelName</code> is empty or is not a valid
	 *   log level name.
	 */
	public static LogLevel getLevel(String levelName)
			throws NullPointerException, IllegalArgumentException {
		
		// Check argument
		if (levelName == null)
			throw new NullPointerException("Cannot get log level for level name: null");
		
		levelName = levelName.trim();
		
		if (levelName.isEmpty())
			throw new IllegalArgumentException("Cannot get log level for level with an empty name");
		
		String orgLevelName = levelName;  // Save for error reporting
		levelName = levelName.toUpperCase();
		
		// Search
		for (LogLevel logLevel : LogLevel.values())
			if (logLevel.name().equals(levelName))
				return logLevel;
		
		throw new IllegalArgumentException(String.format("Unknown log level name: %s: Valid leveles are %s", levelName, levelList(true)));
	}
	
	/**
	 * Get a comma-separated list of the log levels.
	 * 
	 * @param withAnd If <code>true</code>, show &quot;and&quot; before the last level name.
	 * 
	 * @return A list of log levels.
	 */
	public static String levelList(boolean withAnd) {
		StringBuilder builder = new StringBuilder();
		
		LogLevel[] values = LogLevel.values();
		
		builder.append(values[0]);
		
		int lastIndex = values.length;
		
		if (withAnd)
			lastIndex -= 1;
		
		for (int index = 1; index < lastIndex; ++index) {
			builder.append(", ");
			builder.append(values[index]);
		}
		
		if (withAnd) {
			builder.append(", and ");
			builder.append(values[lastIndex]);
		}
		
		return builder.toString();
	}
	
	/**
	 * Get a comma-separated list of the log levels.
	 * @return A list of log levels.
	 */
	public static String levelList() {
		return levelList(false);
	}
}
