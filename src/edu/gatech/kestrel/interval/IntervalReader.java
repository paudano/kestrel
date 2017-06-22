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

package edu.gatech.kestrel.interval;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.lang.reflect.Constructor;
import java.lang.reflect.InvocationTargetException;
import java.net.URL;
import java.net.URLClassLoader;
import java.util.Arrays;
import java.util.Iterator;
import java.util.Set;
import java.util.regex.Pattern;
import java.util.regex.PatternSyntaxException;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import edu.gatech.kanalyze.KAnalyzeConstants;
import edu.gatech.kanalyze.util.StringUtil;
import edu.gatech.kanalyze.util.SystemUtil;

/**
 * Reads a set of intervals from a file.
 */
public abstract class IntervalReader {
	
	/** Logger object. */
	private static final Logger logger = LoggerFactory.getLogger(IntervalReader.class);;
	
	/** Name of this reader. */
	public final String name;
	
	/** Loader for dynamic class loading. */
	protected ClassLoader loader;
	
	/** Name of the package reader implementations are in. */
	private static final String PACKAGE_NAME = "edu.gatech.kestrel.interval";
	
	/**
	 * Create a new interval reader.
	 * 
	 * @param name Name of the reader implementation.
	 * 
	 * @throws NullPointerException If <code>name</code> is <code>null</code>
	 * @throws IllegalArgumentException If <code>name</code> is empty.
	 */
	public IntervalReader(String name)
			throws NullPointerException, IllegalArgumentException {
		
		// Set name
		if (name == null)
			throw new NullPointerException("Reader name is null");
		
		name = name.trim();
		
		if (name.isEmpty())
			throw new IllegalArgumentException("Reader name is empty");
		
		name = name.replaceAll("\\s+", "\t");
		
		this.name = name;
		
		return;
	}
	
	/**
	 * Read an interval file.
	 * 
	 * @param intervalFile File to read.
	 * 
	 * @return An array of intervals found in the file.
	 * 
	 * @throws NullPointerException If <code>intervalFile</code> is <code>null</code>.
	 * @throws FileNotFoundException If <code>intervalFile</code> is not found.
	 * @throws IOException If there is any error reading the file or with the format
	 *   of the file contents.
	 */
	public abstract RegionInterval[] read(File intervalFile)
			throws NullPointerException, FileNotFoundException, IOException;
	
	/**
	 * Get a description of this reader.
	 * 
	 * @return A description of this reader.
	 */
	public abstract String getDescription();
	
	/**
	 * Get a regular expression pattern to match file names that should be sent to this
	 * reader. When the format &quot;auto&quot; is used to automatically resolve the file type, this
	 * this pattern from each reader is used to determine how the file is parsed.
	 * 
	 * @return A regular expression pattern. If &quot;auto&quot; should never resolve a
	 *   file for this reader, <code>null</code> may be returned.
	 * 
	 * @throws IllegalArgumentException If pattern flags are invalid. If this is thrown, all patterns
	 *   are ignored.
	 * @throws PatternSyntaxException If the pattern is invalid. If this is thrown, all patterns are
	 *   ignored. 
	 */
	public Pattern getFilenamePattern()
			throws IllegalArgumentException, PatternSyntaxException {
		
		return null;
	}
	
	/**
	 * Initialize this reader with a set of reader-defined options. An error is thrown if the
	 * reader does not support the options. The  format of the argument string should be a comma-
	 * separated list of attribute-value pairs, but it is up to the implementation.
	 * 
	 * @param readerArgs Arguments to set. This value will never be <code>null</code> and whitespace
	 *   will be trimmed from the ends.
	 * @param loader Loader for dynamic classes. Argument is never <code>null</code>.
	 * 
	 * @throws IllegalArgumentException If an argument is invalid for any reason.
	 * @throws FileNotFoundException If the reader tries to open a file that cannot be found.
	 * @throws IOException If the reader tries to open some resource such as a file, database, or
	 *   network connection and an error occurs.
	 */
	protected void init(String readerArgs, ClassLoader loader)
			throws IllegalArgumentException, FileNotFoundException, IOException {
		
		// Declarations
		String[] tok;  // Attribute and value separated
		
		String attribute;  // Attribute
		String value;      // Value
		
		int avpCount = 0;  // Attribute for error reporting
		
		this.loader = loader;
		
		// Check arguments
		if (readerArgs == null)
			return;
		
		// Process attributes
		ATTRIBUTE_LOOP:
		for (String avpPair : readerArgs.trim().split("\\s*,\\s*")) {
			
			++avpCount;
			
			// Get attribute and value
			tok = avpPair.split("\\s*=\\s*");
			
			attribute = tok[0].toLowerCase();
			
			if (tok.length == 2) {
				value = tok[1].toLowerCase();
				
				if (value.isEmpty())
					value = null;
				
			} else {
				value = null;
			}
			
			if (attribute.isEmpty()) {
				if (value == null)
					continue ATTRIBUTE_LOOP;
				
				throw new IllegalArgumentException(String.format("%s interval reader: Attribute name at positon %d has no name", name, avpCount));
				
			} else if (attribute.matches(".*(\\s|\\p{C}).*")) {
				throw new IllegalArgumentException(String.format("%s interval reader: Attribute name at positon %d contains whitespace or non-printable characters: \"%s\"", name, avpCount, attribute));
			}
			
			// Set argument
			setReaderArg(attribute, value, avpCount);
		}
		
		// Finalize initialization
		readerInit();  // throws IllegalArgumentException, FileNotFoundException, IOException
		
		return;
	}
	
	/**
	 * Set an argument for this reader.
	 * 
	 * @param attribute Attribute. Never <code>null</code> or empty, and never contains whitespace
	 *   or non-printable characters. Whitespace is trimmed, and the name is converted to lower case.
	 * @param value Value. May be <code>null</code> if no value was found or was empty. Whitespace
	 *   is trimmed.
	 * @param avpCount The number of this attribute-value pair in a list of pairs. The first pair is
	 *   <code>1</code>, the second is <code>2</code>, etc. 
	 * 
	 * @throws IllegalArgumentException If the attribute is not supported or there is any error setting
	 *   the value.
	 */
	protected void setReaderArg(String attribute, String value, int avpCount)
			throws IllegalArgumentException {
		
		throw new IllegalArgumentException(String.format("%s interval reader: Unknown attribute (attribute/value pair %d): %d", name, avpCount, attribute));
	}
	
	/**
	 * Called at the end of <code>init()</code> after all options are set. 
	 * 
	 * @throws IllegalArgumentException If the initialization state of this reader is not valid,
	 *   such as a missing required argument.
	 * @throws FileNotFoundException If the reader tries to open a file that cannot be found.
	 * @throws IOException If the reader tries to open some resource such as a file, database, or
	 *   network connection and an error occurs.
	 */
	protected void readerInit()
			throws IllegalArgumentException, FileNotFoundException, IOException {
		
		return;
	}
	
	/**
	 * Read an interval file.
	 * 
	 * @param intervalFile File to read.
	 * @param readerName Name of the interval file reader. If <code>null</code> or empty,
	 *   &quot;auto&quot; is assumed.
	 * @param readerArgs Arguments to be passed to the reader.
	 * @param loader Loader for dynamic classes or <code>null</code> to use the system default
	 *   loader.
	 * 
	 * @return An array of intervals found in the file.
	 * 
	 * @throws NullPointerException If <code>intervalFile</code> is <code>null</code>.
	 * @throws IllegalArgumentException If the reader name is &quot;auto&quot; and the reader
	 *   cannot be resolved by file or if there is an error in <code>readerArgs</code>. 
	 * @throws FileNotFoundException If <code>intervalFile</code> is not found.
	 * @throws IOException If there is any error reading the file or with the format
	 *   of the file contents.
	 * @throws IntervalReaderInitException If any error occurs finding the reader class or creating
	 *   the reader object.
	 */
	public static RegionInterval[] read(File intervalFile, String readerName, String readerArgs, ClassLoader loader)
			throws NullPointerException, FileNotFoundException, IOException, IntervalReaderInitException {
		
		IntervalReader intervalReader;
		
		// Check arguments
		if (intervalFile == null)
			throw new NullPointerException("Cannot read interval file: null");
		
		if (readerName == null)
			readerName = "";
		
		readerName = readerName.trim().toLowerCase();
		
		if (readerName.isEmpty())
			readerName = "auto";
		
		if (loader == null)
			loader = IntervalReader.class.getClassLoader();
		
		// Resolve auto reader
		if (readerName.equals("auto")) {
			readerName = IntervalReader.resolveReaderNameByFile(intervalFile, loader);
			
			if (readerName == null)
				throw new IllegalArgumentException("Reader for file could not be automatically determined by the file name: " + intervalFile.getPath());
		}
		
		// Create reader
		intervalReader = IntervalReader.getReader(readerName, readerArgs, loader);  // throws IllegalArgumentException, FileNotFoundException, IOException, IntervalReaderInitException
		
		// Read intervals
		return intervalReader.read(intervalFile);
	}
	
	/**
	 * Get a reader by its name.
	 *  
	 * @param readerName Reader name.
	 * @param readerArgs A reader-defined set of initialization arguments.
	 * @param loader Class loader or <code>null</code> to use the system default loader.
	 * 
	 * @return An initialized interval reader.
	 * 
	 * @throws NullPointerException If <code>readerName</code> is <code>null</code>.
	 * @throws IllegalArgumentException If <code>readerName</code> is empty or does not contain
	 *   a valid reader name or if <code>readerArgs</code> contains errors.
	 * @throws FileNotFoundException If the reader tries to open a file that cannot be found.
	 * @throws IOException If the reader tries to open some resource such as a file, database, or
	 *   network connection and an error occurs.
	 * @throws IntervalReaderInitException If any error occurs finding the filter class or creating
	 *   the filter object.
	 */
	public static IntervalReader getReader(String readerName, String readerArgs, ClassLoader loader)
			throws NullPointerException, IllegalArgumentException, FileNotFoundException, IOException, IntervalReaderInitException {
		
		return getReader(readerName, readerArgs, loader, true);
	}
	
	/**
	 * Get a reader by its name.
	 *  
	 * @param readerName Reader name.
	 * @param readerArgs A reader-defined set of initialization arguments.
	 * @param loader Class loader or <code>null</code> to use the default loader.
	 * @param initReader Initialize reader if true. Otherwise, an uninitialized reader
	 *   is returned.
	 * 
	 * @return An initialized interval reader.
	 * 
	 * @throws NullPointerException If <code>readerName</code> is <code>null</code>.
	 * @throws IllegalArgumentException If <code>readerName</code> is empty or does not contain
	 *   a valid reader name or if <code>readerArgs</code> contains errors.
	 * @throws FileNotFoundException If the reader tries to open a file that cannot be found.
	 * @throws IOException If the reader tries to open some resource such as a file, database, or
	 *   network connection and an error occurs.
	 * @throws IntervalReaderInitException If any error occurs finding the reader class or creating
	 *   the reader object.
	 */
	private static IntervalReader getReader(String readerName, String readerArgs, ClassLoader loader, boolean initReader)
			throws NullPointerException, IllegalArgumentException, FileNotFoundException, IOException, IntervalReaderInitException {
		
		IntervalReaderClass readerClass;      // Class and ctor for the filter
		IntervalReader intervalReader;  // Instantiated variant filter
		
		if (readerArgs == null)
			readerArgs = "";
		
		readerArgs = readerArgs.trim();
		
		// Check arguments
		if (readerName == null)
			throw new NullPointerException("Cannot get interval reader with name: null");
		
		readerName = readerName.trim();
		
		if (readerName.isEmpty())
			throw new IllegalArgumentException("Cannot get interval reader with an empty name");
		
		if (loader == null)
			loader = IntervalReader.class.getClassLoader();
		
		// Get reader class
		readerClass = getReaderClass(readerName, loader);  // throws IllegalArgumentException, ReaderInitException
		
		// Invoke constructor
		try {
			intervalReader = (IntervalReader) (readerClass.readerCtor.newInstance(new Object[] {}));
			
		} catch (IllegalAccessException ex) {
			throw new IntervalReaderInitException (
					"Illegal access while invoking the variant filter constructor for filter name " + readerName +
					": (" + readerClass.readerClass.getName() + "): " + ex.getMessage(), ex);
			
		} catch (IllegalArgumentException ex) {
			throw new IntervalReaderInitException (
					"Illegal argument while invoking the variant filter constructor for filter name " + readerName +
					": (" + readerClass.readerClass.getName() + "): " + ex.getMessage(), ex);
			
		} catch (InstantiationException ex) {
			throw new IntervalReaderInitException (
					"Instantiation error while invoking the variant filter constructor for filter name " + readerName +
					": (" + readerClass.readerClass.getName() + "): " + ex.getMessage (), ex);
			
		} catch (InvocationTargetException ex) {
			throw new IntervalReaderInitException (
					"Invocation target error while invoking the variant filter constructor for filter name " + readerName +
					": (" + readerClass.readerClass.getName() + "): " + ex.getMessage (), ex);
			
		} catch (Throwable ex) {
			throw new IntervalReaderInitException(
					"Unknown error while invoking the variant filter constructor for filter name " + readerName +
					": (" + readerClass.readerClass.getName() + "): " + ex.getMessage(), ex);
		}
		
		// Initialize this reader
		if (initReader)
			intervalReader.init(readerArgs, loader);  // throws IllegalArgumentException, FileNotFoundException, IOException

		return intervalReader;
	}
	
	/**
	 * Get a reader class by name.
	 * 
	 * @param readerName Interval reader filter name.
	 * @param loader Loader for dynamic classes. If <code>null</code>, the default class
	 *   loader is used. Beware that there are security and stability consequences of
	 *   loading classes from outside the Kestrel project.
	 * 
	 * @return An object with the reader class and constructor for the given format.
	 * 
	 * @throws NullPointerException If <code>readerName</code> is <code>null</code>.
	 * @throws IllegalArgumentException If <code>readerName</code> is not a recognized
	 *   format for a reader name.
	 * @throws IntervalReaderInitException If an error occurs while finding the reader class.
	 */
	public static IntervalReaderClass getReaderClass(String readerName, ClassLoader loader)
			throws NullPointerException, IllegalArgumentException, IntervalReaderInitException {
		
		String className;     // Fully qualified class name
		Class<?> readerClass; // Filter class
		Class<?> superClass;  // Super-class (for testing inheritance)
		
		Constructor<?> readerCtor; // Reader constructor
		
		// Check arguments
		if (readerName == null)
			throw new NullPointerException("Cannot get interval reader with name: null");
		
		if (! readerName.matches(KAnalyzeConstants.FORMAT_TYPE_PATTERN))
			throw new IllegalArgumentException("Interval reader name does not match regular expression \"" + KAnalyzeConstants.FORMAT_TYPE_PATTERN + "\": " + readerName);
		
		if (loader == null)
			loader = IntervalReader.class.getClassLoader();
		
		// Convert format to class name
		className = PACKAGE_NAME + "." + readerName.toLowerCase() + "." + StringUtil.toNameCase(readerName) + "IntervalReader";
		
		// Load class
		try {
			readerClass = Class.forName(className, true, loader);
			
		} catch (ClassNotFoundException ex) {
			throw new IllegalArgumentException("Cannot find class for interval reader: " + readerName + ": Searched " + className);
		}
		
		// Reader must extend this class
		superClass = readerClass.getSuperclass();
		
		while (superClass != null && superClass != IntervalReader.class)
			superClass = superClass.getSuperclass();
		
		if (superClass != IntervalReader.class) {
			throw new IllegalArgumentException(
					"Interval reader class for format " + readerName + " does not extend " +
							IntervalReader.class.getName() + ": " + className);
		}
		
		// Get constructor
		try {
			// Try loading the ctor with a loader argument
			readerCtor = readerClass.getConstructor(); // throws NoSuchMethodException, SecurityException
			
		} catch (NoSuchMethodException ex) {
			throw new IntervalReaderInitException (
					"Cannot find default constructor for the interval reader with name " + readerName +
					": (" + className + "): " + ex.getMessage(), ex);
			
		} catch (SecurityException ex) {
			throw new IntervalReaderInitException (
					"Security error finding the default constructor for the interval reader with name " + readerName +
					": (" + className + "): " + ex.getMessage(), ex);
		}
		
		// Return class
		return new IntervalReaderClass(readerClass, readerCtor);
	}
	
	/**
	 * List any valid readers that can be found.
	 * 
	 * @param loader Loader to search for classes. If <code>null</code>, the system class
	 *   loader is used.
	 * 
	 * @return A sorted list of reader names.
	 */
	public static String[] listReaders(URLClassLoader loader) {
		
		// Create structures
		Set<String> urlSet = SystemUtil.findSubPackages(PACKAGE_NAME, loader.getURLs(), true);
		
		Iterator<String> setIter = urlSet.iterator();
		
		while (setIter.hasNext()) {
			String url = setIter.next();
			
			try {
				getReaderClass(url, loader);
				
			} catch (Throwable ex) {
				setIter.remove();
			}
		}
		
		// Sort and return
		String[] readers = urlSet.toArray(new String[0]);
		Arrays.sort(readers);
		
		return readers;
	}
	
	/**
	 * Get a description for a reader by name.
	 * 
	 * @param readerName Reader name.
	 * @param loader Class loader or <code>null</code> to use the default class loader.
	 * 
	 * @return Description or <code>null</code> if it could not be found for the reader
	 *   with name <code>filterName</code>.
	 */
	public static String getReaderDescription(String readerName, ClassLoader loader) {
		
		IntervalReaderClass readerClass;
		IntervalReader reader;
		
		// Get class
		try {
			readerClass = getReaderClass(readerName, loader);
			
		} catch (Exception ex) {
			return null;
		}
		
		// Get filter
		try {
			reader = (IntervalReader) (readerClass.readerCtor.newInstance(new Object[] {}));
			
		} catch (Exception ex) {
			return null;
		}
		
		return reader.getDescription();
	}
	
	/**
	 * Get the description of the files processed by a reader.
	 * 
	 * @param readerName Reader name.
	 * @param loader Class loader. If <code>null</code>, the system class loader is used. 
	 * 
	 * @return A pattern for the the files processed by the reader for <code>readerName</code>,
	 *   or <code>null</code> if the reader does not exist, does not have a pattern, or a reader
	 *   object could not be created and queried.
	 * 
	 * @throws NullPointerException If <code>readerName</code> is <code>null</code>.
	 * @throws IllegalArgumentException If <code>readerName</code> is not a valid name or the name
	 *   contains an invalid regular expression flag.
	 * @throws PatternSyntaxException If the format pattern is not valid regular expression. 
	 */
	public static Pattern getFormatPattern(String readerName, ClassLoader loader)
			throws NullPointerException, IllegalArgumentException, PatternSyntaxException {
		
		// Get reader
		try {
			return getReader(readerName, null, loader, false).getFilenamePattern();  // throws NullPointerException, IllegalArgumentException, PatternSyntaxException
			
		} catch (Exception ex) {
			return null;
		}
	}
	
	/**
	 * Resolve reader by file name.
	 * 
	 * @param intervalFile File.
	 * @param loader Class loader or <code>null</code> to use the default class loader.
	 * 
	 * @return Reader name or <code>null</code> if the reader could not be resolved.
	 * 
	 * @throws NullPointerException If <code>intervalFile</code> is <code>null</code>.
	 */
	public static String resolveReaderNameByFile(File intervalFile, ClassLoader loader)
			throws NullPointerException {
		
		// Declarations
		URLClassLoader urlLoader;  // URL class loader
		
		// Check arguments
		if (intervalFile == null)
			throw new NullPointerException("Cannot resolve reader name by file: null");
		
		if (loader == null)
			loader = IntervalReader.class.getClassLoader();
		
		// Set local URL loader
		if (loader instanceof URLClassLoader) {
			urlLoader = (URLClassLoader) loader;
		
		} else {
			urlLoader = new URLClassLoader(new URL[0], loader);
		}
		
		// Search each available format
		for (String readerName : IntervalReader.listReaders(urlLoader)) {
			
			try {
				Pattern pattern = IntervalReader.getFormatPattern(readerName, urlLoader);
				
				if (pattern != null && pattern.matcher(intervalFile.getName()).find())
					return readerName;
					
			} catch (PatternSyntaxException ex) {
				logger.warn("Error getting pattern for reader with name {}: {} (PatternSyntaxException)", readerName, ex.getMessage());
			}
		}
		
		// No reader found
		return null;
	}
	
	/**
	 * Represents a reader class and constructor.
	 */
	public static class IntervalReaderClass {
		
		/** Reader class. */
		public final Class<?> readerClass;
		
		/** Reader constructor. */
		public final Constructor<?> readerCtor;
		
		/**
		 * Create a reader class.
		 * 
		 * @param readerClass Reader class.
		 * @param readerCtor Reader constructor.
		 * 
		 * @throws NullPointerException If any argument is <code>null</code>.
		 */
		public IntervalReaderClass(Class<?> readerClass, Constructor<?> readerCtor)
				throws NullPointerException {
			
			// Check arguments
			if (readerClass == null)
				throw new NullPointerException("Reader class is null");
			
			if (readerCtor == null)
				throw new NullPointerException("Reader constructor is null");
			
			// Set fields
			this.readerClass = readerClass;
			this.readerCtor = readerCtor;
			
			return;
		}
	}
}
