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

package edu.gatech.kestrel.hapwriter;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.lang.reflect.Constructor;
import java.lang.reflect.InvocationTargetException;
import java.net.URLClassLoader;
import java.util.Arrays;
import java.util.Iterator;
import java.util.Set;

import edu.gatech.kanalyze.KAnalyzeConstants;
import edu.gatech.kanalyze.util.StringUtil;
import edu.gatech.kanalyze.util.SystemUtil;
import edu.gatech.kestrel.activeregion.Haplotype;
import edu.gatech.kestrel.io.StreamableOutput;
import edu.gatech.kestrel.refreader.ReferenceRegion;
import edu.gatech.kestrel.refreader.ReferenceSequence;


/**
 * Base class of all haplotype writers.
 */
public abstract class HaplotypeWriter {
	
	/** A descriptive name for this writer. */
	public final String name;
	
	/** Output file descriptor. */
	protected StreamableOutput output;
	
	/** Name of the sample currently being processed. */
	protected String sampleName;
	
	/** Current reference region. */
	protected ReferenceRegion referenceRegion;
	
	/** An array of reference sequences. */
	protected ReferenceSequence[] referenceSequenceArray;
	
	/** Class loader. */
	protected ClassLoader loader;
	
	/** Package where writer implementations should be found. */
	private static final String PACKAGE_NAME = "edu.gatech.kestrel.hapwriter";
	
	/** Sample name if it is not set. */
	public static final String DEFAULT_SAMPLE_NAME = "UnknownSample"; 
	
	/**
	 * Create a new haplotype writer.
	 * 
	 * @param name Name of this writer.
	 */
	protected HaplotypeWriter(String name) {
		
		// Normalize name
		if (name == null)
			name = "";
		
		name = name.trim().replaceAll("\\s", "_");
		
		if (name.isEmpty())
			name = "UnknownHaplotypeWriter";
		
		// Set fields
		this.name = name;
		output = StreamableOutput.STDOUT;
		sampleName = DEFAULT_SAMPLE_NAME;
		referenceRegion = null;
		
		return;
	}
	
	/**
	 * Add a resolved haplotype to this writer.
	 * 
	 * @param haplotype Haplotype to add.
	 * 
	 * @throws NullPointerException If <code>haplotype</code> is <code>null</code>.
	 */
	public abstract void add(Haplotype haplotype)
			throws NullPointerException;
	
	/**
	 * Implementation-defined initialization.
	 * 
	 * @param writerSpec Implementation defined specification.
	 * 
	 * @throws IllegalArgumentException If there is a problem with the format of
	 *   <code>writerSpec</code> (as defined by the implementation).
	 * @throws FileNotFoundException If the writer tries to open a file that cannot be found.
	 * @throws IOException If an IO error occurs opening resources.
	 */
	protected abstract void implInit(String writerSpec)
			throws IllegalArgumentException, FileNotFoundException, IOException;
	
	/**
	 * Get a short description of this writer.
	 * 
	 * @return A short description of this writer.
	 */
	public abstract String getDescription();
	
	/**
	 * Called when the sample name changes. The implementation may override this method.
	 */
	protected void sampleChanged() {
		return;
	}
	
	/**
	 * Called when the reference region changes. The implementation may override this method.
	 */
	protected void referenceRegionChanged() {
		return;
	}
	
	/**
	 * Writes all haplotypes to output. This should only be called once after all regions are added to the writer.
	 */
	public void flush() {
		
	}
	
	/**
	 * Initialize this writer.
	 * 
	 * @param writerSpec Writer specification. Implementation-defined string describing the behavior
	 *   of this writer.
	 * @param output Output location. If <code>null</code>, output is written to standard out.
	 * @param referenceSequenceArray Array of reference sequences.
	 * @param loader Class loader or <code>null</code> to use the default loader.
	 * 
	 * @throws NullPointerException If <code>referenceSequenceArray</code> is <code>null</code>.
	 * @throws IllegalArgumentException If there is a problem with the format of
	 *   <code>writerSpec</code> (as defined by the implementation), or if <code>referenceSequenceArray</code> contains
	 *   <code>null</code> references or is empty.
	 * @throws FileNotFoundException If the writer tries to open a file that cannot be found.
	 * @throws IOException If an IO error occurs opening resources.
	 * @throws SecurityException If a security error occurs opening a resource.
	 */
	public final void init(String writerSpec, StreamableOutput output, ReferenceSequence[] referenceSequenceArray, ClassLoader loader)
			throws NullPointerException, IllegalArgumentException, FileNotFoundException, IOException, SecurityException {
		
		ReferenceSequence[] refSeqs;
		
		// Check arguments
		if (writerSpec == null)
			writerSpec = "";
		
		writerSpec = writerSpec.trim();
		
		if (output == null)
			output = StreamableOutput.STDOUT;
		
		if (loader == null)
			loader = HaplotypeWriter.class.getClassLoader();
		
		if (referenceSequenceArray == null)
			throw new NullPointerException("Reference sequence array is null");
		
		if (referenceSequenceArray.length == 0)
			throw new IllegalArgumentException("Reference sequence array is empty");
		
		refSeqs = Arrays.copyOf(referenceSequenceArray, referenceSequenceArray.length);
		
		for (int index = 0; index < refSeqs.length; ++index)
			if (refSeqs[index] == null)
				throw new IllegalArgumentException(String.format("Reference sequence array contains a null reference at index %d", index));
		
		// Set fields
		this.loader = loader;
		this.output = output;
		this.referenceSequenceArray = refSeqs;
		
		referenceRegion = null;
		
		// Initialize implementation
		implInit(writerSpec);  // throws IllegalArgumentException, IOException, SecurityException
		
		return;
	}
	
	/**
	 * Set the name of the sample currently being written. When the sample changes, this method should
	 * be called.
	 * 
	 * @param sampleName Name of the sample. If <code>null</code> or empty, &quot;Unknown&quot;
	 *   is set.
	 */
	public void setSampleName(String sampleName) {
		if (sampleName == null)
			sampleName = "";
		
		sampleName = sampleName.trim();
		
		if (sampleName.isEmpty())
			sampleName = DEFAULT_SAMPLE_NAME;
		
		if (! sampleName.equals(this.sampleName)) {
			this.sampleName = sampleName;
			sampleChanged();
		}
		
		return;
	}
	
	/**
	 * Set the current reference region.
	 * 
	 * @param referenceRegion Reference region.
	 * 
	 * @throws NullPointerException If <code>referenceRegion</code> is <code>null</code>.
	 */
	public void setReferenceRegion(ReferenceRegion referenceRegion)
			throws NullPointerException {
		
		if (referenceRegion == null)
			throw new NullPointerException("Cannot set reference region: null");
		
		if (referenceRegion != this.referenceRegion) {
			this.referenceRegion = referenceRegion;
			
			referenceRegionChanged();
		}
		
		return;
	}
	
	/**
	 * Get a writer by its writer specification. The writer specification is a writer name
	 * followed by a colon and a list of arguments. Whitespace around the colon is
	 * discarded. The arguments are typically a comma-separated list of attribute/value pairs,
	 * however, it is up to the implementation to interpret the arguments.
	 *  
	 * @param writerSpec Writer specification.
	 * @param output Output location. If <code>null</code>, output is written to standard out.
	 * @param referenceSequenceArray Array of reference sequences.
	 * @param loader Class loader or <code>null</code> to use the system default loader.
	 * 
	 * @return An initialized writer.
	 * 
	 * @throws NullPointerException If <code>writerSpec</code> or <code>referenceSequenceArray</code> is <code>null</code>.
	 * @throws IllegalArgumentException If <code>writerSpec</code> is empty, does not contain
	 *   a valid writer name, if the arguments section contains errors, or if <code>referenceSequneceArray</code> is
	 *   empty or contains <code>null</code> references.
	 * @throws FileNotFoundException If the writer tries to open a file that cannot be found.
	 * @throws IOException If the writer tries to open some resource and an error occurs.
	 * @throws HaplotypeWriterInitException If any error occurs finding the writer class or creating
	 *   the object.
	 */
	public static HaplotypeWriter getWriter(String writerSpec, StreamableOutput output, ReferenceSequence[] referenceSequenceArray, ClassLoader loader)
			throws NullPointerException, IllegalArgumentException, FileNotFoundException, IOException, HaplotypeWriterInitException {
		
		String[] writerSpecTok;  // writerSpec split on ":"
		
		WriterClass writerClass;    // Class and ctor for the writer
		HaplotypeWriter hapWriter;  // Instantiated writer
		
		String writerName;
		String writerArgs;
		
		// Check arguments
		if (writerSpec == null)
			throw new NullPointerException("Cannot get haplotype writer with specification: null");
		
		writerSpec = writerSpec.trim();
		
		if (writerSpec.isEmpty())
			throw new IllegalArgumentException("Cannot get haplotype writer with an empty specification");
		
		// Split writer specification into the writer name and arguments
		writerSpecTok = writerSpec.split("\\s*:\\s*", 2);
		
		writerName = writerSpecTok[0];
		
		if (writerSpecTok.length > 1)
			writerArgs = writerSpecTok[1];
		else
			writerArgs = "";
		
		// Get writer class
		writerClass = getWriterClass(writerSpecTok[0], loader);  // throws IllegalArgumentException, ReaderInitException
		
		// Invoke constructor
		try {
			hapWriter = (HaplotypeWriter) (writerClass.writerCtor.newInstance(new Object[] {}));
			
		} catch (IllegalAccessException ex) {
			throw new HaplotypeWriterInitException (
					"Illegal access while invoking the haplotype writer constructor for writer name " + writerName +
					": (" + writerClass.writerClass.getName() + "): " + ex.getMessage(), ex);
			
		} catch (IllegalArgumentException ex) {
			throw new HaplotypeWriterInitException (
					"Illegal argument while invoking the haplotype writer constructor for writer name " + writerName +
					": (" + writerClass.writerClass.getName() + "): " + ex.getMessage(), ex);
			
		} catch (InstantiationException ex) {
			throw new HaplotypeWriterInitException (
					"Instantiation error while invoking the haplotype writer constructor for writer name " + writerName +
					": (" + writerClass.writerClass.getName() + "): " + ex.getMessage (), ex);
			
		} catch (InvocationTargetException ex) {
			throw new HaplotypeWriterInitException (
					"Invocation target error while invoking the haplotype writer constructor for writer name " + writerName +
					": (" + writerClass.writerClass.getName() + "): " + ex.getMessage (), ex);
			
		} catch (Throwable ex) {
			throw new HaplotypeWriterInitException(
					"Unknown error while invoking the haplotype writer constructor for writer name " + writerName +
					": (" + writerClass.writerClass.getName() + "): " + ex.getMessage(), ex);
		}
		
		// Initialize this writer
		hapWriter.init(writerArgs, output, referenceSequenceArray, loader);  // throws IllegalArgumentException, FileNotFoundException, IOException

		return hapWriter;
	}
	
	/**
	 * Get an haplotype writer class by name.
	 * 
	 * @param writerName Haplotype writer name.
	 * @param loader Loader for dynamic classes. If <code>null</code>, the defaulst class
	 *   loader is used. Beware that there are security and stability consequences of
	 *   loading classes from outside the Kestrel project.
	 * 
	 * @return An object with the haplotype writer class and constructor for the given format.
	 * 
	 * @throws NullPointerException If <code>writerName</code> is <code>null</code>.
	 * @throws IllegalArgumentException If <code>writerName</code> is not a recognized
	 *   format for a haplotype writer name.
	 * @throws HaplotypeWriterInitException If an error occurs while finding the writer class.
	 */
	public static WriterClass getWriterClass(String writerName, ClassLoader loader)
			throws NullPointerException, IllegalArgumentException, HaplotypeWriterInitException {
		
		String className;     // Fully qualified class name
		Class<?> writerClass; // Writer class
		Class<?> superClass;  // Super-class (for testing inheritance)
		
		Constructor<?> writerCtor; // Reader constructor
		
		// Check arguments
		if (writerName == null)
			throw new NullPointerException("Cannot get haplotype writer with name: null");
		
		if (! writerName.matches(KAnalyzeConstants.FORMAT_TYPE_PATTERN))
			throw new IllegalArgumentException("Haplotype writer name does not match regular expression \"" + KAnalyzeConstants.FORMAT_TYPE_PATTERN + "\": " + writerName);
		
		if (loader == null)
			loader = HaplotypeWriter.class.getClassLoader();
		
		// Convert format to class name
		className = PACKAGE_NAME + "." + writerName.toLowerCase() + "." + StringUtil.toNameCase(writerName) + "HaplotypeWriter";
		
		// Load class
		try {
			writerClass = Class.forName(className, true, loader);
			
		} catch (ClassNotFoundException ex) {
			throw new IllegalArgumentException("Cannot find class for haplotype writer: " + writerName + ": Searched " + className);
		}
		
		// Reader must extend this class
		superClass = writerClass.getSuperclass();
		
		while (superClass != null && superClass != HaplotypeWriter.class)
			superClass = superClass.getSuperclass();
		
		if (superClass != HaplotypeWriter.class) {
			throw new IllegalArgumentException(
					"Haplotype writer class for format " + writerName + " does not extend " +
							HaplotypeWriter.class.getName() + ": " + className);
		}
		
		// Get constructor
		try {
			// Try loading the ctor with a loader argument
			writerCtor = writerClass.getConstructor(); // throws NoSuchMethodException, SecurityException
			
		} catch (NoSuchMethodException ex) {
			throw new HaplotypeWriterInitException (
					"Cannot find default constructor for the haplotype writer with name " + writerName +
					": (" + className + "): " + ex.getMessage(), ex);
			
		} catch (SecurityException ex) {
			throw new HaplotypeWriterInitException (
					"Security error finding the default constructor for the haplotype writer with name " + writerName +
					": (" + className + "): " + ex.getMessage(), ex);
		}
		
		// Return class
		return new WriterClass(writerClass, writerCtor);
	}
	
	/**
	 * List any valid writers that can be found.
	 * 
	 * @param loader Loader to search for classes. If <code>null</code>, the system class
	 *   loader is used.
	 * 
	 * @return A sorted list of writer names.
	 */
	public static String[] listWriters(URLClassLoader loader) {
		
		// Create structures
		Set<String> urlSet = SystemUtil.findSubPackages(PACKAGE_NAME, loader.getURLs(), true);
		
		Iterator<String> setIter = urlSet.iterator();
		
		while (setIter.hasNext()) {
			String url = setIter.next();
			
			try {
				getWriterClass(url, loader);
				
			} catch (Throwable ex) {
				setIter.remove();
			}
		}
		
		// Sort and return
		String[] writers = urlSet.toArray(new String[0]);
		Arrays.sort(writers);
		
		return writers;
	}
	
	/**
	 * Get a description for a writer by name.
	 * 
	 * @param writerName Writer name.
	 * @param loader Class loader or <code>null</code> to use the default class loader.
	 * 
	 * @return Description or <code>null</code> if it could not be found for the writer
	 *   with name <code>writerName</code>.
	 */
	public static String getWriterDescription(String writerName, ClassLoader loader) {
		
		WriterClass writerClass;
		HaplotypeWriter arWriter;
		
		// Get class
		try {
			writerClass = getWriterClass(writerName, loader);
			
		} catch (Exception ex) {
			return null;
		}
		
		// Get writer
		try {
			arWriter = (HaplotypeWriter) (writerClass.writerCtor.newInstance(new Object[] {}));
			
		} catch (Exception ex) {
			return null;
		}
		
		return arWriter.getDescription();
	}
	
	/**
	 * Represents a writer class and constructor.
	 */
	public static class WriterClass {
		
		/** Writer class. */
		public final Class<?> writerClass;
		
		/** Writer constructor. */
		public final Constructor<?> writerCtor;
		
		/**
		 * Create a writer class.
		 * 
		 * @param writerClass Writer class.
		 * @param writerCtor Writer constructor.
		 * 
		 * @throws NullPointerException If any argument is <code>null</code>.
		 */
		public WriterClass(Class<?> writerClass, Constructor<?> writerCtor)
				throws NullPointerException {
			
			// Check arguments
			if (writerClass == null)
				throw new NullPointerException("Writer class is null");
			
			if (writerCtor == null)
				throw new NullPointerException("Writer constructor is null");
			
			// Set fields
			this.writerClass = writerClass;
			this.writerCtor = writerCtor;
			
			return;
		}
	}
}
