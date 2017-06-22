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

package edu.gatech.kestrel.varfilter;

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
import edu.gatech.kestrel.variant.VariantCall;

/**
 * Super class of all variant callers.
 */
public abstract class VariantFilter {
	
	/** Name of this filter. */
	public final String name;
	
	/** Name of the package filter classes must be in. */
	public static final String PACKAGE_NAME = "edu.gatech.kestrel.varfilter";
	
	/**
	 * Create a new variant call filter.
	 * 
	 * @param name Name of this filter.
	 */
	public VariantFilter(String name) {
		
		// Normalize name
		if (name == null)
			name = "";
		
		name = name.trim().replaceAll("\\s", "_");
		
		if (name.isEmpty())
			name = "UnknownFilter";
		
		// Set fields
		this.name = name;
		
		return;
	}
	
	/**
	 * Initialize this filter. This method must be called before <code>filter()</code>.
	 * 
	 * @param argStr A filter-defined argument list. This argument will never be <code>null</code>. 
	 * 
	 * @throws IllegalArgumentException If <code>argStr</code> is invalid for this filter.
	 * @throws FileNotFoundException If the filter must open a file that cannot be found.
	 * @throws IOException If any IO error occurs, such as reading a file, opening a database, or
	 *   making a network connection.
	 */
	protected abstract void init(String argStr)
			throws IllegalArgumentException, FileNotFoundException, IOException;
	
	/**
	 * Get a description of this filter.
	 * 
	 * @return A description of this filter.
	 */
	public abstract String getDescription();
	
	/**
	 * Perform filtering on a variant.
	 * 
	 * @param variantCall Variant call to filter.
	 * 
	 * @return A variant call object if the call was not filtered out, or <code>null</code> if this
	 *   variant must be discarded. Note that filters may re-write variants, so the returned object
	 *   may not be the same as the variant call argument.
	 */
	public abstract VariantCall filter(VariantCall variantCall);
	
	/**
	 * Get a filter by its filter specification. The filter specification is a filter name
	 * followed by a colon and a list of arguments. Whitespace around the colon is
	 * discarded. The arguments are typically a comma-separated list of attribute/value pairs,
	 * however, it is up to the filter implementation to interpret the arguments.
	 *  
	 * @param filterSpec Filter specification.
	 * @param loader Class loader or <code>null</code> to use the system default loader.
	 * 
	 * @return An initialized filter.
	 * 
	 * @throws NullPointerException If <code>filterSpec</code> is <code>null</code>.
	 * @throws IllegalArgumentException If <code>filterSpec</code> is empty, does not contain
	 *   a valid filter name, or if the arguments section contains errors.
	 * @throws FileNotFoundException If the filter tries to open a file that cannot be found.
	 * @throws IOException If the filter tries to open some resource such as a file, database, or
	 *   network connection and an error occurs.
	 * @throws VariantFilterInitException If any error occurs finding the filter class or creating
	 *   the filter object.
	 */
	public static VariantFilter getFilter(String filterSpec, ClassLoader loader)
			throws NullPointerException, IllegalArgumentException, FileNotFoundException, IOException, VariantFilterInitException {
		
		String[] filterSpecTok;  // filterSpec split on ":"
		
		FilterClass filterClass;      // Class and ctor for the filter
		VariantFilter variantFilter;  // Instantiated variant filter
		
		String filterName;
		String filterArgs;
		
		// Check arguments
		if (filterSpec == null)
			throw new NullPointerException("Cannot get variant filter with specification: null");
		
		filterSpec = filterSpec.trim();
		
		if (filterSpec.isEmpty())
			throw new IllegalArgumentException("Cannot get variant filter with an empty specification");
		
		// Split filter specification into the filter name and arguments
		filterSpecTok = filterSpec.split("\\s*:\\s*", 2);
		
		filterName = filterSpecTok[0];
		
		if (filterSpecTok.length > 1)
			filterArgs = filterSpecTok[1];
		else
			filterArgs = "";
		
		// Get reader class
		filterClass = getFilterClass(filterSpecTok[0], loader);  // throws IllegalArgumentException, ReaderInitException
		
		// Invoke constructor
		try {
			variantFilter = (VariantFilter) (filterClass.filterCtor.newInstance(new Object[] {}));
			
		} catch (IllegalAccessException ex) {
			throw new VariantFilterInitException (
					"Illegal access while invoking the variant filter constructor for filter name " + filterName +
					": (" + filterClass.filterClass.getName() + "): " + ex.getMessage(), ex);
			
		} catch (IllegalArgumentException ex) {
			throw new VariantFilterInitException (
					"Illegal argument while invoking the variant filter constructor for filter name " + filterName +
					": (" + filterClass.filterClass.getName() + "): " + ex.getMessage(), ex);
			
		} catch (InstantiationException ex) {
			throw new VariantFilterInitException (
					"Instantiation error while invoking the variant filter constructor for filter name " + filterName +
					": (" + filterClass.filterClass.getName() + "): " + ex.getMessage (), ex);
			
		} catch (InvocationTargetException ex) {
			throw new VariantFilterInitException (
					"Invocation target error while invoking the variant filter constructor for filter name " + filterName +
					": (" + filterClass.filterClass.getName() + "): " + ex.getMessage (), ex);
			
		} catch (Throwable ex) {
			throw new VariantFilterInitException(
					"Unknown error while invoking the variant filter constructor for filter name " + filterName +
					": (" + filterClass.filterClass.getName() + "): " + ex.getMessage(), ex);
		}
		
		// Initialize this filter
		variantFilter.init(filterArgs);  // throws IllegalArgumentException, FileNotFoundException, IOException

		return variantFilter;
	}
	
	/**
	 * Get a variant filter class by name.
	 * 
	 * @param filterName Variant filter name.
	 * @param loader Loader for dynamic classes. If <code>null</code>, the defaulte class
	 *   loader is used. Beware that there are security and stability consequences of
	 *   loading classes from outside the Kestrel project.
	 * 
	 * @return An object with the variant filter class and constructor for the given format.
	 * 
	 * @throws NullPointerException If <code>filterName</code> is <code>null</code>.
	 * @throws IllegalArgumentException If <code>filterName</code> is not a recognized
	 *   format for a variant filter name.
	 * @throws VariantFilterInitException If an error occurs while finding the filter class.
	 */
	public static FilterClass getFilterClass(String filterName, ClassLoader loader)
			throws NullPointerException, IllegalArgumentException, VariantFilterInitException {
		
		String className;     // Fully qualified class name
		Class<?> filterClass; // Filter class
		Class<?> superClass;  // Super-class (for testing inheritance)
		
		Constructor<?> filterCtor; // Reader constructor
		
		// Check arguments
		if (filterName == null)
			throw new NullPointerException("Cannot get filter with name: null");
		
		if (! filterName.matches(KAnalyzeConstants.FORMAT_TYPE_PATTERN))
			throw new IllegalArgumentException("Filter name does not match regular expression \"" + KAnalyzeConstants.FORMAT_TYPE_PATTERN + "\": " + filterName);
		
		if (loader == null)
			loader = VariantFilter.class.getClassLoader();
		
		// Convert format to class name
		className = PACKAGE_NAME + "." + filterName.toLowerCase() + "." + StringUtil.toNameCase(filterName) + "VariantFilter";
		
		// Load class
		try {
			filterClass = Class.forName(className, true, loader);
			
		} catch (ClassNotFoundException ex) {
			throw new IllegalArgumentException("Cannot find class for variant filter: " + filterName + ": Searched " + className);
		}
		
		// Reader must extend this class
		superClass = filterClass.getSuperclass();
		
		while (superClass != null && superClass != VariantFilter.class)
			superClass = superClass.getSuperclass();
		
		if (superClass != VariantFilter.class) {
			throw new IllegalArgumentException(
					"Variant filter class for format " + filterName + " does not extend " +
							VariantFilter.class.getName() + ": " + className);
		}
		
		// Get constructor
		try {
			// Try loading the ctor with a loader argument
			filterCtor = filterClass.getConstructor(); // throws NoSuchMethodException, SecurityException
			
		} catch (NoSuchMethodException ex) {
			throw new VariantFilterInitException (
					"Cannot find default constructor for the variant filter with name " + filterName +
					": (" + className + "): " + ex.getMessage(), ex);
			
		} catch (SecurityException ex) {
			throw new VariantFilterInitException (
					"Security error finding the default constructor for the variant filter with name " + filterName +
					": (" + className + "): " + ex.getMessage(), ex);
		}
		
		// Return class
		return new FilterClass(filterClass, filterCtor);
	}
	
	/**
	 * List any valid filters that can be found.
	 * 
	 * @param loader Loader to search for classes. If <code>null</code>, the system class
	 *   loader is used.
	 * 
	 * @return A sorted list of filter names.
	 */
	public static String[] listFilters(URLClassLoader loader) {
		
		// Create structures
		Set<String> urlSet = SystemUtil.findSubPackages(PACKAGE_NAME, loader.getURLs(), true);
		
		Iterator<String> setIter = urlSet.iterator();
		
		while (setIter.hasNext()) {
			String url = setIter.next();
			
			try {
				getFilterClass(url, loader);
				
			} catch (Throwable ex) {
				setIter.remove();
			}
		}
		
		// Sort and return
		String[] filters = urlSet.toArray(new String[0]);
		Arrays.sort(filters);
		
		return filters;
	}
	
	/**
	 * Get a description for a filter by name.
	 * 
	 * @param filterName Filter name.
	 * @param loader Class loader or <code>null</code> to use the default class loader.
	 * 
	 * @return Description or <code>null</code> if it could not be found for the filter
	 *   with name <code>filterName</code>.
	 */
	public static String getFilterDescription(String filterName, ClassLoader loader) {
		
		FilterClass filterClass;
		VariantFilter variantFilter;
		
		// Get class
		try {
			filterClass = getFilterClass(filterName, loader);
			
		} catch (Exception ex) {
			return null;
		}
		
		// Get filter
		try {
			variantFilter = (VariantFilter) (filterClass.filterCtor.newInstance(new Object[] {}));
			
		} catch (Exception ex) {
			return null;
		}
		
		return variantFilter.getDescription();
	}
	
	/**
	 * Represents a filter class and constructor.
	 */
	public static class FilterClass {
		
		/** Filter class. */
		public final Class<?> filterClass;
		
		/** Filter constructor. */
		public final Constructor<?> filterCtor;
		
		/**
		 * Create a filter class.
		 * 
		 * @param filterClass Filter class.
		 * @param filterCtor Filter constructor.
		 * 
		 * @throws NullPointerException If any argument is <code>null</code>.
		 */
		public FilterClass(Class<?> filterClass, Constructor<?> filterCtor)
				throws NullPointerException {
			
			// Check arguments
			if (filterClass == null)
				throw new NullPointerException("Filter class is null");
			
			if (filterCtor == null)
				throw new NullPointerException("Filter constructor is null");
			
			// Set fields
			this.filterClass = filterClass;
			this.filterCtor = filterCtor;
			
			return;
		}
	}
}
