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

import edu.gatech.kanalyze.KAnalyzeConstants;

/**
 * Constants for the Kestrel system.
 */
public class KestrelConstants {
	
	//
	// Program basics
	//
	
	/** Name of this program. */
	public static final String PROG_NAME = "kestrel";
	
	//
	// Version constants
	//
	
	/** Major version number. */
	public static final int VERSION_MAJOR = 1;
	
	/** Minor version number. */
	public static final int VERSION_MINOR = 0;
	
	/** Revision version numbers. */
	public static final int VERSION_REVISION = 1;
	
	/** Development version. If not <code>0</code>, this is not a release version. */
	public static final int VERSION_DEV = 0;
	
	/** Full version string. */
	public static final String VERSION = VERSION_MAJOR + "." + VERSION_MINOR + "." + VERSION_REVISION + ((VERSION_DEV == 0) ? "" : "dev" + VERSION_DEV);
	
	
	// Limits
	
	/** Most implementations fail to create arrays larger than this size. */
	public static final int MAX_ARRAY_SIZE = Integer.MAX_VALUE - 8;
	
	/** Minimum k-mer size. */
	public static final int MIN_KMER_SIZE = 4;
	
	/** Dynamic arrays are expanded by this factor by default. */
	public static final float ARRAY_EXPAND_FACTOR = 1.5F;
	
	
	// Return codes
	//
	// Constants starting with "ERR_" will be returned by main () when a
	// program exits. Non-zero return codes indicate an error, and scripts
	// or the command line can use that return code to determine if the
	// command completed successfully.

	/** Error: None. No error occurred. */
	public static final int ERR_NONE = 0;

	/** Indicates a usage error. An argument to the program is incorrect. */
	public static final int ERR_USAGE = 1;

	/** Indicates an I/O error. Reading or writing data failed. */
	public static final int ERR_IO = 2;
	
	/** Indicates a security or permissions error. */
	public static final int ERR_SECURITY = 3;

	/** Indicates that a specified file was not found. */
	public static final int ERR_FILENOTFOUND = 4;

	/**
	 * Indicates that data was malformed. For example, if a data file is corrupt
	 * or invalid, this error would be returned.
	 */
	public static final int ERR_DATAFORMAT = 5;

	/** Analysis was not able to complete for some reason. */
	public static final int ERR_ANALYSIS = 6;
	
	/** Indicates a thread was interrupted before it completed. */
	public static final int ERR_INTERRUPTED = 7;
	
	/** Indicates that a limit was exceeded. */
	public static final int ERR_LIMITS = 8;
	
	/** Termination was requested by an abort. */
	public static final int ERR_ABORT = 98;
	
	/** An unexpected system error occurred. These should be reported as a bug. */
	public static final int ERR_SYSTEM = 99;
	
	
	//
	// Resource Locator Constants
	//
	
	/** Root of all resources. */
	public static final String RESOURCE_ROOT = "edu/gatech/kestrel";
	
	/** Location of test resources and files. */
	public static final String RESOURCE_TEST = RESOURCE_ROOT + "/test";
	
	
	// Formats
	
	/** Format of input file type names. */
	public static final String FORMAT_TYPE_PATTERN = KAnalyzeConstants.FORMAT_TYPE_PATTERN;
}
