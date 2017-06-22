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

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;

import edu.gatech.kanalyze.comp.reader.FileSequenceSource;
import edu.gatech.kanalyze.io.ikc.IkcReader;
import edu.gatech.kanalyze.module.count.CountModule;
import edu.gatech.kanalyze.util.kmer.KmerUtil;

/**
 * A k-mer counter that reads from indexed k-mer counts stored in a file.
 */
public class IkcCountMap extends CountMap {
	
	/** K-mer count reader. */
	private IkcReader reader;
	
	/** Temporary file directory. */
	private File tempDir;
	
	/** Temporary file with k-mer counts. */
	private File tempFile;
	
	/** Remove <code>tempFile</code> if <code>true</code>. */
	private boolean rmLastTemp;
	
	/** Remove the temporary file after each sample if <code>true</code>. */
	private final boolean rmTemp;
	
	/** Suffix for indexed count temporary files. */
	public static final String IKC_TEMP_FILE_SUFFIX = ".ikc";
	
	/** Indexed k-mer count file pattern. Used to determine an IKC file when type is &quot;auto&quot;.*/
	public static final String IKC_FILE_PATTERN = ".*\\.ikc";
	
	/**
	 * Create an in-memory count mapper.
	 * 
	 * @param kUtil K-mer utility.
	 * @param countModule A configured count module.
	 * @param tempDir Temporary directory or <code>null</code> to use the current
	 *   working directory. The k-mer count file will be created in this location.
	 * @param rmTemp If the indexed k-mer count file (ikc) should be removed after each
	 *   sample is processed.
	 * 
	 * @throws NullPointerException If <code>kUtil</code> or <code>countModule</code>
	 *   is <code>null</code>.
	 * @throws IllegalArgumentException If <code>kUtil</code> is not configured with a minimizer size.
	 * @throws IOException If an IO error occurs creating the temporary file.
	 */
	public IkcCountMap(KmerUtil kUtil, CountModule countModule, File tempDir, boolean rmTemp)
			throws NullPointerException, IllegalArgumentException, IOException {
		
		super(kUtil, countModule);  // throws NullPointerException
		
		if (kUtil.kMinSize == 0)
			throw new IllegalArgumentException("K-mer utility was not configured with a minimizer size");
		
		if (tempDir == null)
			tempDir = new File(".");
		
		// Set fields
		this.rmTemp = rmTemp;
		this.tempDir = tempDir;
		
		reader = null;
		tempFile = null;
		
		// Configure module
		countModule.setOutputFormat("ikc");
		
		countModule.setMinimizerSize(kUtil.kMinSize);
		countModule.setMinimizerMask(kUtil.kMinMask);
		
		return;
	}
	
	/**
	 * Setup the count module with the output file.
	 * 
	 * @throws IOException If there is an IO error opening the temporary file.
	 * 
	 * @return <code>true</code> to run the count module, and <code>false</code> if the group
	 *   contains one IKC file that can be opened directly.
	 */
	@Override
	public boolean preModuleRun()
			throws IOException {
		
		// Remove temporary file
		if (tempFile != null && rmLastTemp) {
			tempFile.delete();
			tempFile = null;
		}
		
		// Check for a sample with one file
		if (sample.sources.length == 1 &&
			sample.sources[0] instanceof FileSequenceSource) {
			
			File ikcFile = ((FileSequenceSource) sample.sources[0]).file;
			
			// If the file is an IKC file, read it and skip the k-mer counter
			if (sample.sources[0].formatType.equals("auto") &&
				ikcFile.getName().matches(IKC_FILE_PATTERN) ||
				sample.sources[0].formatType.equals("ikc")) {
				
				tempFile = ikcFile;
				rmLastTemp = false;
				
				return false;
			}
		}
		
		// Create temporary file
		tempFile = File.createTempFile("IKC_" + sample.name + "_", IKC_TEMP_FILE_SUFFIX, tempDir);  // throws IOException
		
		if (rmTemp) {
			tempFile.deleteOnExit();
			rmLastTemp = true;
		}
		
		countModule.setOutput(tempFile);
		
		return true;
	}
	
	/**
	 * Called after the sample is set on the count module, but before it is run.
	 * 
	 * @param onError <code>true</code> if the module or <code>preModuleRun()</code>
	 *   threw an exception.
	 * @param aborted <code>true</code> if the pipeline is stopped because <code>abort()</code>
	 *   was called.
	 * 
	 * @throws FileNotFoundException If the IKC file does not exist.
	 * @throws SecurityException If the IKC file cannot be opened.
	 * @throws IOException If there is an IO error reading the IKC file or
	 *   it does not contain a known or properly formatted header.
	 */
	@Override
	protected void postModuleRun(boolean onError, boolean aborted)
			throws FileNotFoundException, SecurityException, IOException {
		
		// Stop if error or aborted
		if (onError || aborted)
			return;
		
		reader = new IkcReader(kUtil, tempFile);
		
		return;
	}
	
	/**
	 * Get a k-mer from this map.
	 * 
	 * @see edu.gatech.kestrel.counter.CountMap#get(int[])
	 */
	@Override
	public int get(int[] kmer) throws NullPointerException, IllegalArgumentException {
		
		if (reader == null)
			throw new IllegalStateException("Attempted to get a k-mer count before a sample was set for this map");
		
		return reader.get(kmer);
	}
}
