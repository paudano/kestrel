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

package edu.gatech.kestrel.writer.table;

import java.io.FileNotFoundException;
import java.io.PrintStream;

import edu.gatech.kestrel.variant.VariantCall;
import edu.gatech.kestrel.writer.VariantWriter;

/**
 * Write a tab-delimited table of variants.
 */
public class TableVariantWriter extends VariantWriter {
	
	/** Output writer. */
	private PrintStream outStream;
	
	/** Set to <code>true</code> once this writer is initialized. */
	private boolean isInit;
	
	/** Name of the current reference sequence. */
	private String referenceName;
	
	/** Name of the current region. */
	private String regionName;
	
	/** Default reference sequence name if none was set. */
	private static final String DEFAULT_REFERENCE_NAME = "UnknownRef";
	
	/** Default region name if none was set. */
	private static final String DEFAULT_REGION_NAME = "UnknownRegion";
	
	/**
	 * Create writer.
	 */
	public TableVariantWriter() {
		super("table");
		
		referenceName = DEFAULT_REFERENCE_NAME;
		regionName = DEFAULT_REGION_NAME;
		
		return;
	}
	
	/**
	 * Write variant.
	 * 
	 * 
	 */
	@Override
	public void writeVariant(VariantCall variant)
			throws NullPointerException {
		
		// Check arguments
		if (variant == null)
			throw new NullPointerException("Cannot write variant: null");
		
		// Write
		outStream.printf("%s\t%s\t%s\t%d\t%s\t%s\t%d\t%d\n",
				sampleName, referenceName, regionName,
				variant.start,
				variant.ref, variant.alt,
				variant.variantDepth, variant.locusDepth
		);
		
		return;
	}
	
	/**
	 * Initialize.
	 * 
	 * @throws FileNotFoundException If the output file cannot be found.
	 */
	@Override
	protected void implInit(String writerSpec)
			throws FileNotFoundException {
		
		referenceName = DEFAULT_REFERENCE_NAME;
		regionName = DEFAULT_REGION_NAME;
		
		// Close if initialized
		if (isInit) {
			outStream.close();
			outStream = null;
			
			isInit = false;
		}
		
		// Open output file
		outStream = new PrintStream(output.getStream());  // throws FileNotFoundException
		
		// Write header
		outStream.println("sample\treference\tregion\tlocus\tref\talt\tvar_depth\tregion_depth");
		
		// Flag init
		isInit = true;
		
		return;
	}
	
	/**
	 * Get a description of this writer.
	 * 
	 * @return Description of this writer.
	 */
	@Override
	public String getDescription() {
		return "Tab-delimited table of variant information";
	}
	
	/**
	 * Called when the reference region changes. The implementation may override this method.
	 */
	@Override
	protected void referenceRegionChanged() {
		
		referenceName = referenceRegion.referenceSequence.name;
		regionName = referenceRegion.name;
		
		return;
	}
}
