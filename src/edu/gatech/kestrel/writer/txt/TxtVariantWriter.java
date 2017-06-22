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

package edu.gatech.kestrel.writer.txt;

import java.io.FileNotFoundException;
import java.io.PrintStream;

import edu.gatech.kestrel.refreader.ReferenceRegion;
import edu.gatech.kestrel.variant.VariantCall;
import edu.gatech.kestrel.writer.VariantWriter;

/**
 * Outputs variants in a plain-text format.
 */
public class TxtVariantWriter extends VariantWriter {
	
	/** Output writer. */
	private PrintStream outStream;
	
	/** Set to <code>true</code> once this writer is initialized. */
	private boolean isInit;
	
	/** True when the first sample name is written. */
	private boolean firstSample;
	
	/** Last reference region before it was changed. */
	private ReferenceRegion lastRegion;
	
	/**
	 * Create a text variant writer.
	 */
	public TxtVariantWriter() {
		super("txt");
		
		outStream = null;
		isInit = false;
		firstSample = true;
		lastRegion = null;
		
		return;
	}
	
	/**
	 * Write a variant to output.
	 */
	@Override
	public void writeVariant(VariantCall variant) {
		
		// Check arguments
		if (variant == null)
			throw new NullPointerException("Cannot write variant: null");
		
		// Write variant
		outStream.printf("%s (%d/%d)\n", variant.toString(), variant.variantDepth, variant.locusDepth);
		
		return;
	}

	/**
	 * Initialize this writer.
	 * 
	 * @param writerSpec Ignored.
	 * 
	 * @throws FileNotFoundException If the output file cannot be found.
	 */
	@Override
	protected void implInit(String writerSpec)
			throws FileNotFoundException {
		
		// Close if initialized
		if (isInit) {
			outStream.close();  // throws IOException
			outStream = null;
			
			firstSample = true;
			
			isInit = false;
		}
		
		// Open output file
		outStream = new PrintStream(output.getStream());  // throws FileNotFoundException
		
		// Flag init
		isInit = true;
		lastRegion = null;
		
		return;
	}

	@Override
	public String getDescription() {
		
		return "Plain-text list of variants";
	}
	
	/**
	 * Called when the sample name changes. The implementation may override this method.
	 */
	@Override
	protected void sampleChanged() {
		
		// Write newlines separating samples
		if (firstSample) {
			firstSample = false;
			
		} else {
			outStream.println();
			outStream.println();
		}
		
		outStream.println("** Sample: " + sampleName);
		
		return;
	}
	
	@Override
	protected void referenceRegionChanged() {
		
		if (byRegion) {
			outStream.println("\n* Reference region: " + referenceRegion.name);
			
		} else {
			
			if (lastRegion == null || lastRegion.referenceSequence != referenceRegion.referenceSequence)
				outStream.println("\n* Reference: " + referenceRegion.referenceSequence.name);
		}
		
		lastRegion = referenceRegion;
		
		return;
	}
	
	/**
	 * Writes all variants to output.
	 */
	@Override
	public void flush() {
		
		if (isInit) {
			
			outStream.flush();
			outStream.close();
			
			outStream = null;
			lastRegion = null;
			
			isInit = false;
		}
	}
}
