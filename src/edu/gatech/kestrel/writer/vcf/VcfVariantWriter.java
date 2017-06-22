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

package edu.gatech.kestrel.writer.vcf;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintStream;
import java.util.Iterator;

import edu.gatech.kestrel.KestrelConstants;
import edu.gatech.kestrel.refreader.ReferenceSequence;
import edu.gatech.kestrel.variant.VariantCall;
import edu.gatech.kestrel.writer.VariantWriter;

/**
 * Writer for VCF files.
 */
public class VcfVariantWriter extends VariantWriter {
	
	/** A container of VCF records. */
	private VcfRecordContainer recordContainer;
	
	/** Output writer. */
	private PrintStream outStream;
	
	/** Set to <code>true</code> once this writer is initialized. */
	private boolean isInit;
	
	/**
	 * Create a new VCF writer.
	 */
	public VcfVariantWriter() {
		super("VCF");
		
		recordContainer = new VcfRecordContainer();
		
		return;
	}
	
	/**
	 * Called when the sample name changes. The implementation may override this method.
	 */
	@Override
	protected void sampleChanged() {
		recordContainer.newSample(sampleName);
		
		return;
	}
	
	/**
	 * Add a variant to this writer.
	 * 
	 * @param variant Variant to add.
	 * 
	 * @throws NullPointerException If <code>variant</code> is <code>null</code>.
	 */
	@Override
	public void writeVariant(VariantCall variant)
			throws NullPointerException {
		
		recordContainer.addVariant(variant);
		
		return;
	}

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
	@Override
	protected void implInit(String writerSpec) throws IllegalArgumentException, FileNotFoundException, IOException {
		
		recordContainer.clear();
		
		// Close if initialized
		if (isInit) {
			outStream.close();
			outStream = null;
			
			isInit = false;
		}
		
		// Open output file
		outStream = new PrintStream(output.getStream());  // throws FileNotFoundException
		
		// Write header
		outStream.println("##fileformat=VCF4.2");
		outStream.println(String.format("##source=Kestrel%s", KestrelConstants.VERSION));
		
		// Flag init
		isInit = true;
		
		return;
	}
	
	/**
	 * Write all variants to output.
	 */
	@Override
	public void flush() {
		
		Iterator<String> recordIter = recordContainer.vcfLineIterator();
		
		// Write contig headers
		for (ReferenceSequence refSequence : referenceSequenceArray) {
			outStream.printf("##contig=<ID=%s,length=%d,md5=%s>\n", refSequence.name, refSequence.size, refSequence.digest.toString());
		}
		
		// Write FORMAT header lines
		for (String formatLine : recordContainer.getFormatHeaderStringArray()) {
			outStream.print(formatLine);
			outStream.println();
		}
		
		// Write header line
		outStream.print("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");
		
		for (String sampleName : recordContainer.getSampleNameArray()) {
			outStream.print("\t");
			outStream.print(sampleName);
		}
		
		outStream.println();
		
		// Write records
		while (recordIter.hasNext())
			outStream.println(recordIter.next());
		
		return;
	}
	
	/**
	 * Get a short description of this writer.
	 * 
	 * @return A short description of this writer.
	 */
	@Override
	public String getDescription() {
		return "Writes variant call format (VCF) files.";
	}
}
