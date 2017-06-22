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

package edu.gatech.kestrel.hapwriter.sam;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;

import edu.gatech.kestrel.KestrelConstants;
import edu.gatech.kestrel.activeregion.ActiveRegion;
import edu.gatech.kestrel.activeregion.Haplotype;
import edu.gatech.kestrel.hapwriter.HaplotypeWriter;
import edu.gatech.kestrel.refreader.ReferenceRegion;
import edu.gatech.kestrel.refreader.ReferenceSequence;

/**
 * Output resolved haplotype alignments in SAM format.
 */
public class SamHaplotypeWriter extends HaplotypeWriter {
	
	/** A list of SAM records. */
	private final ArrayList<SamRecord> recordList;
	
	/** Output stream. */
	private PrintStream outStream;
	
	/** Init flag. */
	private boolean isInit;
	
	/**
	 * Create an uninitialized writer.
	 */
	public SamHaplotypeWriter() {
		super("SamHaplotypeWriter");
		
		recordList = new ArrayList<SamRecord>();
		
		isInit = false;
		
		return;
	}
	
	/**
	 * Add a resolved haplotype to this writer.
	 * 
	 * @param haplotype Haplotype to add.
	 * 
	 * @throws NullPointerException If <code>haplotype</code> is <code>null</code>.
	 */
	@Override
	public void add(Haplotype haplotype)
			throws NullPointerException {
		
		// Check arguments
		if (haplotype == null)
			throw new NullPointerException("Cannot add haplotype: null");
		
		recordList.add(new SamRecord(haplotype));
		
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
	protected void implInit(String writerSpec)
			throws IllegalArgumentException, FileNotFoundException, IOException {
		
		// Close if initialized
		if (isInit) {
			outStream.close();
			outStream = null;
			
			isInit = false;
		}
		
		// Open output file
		outStream = new PrintStream(output.getStream());  // throws FileNotFoundException
		
		// Reset active regions
		recordList.clear();
		
		// Sort reference sequences
		Arrays.sort(referenceSequenceArray);
		
		// Flag init
		isInit = true;
		
		return;
	}
	
	/**
	 * Get a description of this writer.
	 * 
	 * @return A description of this writer.
	 */
	@Override
	public String getDescription() {
		return "Write resolved haplotype sequences in SAM format.";
	}
	
	/**
	 * Writes all haplotypes to output. This should only be called once after all regions are added to the writer.
	 */
	@Override
	public void flush() {
		
		// Declarations
		SamRecord[] recordArray;
		
		// Get sorted regions
		recordArray = recordList.toArray(new SamRecord[0]);
		Arrays.sort(recordArray);
		
		// Write headers
		outStream.println("@HD\tVN:1.5\tSO:coordinate");  // Version
		
		for (ReferenceSequence ref : referenceSequenceArray) {
			outStream.printf("@SQ\tSN:%s\tLN:%d\tM5:%s\n", ref.name, ref.size, ref.digest.toString());  // Reference sequence
		}
		
		outStream.printf("@PG\tID:Kestrel\tVN:%s\n", KestrelConstants.VERSION);  // Program
		
		// Write records
		for (SamRecord record : recordArray) {
			
			outStream.printf(
					"%s-%s-%d-%d\t0\t%s\t%d\t255\t%s\t*\t0\t0\t%s\t*\tXD:i:%d\tXN:Z:%s\tXL:i:%d\tXR:i:%d\n",
					record.sampleName, record.hap.activeRegion.refRegion.referenceSequence.name, record.hap.activeRegion.startIndex + 1, record.hap.length,  // QNAME
					// FLAG = 0
					record.refName,  // RNAME
					record.pos,  // POS
					// MAPQ = 255 (unassigned)
					record.hap.alignment.getCigarString(),  // CIGAR
					// RNEXT = * (unassigned)
					// PNEXT = 0 (unassigned)
					// TLEN = 0 (unassigned)
					new String(record.hap.sequence),  // SEQ
					// QUAL = * (unassigned)
					record.hap.stats.min,  // XD (min k-mer depth over the haplotype)
					record.hap.activeRegion.refRegion.name,  // XN (reference region name)
					(record.hap.activeRegion.leftEnd ? 0 : 1),  // XL (left end has an anchor)
					(record.hap.activeRegion.rightEnd ? 0 : 1)  // XR (right end has an anchor)
			);
		}
		
		return;
	}
	
	/**
	 * One SAM record.
	 */
	class SamRecord implements Comparable<SamRecord> {
		
		/** Haplotype. */
		public final Haplotype hap;
		
		/** Sample name. */
		public final String sampleName;
		
		/** Reference sequence name. */
		public final String refName;
		
		/** Reference position. */
		public final int pos;
		
		/**
		 * Add a SAM record. The reference and reference region name are read from the writer object.
		 * at the time the record is created.
		 * 
		 * @param haplotype Haplotype for this record.
		 */
		public SamRecord(Haplotype haplotype) {
			
			// Check arguments
			assert(haplotype != null) :
				"Cannot add halpotype to SAM record: null";
			
			// Add
			this.sampleName = SamHaplotypeWriter.this.sampleName;
			this.hap = haplotype;
			this.refName = haplotype.activeRegion.refRegion.referenceSequence.name;
			this.pos = haplotype.activeRegion.startIndex - haplotype.activeRegion.refRegion.leftFlankLength + haplotype.activeRegion.refRegion.interval.start;
			
			return;
		}
		
		/**
		 * Compares this record with another.
		 */
		@Override
		public int compareTo(SamRecord other) {
			
			int cmpVal;
			
			cmpVal = refName.compareTo(other.refName);
			
			if (cmpVal != 0)
				return cmpVal;
			
			cmpVal = pos - other.pos;
			
			if (cmpVal != 0)
				return cmpVal;
			
			return hap.length - other.hap.length;
		}
	}
}
