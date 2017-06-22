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

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.NoSuchElementException;

import edu.gatech.kestrel.variant.VariantCall;

/**
 * Associates the same variant call from multiple samples.
 */
public class VcfRecordContainer {
	
	/** Map variant call to a VCF Record. Records may store several calls from multiple samples. */
	private HashMap<VariantCall, VcfVariantRecord> recordMap;
	
	/** Name of the current sample. */
	private String sampleName;
	
	/** A list of sample names used. */
	private List<String> sampleNameList;
	

	/**
	 * Create a new record container.
	 */
	public VcfRecordContainer() {
		
		recordMap = new HashMap<VariantCall, VcfVariantRecord>();
		sampleName = null;
		
		sampleNameList = new ArrayList<String>();
		
		return;
	}
	
	/**
	 * Start a new sample.
	 * 
	 * @param sampleName Sample name.
	 * 
	 * @throws NullPointerException If <code>sampleName</code> is <code>null</code>.
	 * @throws IllegalArgumentException If <code>sampleName</code> is empty, contains whitespace, or
	 *   is the same as a sample already added.
	 */
	public void newSample(String sampleName)
			throws NullPointerException, IllegalArgumentException {
		
		// Check arguments
		if (sampleName == null)
			throw new NullPointerException("Cannot set new sample by name: null");
		
		sampleName = sampleName.trim();
		
		if (sampleName.matches(".*\\s.*"))
			throw new IllegalArgumentException(String.format("Sample name contains whitespace: \"%s\"", sampleName));
		
		if (sampleName.isEmpty())
			throw new IllegalArgumentException("Sample name is empty");
		
		if (sampleNameList.contains(sampleName))
			throw new IllegalArgumentException(String.format("Sample name already used in this VCF record container: %s", sampleName));
		
		// Set sample
		this.sampleName = sampleName;
		sampleNameList.add(sampleName);
		
		return;
	}
	
	/**
	 * Clear all variants and samples from this record container.
	 */
	public void clear() {
		sampleNameList.clear();
		recordMap.clear();
		sampleName = null;
		
		return;
	}
	
	/**
	 * Add a variant call to this container.
	 * 
	 * @param variantCall Variant call to add.
	 * 
	 * @throws NullPointerException If <code>variantCall</code> is <code>null</code>.
	 * @throws IllegalStateException If <code>newSample()</code> has not been called.
	 */
	public void addVariant(VariantCall variantCall)
			throws NullPointerException, IllegalStateException {
		
		VcfVariantRecord vcfRecord;
		
		// Check arguments
		if (variantCall == null)
			throw new NullPointerException("Cannot add variant call: null");
		
		if (sampleName == null)
			throw new IllegalStateException("Cannot add variants to VCF container before the first sample (newSample()) is declared");
		
		// Get existing VCF record
		vcfRecord = recordMap.get(variantCall);
		
		if (vcfRecord == null) {
			vcfRecord = new VcfVariantRecord(sampleName, variantCall);
			
			recordMap.put(variantCall, vcfRecord);
			
		} else {
			vcfRecord.addVariant(sampleName, variantCall);
		}
		
		return;
	}
	
	/**
	 * Get an iterator to generate a string for each VCF record, and return records
	 * in VCF order.
	 * 
	 * @return VCF formatted record iterator.
	 */
	public Iterator<String> vcfLineIterator() {
		return new VcfStringRecordIterator();
	}
	
	/**
	 * Get an array containing the names of samples added to this container.
	 * 
	 * @return An array of sample names.
	 */
	public String[] getSampleNameArray() {
		return sampleNameList.toArray(new String[0]);
	}
	
	/**
	 * Get an array of format headers for the VCF file.
	 * 
	 * @return Array of format headers for the VCF file.
	 */
	public String[] getFormatHeaderStringArray() {
		return new String[] {
				"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
				"##FORMAT=<ID=GDP,Number=A,Type=Integer,Description=\"Estimated depth of all haplotypes supporting the alternate variant\">",
				"##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Estimated depth of all haplotypes in the variant active region\">"
		};
	}
	
	/**
	 * Represents one variant that was found in at least one sample. Variants are maintained in the order in
	 * which they are added.
	 */
	private class VcfVariantRecord implements Comparable<VcfVariantRecord> {
		
		/** First sample info element. */
		public VcfVariantSampleInfo sampleInfoHead;
		
		/** Last sample info element (new samples are appended here). */
		public VcfVariantSampleInfo sampleInfoTail;
		
		/** Chromosome. */
		public final String chrom;
		
		/** Position where this variant starts (1-based). */
		public final int pos;
		
		/** Accession (always missing). */
		public static final String id = ".";
		
		/** Reference sequence string. */
		public final String ref;
		
		/** Alternate sequence string. */
		public final String alt;
		
		/** Variant call quality. */
		public static final String qual = ".";
		
		/** Filtered reason. */
		public static final String filter = ".";
		
		/** Info field. */
		public static final String info = ".";
		
		/** Format of the sample field. */
		public static final String format = "GT:GDP:DP";
		
		/**
		 * Create a container of equivalent variant calls (from multiple samples).
		 * 
		 * @param sampleName Name of the sample <code>initCall</code> is from.
		 * @param variantCall Initial variant call. Information for the container is extracted from this call.
		 * 
		 * @throws NullPointerException If any argument is <code>null</code>.
		 * @throws IllegalArgumentException If <code>sampleName</code> is <code>null</code>.
		 */
		public VcfVariantRecord(String sampleName, VariantCall variantCall)
				throws NullPointerException, IllegalArgumentException {
			
			// Check arguments
			if (sampleName == null)
				throw new NullPointerException("sampleName is null");
			
			if (sampleName.isEmpty())
				throw new IllegalArgumentException("sampleName is empty");
			
			if (variantCall == null)
				throw new NullPointerException("variantCall is null");
			
			// Get fields
			chrom = variantCall.activeRegion.refRegion.referenceSequence.name;
			
			pos = variantCall.getVcfPos();
			ref = variantCall.getVcfRef();
			alt = variantCall.getVcfAlt();
			
			sampleInfoHead = sampleInfoTail = new VcfVariantSampleInfo(sampleName, variantCall);
			
			return;
		}
		
		/**
		 * Add a variant to this record.
		 * 
		 * @param sampleName Name of the sample the variant was called in.
		 * @param variantCall Variant call.
		 * 
		 * @throws NullPointerException If any argument is <code>null</code>.
		 * @throws IllegalArgumentException If <code>sampleName</code> is <code>null</code>.
		 */
		public void addVariant(String sampleName, VariantCall variantCall)
				throws NullPointerException, IllegalArgumentException {
			
			// Check arguments
			if (sampleName == null)
				throw new NullPointerException("sampleName is null");
			
			if (sampleName.isEmpty())
				throw new IllegalArgumentException("sampleName is empty");
			
			if (variantCall == null)
				throw new NullPointerException("variantCall is null");
			
			// Append record
			sampleInfoTail.nextSampleInfo = new VcfVariantSampleInfo(sampleName, variantCall);
			sampleInfoTail = sampleInfoTail.nextSampleInfo;
		}
		
		/**
		 * Compare VCF records by the variant they represent. This does not include
		 * the samples in the VCF record, only the pos, ref, and alt fields. 
		 */
		@Override
		public int compareTo(VcfVariantRecord vcfRecord)
				throws NullPointerException {
			
			int cmpVal;
			
			// Check arguments
			if (vcfRecord == null)
				throw new NullPointerException("Cannot compare VCF record to other VCF record: null");
			
			// Compare
			if (pos != vcfRecord.pos)
				return pos - vcfRecord.pos;
			
			cmpVal = ref.compareTo(vcfRecord.ref);
			
			if (cmpVal != 0)
				return cmpVal;
			
			cmpVal = alt.compareTo(vcfRecord.alt);
			
			if (cmpVal != 0)
				return cmpVal;
			
			return 0;
		}
	}
	
	/**
	 * One variant contained within a <code>VcfVariantRecord</code>.
	 */
	private class VcfVariantSampleInfo {
		
		/**
		 * Point to the next sample info in this record or <code>null</code> if this is the last
		 * sample info element.
		 */
		public VcfVariantSampleInfo nextSampleInfo;
		
		/** Sample name */
		public final String sampleName;
		
		/** Variant call */
		public final VariantCall variantCall;
		
		/**
		 * Create a record for one variant.
		 * 
		 * @param sampleName Name of the sample.
		 * @param variantCall Variant call.
		 */
		public VcfVariantSampleInfo(String sampleName, VariantCall variantCall) {
			
			// Check arguments
			assert (sampleName != null) :
				"sampleName is null";
			
			assert (! sampleName.isEmpty()) :
				"sampleName is empty";
			
			assert (variantCall != null) :
				"variantCall is null";
			
			// Save fields
			nextSampleInfo = null;
			this.sampleName = sampleName;
			this.variantCall = variantCall;
			
			return;
		}
	}
	
	/**
	 * Iterate over VCF records and output a formatted string in VCF format.
	 */
	private class VcfStringRecordIterator implements Iterator<String> {
		
		/** Next variant record. */
		private VcfVariantRecord[] vcfRecordList;
		
		/** Index of the next record in <code>vcfRecordList</code>. */
		private int vcfRecordIndex;
		
		/** Size of <code>vcfRecordList</code>. */
		private final int vcfRecordSize;
		
		/**
		 * Create a new record iterator.
		 */
		public VcfStringRecordIterator() {
			
			Collection<VcfVariantRecord> recordCollection = recordMap.values();
			
			vcfRecordSize = recordCollection.size();
			vcfRecordList = recordCollection.toArray(new VcfVariantRecord[vcfRecordSize]);
			vcfRecordIndex = 0;
			
			Arrays.sort(vcfRecordList);
			
			return;
		}
		
		/**
		 * Determine if this iterator has another record.
		 * 
		 * @return <code>true</code> if this variant has another record.
		 */
		@Override
		public boolean hasNext() {
			return vcfRecordIndex < vcfRecordSize;
		}
		
		@Override
		public String next()
				throws NoSuchElementException {
			
			StringBuilder record = new StringBuilder();  // Build VCF record string
			
			VcfVariantRecord nextRecord;  // Next VCF record
			
			VcfVariantSampleInfo sampleInfo;  // For traversing sample info
			
			// Check state
			if (vcfRecordIndex >= vcfRecordSize)
				throw new NoSuchElementException("No more VCF string records");
			
			// Add record fields up to the samples
			nextRecord = vcfRecordList[vcfRecordIndex];
			
			record.append(nextRecord.chrom);
			record.append("\t");
			
			record.append(nextRecord.pos);
			record.append("\t");
			
			record.append(VcfVariantRecord.id);
			record.append("\t");
			
			record.append(nextRecord.ref);
			record.append("\t");
			
			record.append(nextRecord.alt);
			record.append("\t");
			
			record.append(VcfVariantRecord.qual);
			record.append("\t");
			
			record.append(VcfVariantRecord.filter);
			record.append("\t");
			
			record.append(VcfVariantRecord.info);
			record.append("\t");
			
			record.append(VcfVariantRecord.format);
			
			// Add samples
			sampleInfo = nextRecord.sampleInfoHead;
			
			for (String sampleName : sampleNameList) {
				
				if (sampleName.equals(sampleInfo.sampleName)) {
					// Sample has an info record
					record.append("\t");
					record.append("1:");
					record.append(sampleInfo.variantCall.variantDepth);
					record.append(":");
					record.append(sampleInfo.variantCall.locusDepth);
					
				} else {
					// Sample does not have an info record
					record.append("\t");
					record.append("0:.:.");
				}
				
			}
			
			// Advance to the next record
			++vcfRecordIndex;
			
			// Return record string
			return record.toString();
		}

		@Override
		public void remove()
				throws UnsupportedOperationException {
			
			throw new UnsupportedOperationException("VCF string iterator does not support remove()");
		}
	}
}
