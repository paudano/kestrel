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

package edu.gatech.kestrel.runner;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.Iterator;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import edu.gatech.kanalyze.comp.reader.SequenceSource;
import edu.gatech.kanalyze.util.kmer.KmerUtil;
import edu.gatech.kestrel.activeregion.ActiveRegionContainer;
import edu.gatech.kestrel.activeregion.ActiveRegionDetector;
import edu.gatech.kestrel.activeregion.Haplotype;
import edu.gatech.kestrel.activeregion.RegionHaplotype;
import edu.gatech.kestrel.counter.CountMap;
import edu.gatech.kestrel.counter.IkcCountMap;
import edu.gatech.kestrel.counter.MemoryCountMap;
import edu.gatech.kestrel.hapwriter.HaplotypeWriter;
import edu.gatech.kestrel.hapwriter.HaplotypeWriterInitException;
import edu.gatech.kestrel.hapwriter.NullHaplotypeWriter;
import edu.gatech.kestrel.io.InputSample;
import edu.gatech.kestrel.refreader.ReferenceReader;
import edu.gatech.kestrel.refreader.ReferenceRegion;
import edu.gatech.kestrel.refreader.ReferenceRegionContainer;
import edu.gatech.kestrel.refreader.ReferenceSequence;
import edu.gatech.kestrel.varfilter.VariantFilterRunner;
import edu.gatech.kestrel.variant.VariantCall;
import edu.gatech.kestrel.variant.VariantCaller;
import edu.gatech.kestrel.writer.VariantWriter;
import edu.gatech.kestrel.writer.VariantWriterInitException;

/**
 * Configures and runs Kestrel. This is the key API class for defining a Kestrel run, executing it, and monitoring
 * its progress.
 * </p>
 * The configuration logic for this class is implemented in <code>KestrelRunnerBase</code>, which this class
 * extends. All the runtime logic for executing and monitoring a Kestrel task is defined in this class.
 */
public class KestrelRunner extends KestrelRunnerBase implements Runnable {
	
	/** Logger object. */
	private final Logger logger;
	
	/** Set to a throwable object if <code>exec()</code> fails for any reason. */
	private Throwable errThrowable;
	
	/**
	 * Create a Kestrel runner.
	 */
	public KestrelRunner() {
		logger = LoggerFactory.getLogger(KestrelRunner.class);
		
		errThrowable = null;
		
		return;
	}
	
	/**
	 * Execute this runner.
	 */
	@Override
	public void exec() {
		
		// Declarations
		KmerUtil kUtil;  // k-mer utility
		
		ReferenceRegionContainer refRegionContainer;            // Reference sequences regions
		Iterator<ReferenceSequence> refRegionSequenceIterator;  // Iterator for references in refRegionContainer
		
		ReferenceSequence refSequence;     // A reference sequence from refRegionSequenceIterator
		ReferenceRegion[] refRegionArray;  // Reference regions associated with refSequence
		
		CountMap counter;             // K-mer counter
		
		ActiveRegionDetector arDetector;        // Detects active regions
		ActiveRegionContainer regionContainer;  // A container of active regions
		
		VariantCaller varCaller;  // Calls variants from alignments
		VariantFilterRunner variantFilterRunner;
		
		VariantWriter variantWriter;      // Writes variants to output
		HaplotypeWriter haplotypeWriter;  // Writes resolved haplotypes to output
		
		int flankLength;  // Length of flanks to add to regions extracted from reference sequences
		
		// Init
		logger.trace("exec(): Started");
		errThrowable = null;
		
		try {
			
			// Set flank length
			flankLength = this.flankLength;
			
			if (flankLength < 0) {
				flankLength = (int) (kSize * DEFAULT_FLANK_LENGTH_MULTIPLIER);
			
				if (flankLength < 0)
					flankLength = Integer.MAX_VALUE;
			}
			
			// Normalize temporary directory
			if (tempDirName == null)
				tempDirName = "";
			
			tempDirName = tempDirName.trim();
			
			if (tempDirName.isEmpty())
				tempDirName = ".";
			
			// Setup counter and kUtil
			if (kmerCountInMemory) {
				kUtil = KmerUtil.get(kSize);
				
				logger.info("Counting k-mers in memory");
				counter = new MemoryCountMap(kUtil, getCountModule());
				
			} else {
				kUtil = KmerUtil.get(kSize, kMinSize, kMinMask);
				
				logger.info("Counting k-mers from file");
				counter = new IkcCountMap(kUtil, getCountModule(), new File(tempDirName), removeIkc);
			}
			
			// Create filter runner
			variantFilterRunner = new VariantFilterRunner();
			variantFilterRunner.addFilter(variantFilterList);
			
			// Create aligner and variant caller
			varCaller = new VariantCaller();
			varCaller.setCallAmbiguousVariant(callAmbiguousVariant);
			
			if (variantCallByRegion)
				varCaller.setVariantCallByRegion();
			else
				varCaller.setVariantCallByReference();
			
			// Read reference sequences
			logger.info("Reading references");
			
			ReferenceReader referenceReader = new ReferenceReader(kUtil, loader);
			referenceReader.setFlankLength(flankLength);
			referenceReader.setRemoveDescription(removeReferenceSequenceDescription);
			referenceReader.setRevComplementNegStrand(reverseComplementNegativeStrand);
			
			try {
				refRegionContainer = referenceReader.read(
						referenceList.toArray(new SequenceSource[0]),
						((intervalContainer.isEmpty()) ? null : intervalContainer.getMap())
				);
				
			} catch (IOException ex) {
				err("Error reading reference sequence(s)", ex);
				
				return;
			}
			
			// Check reference sequences
			if (refRegionContainer.isEmpty()) {
				err("No reference sequences (see -r option)", null);
				
				return;
			}
			
			// Sort references
			refRegionContainer.sortReferences();
			
			// Get k-mer counter and active region detector
			arDetector = new ActiveRegionDetector(counter);
			
			arDetector.setAlignmentWeight(alignmentWeight);
			arDetector.setMinimumDifference(minimumDifference);
			arDetector.setDifferenceQuantile(differenceQuantile);
			arDetector.setAnchorBothEnds(anchorBothEnds);
			arDetector.setCountReverseKmers(countReverseKmers);
			arDetector.setPeakScanLength(peakScanLength);
			arDetector.setScanLimitFactor(scanLimitFactor);
			arDetector.setCallAmbiguousRegions(callAmbiguousRegions);
			arDetector.setDecayMinimum(expDecayMin);
			arDetector.setDecayAlpha(expDecayAlpha);
			arDetector.setMaxAlignerState(maxAlignerState);
			arDetector.setMaxHaplotypes(maxHaplotypes);
			arDetector.setMaxRepeatCount(maxRepeatCount);
			
			// Open output
			try {
				variantWriter = VariantWriter.getWriter(outputFormat, outputFile, refRegionContainer.getReferenceSequenceArray(), variantCallByRegion, loader);
				
			} catch (IllegalArgumentException ex) {
				err("Error opening variant writer: " + ex.getMessage(), ex);
				return;
				
			} catch (FileNotFoundException ex) {
				err("File not found while variant writer: " + ex.getMessage(), ex);
				return;
				
			} catch (IOException ex) {
				err("IO error while variant writer: " + ex.getMessage(), ex);
				return;
				
			} catch (VariantWriterInitException ex) {
				err("Error loading variant writer: " + ex.getMessage(), ex);
				return;
			}
			
			// Open haplotype output
			try {
				if (haplotypeOutputFile != null)
					haplotypeWriter = HaplotypeWriter.getWriter(haplotypeOutputFormat, haplotypeOutputFile, refRegionContainer.getReferenceSequenceArray(), loader);
				else
					haplotypeWriter = new NullHaplotypeWriter();
				
			} catch (IllegalArgumentException ex) {
				err("Error opening haplotype writer: " + ex.getMessage(), ex);
				return;
				
			} catch (FileNotFoundException ex) {
				err("File not found while haplotype writer: " + ex.getMessage(), ex);
				return;
				
			} catch (IOException ex) {
				err("IO error while haplotype writer: " + ex.getMessage(), ex);
				return;
				
			} catch (HaplotypeWriterInitException ex) {
				err("Error loading haplotype writer: " + ex.getMessage(), ex);
				return;
			}
			
			// Process samples
			for (InputSample sample : sampleList) {
				logger.info("Processing sample: {}", sample.name);
				
				variantWriter.setSampleName(sample.name);
				haplotypeWriter.setSampleName(sample.name);
				
				// Get counts
				try {
					counter.set(sample);
					
				} catch (FileNotFoundException ex) {
					err(null, ex);
					return;
				
				} catch (IOException ex) {
					err(null, ex);
					return;
				}
				
				// Find active regions
				logger.trace("Getting active regions: {}", sample.name);
				
				// Iterate over reference sequences
				refRegionSequenceIterator = refRegionContainer.refSequenceIterator();
				
				while (refRegionSequenceIterator.hasNext()) {
					
					refSequence = refRegionSequenceIterator.next();
					
					refRegionArray = refRegionContainer.get(refSequence);
					
					// Iterate over reference regions
					for (ReferenceRegion refRegion : refRegionArray) {
						logger.trace("Searching for active regions in {}", refRegion);
						
						variantWriter.setReferenceRegion(refRegion);
						
						regionContainer = arDetector.getActiveRegions(refRegion);
						
						// Find variants in active regions
						for (RegionHaplotype thisRegionHaplotype : regionContainer.haplotypes) {
							
							logger.trace("Found: {} ({} haplotypes)", thisRegionHaplotype.toString(), thisRegionHaplotype.haplotype.length);
							
							varCaller.init(thisRegionHaplotype.activeRegion);
							
							// Find variants in each haplotype
							for (Haplotype haplotype : thisRegionHaplotype.haplotype) {
								
								logger.trace("Found: {}", haplotype.toString());
								
								varCaller.add(haplotype);
								haplotypeWriter.add(haplotype);
							}
							
							// Output variant
							for (VariantCall var : varCaller.getVariants()) {
								var = variantFilterRunner.filter(var);
								
								if (var != null)
									variantWriter.writeVariant(var);
							}
						}
					}
				}
			}
			
			/** Flush variants and haplotypes. */
			logger.trace("exec(): Flushing output");
			variantWriter.flush();
			haplotypeWriter.flush();
						
			/** Complete. */
			logger.trace("exec(): Complete");
			
		} catch (Exception ex) {
			logger.error("Unexpected exception in run(): {} ({})", ex.getMessage(), ex.getClass().getSimpleName());
			errThrowable = ex;
			
			ex.printStackTrace();
		}
		
		return;
	}
	
	/**
	 * If <code>run()</code> fails for any reason, the reason is logged and a <code>Throwable</code> of
	 * the first error is saved. This method retrieves that first error if it occurs. If this method returns
	 * <code>null</code> after <code>run()</code> has completed, then <code>run()</code> was successful.
	 * This object is cleared next time <code>run()</code> is invoked.
	 * 
	 * @return A <code>Throwable</code> object if <code>run()</code> failed, or <code>null</code> if
	 *   <code>run()</code> was successful.
	 */
	public Throwable getThrowable() {
		return errThrowable;
	}
	
	/**
	 * The implementation may reset local values by overriding this method. It is called
	 * at the end of <code>configure</code>.
	 */
	@Override
	protected void implReset() {
		errThrowable = null;
		
		return;
	}
	
	/**
	 * Report an error in <code>run()</code>. This method logs the error and sets
	 * <code>errThrowable</code>.
	 * 
	 * @param message Error message.
	 * @param ex Cause or <code>null</code>.
	 */
	private void err(String message, Throwable ex) {
		
		// Log
		if (ex != null) {
			logger.error("{}: {} ({})", message, ex.getMessage(), ex.getClass().getSimpleName());
			
		} else {
			logger.error(message);
			ex = new KestrelRunnerException(message, ex);
		}
		
		// Save if this is the first error
		if (errThrowable == null)
			errThrowable = ex;
		
		return;
	}
	
	/**
	 * Set to <code>errThrowable</code> by <code>err</code> when the throwable cause is
	 * <code>null</code>.
	 */
	private class KestrelRunnerException extends RuntimeException {
		
		/** Serial version UID. */
		private static final long serialVersionUID = 1L;

		/**
		 * Create a runner exception.
		 * 
		 * @param message Message.
		 * @param cause Cause.
		 */
		public KestrelRunnerException(String message, Throwable cause) {
			super(message, cause);
			
			return;
		}
	}
}
