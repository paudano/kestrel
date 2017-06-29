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
import java.io.FileDescriptor;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.net.MalformedURLException;
import java.net.URL;
import java.net.URLClassLoader;
import java.nio.charset.Charset;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.locks.ReentrantLock;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import edu.gatech.kanalyze.KAnalyzeConstants;
import edu.gatech.kanalyze.KAnalyzeRunnable;
import edu.gatech.kanalyze.comp.reader.SequenceSource;
import edu.gatech.kanalyze.condition.ConditionEvent;
import edu.gatech.kanalyze.condition.StreamConditionListener;
import edu.gatech.kanalyze.module.count.CountModule;
import edu.gatech.kanalyze.util.argparse.HelpRunException;
import edu.gatech.kanalyze.util.kmer.KmerUtil;
import edu.gatech.kestrel.KestrelConstants;
import edu.gatech.kestrel.LogLevel;
import edu.gatech.kestrel.activeregion.ActiveRegionDetector;
import edu.gatech.kestrel.align.AlignmentWeight;
import edu.gatech.kestrel.align.KmerAligner;
import edu.gatech.kestrel.align.KmerAlignmentBuilder;
import edu.gatech.kestrel.interval.IntervalReader;
import edu.gatech.kestrel.interval.IntervalReaderInitException;
import edu.gatech.kestrel.interval.RegionInterval;
import edu.gatech.kestrel.interval.RegionIntervalContainer;
import edu.gatech.kestrel.io.InputSample;
import edu.gatech.kestrel.io.StreamableOutput;
import edu.gatech.kestrel.refreader.ReferenceReader;
import edu.gatech.kestrel.varfilter.VariantFilter;
import edu.gatech.kestrel.varfilter.VariantFilterInitException;

/**
 * Configuration logic for <code>KestrelRunner</code>.
 */
public abstract class KestrelRunnerBase extends KAnalyzeRunnable {
	
	//
	// Private fields
	//
	
	/** Lock for making configuration and runtime tasks atomic (thread-safe). */
	protected ReentrantLock runnerLock;
	
	/** List of input sources. */
	protected List<InputSample> sampleList;
	
	/** List of reference sources. */
	protected List<SequenceSource> referenceList;
	
	
	//
	// Options
	//
	
	/** Size of k-mers to use for analysis. */
	protected int kSize;
	
	/** K-mer minimizer size or <code>0</code> if a minimizer is not applied. */
	protected int kMinSize;
	
	/** K-mer minimizer mask or <code>0x00000000</code> if a minimizer mask is not applied. */
	protected int kMinMask;
	
	/**
	 * Output file.
	 */
	protected StreamableOutput outputFile;
	
	/**
	 * Output format specifier. The output format name may be followed by a colon and an implementation-defined
	 * string that modifies the format.
	 */
	protected String outputFormat;
	
	/**
	 * Suggested log file. The Kestrel API does not configure logging, but this may be set by the command-line.
	 * The CLUI will log to this file. 
	 */
	protected StreamableOutput logFile;
	
	/**
	 * Log level. Controls which messages are sent to <code>logOutputObj</code>. The Kestrel API
	 * does not configure logging, but this may be set by the command-line. The CLUI will log to
	 * <code>logFile</code> at this level.
	 */
	protected LogLevel logLevel;
	
	/** Loader for the classpath and elements in <code>libraryUrlList</code>. */
	protected URLClassLoader loader;
	
	/** List of locations to be searched when loading dynamic classes. */
	private List<URL> libraryUrlList;
	
	/**
	 * Handles argument parsing and setting default values.
	 */
	protected KestrelArgumentParser argumentParser;
	
	/** Structures of weights for the alignment. */
	protected AlignmentWeight alignmentWeight;
	
	/** Location of temporary files. */
	protected String tempDirName;
	
	/** Minimum k-mer count. K-mers with lower counts are filtered out. */
	protected int minKmerCount;
	
	/**
	 * If <code>true</code>, Keep k-mer counts from each sample in memory while it is being processed.
	 * If <code>false</code>, k-mer counts are offloaded to an indexed k-mer count file.
	 */
	protected boolean kmerCountInMemory;
	
	/**
	 * Free resources as soon as possible if set to <code>true</code>. This may use less memory, but
	 * some expensive resources may have to be recreated for each sample.
	 */
	protected boolean freeResources;
	
	/** If <code>true</code>, delete the indexed k-mer count (IKC) file after each sample finishes. */
	protected boolean removeIkc;
	
	/**
	 * If <code>true</code>, active regions must be bordered on both sides by unaltered k-mers. The
	 * variants called in these active regions are supported by better evidence than those called
	 * in active regions that reach an end, however, variants close to the end of a reference may
	 * be missed.
	 */
	protected boolean anchorBothEnds;
	
	/** If <code>true</code>, count reverse k-mers in region statistics. */
	protected boolean countReverseKmers;
	
	/** Minimum k-mer difference to flag a potential active region. */
	protected int minimumDifference;
	
	/**
	 * K-mer count differences must be at least this quantile of total difference to trigger an
	 * active region scan over the reference sequence. <code>0.0</code> disables picking the
	 * threshold by quantile and <code>minimumDifference</code> is used as the threshold. If
	 * the quantile is less than <code>minimumDifference</code>, then the minimum is used. The
	 * difference of neighboring k-mers is defined as their distance apart (the absolute value
	 * of one k-mer count subtracted by the count of its neighbor).
	 */
	protected double differenceQuantile;
	
	/**
	 * For peak detection, scan forward over a suspected peak this number of k-mers. If
	 * <code>0</code>, disable peak detection.
	 */
	protected int peakScanLength;
	
	/**
	 * Limit the length of an active region scan. This is computed by
	 * taking maximum length of a gap (assuming no other variants), adding one, and then multiplying
	 * by this factor.
	 */
	protected double scanLimitFactor;
	
	/**
	 * The alpha value used to calculate lambda (log(alpha)/k) for the exponential decay function
	 * in the active region detector.
	 */
	protected double expDecayAlpha;
	
	/**
	 * Set the minimum (lower asymptotic bound) of the exponential decay function as a proportion
	 * of the anchor k-mer count for active region detection.
	 */
	protected double expDecayMin;
	
	/**
	 * If <code>true</code>, then allow haplotypes and variants to occur over ambiguous regions.
	 */
	protected boolean callAmbiguousRegions;
	
	/** If <code>true</code>, allow variants that span ambiguous regions of the reference. */
	protected boolean callAmbiguousVariant;
	
	/** List of variant filters. */
	protected List<VariantFilter> variantFilterList;
	
	/** Manages region intervals variants are called on. */
	protected RegionIntervalContainer intervalContainer;
	
	/** Length of flanks to add to intervals over the reference sequence. */
	protected int flankLength;
	
	/** Output variant calls relative to the reference regions when <code>true</code>. */
	protected boolean variantCallByRegion;
	
	/** The maximum number of alignment states that may be saved by the aligner. */
	protected int maxAlignerState;
	
	/** The maximum number of haplotypes that may be saved by the aligner. */
	protected int maxHaplotypes;
	
	/** <code>true</code> if the sequence description should be removed from each reference sequence name. */
	protected boolean removeReferenceSequenceDescription;
	
	/** <code>true</code> if reference regions on the negative strand should be reverse complemented before variant calling. */
	protected boolean reverseComplementNegativeStrand;
	
	/** Haplotype output file or <code>null</code> if haplotypes are not output. */
	protected StreamableOutput haplotypeOutputFile;
	
	/** Haplotype output format. Determines the format of the output file for <code>haplotypeOutputFile</code>. */
	protected String haplotypeOutputFormat;
	
	
	//
	// Limits and defaults
	//
	
	/** Default input file format. */
	public static final String DEFAULT_FORMAT = "auto";
	
	/** Default k-mer size. */
	public static final int DEFAULT_KSIZE = 31;
	
	/** Default output file name. */
	public static final StreamableOutput DEFAULT_OUTPUT_FILE = StreamableOutput.STDOUT;
	
	/** Default output format. */
	public static final String DEFAULT_OUTPUT_FORMAT = "vcf";
	
	/** Default log file location. */
	public static final StreamableOutput DEFAULT_LOG_FILE = StreamableOutput.get(FileDescriptor.err);
	
	/** Default log level. */
	public static final LogLevel DEFAULT_LOG_LEVEL = LogLevel.WARN;
	
	/** Default character-set for input files. */
	public static final Charset DEFAULT_CHARSET = Charset.forName("UTF-8");
	
	/** Default minimizer size (0 disables minimizers). */
	public static final int DEFAULT_MINIMIZER_SIZE = 15;
	
	/** Default minimizer mask. */
	public static final int DEFAULT_MINIMIZER_MASK = 0x00000000;
	
	/** Default minimum k-mer count. */
	public static final int DEFAULT_MIN_KMER_COUNT = 5;
	
	/** Default size of the sequence buffer for reading reference sequences. */ 
	public static final int DEFAULT_READER_SEQUENCE_BUFFER_SIZE = 1024;
	
	/** Default size of sequence queues. */
	public static final int DEFAULT_SEQUENCE_QUEUE_SIZE = KAnalyzeConstants.DEFAULT_SEQUENCE_QUEUE_SIZE;
	
	/** Default sequence batch size. */
	public static final int DEFAULT_SEQUENCE_BATCH_SIZE = KAnalyzeConstants.DEFAULT_SEQUENCE_BATCH_SIZE;
	
	/** Default size of k-mer queues. */
	public static final int DEFAULT_KMER_QUEUE_SIZE = KAnalyzeConstants.DEFAULT_KMER_QUEUE_SIZE;
	
	/** Default k-mer batch size. */
	public static final int DEFAULT_KMER_BATCH_SIZE = KAnalyzeConstants.DEFAULT_KMER_BATCH_SIZE;
	
	/** Default size of batch caches. */
	public static final int DEFAULT_CACHE_SIZE = 100;
	
	/** Default number of KAnalyze k-mer component threads. */
	public static final int DEFAULT_KMER_THREAD_COUNT = KAnalyzeConstants.DEFAULT_KMER_THREAD_COUNT;
	
	/** Default number of KAnalyze k-mer component threads. */
	public static final int DEFAULT_SPLIT_THREAD_COUNT = KAnalyzeConstants.DEFAULT_SPLIT_THREAD_COUNT;
	
	/** Default option for storing k-mer counts in memory. */
	public static final boolean DEFAULT_KMER_COUNT_IN_MEMORY = false;
	
	/** Free resources as soon as possible if set. */
	public static final boolean DEFAULT_FREE_RESOURCES = false;
	
	/** Remove the IKC (indexed k-mer count) file after running. */
	public static final boolean DEFAULT_REMOVE_IKC = true;
	
	/** Count reverse complement k-mers in region statistics. */
	public static final boolean DEFAULT_COUNT_REV_KMER = true;
	
	/** This value is multiplied by the k-mer size to determine the default flank length. */
	public static final double DEFAULT_FLANK_LENGTH_MULTIPLIER = 3.5;
	
	/** Default reverse-complement negative strand reference regions. */
	public static final boolean DEFAULT_REVERSE_COMPLEMENT_NEGATIVE_STRAND = ReferenceReader.DEFAULT_REVERSE_COMPLEMENT_NEGATIVE_STRAND;
	
	/** Default haplotype output format. */
	public static final String DEFAULT_HAPLOTYPE_OUTPUT_FORMAT = "sam";
	
	
	//
	// Other constants
	//
	
	/** Root logger for all Kestrel classes. */
	private static final Logger KESTREL_ROOT_LOGGER = LoggerFactory.getLogger("edu.gatech.kestrel");
	
	
	/**
	 * Create a runner base and configure it to the default state.
	 */
	public KestrelRunnerBase() {
		
		// Initialize lock
		runnerLock = new ReentrantLock();
		
		// Initialize data structures
		referenceList = new ArrayList<SequenceSource>();
		sampleList = new ArrayList<InputSample>();
		variantFilterList = new ArrayList<VariantFilter>();
		intervalContainer = new RegionIntervalContainer();
		
		// Initialize loader and URLs
		libraryUrlList = new ArrayList<URL>();
		loader = new URLClassLoader(new URL[0]);
		
		// Create argument parser
		argumentParser = new KestrelArgumentParser(this);
		
		// Set default values;
		argumentParser.init();
		
		return;
	}
	
	/**
	 * Configure this runner given a set of command-line arguments.
	 * 
	 * @param args Command-line arguments. If <code>null</code>, an empty set of arguments is used.
	 * 
	 * @throws HelpRunException If a help option was run and this program should terminate. When this
	 *   exception is thrown, the runner may not be left in a runnable state.
	 * @throws ConfigurationErrorException If an error occurs while configuring this runner.
	 */
	public void configure(String[] args)
			throws HelpRunException, ConfigurationErrorException {
		
		// Check arguments
		if (args == null)
			args = new String[0];
		
		// Parse arguments
		argumentParser.parse(args);  // throws HelpRunException
		
		// Set warning
		argumentParser.throwConfigurationError();  // throws ConfigurationErrorException
		
		// Reset implementation
		implReset();
		
		return;
	}
	
	/**
	 * The implementation may reset local values by overriding this method. It is called
	 * at the end of <code>configure</code>.
	 */
	protected void implReset() {
		return;
	}
	
	/**
	 * Add an input sample.
	 * 
	 * @param sample Input sample.
	 * 
	 * @throws NullPointerException If <code>sample</code> is <code>null</code>.
	 */
	public void addSample(InputSample sample)
			throws NullPointerException {
		
		if (sample == null)
			throw new NullPointerException("Cannot add input sample: " + null);
		
		sampleList.add(sample);
		
		return;
	}
	
	/**
	 * Add an array of input samples.
	 * 
	 * @param sampleArray Input sources.
	 * 
	 * @throws NullPointerException If <code>sampleArray</code> is <code>null</code>.
	 * @throws IllegalArgumentException If <code>sampleArray</code> contains
	 *   <code>null</code> references.
	 */
	public void addSource(InputSample[] sampleArray)
			throws NullPointerException {
		
		if (sampleArray == null)
			throw new NullPointerException("Cannot add sample sources from array: " + null);
		
		// Check for null references before adding sources
		for (int index = 0; index < sampleArray.length; ++index)
			if (sampleArray[index] == null)
				throw new IllegalArgumentException("Cannot add null sample source in array at index " + index);
		
		// Add sources
		for (int index = 0; index < sampleArray.length; ++index)
			sampleList.add(sampleArray[index]);
		
		return;
	}
	
	/**
	 * Clear all sequence sources and reset the sequence source index.
	 */
	public void clearSamples() {
		sampleList.clear();
		
		return;
	}
	
	/**
	 * Add a sequence source containing references variants will be called against.
	 * 
	 * @param referenceSource Reference sequence source.
	 * 
	 * @throws NullPointerException If <code>referenceSource</code> is <code>null</code>.
	 */
	public void addReference(SequenceSource referenceSource)
			throws NullPointerException {
		
		if (referenceSource == null)
			throw new NullPointerException("Cannot add reference: null");
		
		referenceList.add(referenceSource);
		
		return;
	}
	
	/**
	 * Clear reference sequences.
	 */
	public void clearReference() {
		referenceList.clear();
		
		return;
	}
	
	/**
	 * Set the k-mer size.
	 * 
	 * @param kSize K-mer size.
	 * 
	 * @throws IllegalArgumentException If <code>kSize</code> is less than <code>Constants.MIN_KMER_SIZE</code>.
	 */
	public void setKSize(int kSize)
			throws IllegalArgumentException {
		
		if (kSize < KestrelConstants.MIN_KMER_SIZE)
			throw new IllegalArgumentException("K-mer size must not be less than " + KestrelConstants.MIN_KMER_SIZE + ": " + kSize);
		
		this.kSize = kSize;
	}
	
	/**
	 * Get k-mer size.
	 * 
	 * @return K-mer size.
	 */
	public int getKSize() {
		return kSize;
	}
	
	/**
	 * Set the output file. 
	 * 
	 * @param outputObj Output object identifying where the output is
	 *   written. This object must be of class <code>String</code>,
	 *   <code>File</code>, or <code>FileDescriptor</code>, or
	 *   <code>StreamableOutput</code>. If <code>null</code> or an empty string,
	 *   the default output file, <code>DEFAULT_OUTPUT_FILE</code>, is used.
	 * 
	 * @throws IllegalArgumentException If <code>outputObj</code> is not
	 *   <code>null</code> and is not an instance of <code>String</code>,
	 *   <code>File</code>, <code>FileDescriptor</code>, or <code>StreamableOutput</code>.
	 * 
	 * @see #DEFAULT_OUTPUT_FILE
	 */
	public void setOutputFile(Object outputObj)
			throws IllegalArgumentException {

		// Check null
		if (outputObj == null) {
			outputObj = DEFAULT_OUTPUT_FILE;
			
			return;
		}
		
		// Check empty string
		if (outputObj instanceof String) {
			if (((String) outputObj).trim().isEmpty()) {
				outputObj = DEFAULT_OUTPUT_FILE;
				
				return;
			}
		}
		
		// Set output file
		outputFile = StreamableOutput.get(outputObj);  // throws IllegalArgumentException 
		
		return;
	}
	
	/**
	 * Get the output file.
	 * 
	 * @return Output file.
	 */
	public StreamableOutput getOutputFile() {
		return outputFile;
	}
	
	/**
	 * Set the output format for this runner.
	 * 
	 * @param outputFormat Output format.
	 * 
	 * @throws NullPointerException If <code>outputFormat</code> is <code>null</code>.
	 * @throws IllegalArgumentException If <code>outputFormat</code> is empty.
	 */
	public void setOutputFormat(String outputFormat)
			throws NullPointerException, IllegalArgumentException {
		
		if (outputFormat == null)
			throw new NullPointerException("Cannot set output format: null");
		
		outputFormat = outputFormat.trim();
		
		if (outputFormat.isEmpty())
			throw new IllegalArgumentException("Cannot set output format to an empty string");
		
		this.outputFormat = outputFormat;
		
		return;
	}
	
	/**
	 * Get the output format set on this runner.
	 * 
	 * @return Output format.
	 */
	public String getOutputFormat() {
		return outputFormat;
	}
	
	/**
	 * Set the log file. 
	 * <p/>
	 * It is up to the software using the Kestrel API to configure the logger. Kestrel
	 * logs messages with SLF4J, and this can be sent to a number of frameworks.
	 * 
	 * @param logObj Output object identifying where the output is
	 *   written. This object must be of class <code>String</code>,
	 *   <code>File</code>, or <code>FileDescriptor</code>, or
	 *   <code>StreamableOutput</code>. If <code>null</code> or an empty string,
	 *   the default output file, <code>DEFAULT_LOG_FILE</code>, is used.
	 * 
	 * @throws IllegalArgumentException If <code>logObj</code> is not
	 *   <code>null</code> and is not an instance of <code>String</code>,
	 *   <code>File</code>, <code>FileDescriptor</code>, or <code>StreamableOutput</code>.
	 * 
	 * @see #DEFAULT_LOG_FILE
	 * @see <a href="http://www.slf4j.org/">http://www.slf4j.org/</a>
	 */
	public void setLogFile(Object logObj)
			throws IllegalArgumentException {

		// Check null
		if (logObj == null) {
			logObj = DEFAULT_LOG_FILE;
			
			return;
		}
		
		// Check empty string
		if (logObj instanceof String) {
			if (((String) logObj).trim().isEmpty()) {
				logObj = DEFAULT_LOG_FILE;
				
				return;
			}
		}
		
		// Set output file
		logFile = StreamableOutput.get(logObj);  // throws IllegalArgumentException 
		
		return;
	}
	
	/**
	 * Get the log file.
	 * <p/>
	 * It is up to the software using the Kestrel API to configure the logger. Kestrel
	 * logs messages with SLF4J, and this can be sent to a number of frameworks.
	 * 
	 * @return Log file.
	 * 
	 * @see <a href="http://www.slf4j.org/">http://www.slf4j.org/</a>
	 */
	public StreamableOutput getLogFile() {
		return logFile;
	}
	
	/**
	 * Set the log level returned by <code>getLogLevel()</code>. Since different
	 * implementations using Kestrel as an API may have different logging frameworks,
	 * the Kestrel runner does not attempt to configure logging. This method is
	 * designed to capture the suggested logging level and allow the implementation
	 * to decide how to use it. The command-line interface uses these settings, which are
	 * configured from command-line options, to implement logging in LOGBack framework.
	 * <p/>
	 * Kestrel itself uses SLF4J to log messages, and these messages can be sent to many
	 * logging frameworks including JUL (Java Utility Logging), Log4j, and LOGBack. See the
	 * documentation for SLF4J for details on how to send log messages to your framework, and
	 * keep in mind that LOGBack libraries are in the classpath of the Kestrel JAR file.  
	 * 
	 * @param levelName Level of this runner. If <code>null</code>, <code>LogLevel.OFF</code> is
	 *   set.
	 * 
	 * @throws IllegalArgumentException If <code>levelName</code> is not a
	 *   recognized level name.
	 * 
	 * @see <a href="http://www.slf4j.org/">http://www.slf4j.org/</a>
	 */
	public void setLogLevel(String levelName)
			throws IllegalArgumentException {
		
		if (levelName == null)
			logLevel = LogLevel.OFF;
		
		setLogLevel(LogLevel.getLevel(levelName));
		
		return;
	}
	
	/**
	 * Set the log level returned by <code>getLogLevel()</code>. Since different
	 * implementations using Kestrel as an API may have different logging frameworks,
	 * the Kestrel runner does not attempt to configure logging. This method is
	 * designed to capture the suggested logging level and allow the implementation
	 * to decide how to use it. The command-line interface uses these settings, which are
	 * configured from command-line options, to implement logging in LOGBack framework.
	 * <p/>
	 * Kestrel itself uses SLF4J to log messages, and these messages can be sent to many
	 * logging frameworks including JUL (Java Utility Logging), Log4j, and LOGBack. See the
	 * documentation for SLF4J for details on how to send log messages to your framework, and
	 * keep in mind that LOGBack libraries are in the classpath of the Kestrel JAR file.  
	 * 
	 * @param logLevel Level of this runner. If <code>null</code>, <code>DEFAULT_LOG_LEVEL</code> is
	 *   set.
	 * 
	 * @throws IllegalArgumentException If <code>levelName</code> is not a
	 *   recognized level name.
	 * 
	 * @see <a href="http://www.slf4j.org/">http://www.slf4j.org/</a>
	 */
	public void setLogLevel(LogLevel logLevel)
			throws IllegalArgumentException {
		
		if (logLevel == null)
			logLevel = DEFAULT_LOG_LEVEL;
		
		this.logLevel = logLevel;
		
		if (KESTREL_ROOT_LOGGER instanceof ch.qos.logback.classic.Logger)
			((ch.qos.logback.classic.Logger) KESTREL_ROOT_LOGGER).setLevel(logLevel.level);
		
		return;
	}
	
	/**
	 * Get the log level suggested by this runner. It is up to the software using the Kestrel
	 * API to configure the logger. Kestrel logs messages with SLF4J, and this can be sent to
	 * a number of frameworks.
	 * 
	 * @return The suggested log level of this runner.
	 * 
	 * @see <a href="http://www.slf4j.org/">http://www.slf4j.org/</a>
	 */
	public LogLevel getLogLevel() {
		return logLevel;
	}
	
	/**
	 * Set the k-mer minimizer size for k-mers stored in indexed count files.
	 * 
	 * @param kMinSize Minimizer size.
	 * 
	 * @throws IllegalArgumentException If <code>kMinSize</code> is not a valid 
	 *   minimizer size.
	 */
	public void setMinimizerSize(int kMinSize)
			throws IllegalArgumentException {
		
		KmerUtil.validateMinimizerSize(kMinSize);
		
		this.kMinSize = kMinSize;
		
		return;
	}
	
	/**
	 * Get the k-mer minimizer size.
	 * 
	 * @return Get the minimizer size.
	 */
	public int getMinimizerSize() {
		return kMinSize;
	}
	
	/**
	 * Set the minimizer mask, which is an XOR mask applied to sub-k-mers while choosing the k-mer
	 * minimizer. If <code>0x00000000</code>, the mask is essentially disabled. This mask may be used
	 * to break large minimizer groups due to low-complexity k-mers. If the k-mer minimizer is not used,
	 * then this mask has no effect.
	 * 
	 * @param kMinMask Minimizer mask.
	 */
	public void setMinimizerMask(int kMinMask) {
		
		this.kMinMask = kMinMask;
		
		return;
	}
	
	/**
	 * Get the k-mer minimizer mask.
	 * 
	 * @return K-mer minimizer mask.
	 */
	public int getMinimizerMask() {
		return kMinMask;
	}
	
	/**
	 * Add a library URL.
	 * 
	 * @param url URL to add. If this URL is already in the list, it is not added.
	 * 
	 * @throws NullPointerException If <code>url</code> is <code>null</code>.
	 */
	public void addLibraryURL(URL url)
			throws NullPointerException {
		
		if (url == null)
			throw new NullPointerException("Cannot add library URL: null");
		
		if (libraryUrlList.contains(url))
			return;
		
		libraryUrlList.add(url);
		loader = new URLClassLoader(libraryUrlList.toArray(new URL[0]));
		
		return;
	}
	
	/**
	 * Add a library file to the dynamic class loader.
	 * 
	 * @param file File to add.
	 * 
	 * @throws NullPointerException If <code>file</code> is <code>null</code>.
	 * @throws MalformedURLException If there is a bug causing <code>File.toURI()</code> to
	 *   output a URI that cannot be translated to a URL. Normally, this exception means
	 *   that no legal protocol could be found in a specification string, or the string could
	 *   not be parsed (text from <code>MalformedURLException</code>).
	 * @throws FileNotFoundException If <code>file</code> does not exist.
	 * @throws IOException If <code>file</code> is not a regular file or directory.
	 */
	public void addLibraryFile(File file)
			throws NullPointerException, MalformedURLException, FileNotFoundException, IOException {
		
		// Check arguments
		if (file == null)
			throw new NullPointerException("Cannot add library file to the dynamic loader: null");
		
		if (! file.exists())
			throw new FileNotFoundException("Library file was not found: " + file.getPath());
		
		if (file.isFile()) {
			addLibraryURL(new URL("file:" + file.getAbsolutePath()));  // throws MalformedURLException
			
		} else if (file.isDirectory()) {
			addLibraryURL(new URL("file:" + file.getAbsolutePath() + "/"));  // throws MalformedURLException
			
		} else {
			throw new IOException("Library file was is not a regular file or directory: " + file.getPath());
		}
		
		return;
	}
	
	/**
	 * Clear all loaders configured on this module. These loaders are used to locate dynamic
	 * classes, and they may be removed to avoid conflicts or for security reasons. Generally,
	 * it is desirable to keep the loaders.
	 */
	public void clearLibraries() {
		libraryUrlList.clear();
		loader = new URLClassLoader(new URL[0]);
		
		return;
	}
	
	/**
	 * Set the alignment weight from a string representing the weights as a
	 * comma-separated vector of values. The order of the values is match, mismatch, gap-open, gap-extend, and
	 * initial score. If less than five elements are defined, then they are assigned in this order and missing elements are
	 * assigned their default values, as defined by <code>SequenceWeight</code> (see &quot;DEFAULT_&quot; constants).
	 * <p/>
	 * The vector may be surrounded by matching parenthesis, angle-braces, square-braces, curly-braces, or
	 * it may not be surrounded at all. Each value may be any valid floating-point number, including ones
	 * in scientific notation (e.g. 2.1e2), as well as hexadecimal (0xD2) and octal (0322).
	 * 
	 * @param weightString A string of sequence weights.
	 * 
	 * @throws IllegalArgumentException If <code>weightString</code> is not properly formatted or
	 *   contains values that cannot be converted to a floating-point number.
	 * 
	 * @see AlignmentWeight#DEFAULT_MATCH
	 * @see AlignmentWeight#DEFAULT_MISMATCH
	 * @see AlignmentWeight#DEFAULT_GAP_OPEN
	 * @see AlignmentWeight#DEFAULT_GAP_EXTEND
	 * @see AlignmentWeight#DEFAULT_INIT_SCORE
	 */
	public void setAlignmentWeight(String weightString)
			throws IllegalArgumentException {
		
		alignmentWeight = AlignmentWeight.get(weightString);
	}
	
	/**
	 * Set the weight of aligned sequences that match.
	 * 
	 * @param match Weight for aligned bases that match. If negative, this weight is automatically
	 *   converted to a positive number.
	 * 
	 * @throws IllegalArgumentException If <code>match</code> is <code>0</code>.
	 * 
	 * @see AlignmentWeight#getWithMatch(float)
	 */
	public void setAlignmentWeightMatch(float match)
			throws IllegalArgumentException {
		
		alignmentWeight = alignmentWeight.getWithMatch(match);
	}
	
	/**
	 * Set the weight of aligned sequences that do not match.
	 * 
	 * @param mismatch Weight for aligned bases that do not match. If positive, this weight is
	 *   automatically converted to a negative number.
	 * 
	 * @throws IllegalArgumentException If <code>mismatch</code> is <code>0</code>.
	 * 
	 * @see AlignmentWeight#getWithMismatch(float)
	 */
	public void setAlignmentWeightMismatch(float mismatch)
			throws IllegalArgumentException {
		
		alignmentWeight = alignmentWeight.getWithMismatch(mismatch);
	}
	
	/**
	 * Set the weight of opening a gap in the alignment.
	 * 
	 * @param gapOpen Weight of opening a gap in the alignment. If positive, this weight is
	 *   automatically converted to a negative number.
	 * 
	 * @see AlignmentWeight#getWithGapOpen(float)
	 */
	public void setAlignmentWeightGapOpen(float gapOpen) {
		
		alignmentWeight = alignmentWeight.getWithGapOpen(gapOpen);
	}
	
	/**
	 * Set the weight of extending a gap in the alignment.
	 * 
	 * @param gapExtend Weight of extending a gap in the alignment. If positive, this weight
	 *   is automatically converted to a negative number.
	 * 
	 * @throws IllegalArgumentException If <code>gapExtend</code> is <code>0</code>.
	 * 
	 * @see AlignmentWeight#getWithGapExtend(float)
	 */
	public void setAlignmentWeightGapExtend(float gapExtend)
			throws IllegalArgumentException {
		
		alignmentWeight = alignmentWeight.getWithGapExtend(gapExtend);
	}
	
	/**
	 * Get a structure containing sequence alignment weights.
	 * 
	 * @return Sequence alignment weights.
	 */
	public AlignmentWeight getAlignmentWeight() {
		return alignmentWeight;
	}
	
	/**
	 * Get the class loader for this runner. If any external libraries were added, this loader will search them.
	 * 
	 * @return This runner&apos;s class loader.
	 */
	public URLClassLoader getLoader() {
		return loader;
	}
	
	/**
	 * Set the temporary directory (location where segment files are offloaded).
	 * If <code>null</code>, the default temporary directory (current directory)
	 * is used.
	 * 
	 * @param tempDirName Name of the temporary directory.
	 */
	public void setTempDirName(String tempDirName) {
		this.tempDirName = tempDirName;
		
		return;
	}
	
	/**
	 * Get temporary directory name.
	 * 
	 * @return Temporary directory name. If an empty string or <code>null</code>,
	 *   the default temporary directory is used.
	 */
	public String getTempDirName() {
		return tempDirName;
	}
	
	/**
	 * Set the minimum k-mer count. K-mers with a lower count are filtered out while processing samples.
	 * 
	 * @param minKmerCount Minimum k-mer count.
	 * 
	 * @throws IllegalArgumentException If <code>minKmerCount</code> is negative.
	 */
	public void setMinKmerCount(int minKmerCount) {
		
		if (minKmerCount < 0)
			throw new IllegalArgumentException("K-mer count must not be less than 0: " + minKmerCount);
		
		this.minKmerCount = minKmerCount;
		
		return;
	}
	
	/**
	 * Get the minimum k-mer count. K-mers with a lower count are filtered out while processing samples.
	 * 
	 * @return Minimum k-mer count.
	 */
	public int getMinKmerCount() {
		return minKmerCount;
	}
	
	/**
	 * Set the minimum difference between neighboring k-mers to trigger a scan for an
	 * active region.
	 *  
	 * @param minimumDifference Minimum difference between two k-mers to trigger a scan
	 *   for an active region. 
	 * 
	 * @throws IllegalArgumentException If <code>minimumDifference</code> is less than
	 *   <code>1</code>.
	 */
	public void setMinimumDifference(int minimumDifference)
			throws IllegalArgumentException {
		
		if (minimumDifference < 1)
			throw new IllegalArgumentException("Cannotset minimum difference (to trigger a correction attempt) to a value less than 1: " + minimumDifference);
		
		this.minimumDifference = minimumDifference;
		
		return;
	}
	
	/**
	 * Get the minimum difference between neighboring k-mers to trigger a scan for an
	 * active region.
	 * 
	 * @return Minimum difference between neighboring k-mers to trigger a scan for an
	 *   active region.
	 * 
	 * @see #setMinimumDifference(int)
	 */
	public int getMinimumDifference() {
		return minimumDifference;
	}
	
	/**
	 * Set the quantile used to dynamically dermine k-mer count difference thresholds that trigger
	 * a correction attempt. If greater than <code>0.0</code>, determine the count difference
	 * threshold that by taking the difference between each neighboring k-mer in the reference
	 * sequence and finding this quantile over that set of differences. For example, if set to
	 * <code>0.95</code>, then at most 5% (100% - 95%) differences will trigger a correction attempt.
	 * If this value, once computed over the reference, is less than the minimum k-mer count difference,
	 * then it is ignored and that minimum is the count difference threshold. If <code>0.0</code>,
	 * then no quantile is performed and the minimum is the threshold. Must not be less than
	 * <code>0.0</code> or greater than (or equal to) <code>1.0</code>.
	 * 
	 * @param kmerDiffQuantile K-mer count difference quantile.
	 * 
	 * @throws IllegalArgumentException If <code>kmerDiffQuantile</code> is less than
	 *   <code>0.0</code> or is greater than (or equal to) <code>1.0</code>.
	 */
	public void setDifferenceQuantile(double kmerDiffQuantile)
			throws IllegalArgumentException {
		
		if (kmerDiffQuantile < 0.0F || kmerDiffQuantile >= 1.0F)
			throw new IllegalArgumentException("K-mer difference threshold is out of range [0.0, 1.0): " + kmerDiffQuantile);
		
		this.differenceQuantile = kmerDiffQuantile;
		
		return;
	}
	
	/**
	 * Get the k-mer count difference quantile used to set the difference threshold for triggerring
	 * correction attempts.
	 * 
	 * @return K-mer count difference quantile.
	 */
	public double getDifferenceQuantile() {
		return differenceQuantile;
	}
	
	/**
	 * Set the k-mer anchoring flag. If this value is set to <code>true</code>,
	 * an active region must be bordered on each end by an unaltered k-mer. If
	 * <code>false</code>, an active region may extend to the end of the sequence.
	 * 
	 * @param anchorBothEnds Anchor flag.
	 */
	public void setAnchorBothEnds(boolean anchorBothEnds) {
		this.anchorBothEnds = anchorBothEnds;
		
		return;
	}
	
	/**
	 * Get the k-mer anchoring flag. If this value is set to <code>true</code>,
	 * an active region must be bordered on each end by an unaltered k-mer. If
	 * <code>false</code>, an active region may extend to the end of the sequence.
	 * 
	 * @return Anchor flag.
	 */
	public boolean getAnchorBothEnds() {
		return anchorBothEnds;
	}
	
	/**
	 * Set the property to keep k-mer counts in memory. If set to <code>true</code>,
	 * k-mer counts from each sample are not offloaded to an indexed k-mer count file.
	 * Setting this option to <code>true</code> should only be done when all samples
	 * are relatively small or the system has sufficient memory to handle them.
	 * 
	 * @param kmerCountInMemory The property to keep k-mer counts in memory.
	 */
	public void setKmerCountInMemory(boolean kmerCountInMemory) {
		this.kmerCountInMemory = kmerCountInMemory;
		
		return;
	}
	
	/**
	 * Get the property to keep k-mer counts in memory.
	 * 
	 * @return <code>true</code> if k-mer counts from each sample are stored in memory.
	 */
	public boolean isKmerCountInMemory() {
		return kmerCountInMemory;
	}
	
	/**
	 * Free resources as soon as possible is set to <code>true</code>. This may force expensive
	 * resource allocation between samples.
	 * 
	 * @param freeResources Free-resources property.
	 */
	public void setFreeResources(boolean freeResources) {
		this.freeResources = freeResources;
		
		return;
	}
	
	/**
	 * Free resources as soon as possible is set to <code>true</code>. This may force expensive
	 * resource allocation between samples.
	 * 
	 * @return Free-resources property.
	 */
	public boolean isFreeResources() {
		return freeResources;
	}
	
	/**
	 * Set the property to remove indexed k-mer count (IKC) files for each sample as Kestrel
	 * finishes.
	 * 
	 * @param removeIkc <code>true</code> if IKC files should be removed.
	 */
	public void setRemoveIkc(boolean removeIkc) {
		this.removeIkc = removeIkc;
		
		return;
	}
	
	/**
	 * Get the property to remove indexed k-mer count (IKC) files for each sample as Kestrel
	 * finishes.
	 * 
	 * @return <code>true</code> if IKC files should be removed.
	 */
	public boolean isRemoveIck() {
		return removeIkc;
	}
	
	/**
	 * Set the property to count reverse complement k-mers in sequencing sets. This
	 * should be set whenever the sequence data may be in the direction of the reference
	 * or its reverse complement. If the sequencing data should not contain reverse
	 * complement data, which is true of some RNASeq protocols, then this should be
	 * set to <code>false</code>.
	 * <p/>
	 * This option affects all read-depth estimates based on k-mer counts.
	 * 
	 * @param countReverseKmers Count k-mers and their reverse complements in region
	 *   statistics.
	 */
	public void setCountReverseKmers(boolean countReverseKmers) {
		this.countReverseKmers = countReverseKmers;
	}
	
	/**
	 * Get the property to count reverse complement k-mers while generating read depth
	 * estimates for regions and variants.
	 * 
	 * @return <code>true</code> if reverse complement k-mers are counted in read depth
	 *   estimates.
	 */
	public boolean getCountReverseKmers() {
		return countReverseKmers;
	}
	
	/**
	 * If set to <code>true</code>, allow active regions to include ambiguous bases.
	 * 
	 * @param callAmbiguousRegions Allow ambiguous bases in active regions if
	 *   <code>true</code>.
	 */
	public void setCallAmbiguousRegions(boolean callAmbiguousRegions) {
		this.callAmbiguousRegions = callAmbiguousRegions;
		
		return;
	}
	
	/**
	 * Get the property allow ambiguous bases in active regions.
	 *  
	 * @return <code>true</code> if active regions may contain ambiguous bases.
	 * 
	 * @see #setCallAmbiguousRegions(boolean)
	 */
	public boolean getCallAmbiguousRegions() {
		return callAmbiguousRegions;
	}
	
	/**
	 * Set the property to allow variants over ambiguous bases in the reference.
	 * 
	 * @param allowAmbiguousReference <code>true</code> to allow variants over ambiguous
	 *   bases in the reference.
	 */
	public void setCallAmbiguousVariant(boolean allowAmbiguousReference) {
		this.callAmbiguousVariant = allowAmbiguousReference;
		
		return;
	}
	
	/**
	 * Get the property to allow variants over ambiguous bases in the reference.
	 * 
	 * @return The property to allow variants over ambiguous bases in the reference.
	 * 
	 * @see #setCallAmbiguousVariant(boolean)
	 */
	public boolean getCallAmbiguousVariant() {
		return callAmbiguousVariant;
	}
	
	/**
	 * Set the number of k-mers to scan when detecting and stepping over a suspected peak.
	 * A peak can occur when a k-mer matches another region of the genome and the k-mer counts
	 * from that region are added. If set to <code>0</code>, peak detection is disabled.
	 * 
	 * @param peakScanLength Number of k-mers to scan when detecting a peak or <code>0</code>
	 *   to disable peak detection.
	 * 
	 * @throws IllegalArgumentException If <code>peakScanLength</code> is negative.
	 */
	public void setPeakScanLength(int peakScanLength)
			throws IllegalArgumentException {
		
		if (peakScanLength < 0)
			throw new IllegalArgumentException("Peak scan length is negative: " + peakScanLength);
		
		this.peakScanLength = peakScanLength;
		
		return;
	}
	
	/**
	 * Get the number of k-mers to scan during peak detection.
	 * 
	 * @return The number of k-mers to scan during peak detection.
	 * 
	 * @see #setPeakScanLength(int)
	 */
	public int getPeakScanLength() {
		return peakScanLength;
	}
	
	/**
	 * Limit the length of an active region scan. This is computed by taking maximum length of a
	 * gap (assuming no other variants), adding one, and then multiplying by this factor. If the
	 * limit is less than the k-mer size, then the k-mer size is the limit. Setting a large value
	 * does not limit the length of active regions.
	 * 
	 * @param scanLimitFactor The end-scan limit factor.
	 * 
	 * @throws IllegalArgumentException If <code>endScanLimit</code> is negative.
	 */
	public void setScanLimitFactor(double scanLimitFactor)
			throws IllegalArgumentException {
		
		if (scanLimitFactor < 0.0)
			throw new IllegalArgumentException("Scan limit factor may not be negative: " + scanLimitFactor);
		
		this.scanLimitFactor = scanLimitFactor;
		
		return;
	}
	
	/**
	 * Get the end-scan limit.
	 * 
	 * @return End-scan limit.
	 * 
	 * @see #setScanLimitFactor(double)
	 */
	public double getScanLimitFactor() {
		return scanLimitFactor;
	}
	

	/**
	 * Set the exponential decay minimum. This is the minimum value (asymptotic lower bound)
	 * of the exponential decay function as a proportion of the anchor k-mer count. If this
	 * value is <code>0.0</code>, k-mer count recovery threshold may decline to <code>0</code>.
	 * If this value is <code>1.0</code>, the decay function is not used and the detector falls
	 * back to finding a k-mer with a count within the difference threshold of the anchor k-mer
	 * count.
	 * 
	 * @param expDecayMin Exponential decay minimum. Must be between <code>0.0</code>
	 *   and <code>1.0</code> (inclusive).
	 * 
	 * @throws IllegalArgumentException If <code>expDecayMin</code> is less than <code>0.0</code>
	 *   or greater than <code>1.0</code>.
	 */
	public void setDecayMinimum(double expDecayMin)
			throws IllegalArgumentException {
		
		if (expDecayMin < 0.0 || expDecayMin > 1.0)
			throw new IllegalArgumentException("Exponential decay minimum must not be negative or greater than 1.0: " + expDecayMin);
		
		this.expDecayMin = expDecayMin;
		
		return;
	}
	
	/**
	 * Get the exponential decay minimum as a proportion of the anchor k-mer count.
	 * 
	 * @return Exponential decay minimum.
	 * 
	 * @see #setDecayMinimum(double)
	 */
	public double getDecayMinimum() {
		return expDecayMin;
	}
	
	/**
	 * Set the exponential decay alpha, which controls how quickly the recovery threshold
	 * declines. Alpha is defined as the proportion of decay from the lower bound to the
	 * minimum count at <code>k</code> k-mers from the anchor. A lower value causes a steeper
	 * decline toward the minimum threshold value.
	 * 
	 * @param expDecayAlpha Exponential decay alpha. Must be between <code>0.0</code>
	 *   and <code>1.0</code> (exclusive).
	 * 
	 * @throws IllegalArgumentException If <code>expDecayMin</code> is not between <code>0.0</code>
	 *   and <code>1.0</code> (exclusive).
	 * 
	 * @see ActiveRegionDetector#DEFAULT_EXP_ALPHA
	 */
	public void setDecayAlpha(double expDecayAlpha)
			throws IllegalArgumentException {
		
		if (expDecayAlpha <= 0.0 || expDecayAlpha >= 1.0)
			throw new IllegalArgumentException("Exponential decay alpha must be between 0.0 and 1.0 (exclusive): " + expDecayAlpha);
		
		this.expDecayAlpha = expDecayAlpha;
		
		return;
	}
	
	/**
	 * Get the exponential decay alpha.
	 * 
	 * @return Exponential decay alpha.
	 * 
	 * @see #setDecayAlpha(double)
	 */
	public double getDecayAlpha() {
		return expDecayAlpha;
	}
	
	/**
	 * Add a variant filter by its string specification. The filter specification is a filter
	 * name followed by a list of arguments separated by a colon. Whitespace around the colon
	 * is discarded. The arguments are typically a comma-separated list of attribute/value
	 * pairs, however, it is up to the filter implementation to interpret the arguments.
	 *  
	 * @param filterSpec Filter specification.
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
	public void addVariantFilter(String filterSpec)
			throws NullPointerException, IllegalArgumentException, FileNotFoundException, IOException, VariantFilterInitException {
		
		variantFilterList.add(VariantFilter.getFilter(filterSpec, loader));  // throws NullPointerException, IllegalArgumentException, FileNotFoundException, IOException, VariantFilterInitException
		
		return;
	}
	
	/**
	 * Remove variant filters from this runner.
	 */
	public void clearVariantFilters() {
		variantFilterList.clear();
	}
	
	/**
	 * Add an interval to call variants on. If no intervals are added, variants are called
	 * over the whole of each reference sequence.
	 * 
	 * @param interval Interval to add.
	 * 
	 * @throws NullPointerException If <code>interval</code> is <code>null</code>.
	 * @throws IllegalArgumentException If <code>interval</code> overlaps with another interval.
	 * @throws IllegalStateException If the maximum number of intervals that can be associated with
	 *   a reference sequence has been reached. This is the maximum array size.
	 */
	public void addInterval(RegionInterval interval)
			throws NullPointerException, IllegalArgumentException, IllegalStateException {
		
		intervalContainer.add(interval);
		
		return;
	}
	
	/**
	 * Read intervals from a file.
	 * 
	 * @param intervalFile Interval file.
	 * @param format Format of the file or &quot;auto&quot; to automatically detect the file
	 *   type using the file name. If <code>null</code> or empty, &quot;auto&quot; is assumed.
	 *   this type string is not case sensitive.
	 * @param args Arguments for the interval file reader. These are defined by the interval reader,
	 *   and they are normally key=value pairs (separated by an equal sign) with pairs separated
	 *   by commas. May be <code>null</code> or empty to indicate that there are no options for
	 *   the reader.
	 * 
	 * @return The number of intervals added.
	 * 
	 * @throws NullPointerException If <code>intervalFile</code> is <code>null</code>.
	 * @throws FileNotFoundException If <code>intervalFile</code> cannot be found.
	 * @throws IOException If there is any error reading <code>intervalFile</code> including
	 *   errors in the file contents.
	 * @throws IntervalReaderInitException If the interval reader for this file cannot be found
	 *   or if there is an error initializing it.
	 */
	public int readIntervals(File intervalFile, String format, String args)
			throws NullPointerException, FileNotFoundException, IOException, IntervalReaderInitException {
		
		RegionInterval[] intervals = IntervalReader.read(intervalFile, format, args, loader);
		
		addInterval(intervals);
		
		return intervals.length;
	}
	
	/**
	 * Add an array of intervals to call variants on. If no intervals are added, variants are
	 * called over the whole of each reference sequence.
	 * 
	 * @param intervals Interval array to add.
	 * 
	 * @throws NullPointerException If <code>intervals</code> is <code>null</code>.
	 * @throws IllegalArgumentException If an interval in the array overlaps with another interval.
	 * @throws IllegalStateException If the maximum number of intervals that can be associated with
	 *   a reference sequence has been reached. This is the maximum array size.
	 */
	public void addInterval(RegionInterval[] intervals)
			throws NullPointerException, IllegalArgumentException, IllegalStateException {
		
		if (intervals == null)
			throw new NullPointerException("Cannot add intervals from array: null");
		
		for (int index = 0; index < intervals.length; ++index)
			if (intervals[index] != null)
				intervalContainer.add(intervals[index]);
		
		return;
	}
	
	/**
	 * Clear all intervals.
	 */
	public void clearIntervals() {
		intervalContainer.clear();
		
		return;
	}
	
	/**
	 * Clear all intervals associated with a reference sequence.
	 * 
	 * @param referenceName Reference sequence name. If <code>null</code>, no action
	 *   is taken.
	 */
	public void clearIntervals(String referenceName) {
		intervalContainer.clear(referenceName);
		
		return;
	}
	
	/**
	 * Set the length of flanks to add to intervals where variants are called. Adding flanks
	 * improves active region detection, and therefore variant recovery, but it does not otherwise
	 * affect variants.
	 * <p/>
	 * If a flank length is not set, Kestrel determines the length by multiplying the k-mer size
	 * by <code>DEFAULT_FLANK_LENGTH_MULTIPLIER</code>.
	 * 
	 * @param flankLength Length of flanks. Set <code>0</code> to disable flanks and <code>-1</code>
	 *   for Kestrel to determine the appropriate length when it runs.
	 * @throws IllegalArgumentException If <code>flankLength</code> is less than <code>-1</code>.
	 * 
	 * @see #DEFAULT_FLANK_LENGTH_MULTIPLIER
	 */
	public void setFlankLength(int flankLength)
			throws IllegalArgumentException {
		
		if (flankLength < -1)
			throw new IllegalArgumentException("Flank length must be 0 or greater (or -1 to automatically determine the length): " + flankLength);
		
		this.flankLength = flankLength;
		
		return;
	}
	
	/**
	 * Kestrel will automatically determine the correct flank length. This is the default behavior.
	 * This function can undo calls to <code>setFlankLength()</code>.
	 * 
	 * @see #setFlankLength(int)
	 * @see #DEFAULT_FLANK_LENGTH_MULTIPLIER
	 */
	public void setDefaultFlankLength() {
		flankLength = -1;
		
		return;
	}
	
	/**
	 * When set, variant calls are relative to the reference sequence. This
	 * will undo a call to <code>setVariantCallByRegion()</code>.
	 */
	public void setVariantCallByReference() {
		
		this.variantCallByRegion = false;
		
		return;
	}
	
	/**
	 * Determine if variant calls are relative to the reference sequence.
	 * 
	 * @return <code>true</code> if variant calls are relative to the reference
	 *   sequence.
	 */
	public boolean isVariantCallByReference() {
		return ! variantCallByRegion;
	}
	
	/**
	 * When set, variant calls are relative to the reference region. This
	 * will undo a call to <code>setVariantCallByReference()</code>.
	 */
	public void setVariantCallByRegion() {
		
		this.variantCallByRegion = true;
		
		return;
	}
	
	/**
	 * Determine if variant calls are relative to reference regions.
	 * 
	 * @return <code>true</code> if variant calls are relative to reference
	 *   regions.
	 */
	public boolean isVariantCallByRegion() {
		return variantCallByRegion;
	}
	
	/**
	 * Set the maximum number of alignment states the aligner may save. The least-likely states
	 * are trimmed when this threshold is reached.
	 * 
	 * @param maxAlignerState Maximum number of states.
	 * 
	 * @throws IllegalArgumentException If <code>maxAlignerState</code> is less than <code>1</code>.
	 * 
	 * @see KmerAligner#DEFAULT_MAX_STATE
	 */
	public void setMaxAlignerState(int maxAlignerState)
		throws IllegalArgumentException {
		
		if (maxAlignerState < 1)
			throw new IllegalArgumentException("Maximum number of aligner states must not be less than 1: " + maxAlignerState);
		
		this.maxAlignerState = maxAlignerState;
		
		return;
	}
	
	/**
	 * Get the maximum number of states that may be saved.
	 * 
	 * @return The maximum number of states that may be saved.
	 * 
	 * @see #setMaxAlignerState(int)
	 */
	public int getMaxAlignerState() {
		return maxAlignerState;
	}
	
	/**
	 * Set the maximum number of haplotypes the aligner may save. The least-likely states
	 * are trimmed when this threshold is reached.
	 * 
	 * @param maxHaplotypes Maximum number of haplotypes.
	 * 
	 * @throws IllegalArgumentException If <code>maxHaplotypes</code> is less than <code>1</code>.
	 * 
	 * @see KmerAlignmentBuilder#DEFAULT_MAX_HAPLOTYPES
	 */
	public void setMaxHaplotypes(int maxHaplotypes)
		throws IllegalArgumentException {
		
		if (maxHaplotypes < 1)
			throw new IllegalArgumentException("Maximum number of haplotypes not be less than 1: " + maxAlignerState);
		
		this.maxHaplotypes = maxHaplotypes;
		
		return;
	}
	
	/**
	 * Get the maximum number of states that may be saved.
	 * 
	 * @return The maximum number of states that may be saved.
	 * 
	 * @see #setMaxAlignerState(int)
	 */
	public int getMaxHaplotypes() {
		return maxHaplotypes;
	}
	
	/**
	 * Set the property to remove the  reference sequence description from each sequence name.
	 * The sequence name is defined as everything up to the first whitespace character, and the
	 * description follows. When this option is enabled, the whitespace and description are
	 * removed.
	 * 
	 * @param removeReferenceSequenceDescription Remove reference sequence description parameter.
	 * 
	 * @see ReferenceReader#DEFAULT_REMOVE_SEQUENCE_DESCRIPTION
	 */
	public void setRemoveReferencDescription(boolean removeReferenceSequenceDescription) {
		this.removeReferenceSequenceDescription = removeReferenceSequenceDescription;
		
		return;
	}
	
	/**
	 * Get the property to remove reference sequence descriptions.
	 * 
	 * @return The remove sequence description property.
	 * 
	 * @see #setRemoveReferencDescription(boolean)
	 */
	public boolean getRemoveReferenceDescription() {
		return removeReferenceSequenceDescription;
	}
	
	/**
	 * When set to <code>true</code>, the &quot;reverse complement negative strand&quot; property reverse
	 * complements any reference regions defined on a negative strand interval. The interval
	 * object has an attribute delineating strandedness, and this is usually set by reading
	 * a file with this information (e.g. BED file).
	 * 
	 * @param reverseComplementNegativeStrand Reverse complement negative strand regions if <code>true</code>.
	 *   Otherwise, regions are left in the orientation found in the reference sequence.
	 * 
	 * @see #DEFAULT_REVERSE_COMPLEMENT_NEGATIVE_STRAND
	 */
	public void setRevComplementNegReferenceStrand(boolean reverseComplementNegativeStrand) {
		this.reverseComplementNegativeStrand = reverseComplementNegativeStrand;
	}
	
	/**
	 * Get the property to reverse complement negative strand regions.
	 * 
	 * @return The &quot;Reverse complement negative strand&quot; property.
	 * 
	 * @see #setRevComplementNegReferenceStrand(boolean)
	 */
	public boolean getRevComplementNegReferenceStrand() {
		return reverseComplementNegativeStrand;
	}
	
	/**
	 * Set the file haplotypes will be output to. 
	 * 
	 * @param outputObj Output object identifying where the output is
	 *   written. This object must be of class <code>String</code>,
	 *   <code>File</code>, or <code>FileDescriptor</code>, or
	 *   <code>StreamableOutput</code>. If <code>null</code> or an empty string,
	 *   the default output file, <code>DEFAULT_OUTPUT_FILE</code>, is used.
	 * 
	 * @throws IllegalArgumentException If <code>outputObj</code> is not
	 *   <code>null</code> and is not an instance of <code>String</code>,
	 *   <code>File</code>, <code>FileDescriptor</code>, or <code>StreamableOutput</code>.
	 * 
	 * @see #DEFAULT_OUTPUT_FILE
	 */
	public void setHaplotypeOutputFile(Object outputObj)
			throws IllegalArgumentException {

		// Check null
		if (outputObj == null) {
			haplotypeOutputFile = null;
			
			return;
		}
		
		// Check empty string
		if (outputObj instanceof String) {
			if (((String) outputObj).trim().isEmpty()) {
				haplotypeOutputFile = null;
				
				return;
			}
		}
		
		// Set output file
		haplotypeOutputFile = StreamableOutput.get(outputObj);  // throws IllegalArgumentException 
		
		return;
	}
	
	/**
	 * Get the output object for resolved haplotypes or <code>null</code> if haplotypes
	 * are not output.
	 * 
	 * @return Haplotype output object.
	 */
	public StreamableOutput getHaplotypeOutputFile() {
		return haplotypeOutputFile;
	}
	
	/**
	 * Set the output format for haplotypes if they are written. If the haplotype output object
	 * is <code>null</code>, this parameter has no effect.
	 * 
	 * @param haplotypeOutputFormat Output format.
	 * 
	 * @throws NullPointerException If <code>haplotypeOutputFormat</code> is <code>null</code>.
	 * @throws IllegalArgumentException If <code>haplotypeOutputFormat</code> is empty.
	 */
	public void setHaplotypeOutputFormat(String haplotypeOutputFormat)
			throws NullPointerException, IllegalArgumentException {
		
		// Check arguments
		if (haplotypeOutputFormat == null)
			throw new NullPointerException("Cannot set haplotype output format: null");
		
		haplotypeOutputFormat = haplotypeOutputFormat.trim();
		
		if (haplotypeOutputFormat.isEmpty()) {
			throw new IllegalArgumentException("Cannot set haplotype output format to an empty string");
		}
		
		this.haplotypeOutputFormat = haplotypeOutputFormat;
	}
	
	public String getHaplotypeOutputFormat() {
		return haplotypeOutputFormat;
	}
	
	/**
	 * Create the count module and configure parameters from this runner. The caller
	 * will need to set the output location and type.
	 * 
	 * @return A configured module (except output location and type).
	 */
	protected CountModule getCountModule() {
		CountModule countModule = new CountModule();
		
		countModule.configure(null);
		
		// Add libraries to KAnalyze
		for (URL libUrl : libraryUrlList)
			countModule.addLibraryURL(libUrl);
		
		// Set k-mer size
		countModule.setKSize(kSize);
		
		// Set options
		countModule.setTempDirName(tempDirName);
		
		if (minKmerCount > 0)
			countModule.addPostCountFilterDefinition("kmercount:" + minKmerCount);
		
		countModule.setFreeSegment(freeResources);
		
		countModule.addListener(new StreamConditionListener("Kestrel", System.err, null, true));  // TODO: Replace with a Kestrel condition listener
		
		return countModule;
	}
	
	/**
	 * Get a list of all errors that have occurred when running the argument parser.
	 * Calling this method does not reset the error condition.
	 * 
	 * @return A list of all errors that have occurred, or an empty list if there are no errors to report.
	 */
	public ConditionEvent[] getErrorConditions() {
		return argumentParser.getErrorConditions();
	}
	
	/**
	 * Get a list of all warnings that have occurred when running the argument parser.
	 * Calling this method does not reset the warning condition.
	 * 
	 * @return A list of all warnings that have occurred, or an empty list if there are no warnings to report.
	 */
	public ConditionEvent[] getWarningConditions() {
		return argumentParser.getWarningConditions();
	}
}
