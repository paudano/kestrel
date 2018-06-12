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
import java.nio.charset.IllegalCharsetNameException;
import java.nio.charset.UnsupportedCharsetException;
import java.util.ArrayList;
import java.util.List;

import edu.gatech.kanalyze.KAnalyzeConstants;
import edu.gatech.kanalyze.comp.reader.FileSequenceSource;
import edu.gatech.kanalyze.comp.reader.SequenceReader;
import edu.gatech.kanalyze.comp.reader.SequenceSource;
import edu.gatech.kanalyze.condition.ConditionEvent;
import edu.gatech.kanalyze.condition.ConditionListener;
import edu.gatech.kanalyze.util.StringUtil;
import edu.gatech.kanalyze.util.argparse.ArgumentParser;
import edu.gatech.kanalyze.util.argparse.HelpTopicElement;
import edu.gatech.kanalyze.util.argparse.NonOptionElement;
import edu.gatech.kanalyze.util.argparse.OptionArgumentType;
import edu.gatech.kanalyze.util.argparse.OptionSpecElement;
import edu.gatech.kestrel.KestrelConstants;
import edu.gatech.kestrel.LogLevel;
import edu.gatech.kestrel.activeregion.ActiveRegionDetector;
import edu.gatech.kestrel.align.AlignmentWeight;
import edu.gatech.kestrel.align.KmerAligner;
import edu.gatech.kestrel.align.KmerAlignmentBuilder;
import edu.gatech.kestrel.interval.IntervalReaderInitException;
import edu.gatech.kestrel.io.InputSample;
import edu.gatech.kestrel.io.StreamableOutput;
import edu.gatech.kestrel.refreader.ReferenceReader;
import edu.gatech.kestrel.variant.VariantCaller;

/**
 * Parses arguments for the Kestrel command-line interface.
 */
public class KestrelArgumentParser extends ArgumentParser {
	
	/** The Kestrel runner object this parser configures. */
	protected KestrelRunnerBase runnerBase;
	
	/** Currently set format type for input files. */
	private String formatType;
	
	/** Currently set character-set for input files. */
	private Charset charset;
	
	/** Currently set filter specification for input files. */
	private String filterSpec;
	
	/** Sequence source index. Automatically incremented as sources are added. */
	private int seqSourceIndex;
	
	/** Reference source index. Automatically incremented as reference sequences are added. */
	private int refSourceIndex;
	
	/** List of errors generated while processing command-line arguments. */
	public List<ConditionEvent> errorList;
	
	/** List of warnings generated while processing command-line arguments. */
	private List<ConditionEvent> warningList;
	
	/**
	 * List of queued sequence sources. These are accumulated until they can be grouped
	 * into an <code>InputSample</code>.
	 */
	private List<SequenceSource> inputSourceQueue;
	
	/**
	 * If a name is explicitly assigned to the sample being accumulated in <code>inputSourceQueue</code>,
	 * it is assigned here. If there is no name, then this value is <code>null</code> and the name will be
	 * automatically assigned by name of the first file assigned to the sample.
	 */
	private String sampleName;
	
	/**
	 * If greater than <code>0</code>, then this many input files automatically becomes one input sample.
	 * Setting this to <code>2</code> for paired-end reads will automatically assign every 2 files to
	 * one input sample without using multiple -s arguments on the command-line.
	 */
	private int filesPerSample;
	
	
	/**
	 * Create a new Kestrel command-line argument parser.
	 * 
	 * @param runnerBase Kestrel runner this parser configures.
	 * 
	 * @throws NullPointerException If <code>runner</code> is <code>null</code>.
	 */
	public KestrelArgumentParser(KestrelRunnerBase runnerBase)
			throws NullPointerException {
		
		super(KestrelConstants.PROG_NAME, KestrelConstants.VERSION);
		
		// Check arguments
		if (runnerBase == null)
			throw new NullPointerException("Parser cannot configure kestrel runner: null");
		
		// Set fields
		this.runnerBase = runnerBase;
		formatType = "auto";
		inputSourceQueue = new ArrayList<SequenceSource>();
		
		// Set condition listener and condition lists
		errorList = new ArrayList<ConditionEvent>();
		warningList = new ArrayList<ConditionEvent>();
		
		clearListeners();
		addListener(new ConfigurationConditionListener());
		
		// Add options
		addSpecification(new OptAddVariantFilter());
		addSpecification(new OptAutoFlankLength());
		addSpecification(new OptAlignmentWeightVector());
		addSpecification(new OptAnchorBothEnds());
		addSpecification(new OptCallAmbiguousRegions());
		addSpecification(new OptCallAmbiguousVariant());
		addSpecification(new OptVarCallRelativeReference());
		addSpecification(new OptVarCallRelativeRegion());
		addSpecification(new OptCharset());
		addSpecification(new OptCountReverseKmers());
		addSpecification(new OptDecayAlpha());
		addSpecification(new OptDecayMinimum());
		addSpecification(new OptScanLimitFactor());
		addSpecification(new OptFilesPerSample());
		addSpecification(new OptFormat());
		addSpecification(new OptFreeResources());
		addSpecification(new OptHaplotypeOutputFileName());
		addSpecification(new OptHaplotypeOutputFormat());
		addSpecification(new OptKmerCountDiffQuantile());
		addSpecification(new OptKmerCountInMemory());
		addSpecification(new OptKmerSize());
		addSpecification(new OptLoadLibraryFile());
		addSpecification(new OptLoadLibraryUrl());
		addSpecification(new OptLogFileName());
		addSpecification(new OptLogLevel());
		addSpecification(new OptLogStderr());
		addSpecification(new OptLogStdout());
		addSpecification(new OptMaxAlignStates());
		addSpecification(new OptMaxHaplotypeStates());
		addSpecification(new OptMaxRepeatCount());
		addSpecification(new OptMinimizerMask());
		addSpecification(new OptMinimizerSize());
		addSpecification(new OptMinKmerCount());
		addSpecification(new OptMinKmerCountDiff());
		addSpecification(new OptNoAnchorBothEnds());
		addSpecification(new OptNoCallAmbiguousRegions());
		addSpecification(new OptNoCallAmbiguousVariant());
		addSpecification(new OptNoCountReverseKmers());
		addSpecification(new OptNoFreeResources());
		addSpecification(new OptNoKmerCountInMemory());
		addSpecification(new OptNoRemoveRefDescription());
		addSpecification(new OptNoRevComplNegRegStrand());
		addSpecification(new OptNoRmIkc());
		addSpecification(new OptNoSequenceFilter());
		addSpecification(new OptOutFileName());
		addSpecification(new OptOutFormat());
		addSpecification(new OptPeakScanLength());
		addSpecification(new OptQuality());
		addSpecification(new OptReadIntervalFile());
		addSpecification(new OptReference());
		addSpecification(new OptRemoveRefDescription());
		addSpecification(new OptRevComplNegRegStrand());
		addSpecification(new OptRmIkc());
		addSpecification(new OptSampleSet());
		addSpecification(new OptSequenceFilter());
		addSpecification(new OptSetFlankLength());
		addSpecification(new OptTempFileLocation());
		addSpecification(new OptWriteStdout());
		
		setNonOptionElement(new InputSourceFile());
		
		addHelpTopic(new HelpFormat());
		addHelpTopic(new HelpReader());
		addHelpTopic(new HelpCite());
		addHelpTopic(new HelpBib());
		addHelpTopic(new HelpRefReader());
		
		return;
	}
	
	/**
	 * Reset fields not managed by options.
	 */
	@Override
	protected void initImplementation() {
		
		errorList.clear();
		warningList.clear();
		
		return;
	}
	
	/**
	 * Post-parsing tasks. Flushes the sample queue to ensure all sample files are read.
	 */
	@Override
	public void postParse() {
		
		if (filesPerSample > 0 && filesPerSample != inputSourceQueue.size()) {
		
			warn(
				String.format(
					"Files in sample is too small: Last sample contains %d files and the automatic group size (--filespersample) was set to %d",
					inputSourceQueue.size(),
					filesPerSample
				),
				KestrelConstants.ERR_USAGE
				);
		}
		
		flushSample();
		
		return;
	}
	
	//
	// Options
	//
	
	/**
	 * Option: Format type
	 */
	protected class OptFormat extends OptionSpecElement {
		
		/**
		 * Create option.
		 */
		public OptFormat() {
			super('f', "format",
					OptionArgumentType.REQUIRED,
					"INPUT_FORMAT", KestrelRunnerBase.DEFAULT_FORMAT,
					"Set the input sequence format type. This option determines how the " +
					"format files are read. This option may be set multiple times when reading " +
					"files with different formats. See \"count -hreader\" for a full list of " +
					"readers."
					);
		}
		
		/**
		 * Invoke this option
		 * 
		 * @param option Option.
		 * @param argument Option argument.
		 */
		@Override
		public boolean invoke(String option, String argument) {
			
			argument = argument.trim();
			
			if (argument.isEmpty()) {
				error("Cannot set format type (" + option + "): Type name is empty", KestrelConstants.ERR_USAGE);
				return false;
			}
			
			if (! argument.matches(KAnalyzeConstants.FORMAT_TYPE_PATTERN)) {
				error("Cannot set format type (" + option + "): Type is not alpha-numeric: " + argument, KestrelConstants.ERR_USAGE);
				return false;
			}
			
			formatType = argument;
			
			return true;
		}
		
		/**
		 * Initialize this option.
		 */
		@Override
		public void init() {
			formatType = KestrelRunnerBase.DEFAULT_FORMAT;
		}
	}
	
	/**
	 * Option: K-mer size
	 */
	protected class OptKmerSize extends OptionSpecElement {
		
		/**
		 * Create option.
		 */
		public OptKmerSize() {
			super('k', "ksize",
					OptionArgumentType.REQUIRED,
					"KSIZE", Integer.toString(KestrelRunnerBase.DEFAULT_KSIZE),
					"Size of k-mers sequence data is translated to during analysis. If unsure, use the " +
					"default value. If the sequencing error rate is very high, or if the reference is very " +
					"short, a small (e.g. a single short gene), then a smaller k-mer size, such as 21, may be " +
					"useful if the defalt value does not produce meaningful results."
					);
		}
		
		/**
		 * Invoke this option
		 * 
		 * @param option Option.
		 * @param argument Option argument.
		 */
		@Override
		public boolean invoke(String option, String argument) {
			
			argument = argument.trim();
			
			if (argument.isEmpty()) {
				error("Cannot set k-mer size (" + option + "): Size is empty", KestrelConstants.ERR_USAGE);
				return false;
			}
			
			try {
				runnerBase.setKSize(Integer.parseInt(argument));
				
			} catch (NumberFormatException ex) {
				error("Cannot set k-mer size (" + option + "): Size is not an integer: " + argument, KestrelConstants.ERR_USAGE);
				return false;
				
			} catch (IllegalArgumentException ex) {
				error("Cannot set k-mer size (" + option + "): " + ex.getMessage(), KestrelConstants.ERR_USAGE, ex);
				return false;
			}
			
			return true;
		}
		
		/**
		 * Initialize this option.
		 */
		@Override
		public void init() {
			runnerBase.setKSize(KestrelRunnerBase.DEFAULT_KSIZE);
		}
	}
	
	/**
	 * Option: Output file name
	 */
	protected class OptOutFileName extends OptionSpecElement {
		
		/**
		 * Create option.
		 */
		public OptOutFileName() {
			super('o', "out",
					OptionArgumentType.REQUIRED,
					"OUT_FILE", (KestrelRunnerBase.DEFAULT_OUTPUT_FILE != StreamableOutput.STDOUT ? KestrelRunnerBase.DEFAULT_OUTPUT_FILE.name : null),
					"Set output file name."
					);
		}
		
		/**
		 * Invoke this option
		 * 
		 * @param option Option.
		 * @param argument Option argument.
		 */
		@Override
		public boolean invoke(String option, String argument) {
			
			argument = argument.trim();
			
			if (argument.isEmpty()) {
				error("Cannot set output file name (" + option + "): Argument is empty", KestrelConstants.ERR_USAGE);
				return false;
			}
			
			try {
				runnerBase.setOutputFile(argument);
				
			} catch (IllegalArgumentException ex) {
				error("Cannot set output file name (" + option + "): " + ex.getMessage(), KestrelConstants.ERR_SYSTEM, ex);
				return false;
			}
			
			return true;
		}
		
		/**
		 * Initialize this option.
		 */
		@Override
		public void init() {
			runnerBase.setOutputFile(KestrelRunnerBase.DEFAULT_OUTPUT_FILE);
		}
	}
	
	/**
	 * Option: Write to STDOUT
	 */
	protected class OptWriteStdout extends OptionSpecElement {
		
		/**
		 * Create option.
		 */
		public OptWriteStdout() {
			super('\0', "stdout",
					OptionArgumentType.NONE,
					null, (KestrelRunnerBase.DEFAULT_OUTPUT_FILE == StreamableOutput.STDOUT ? "" : null),
					"Write output to standard output instead of a file. Unless redirected, this output " +
					"is written to the the screen."
					);
		}
		
		/**
		 * Invoke this option
		 * 
		 * @param option Option.
		 * @param argument Option argument.
		 */
		@Override
		public boolean invoke(String option, String argument) {
			
			runnerBase.setOutputFile(FileDescriptor.out);
			
			return true;
		}
		
		// Init by OptOutFileName
	}
	
	/**
	 * Option: Output format
	 */
	protected class OptOutFormat extends OptionSpecElement {
		
		/**
		 * Create option.
		 */
		public OptOutFormat() {
			super('m', "outfmt",
					OptionArgumentType.REQUIRED,
					"OUT_FORMAT", KestrelRunnerBase.DEFAULT_OUTPUT_FORMAT,
					"Set output format."
					);
		}
		
		/**
		 * Invoke this option
		 * 
		 * @param option Option.
		 * @param argument Option argument.
		 */
		@Override
		public boolean invoke(String option, String argument) {
			
			if (argument.isEmpty()) {
				error("Cannot set output format (" + option + "): Argument is empty", KestrelConstants.ERR_USAGE);
				return false;
			}
			
			runnerBase.setOutputFormat(argument);
			
			return true;
		}
		
		/**
		 * Initialize this option.
		 */
		@Override
		public void init() {
			runnerBase.setOutputFormat(KestrelRunnerBase.DEFAULT_OUTPUT_FORMAT);
		}
	}
	
	/**
	 * Option: Output file name
	 */
	protected class OptLogFileName extends OptionSpecElement {
		
		/**
		 * Create option.
		 */
		public OptLogFileName() {
			super('\0', "logfile",
					OptionArgumentType.REQUIRED,
					"LOG_FILE", KestrelRunnerBase.DEFAULT_LOG_FILE.name,
					"Set log file name."
					);
		}
		
		/**
		 * Invoke this option
		 * 
		 * @param option Option.
		 * @param argument Option argument.
		 */
		@Override
		public boolean invoke(String option, String argument) {
			
			argument = argument.trim();
			
			if (argument.isEmpty()) {
				error("Cannot set log file name (" + option + "): Argument is empty", KestrelConstants.ERR_USAGE);
				return false;
			}
			
			try {
				runnerBase.setLogFile(argument);
				
			} catch (IllegalArgumentException ex) {
				error("Cannot set log file name (" + option + "): " + ex.getMessage(), KestrelConstants.ERR_SYSTEM, ex);
				return false;
			}
			
			return true;
		}
		
		/**
		 * Initialize this option.
		 */
		@Override
		public void init() {
			runnerBase.setLogFile(KestrelRunnerBase.DEFAULT_LOG_FILE);
		}
	}
	
	/**
	 * Option: Write log messages to STDOUT
	 */
	protected class OptLogStdout extends OptionSpecElement {
		
		/**
		 * Create option.
		 */
		public OptLogStdout() {
			super('\0', "logstdout",
					OptionArgumentType.NONE,
					null, null,
					"Write log messages to standard output instead of a file. Unless redirected, this output " +
					"is written to the the screen."
					);
		}
		
		/**
		 * Invoke this option
		 * 
		 * @param option Option.
		 * @param argument Option argument.
		 */
		@Override
		public boolean invoke(String option, String argument) {
			
			runnerBase.setLogFile(FileDescriptor.out);
			
			return true;
		}
		
		// Init by OptLogFileName
	}
	
	/**
	 * Option: Write log messages to STDERR
	 */
	protected class OptLogStderr extends OptionSpecElement {
		
		/**
		 * Create option.
		 */
		public OptLogStderr() {
			super('\0', "logstderr",
					OptionArgumentType.NONE,
					null, null,
					"Write log messages to standard error instead of a file. Unless redirected, this output " +
					"is written to the the screen."
					);
		}
		
		/**
		 * Invoke this option
		 * 
		 * @param option Option.
		 * @param argument Option argument.
		 */
		@Override
		public boolean invoke(String option, String argument) {
			
			runnerBase.setLogFile(FileDescriptor.err);
			
			return true;
		}
		
		// Init by OptLogFileName
	}
	
	/**
	 * Option: Set log level
	 */
	protected class OptLogLevel extends OptionSpecElement {
		
		/**
		 * Create option.
		 */
		public OptLogLevel() {
			super('\0', "loglevel",
					OptionArgumentType.REQUIRED,
					"LOG_LEVEL", KestrelRunnerBase.DEFAULT_LOG_LEVEL.name(),
					"Set the log level. Valid levels are " + LogLevel.levelList(true) + "."
					);
		}
		
		/**
		 * Invoke this option
		 * 
		 * @param option Option.
		 * @param argument Option argument.
		 */
		@Override
		public boolean invoke(String option, String argument) {
			
			runnerBase.setLogLevel(LogLevel.getLevel(argument));
			
			return true;
		}
		
		/**
		 * Initialize this option.
		 */
		@Override
		public void init() {
			runnerBase.setLogLevel(KestrelRunnerBase.DEFAULT_LOG_LEVEL);
		}
	}
	
	/**
	 * Option: Load library file.
	 */
	protected class OptLoadLibraryFile extends OptionSpecElement {
		
		/**
		 * Create option.
		 */
		public OptLoadLibraryFile() {
			super('\0', "lib",
					OptionArgumentType.REQUIRED,
					"LIB_FILE", null,
					"Load a library file. Kestrel can accept external components, and they must be " +
					"packaged on a JAR file."
					);
		}
		
		/**
		 * Invoke this option
		 * 
		 * @param option Option.
		 * @param argument Option argument.
		 */
		@Override
		public boolean invoke(String option, String argument) {
			
			argument = argument.trim();
			
			if (argument.isEmpty()) {
				error("Cannot load library file/directory (" + option + "): File name is empty", KestrelConstants.ERR_USAGE);
				return false;
			}
			
			File libFile = new File(argument);
			
			try {
				runnerBase.addLibraryFile(libFile);
				
			} catch (MalformedURLException ex) {
				error("Cannot load library file/directory (" + option + "): Bad URL: " + argument + ": " + ex.getMessage(), KestrelConstants.ERR_USAGE);
				return false;
				
			} catch (SecurityException ex) {
				error("Security error loading library file/directory (" + option + "): " + argument + ": " + ex.getMessage(), KestrelConstants.ERR_SECURITY);
				return false;
				
			} catch (FileNotFoundException ex) {
				error("File not found while loading library file/directory (" + option + "): " + argument + ": " + ex.getMessage(), KestrelConstants.ERR_FILENOTFOUND);
				return false;
				
			} catch (IOException ex) {
				error("IO error while loading library file/directory (" + option + "): " + argument + ": " + ex.getMessage(), KestrelConstants.ERR_IO);
				return false;
			}
			
			return true;
		}
		
		/**
		 * Initialize this option.
		 */
		@Override
		public void init() {
			runnerBase.clearLibraries();
		}
	}
	
	/**
	 * Option: Load library by URL.
	 */
	protected class OptLoadLibraryUrl extends OptionSpecElement {
		
		/**
		 * Create option.
		 */
		public OptLoadLibraryUrl() {
			super('\0', "liburl",
					OptionArgumentType.REQUIRED,
					"LIB_URL", null,
					"Load a library by its URL. Kestrel can accept external components, and they must " +
					"be packaged on a JAR file. This option can access JAR files on the local system or " +
					"stored anywhere the program can access and that can be represented as a URL."
					);
		}
		
		/**
		 * Invoke this option
		 * 
		 * @param option Option.
		 * @param argument Option argument.
		 */
		@Override
		public boolean invoke(String option, String argument) {
			
			argument = argument.trim();
			
			if (argument.isEmpty()) {
				error("Cannot load library URL (" + option + "): URL is empty", KestrelConstants.ERR_USAGE);
				return false;
			}
			
			try {
				runnerBase.addLibraryURL(new URL(argument));
				
			} catch (MalformedURLException ex) {
				error("Cannot load library file/directory (" + option + "): Bad URL: " + argument + ": " + ex.getMessage(), KestrelConstants.ERR_USAGE);
				return false;
			}
			
			return true;
		}
		
		// init(): Performed by OptLoadLibraryFile
	}
	
	/**
	 * Option: Number of files to automatically group into input samples.
	 */
	public class OptFilesPerSample extends OptionSpecElement {
		
		/**
		 * Create option.
		 */
		public OptFilesPerSample() {
			super('\0', "filespersample",
					OptionArgumentType.REQUIRED,
					"N_FILES", "0",
					"Set the number of input files per sample. For example, reading paired-end FASTQ files " +
					"(2 files per sample) can be simplified by setting this value to 2. Alternatively, samples " +
					"can be separated by multiple -s (--sample) arguments. The default value, 0, will not automatically " +
					"group input files. If -s (--sample) is read on the command-line, this value is set back to 0. Any " +
					"sequence files found on the command-line before this option are assigned to a group."
					);
		}
		
		/**
		 * Invoke this option
		 * 
		 * @param option Option.
		 * @param argument Option argument.
		 */
		@Override
		public boolean invoke(String option, String argument) {
			
			int filesPerSample;
			
			argument = argument.trim();
			
			if (argument.isEmpty()) {
				error("Cannot set files-per-sample (" + option + "): Argument is empty", KestrelConstants.ERR_USAGE);
				return false;
			}
				
			
			try {
				filesPerSample = Integer.parseInt(argument);
				
			} catch (NumberFormatException ex) {
				error("Cannot set files-per-sample (" + option + "): Argument is not an integer: " + argument, KestrelConstants.ERR_USAGE);
				return false;
			}
			
			if (filesPerSample < 0) {
				error("Cannot set files-per-sample (" + option + "): Argument is a negative integer: " + filesPerSample, KestrelConstants.ERR_USAGE);
				return false;
			}
			
			// Check for an existing automatic grouping
			if (KestrelArgumentParser.this.filesPerSample > 0 && ! inputSourceQueue.isEmpty()) {
				
				warn(
						String.format(
							"Reassigning automatic grouping (%s) from %d to %d leaves %d input files assigned to the last group instead of %d files",
							option,
							KestrelArgumentParser.this.filesPerSample,
							filesPerSample,
							inputSourceQueue.size(),
							KestrelArgumentParser.this.filesPerSample
						),
						KestrelConstants.ERR_USAGE
				);
			}
			
			KestrelArgumentParser.this.filesPerSample = filesPerSample;
			
			// Clear queued samples
			flushSample();
			
			return true;
		}
		
		/**
		 * Initialize this option.
		 */
		@Override
		public void init() {
			filesPerSample = 0;
			
			return;
		}
	}
	
	/**
	 * Option: Sample set
	 */
	public class OptSampleSet extends OptionSpecElement {
		
		/**
		 * Create option.
		 */
		public OptSampleSet() {
			super('s', "sample",
					OptionArgumentType.OPTIONAL,
					"SAMPLE_NAME", null,
					"Set the name of the sample that the next sample files are assigned to. If the argument " +
					"(SAMPLE_NAME) is given, the name of the sample is set to this name. If the argument is " +
					"not given, then the sample name is assigned from the name of the first file after this " +
					"option. Any files on the command-line appearing before this option are assigned to a " +
					"sample and will not be part of this sample. If --filespersample was used on the command-" +
					"line before this option, it is reset and files are no longer automatically grouped."
					);
		}
		
		/**
		 * Invoke this option
		 * 
		 * @param option Option.
		 * @param argument Option argument.
		 */
		@Override
		public boolean invoke(String option, String argument) {
			
			// Null empty argument string
			if (argument != null) {
				argument = argument.trim();
				
				if (argument.isEmpty())
					argument = null;
			}
			
			// Check for a set that is not full
			if (filesPerSample > 0 && inputSourceQueue.size() < filesPerSample) {
				warn(
					String.format(
						"Setting group (%s) when automatic grouping was set to %d leaves %d input files assigned to the last group instead of %d files",
						option,
						KestrelArgumentParser.this.filesPerSample,
						inputSourceQueue.size(),
						KestrelArgumentParser.this.filesPerSample
					),
					KestrelConstants.ERR_USAGE
				);
			}
			
			// Clear queued samples
			flushSample();
			
			// Assign sample name
			sampleName = argument;
			filesPerSample = 0;
			
			return true;
		}
		
		/**
		 * Initialize this option.
		 */
		@Override
		public void init() {
			sampleName = null;
			
			return;
		}
	}
	
	/**
	 * Option: Add a reference sequence
	 */
	protected class OptReference extends OptionSpecElement {
		
		/**
		 * Create option.
		 */
		public OptReference() {
			super('r', "ref",
					OptionArgumentType.REQUIRED,
					"REF_SEQUENCE", null,
					"Add reference sequences variants will be called against. This can be any file that Kestrel " +
					"can read. The format and character-set options apply to reference sequences, but not filters.");
		}
		
		/**
		 * Invoke this option
		 * 
		 * @param option Option.
		 * @param argument Option argument.
		 */
		@Override
		public boolean invoke(String option, String argument) {
			
			// Check file name
			if (argument.isEmpty()) {
				error("Cannot add reference sequence file (" + option + "): File name is empty", KestrelConstants.ERR_USAGE);
				return false;
			}
			
			// Add reference
			runnerBase.addReference(new FileSequenceSource(new File(argument), formatType, charset, ++refSourceIndex, ""));
			
			return true;
		}
		
		/**
		 * Initialize this option.
		 */
		@Override
		public void init() {
			refSourceIndex = 0;
			runnerBase.clearReference();
		}
	}
	
	/**
	 * Option: Alignment weights as a vector
	 */
	protected class OptAlignmentWeightVector extends OptionSpecElement {
		
		/**
		 * Create option.
		 */
		public OptAlignmentWeightVector() {
			super('w', "weight",
					OptionArgumentType.REQUIRED,
					"WEIGHT_VEC", AlignmentWeight.get().toString(1),
					"Set the alignment weights as a comma-separated list of values. The order of weights " +
					"is match, mismatch, gap-open, gap-extend, and initial score. If values are blank or the list has " +
					"fewer than 5 elements, the missing values are assigned their default weight. Each " +
					"value is a floating-point number, and it may be represented in exponential form (e.g. 1.0e2) " +
					"or as an integer in hexadecimal or octal format. Optionally, the list may be surrounded by " +
					"parenthesis or braces (angle, square, or curly)."
			);
		}
		
		/**
		 * Invoke this option
		 * 
		 * @param option Option.
		 * @param argument Option argument.
		 */
		@Override
		public boolean invoke(String option, String argument) {
			
			try {
				runnerBase.setAlignmentWeight(argument);
				
			} catch (IllegalArgumentException ex) {
				error("Invalid vector of alignment weights (" + option + "): " + ex.getMessage(), KestrelConstants.ERR_USAGE, ex);
				return false;
			}
			
			return true;
		}
		
		/**
		 * Initialize this option.
		 */
		@Override
		public void init() {
			runnerBase.setAlignmentWeight(null);
			
			return;
		}
	}
	
	/**
	 * Option: Set character-set encoding of files
	 */
	protected class OptCharset extends OptionSpecElement {
		
		/**
		 * Create option.
		 */
		public OptCharset() {
			super ('\0', "charset",
					OptionArgumentType.REQUIRED,
					"CHARSET", KestrelRunnerBase.DEFAULT_CHARSET.toString(),
					"Character set encoding of input files. This option specifies the character set " +
					"of all files following it. The default, \"UTF-8\", properly handles ASCII files, " +
					"which is a safe assumption for most files. Latin-1 files with values greater than " +
					"127 will not be properly parsed."
					);
		}
		
		/**
		 * Invoke this option
		 * 
		 * @param option Option.
		 * @param argument Option argument.
		 */
		@Override
		public boolean invoke(String option, String argument) {
			
			if (argument.isEmpty()) {
				error("Cannot set input file character-set (" + option + "): Character-set name is empty", KestrelConstants.ERR_USAGE);
				return false;
			}
			
			try {
				charset = Charset.forName(argument);
				
			} catch (IllegalCharsetNameException ex) {
				error("Cannot set input file character-set (" + option + "): Character-set name is illegal: " + argument, KestrelConstants.ERR_USAGE);
				return false;
				
			} catch (UnsupportedCharsetException ex) {
				error("Cannot set input file character-set (" + option + "): Character-set name is unsupported: " + argument, KestrelConstants.ERR_USAGE);
				return false;
			}
			
			return true;
		}
		
		/**
		 * Initialize this option.
		 */
		@Override
		public void init() {
			charset = KestrelRunnerBase.DEFAULT_CHARSET;
		}
	}
	
	/**
	 * Base class for sequence filter and quality options (aliases of each other).
	 */
	protected class OptSequenceFilterBase extends OptionSpecElement {
		
		/**
		 * Create a sequence filter option.
		 * 
		 * @param opt Option character.
		 * @param longOpt Long option string.
		 * @param argType Type of argument.
		 * @param argName Name of argument.
		 * @param defaultValue Default value.
		 * @param helpText Help text.
		 */
		public OptSequenceFilterBase(char opt, String longOpt, OptionArgumentType argType, String argName, String defaultValue, String helpText) {
			super(opt, longOpt, argType, argName, defaultValue, helpText);
		}
		
		/**
		 * Invoke this option
		 * 
		 * @param option Option.
		 * @param argument Option argument.
		 */
		@Override
		public boolean invoke(String option, String argument) {
			
			if (argument.isEmpty()) {
				filterSpec = "";
				return true;
			}
			
			String[] tok = argument.split("\\s*:\\s*", 2);
			
			if (tok.length == 1) {
				filterSpec = "sanger:" + tok[0];
				
			} else {
				
				if (tok[0].isEmpty())
					tok[0] = "sanger";
				
				filterSpec = tok[0] + ":" + tok[1];
			}
			
			return true;
		}
		
		/**
		 * Initialize this option.
		 */
		@Override
		public void init() {
			filterSpec = "";
		}
	}
	
	/**
	 * Option: Sequence filter.
	 */
	protected class OptSequenceFilter extends OptSequenceFilterBase {
		
		/**
		 * Create option.
		 */
		public OptSequenceFilter() {
			super('\0', "seqfilter",
					OptionArgumentType.REQUIRED,
					"SEQ_FILTER", null,
					"Filter sequences as they are read and before k-mers are extracted from them. " +
					"Some sequence readers can filter or alter reads at runtime. The most common " +
					"filter is a quality filter where low-quality bases are removed. The filter " +
					"specification is a filter name followed by a colon (:) and arguments to the " +
					"filter. If a filter name is not specified, then the \"sanger\" quality filter " +
					"is assumed. For example, \"sanger:10\" and \"10\" will filter k-mers with any " +
					"base quality score less than 10. The sequence filter specification is set for " +
					"all files appearing on the command-line after this option. To turn off filtering " +
					"once it has been set, files following --noseqfilter will have no filter specification."
					);
		}
	}
	
	/**
	 * Option: Quality score filtering. Alias for <code>OptSequenceFilter</code>.
	 */
	protected class OptQuality extends OptSequenceFilterBase {
		
		/**
		 * Create option.
		 */
		public OptQuality() {
			super('\0', "quality",
					OptionArgumentType.REQUIRED,
					"SEQ_FILTER", null,
					"This option is an alias for \"seqfilter\""
					);
		}
	}
	
	/**
	 * Option: No sequence (or quality) filtering.
	 */
	protected class OptNoSequenceFilter extends OptionSpecElement {
		
		/**
		 * Create option.
		 */
		public OptNoSequenceFilter() {
			super('\0', "noseqfilter",
					OptionArgumentType.NONE,
					null, "",
					"Turn off sequence filtering for all files following this option. If --seqfilter or " +
					"--quality was specified, this option disables sequence filtering. These options together " +
					"make it possible to specify filtering for some files and disable filtering for others."
					);
		}
		
		/**
		 * Invoke this option
		 * 
		 * @param option Option.
		 * @param argument Option argument.
		 */
		@Override
		public boolean invoke(String option, String argument) {
			
			filterSpec = "";
			
			return true;
		}
	}
	
	/**
	 * Option: Temp file location
	 */
	protected class OptTempFileLocation extends OptionSpecElement {
		
		/**
		 * Create option.
		 */
		public OptTempFileLocation() {
			super('\0', "temploc",
					OptionArgumentType.REQUIRED,
					"TEMP", "Output location",
					"The location where segments are offloaded. This argument must be a " +
					"directory or the location for a new directory. Parent directories will " +
					"be created as needed."
					);
		}
		
		/**
		 * Invoke this option
		 * 
		 * @param option Option.
		 * @param argument Option argument.
		 */
		@Override
		public boolean invoke(String option, String argument) {
			
			if (argument.isEmpty()) {
				error("Cannot set temp directory (" + option + "): Directory name is empty", KestrelConstants.ERR_USAGE);
				return false;
			}
			
			runnerBase.setTempDirName(argument);
			
			return true;
		}
	}
	
	/**
	 * Option: Minimum k-mer count.
	 */
	protected class OptMinKmerCount extends OptionSpecElement {
		
		/**
		 * Create option.
		 */
		public OptMinKmerCount() {
			super('\0', "mincount",
					OptionArgumentType.REQUIRED,
					"MIN_COUNT", Integer.toString(KestrelRunnerBase.DEFAULT_MIN_KMER_COUNT),
					"Set the minimum k-mer count for processing samples. K-mers with a count less than " +
					"this value will be discarded. Sequence read errors produce many erroneous k-mers, " +
					"and this slows the process of variant calling significantly."
					);
		}
		
		/**
		 * Invoke this option
		 * 
		 * @param option Option.
		 * @param argument Option argument.
		 */
		@Override
		public boolean invoke(String option, String argument) {
			
			// Check argument
			if (argument.isEmpty()) {
				error("Cannot set minimum k-mer count (" + option + "): Count is empty", KestrelConstants.ERR_USAGE);
				return false;
			}
			
			try {
				runnerBase.setMinKmerCount(Integer.parseInt(argument));
				
			} catch (NumberFormatException ex) {
				error("Cannot set minimum k-mer count (" + option + "): Count is not an integer: " + argument, KestrelConstants.ERR_USAGE);
				return false;
				
			} catch (IllegalArgumentException ex) {
				error("Cannot set minimum k-mer count (" + option + "): " + ex.getMessage(), KestrelConstants.ERR_USAGE, ex);
				return false;
			}
			
			return true;
		}
		
		/**
		 * Initialize this option.
		 */
		@Override
		public void init() {
			runnerBase.setMinKmerCount(KestrelRunnerBase.DEFAULT_MIN_KMER_COUNT);
		}
	}
	
	/**
	 * Option: Minimum k-mer count difference for detecting active regions.
	 */
	protected class OptMinKmerCountDiff extends OptionSpecElement {
		
		/**
		 * Create option.
		 */
		public OptMinKmerCountDiff() {
			super('\0', "mindiff",
					OptionArgumentType.REQUIRED,
					"COUNT_DIFF", Integer.toString(ActiveRegionDetector.DEFAULT_MINIMUM_DIFFERENCE),
					"Set the minimum k-mer count difference for identifying active regions. When the " +
					"count between neighboring k-mer counts is this or greater, Kestrel will treat it as " +
					"a region where a variant may occur."
					);
		}
		
		/**
		 * Invoke this option
		 * 
		 * @param option Option.
		 * @param argument Option argument.
		 */
		@Override
		public boolean invoke(String option, String argument) {
			
			// Check argument
			if (argument.isEmpty()) {
				error("Cannot set minimum k-mer count difference (" + option + "): Count is empty", KestrelConstants.ERR_USAGE);
				return false;
			}
			
			try {
				runnerBase.setMinimumDifference(Integer.parseInt(argument));
				
			} catch (NumberFormatException ex) {
				error("Cannot set minimum k-mer count difference (" + option + "): Count is not an integer: " + argument, KestrelConstants.ERR_USAGE);
				return false;
				
			} catch (IllegalArgumentException ex) {
				error("Cannot set minimum k-mer count difference (" + option + "): " + ex.getMessage(), KestrelConstants.ERR_USAGE, ex);
				return false;
			}
			
			return true;
		}
		
		/**
		 * Initialize this option.
		 */
		@Override
		public void init() {
			runnerBase.setMinimumDifference(ActiveRegionDetector.DEFAULT_MINIMUM_DIFFERENCE);
		}
	}
	
	/**
	 * Option: K-mer count difference quantile.
	 */
	protected class OptKmerCountDiffQuantile extends OptionSpecElement {
		
		/**
		 * Create option.
		 */
		public OptKmerCountDiffQuantile() {
			super('\0', "diffq",
					OptionArgumentType.REQUIRED,
					"QUANTILE", String.format("%.04f", ActiveRegionDetector.DEFAULT_DIFFERENCE_QUANTILE),
					"If set to a value greater than 0.0, then the k-mer count difference between two k-mers " +
					"that triggers a correction attempt is found dynamically. The difference in k-mer counts " +
					"between each pair of neighboring k-mers over an uncorrected reference region is found, and " +
					"this quantile of is computed over those differences. For example, a value of 0.85 means " +
					"that at most 15% (100% - 85%) of the k-mer count differences will be high enough. If this " +
					"computed value is less than the minimum k-mer count difference (--mindiff), then that " +
					"minimum is the difference threshold. This value may not be 1.0 or greater, and it may not be " +
					"negative. If 0.0, the minimum count difference is always the minimum threshold (--mindiff)."
					);
		}
		
		/**
		 * Invoke this option
		 * 
		 * @param option Option.
		 * @param argument Option argument.
		 */
		@Override
		public boolean invoke(String option, String argument) {
			
			// Check argument
			if (argument.isEmpty()) {
				error("Cannot set k-mer count quantile (" + option + "): Quantile is empty", KestrelConstants.ERR_USAGE);
				return false;
			}
			
			try {
				runnerBase.setDifferenceQuantile(Double.parseDouble(argument));
				
			} catch (NumberFormatException ex) {
				error("Cannot set k-mer count quantile (" + option + "): Quantile is not a number: " + argument, KestrelConstants.ERR_USAGE);
				return false;
				
			} catch (IllegalArgumentException ex) {
				error("Cannot set k-mer count quantile (" + option + "): " + ex.getMessage(), KestrelConstants.ERR_USAGE, ex);
				return false;
			}
			
			return true;
		}
		
		/**
		 * Initialize this option.
		 */
		@Override
		public void init() {
			runnerBase.setDifferenceQuantile(ActiveRegionDetector.DEFAULT_DIFFERENCE_QUANTILE);
		}
	}
	
	/**
	 * Option: Peak scan length
	 */
	protected class OptPeakScanLength extends OptionSpecElement {
		
		/**
		 * Create option.
		 */
		public OptPeakScanLength() {
			super('\0', "peakscan",
					OptionArgumentType.REQUIRED,
					"LENGTH", "" + ActiveRegionDetector.DEFAULT_PEAK_SCAN_LENGTH,
					"Reference regions with sequence homology in other regions of the genome may " +
					"contain k-mers with artificially high frequencies from adding counts for k-mers that " +
					"appear in both regions. This causes a peak in the k-mer frequencies over the reference, " +
					"and it can trigger an erroneous active-region scan for variants. When encountering " +
					"a difference, Kestrel will scan forward this number of k-mers looking for a peak in " +
					"the k-mer frequencies. If the frequencies drop back down to the original range, the " +
					"active-region scan is not performed. This keeps Kestrel from erroneously searching large " +
					"regions of the reference. Setting this value to 0 disables peak detection."
					);
		}
		
		/**
		 * Invoke this option
		 * 
		 * @param option Option.
		 * @param argument Option argument.
		 */
		@Override
		public boolean invoke(String option, String argument) {
			
			// Check argument
			if (argument.isEmpty()) {
				error("Cannot set peak scan length (" + option + "): Length is empty", KestrelConstants.ERR_USAGE);
				return false;
			}
			
			try {
				runnerBase.setPeakScanLength(Integer.parseInt(argument));
				
			} catch (NumberFormatException ex) {
				error("Cannot set peak scan length (" + option + "): Length is not an integer: " + argument, KestrelConstants.ERR_USAGE);
				return false;
				
			} catch (IllegalArgumentException ex) {
				error("Cannot set peak scan length (" + option + "): " + ex.getMessage(), KestrelConstants.ERR_USAGE, ex);
				return false;
			}
			
			return true;
		}
		
		/**
		 * Initialize this option.
		 */
		@Override
		public void init() {
			runnerBase.setPeakScanLength(ActiveRegionDetector.DEFAULT_PEAK_SCAN_LENGTH);
		}
	}
	
	/**
	 * Option: Scan limit factor
	 */
	protected class OptScanLimitFactor extends OptionSpecElement {
		
		/**
		 * Create option.
		 */
		public OptScanLimitFactor() {
			super('\0', "scanlimitfactor",
					OptionArgumentType.REQUIRED,
					"FACTOR", "" + String.format("%.01f", ActiveRegionDetector.DEFAULT_SCAN_LIMIT_FACTOR),
					"Set a limit on how long an active region may be. This is computed " +
					"by multiplying the k-mer size by this factor and adding the maximum length of a gap. The " +
					"computed limit will be adjusted so that active regions are at least large enough to capture " +
					"a SNP in cases where the maximum gap length is 0. Setting this to a low value or \"min\" will " +
					"set the limit so that it is just large enough to catch SNPs and deletions, but it will miss large " +
					"deletions if another variant is within the k-mer size window. Setting this to a high value or " +
					"\"max\" lifts the restrictions on active region lengths, and this may cause the program to take an " +
					"excessive amount of time and memory trying to solve arbitrarily long active regions."
					);
		}
		
		/**
		 * Invoke this option
		 * 
		 * @param option Option.
		 * @param argument Option argument.
		 */
		@Override
		public boolean invoke(String option, String argument) {
			
			// Check argument
			if (argument.isEmpty()) {
				error("Cannot set end scan limit factor (" + option + "): Limit is empty", KestrelConstants.ERR_USAGE);
				return false;
			}
			
			if (argument.equalsIgnoreCase("min")) {
				runnerBase.setScanLimitFactor(0.0);
				
			} if (argument.equalsIgnoreCase("max")) {
				runnerBase.setScanLimitFactor(Integer.MAX_VALUE);
				
			} else {
				
				try {
					runnerBase.setScanLimitFactor(Double.parseDouble(argument));
					
				} catch (NumberFormatException ex) {
					error("Cannot set end scan limit factor (" + option + "): Limit is not a number, \"min\", or \"max\": " + argument, KestrelConstants.ERR_USAGE);
					return false;
					
				} catch (IllegalArgumentException ex) {
					error("Cannot set end scan limit factor (" + option + "): " + ex.getMessage(), KestrelConstants.ERR_USAGE, ex);
					return false;
				}
			}
			
			return true;
		}
		
		/**
		 * Initialize this option.
		 */
		@Override
		public void init() {
			runnerBase.setScanLimitFactor(ActiveRegionDetector.DEFAULT_SCAN_LIMIT_FACTOR);
		}
	}
	
	/**
	 * Option: Set decay min
	 */
	protected class OptDecayMinimum extends OptionSpecElement {
		
		/**
		 * Create option.
		 */
		public OptDecayMinimum() {
			super('\0', "decaymin",
					OptionArgumentType.REQUIRED,
					"MIN_PROPORTION", "" + String.format("%.02f", ActiveRegionDetector.DEFAULT_EXP_MIN),
					"Set the minimum value (asymptotic lower bound) of the exponential decay function " +
					"used in active region detection as a proportion of the anchor k-mer count. If this " +
					"value is 0.0, k-mer count recovery threshold may decline to 1. If this value is " +
					"1.0, the decay function is not used and the detector falls back to finding a k-mer " +
					"with a count within the difference threshold of the anchor k-mer count."
					);
		}
		
		/**
		 * Invoke this option
		 * 
		 * @param option Option.
		 * @param argument Option argument.
		 */
		@Override
		public boolean invoke(String option, String argument) {
			
			// Check argument
			if (argument.isEmpty()) {
				error("Cannot set exponential decay minimum (" + option + "): Minimum is empty", KestrelConstants.ERR_USAGE);
				return false;
			}
			
			try {
				runnerBase.setDecayMinimum(Double.parseDouble(argument));
				
			} catch (NumberFormatException ex) {
				error("Cannot set exponential decay minimum (" + option + "): Minimum is not a number: " + argument, KestrelConstants.ERR_USAGE);
				return false;
				
			} catch (IllegalArgumentException ex) {
				error("Cannot set exponential decay minimum (" + option + "): " + ex.getMessage(), KestrelConstants.ERR_USAGE, ex);
				return false;
			}
			
			return true;
		}
		
		/**
		 * Initialize this option.
		 */
		@Override
		public void init() {
			runnerBase.setDecayMinimum(ActiveRegionDetector.DEFAULT_EXP_MIN);
		}
	}
	
	/**
	 * Option: Set decay alpha
	 */
	protected class OptDecayAlpha extends OptionSpecElement {
		
		/**
		 * Create option.
		 */
		public OptDecayAlpha() {
			super('\0', "alpha",
					OptionArgumentType.REQUIRED,
					"ALPHA", "" + String.format("%.02f", ActiveRegionDetector.DEFAULT_EXP_ALPHA),
					"Set the exponential decay alpha, which controls how quickly the recovery threshold " +
					"declines to its minimum value (see --decaymin) in an active region search. Alpha is " +
					"defined as the rate of decay for every k bases. At k bases from the left anchor, the " +
					"threshold will have declined to alpha * range. At every k bases, the threshold will " +
					"continue to decline at this rate."
					);
		}
		
		/**
		 * Invoke this option
		 * 
		 * @param option Option.
		 * @param argument Option argument.
		 */
		@Override
		public boolean invoke(String option, String argument) {
			
			// Check argument
			if (argument.isEmpty()) {
				error("Cannot set exponential decay alpha (" + option + "): Alpha is empty", KestrelConstants.ERR_USAGE);
				return false;
			}
			
			try {
				runnerBase.setDecayAlpha(Double.parseDouble(argument));
				
			} catch (NumberFormatException ex) {
				error("Cannot set exponential decay alpha (" + option + "): Alpha is not a number: " + argument, KestrelConstants.ERR_USAGE);
				return false;
				
			} catch (IllegalArgumentException ex) {
				error("Cannot set exponential decay alpha (" + option + "): " + ex.getMessage(), KestrelConstants.ERR_USAGE, ex);
				return false;
			}
			
			return true;
		}
		
		/**
		 * Initialize this option.
		 */
		@Override
		public void init() {
			runnerBase.setDecayAlpha(ActiveRegionDetector.DEFAULT_EXP_ALPHA);
		}
	}
	
	/**
	 * Option: Haplotypes must be anchored to unaltered k-mers on both ends.
	 */
	protected class OptAnchorBothEnds extends OptionSpecElement {
		
		/**
		 * Create option.
		 */
		public OptAnchorBothEnds() {
			super('\0', "anchorboth",
					OptionArgumentType.NONE,
					null, (ActiveRegionDetector.DEFAULT_ANCHOR_BOTH_ENDS ? "" : null),
					"Both ends of an active region (region with variants) must be bordered by unaltered " +
					"k-mers or variants will not be called in it. This option may miss variants near the ends " +
					"of a reference sequence, but it forces stronger evidence for the variants that are called."
					);
		}
		
		/**
		 * Invoke this option
		 * 
		 * @param option Option.
		 * @param argument Option argument.
		 */
		@Override
		public boolean invoke(String option, String argument) {
			runnerBase.setAnchorBothEnds(true);
			
			return true;
		}
		
		/**
		 * Initialize this option.
		 */
		@Override
		public void init() {
			runnerBase.setAnchorBothEnds(ActiveRegionDetector.DEFAULT_ANCHOR_BOTH_ENDS);
		}
	}
	
	/**
	 * Option: Haplotypes must be anchored to on at least one end, and possibly both.
	 */
	protected class OptNoAnchorBothEnds extends OptionSpecElement {
		
		/**
		 * Create option.
		 */
		public OptNoAnchorBothEnds() {
			super('\0', "noanchorboth",
					OptionArgumentType.NONE,
					null, (ActiveRegionDetector.DEFAULT_ANCHOR_BOTH_ENDS ? null : ""),
					"An active region (region with variants) must be bordered on at least one side by an unaltered " +
					"k-mers, but it may extend to the end of the sequence. This will allow Kestrel to find variants " +
					"less than a k-mer from the ends, but the evidence supporting these variants is weaker."
					);
		}
		
		/**
		 * Invoke this option
		 * 
		 * @param option Option.
		 * @param argument Option argument.
		 */
		@Override
		public boolean invoke(String option, String argument) {
			runnerBase.setAnchorBothEnds(false);
			
			return true;
		}
		
		// Init by OptAnchorBothEnds
	}
	
	/**
	 * Option: Allow active regions with ambiguous bases.
	 */
	protected class OptCallAmbiguousRegions extends OptionSpecElement {
		
		/**
		 * Create option.
		 */
		public OptCallAmbiguousRegions() {
			super('\0', "ambiregions",
					OptionArgumentType.NONE,
					null, (ActiveRegionDetector.DEFAULT_CALL_AMBIGUOUS_REGIONS ? "" : null),
					"Allow active regions to cover ambiguous bases, such as N."
					);
		}
		
		/**
		 * Invoke this option
		 * 
		 * @param option Option.
		 * @param argument Option argument.
		 */
		@Override
		public boolean invoke(String option, String argument) {
			runnerBase.setCallAmbiguousRegions(true);
			
			return true;
		}
		
		/**
		 * Initialize this option.
		 */
		@Override
		public void init() {
			runnerBase.setCallAmbiguousRegions(ActiveRegionDetector.DEFAULT_CALL_AMBIGUOUS_REGIONS);
		}
	}
	
	/**
	 * Option: Do not allow active regions with ambiguous bases.
	 */
	protected class OptNoCallAmbiguousRegions extends OptionSpecElement {
		
		/**
		 * Create option.
		 */
		public OptNoCallAmbiguousRegions() {
			super('\0', "noambigregions",
					OptionArgumentType.NONE,
					null, (ActiveRegionDetector.DEFAULT_CALL_AMBIGUOUS_REGIONS ? null : ""),
					"An active region may not span any base that is not A, C, G, T, or U."
					);
		}
		
		/**
		 * Invoke this option
		 * 
		 * @param option Option.
		 * @param argument Option argument.
		 */
		@Override
		public boolean invoke(String option, String argument) {
			runnerBase.setCallAmbiguousRegions(false);
			
			return true;
		}
		
		// Init by OptCallAmbiguousRegions
	}
	
	/**
	 * Option: Allow variants over ambiguous bases.
	 */
	protected class OptCallAmbiguousVariant extends OptionSpecElement {
		
		/**
		 * Create option.
		 */
		public OptCallAmbiguousVariant() {
			super('\0', "ambivar",
					OptionArgumentType.NONE,
					null, (VariantCaller.DEFAULT_CALL_AMBIGUOUS_VARIANT ? "" : null),
					"Allow variants over ambiguous bases, such as N."
					);
		}
		
		/**
		 * Invoke this option
		 * 
		 * @param option Option.
		 * @param argument Option argument.
		 */
		@Override
		public boolean invoke(String option, String argument) {
			runnerBase.setCallAmbiguousVariant(true);
			
			return true;
		}
		
		/**
		 * Initialize this option.
		 */
		@Override
		public void init() {
			runnerBase.setCallAmbiguousVariant(VariantCaller.DEFAULT_CALL_AMBIGUOUS_VARIANT);
		}
	}
	
	/**
	 * Option: Do not allow variants over ambiguous bases.
	 */
	protected class OptNoCallAmbiguousVariant extends OptionSpecElement {
		
		/**
		 * Create option.
		 */
		public OptNoCallAmbiguousVariant() {
			super('\0', "noambivar",
					OptionArgumentType.NONE,
					null, (VariantCaller.DEFAULT_CALL_AMBIGUOUS_VARIANT ? null : ""),
					"A variant may not span any base that is not A, C, G, T, or U."
					);
		}
		
		/**
		 * Invoke this option
		 * 
		 * @param option Option.
		 * @param argument Option argument.
		 */
		@Override
		public boolean invoke(String option, String argument) {
			runnerBase.setCallAmbiguousVariant(false);
			
			return true;
		}
		
		// Init by OptCallAmbiguousVariant
	}
	
	/**
	 * Option: Size of k-mer minimizers.
	 */
	protected class OptMinimizerSize extends OptionSpecElement {
		
		/**
		 * Create option.
		 */
		public OptMinimizerSize() {
			super ('\0', "minsize",
					OptionArgumentType.REQUIRED,
					"MIN_SIZE", "" + KestrelRunnerBase.DEFAULT_MINIMIZER_SIZE,
					"Minimizers group k-mers in the indexed k-mer count (IKC) file generated by " +
					"Kestrel when reading sequences, and this parameter controls the size of the " +
					"minimizer."
					);
		}
		
		/**
		 * Invoke this option
		 * 
		 * @param option Option.
		 * @param argument Option argument.
		 */
		@Override
		public boolean invoke(String option, String argument) {
			
			int val;
			
			// Check argument
			if (argument.isEmpty()) {
				error("Cannot set minimizer size (" + option + "): Arugment is empty", KAnalyzeConstants.ERR_USAGE);
			}
			
			// Get integer
			try {
				val = Integer.parseInt(argument);
						
			} catch (NumberFormatException ex) {
				error("Cannot set minimizer size (" + option + "): Arugment is not an integer: " + argument, KAnalyzeConstants.ERR_USAGE);
				return false;
			}
			
			// Set
			try {
				runnerBase.setMinimizerSize(val);
				
			} catch (IllegalArgumentException ex) {
				error("Cannot set minimizer size (" + option + "): " + ex.getMessage(), KAnalyzeConstants.ERR_USAGE);
			}
			
			return true;
		}
		
		/**
		 * Initialize this option.
		 */
		@Override
		public void init() {
			runnerBase.setMinimizerSize(KestrelRunnerBase.DEFAULT_MINIMIZER_SIZE);
		}
	}
	
	/**
	 * Option: Minimizer mask.
	 */
	protected class OptMinimizerMask extends OptionSpecElement {
		
		/**
		 * Create option.
		 */
		public OptMinimizerMask() {
			super ('\0', "minmask",
					OptionArgumentType.REQUIRED,
					"MIN_MASK", "" + KestrelRunnerBase.DEFAULT_MINIMIZER_MASK,
					"Size of k-mer minimizers or 0 to disable processing by minimizers. The minimizer of a k-mer is " +
					"determined by taking all sub-k-mers of a given size (set by this option) from a k-mer and its " +
					"reverse complement and choosing the lesser of the sub-k-mers. Sub-k-mers are XORed with this " +
					"mask while comparing them, but the minimizer is not XORed (it is still a sub-k-mer of the original " +
					"k-mer). This option can be used to break up large minimizer groups due to low-complexity k-mers when " +
					"minimizers are used."
					);
		}
		
		/**
		 * Invoke this option
		 * 
		 * @param option Option.
		 * @param argument Option argument.
		 */
		@Override
		public boolean invoke(String option, String argument) {
			
			int val;
			
			// Check argument
			if (argument.isEmpty()) {
				error("Cannot set minimizer size (" + option + "): Arugment is empty", KAnalyzeConstants.ERR_USAGE);
			}
			
			// Get integer
			try {
				val = StringUtil.toInt(argument, false);
						
			} catch (NumberFormatException ex) {
				error("Cannot set minimizer size (" + option + "): Arugment is not an integer: " + argument + ": " + ex.getMessage(), KAnalyzeConstants.ERR_USAGE);
				return false;
			}
			
			// Set
			try {
				runnerBase.setMinimizerMask(val);
				
			} catch (IllegalArgumentException ex) {
				error("Cannot set minimizer size (" + option + "): " + ex.getMessage(), KAnalyzeConstants.ERR_USAGE);
			}
			
			return true;
		}
		
		/**
		 * Initialize this option.
		 */
		@Override
		public void init() {
			runnerBase.setMinimizerMask(KestrelRunnerBase.DEFAULT_MINIMIZER_MASK);
		}
	}
	
	/**
	 * Option: Haplotypes must be anchored to unaltered k-mers on both ends.
	 */
	protected class OptKmerCountInMemory extends OptionSpecElement {
		
		/**
		 * Create option.
		 */
		public OptKmerCountInMemory() {
			super('\0', "memcount",
					OptionArgumentType.NONE,
					null, (KestrelRunnerBase.DEFAULT_KMER_COUNT_IN_MEMORY ? "" : null),
					"K-mer counts from each sample will be stored in memory. This option assumes that " +
					"samples are relatively small or the machine has enough memory to handle the counts. " +
					"Note that the JVM might need to be run with additional memory (-Xmx option) to support " +
					"this option."
					);
		}
		
		/**
		 * Invoke this option
		 * 
		 * @param option Option.
		 * @param argument Option argument.
		 */
		@Override
		public boolean invoke(String option, String argument) {
			runnerBase.setKmerCountInMemory(true);
			
			return true;
		}
		
		/**
		 * Initialize this option.
		 */
		@Override
		public void init() {
			runnerBase.setKmerCountInMemory(KestrelRunnerBase.DEFAULT_KMER_COUNT_IN_MEMORY);
		}
	}
	
	/**
	 * Option: Haplotypes must be anchored to on at least one end, and possibly both.
	 */
	protected class OptNoKmerCountInMemory extends OptionSpecElement {
		
		/**
		 * Create option.
		 */
		public OptNoKmerCountInMemory() {
			super('\0', "nomemcount",
					OptionArgumentType.NONE,
					null, (KestrelRunnerBase.DEFAULT_KMER_COUNT_IN_MEMORY ? null : ""),
					"K-mer counts for each sample are offloaded to an indexed k-mer count file. This " +
					"option reduces the memory demand of Kestrel."
					);
		}
		
		/**
		 * Invoke this option
		 * 
		 * @param option Option.
		 * @param argument Option argument.
		 */
		@Override
		public boolean invoke(String option, String argument) {
			runnerBase.setKmerCountInMemory(false);
			
			return true;
		}
		
		// Init by OptAnchorBothEnds
	}
	
	/**
	 * Option: Free resources as soon as possible.
	 */
	protected class OptFreeResources extends OptionSpecElement {
		
		/**
		 * Create option.
		 */
		public OptFreeResources() {
			super('\0', "free",
					OptionArgumentType.NONE,
					null, (KestrelRunnerBase.DEFAULT_FREE_RESOURCES ? "" : null),
					"Free resources between processing samples. This may reduce the memory footprint of " +
					"Kestrel, but it may force expensive resources to be recreated and impact performance."
					);
		}
		
		/**
		 * Invoke this option
		 * 
		 * @param option Option.
		 * @param argument Option argument.
		 */
		@Override
		public boolean invoke(String option, String argument) {
			runnerBase.setFreeResources(true);
			
			return true;
		}
		
		/**
		 * Initialize this option.
		 */
		@Override
		public void init() {
			runnerBase.setFreeResources(KestrelRunnerBase.DEFAULT_FREE_RESOURCES);
		}
	}
	
	/**
	 * Option: Retain resources and avoid re-creating them.
	 */
	protected class OptNoFreeResources extends OptionSpecElement {
		
		/**
		 * Create option.
		 */
		public OptNoFreeResources() {
			super('\0', "nofree",
					OptionArgumentType.NONE,
					null, (KestrelRunnerBase.DEFAULT_FREE_RESOURCES ? null : ""),
					"Retain resources between samples. This may use more memory, but it will avoid " +
					"re-creating expensive resources between samples."
					);
		}
		
		/**
		 * Invoke this option
		 * 
		 * @param option Option.
		 * @param argument Option argument.
		 */
		@Override
		public boolean invoke(String option, String argument) {
			runnerBase.setFreeResources(false);
			
			return true;
		}
		
		// Init by OptFreeResources
	}
	
	/**
	 * Option: Remove generated indexed k-mer count (IKC) files.
	 */
	protected class OptRmIkc extends OptionSpecElement {
		
		/**
		 * Create option.
		 */
		public OptRmIkc() {
			super('\0', "rmikc",
					OptionArgumentType.NONE,
					null, (KestrelRunnerBase.DEFAULT_REMOVE_IKC ? "" : null),
					"Remove the indexed k-mer count (IKC) for each sample after kestrel runs."
					);
		}
		
		/**
		 * Invoke this option
		 * 
		 * @param option Option.
		 * @param argument Option argument.
		 */
		@Override
		public boolean invoke(String option, String argument) {
			runnerBase.setRemoveIkc(true);
			
			return true;
		}
		
		/**
		 * Initialize this option.
		 */
		@Override
		public void init() {
			runnerBase.setRemoveIkc(KestrelRunnerBase.DEFAULT_REMOVE_IKC);
		}
	}
	
	/**
	 * Option: Do not delete generated indexed k-mer count (IKC) files.
	 */
	protected class OptNoRmIkc extends OptionSpecElement {
		
		/**
		 * Create option.
		 */
		public OptNoRmIkc() {
			super('\0', "normikc",
					OptionArgumentType.NONE,
					null, (KestrelRunnerBase.DEFAULT_REMOVE_IKC ? null : ""),
					"Do not remove the indexed k-mer count (IKC) file for each sample."
					);
		}
		
		/**
		 * Invoke this option
		 * 
		 * @param option Option.
		 * @param argument Option argument.
		 */
		@Override
		public boolean invoke(String option, String argument) {
			runnerBase.setRemoveIkc(false);
			
			return true;
		}
		
		// Init by OptRmIkc
	}
	
	/**
	 * Option: Count reverse complement k-mers in region statistics.
	 */
	protected class OptCountReverseKmers extends OptionSpecElement {
		
		/**
		 * Create option.
		 */
		public OptCountReverseKmers() {
			super('\0', "countrev",
					OptionArgumentType.NONE,
					null, (KestrelRunnerBase.DEFAULT_COUNT_REV_KMER ? "" : null),
					"Count reverse complement k-mers in region statistics. This should be set " +
					"for most sequencing protocols."
					);
		}
		
		/**
		 * Invoke this option
		 * 
		 * @param option Option.
		 * @param argument Option argument.
		 */
		@Override
		public boolean invoke(String option, String argument) {
			runnerBase.setCountReverseKmers(true);
			
			return true;
		}
		
		/**
		 * Initialize this option.
		 */
		@Override
		public void init() {
			runnerBase.setCountReverseKmers(KestrelRunnerBase.DEFAULT_COUNT_REV_KMER);
		}
	}
	
	/**
	 * Option: Do not delete generated indexed k-mer count (IKC) files.
	 */
	protected class OptNoCountReverseKmers extends OptionSpecElement {
		
		/**
		 * Create option.
		 */
		public OptNoCountReverseKmers() {
			super('\0', "nocountrev",
					OptionArgumentType.NONE,
					null, (KestrelRunnerBase.DEFAULT_COUNT_REV_KMER ? null : ""),
					"Do not include the reverse complement of k-mers in read depth estimates. " +
					"If all sequence reads are in the same orientation as the reference, then " +
					"this option should be used."
					);
		}
		
		/**
		 * Invoke this option
		 * 
		 * @param option Option.
		 * @param argument Option argument.
		 */
		@Override
		public boolean invoke(String option, String argument) {
			runnerBase.setCountReverseKmers(false);
			
			return true;
		}
		
		// Init by OptCountReverseKmers
	}
	
	/**
	 * Option: Add a variant filter.
	 */
	protected class OptAddVariantFilter extends OptionSpecElement {
		
		/**
		 * Create option.
		 */
		public OptAddVariantFilter() {
			super('\0', "varfilter",
					OptionArgumentType.REQUIRED,
					"FILTER_SPEC", null,
					"Add a variant filter specification. The argument should be the name of the " +
					"filter, a colon, and the filter arguments. The correct filter is loaded by " +
					"name and the filter arguments are passed to it."
			);
		}
		
		/**
		 * Invoke this option
		 * 
		 * @param option Option.
		 * @param argument Option argument.
		 */
		@Override
		public boolean invoke(String option, String argument) {
			
			try {
				runnerBase.addVariantFilter(argument);
				
			} catch (Exception ex) {
				error("Error setting filter (" + option + "): " + ex.getMessage(), KAnalyzeConstants.ERR_USAGE);
				return false;
			}
			
			return true;
		}
		
		/**
		 * Initialize this option.
		 */
		@Override
		public void init() {
			runnerBase.clearVariantFilters();
		}
	}
	
	/**
	 * Option: Add intervals
	 */
	protected class OptReadIntervalFile extends OptionSpecElement {
		
		/**
		 * Create option.
		 */
		public OptReadIntervalFile() {
			super('i', "interval",
					OptionArgumentType.REQUIRED,
					"INTERVAL_FILE", null,
					"Reads a file of intervals defining the regions over the reference sequences where " +
					"variants should be detected. If no intervals are specified, variants are detected " +
					"over the full length of each reference sequence. The file type is determined by the " +
					"file name, such as \"intervals.bed\". BED files are supported by Kestrel, and others " +
					"may be added."
			);
		}
		
		/**
		 * Invoke this option
		 * 
		 * @param option Option.
		 * @param argument Option argument.
		 */
		@Override
		public boolean invoke(String option, String argument) {
			
			File intervalFile = new File(argument);
			
			try {
				if (runnerBase.readIntervals(intervalFile, null, null) == 0) {
					error("Error reading interval file (" + option + "): No intervals found in file: " + intervalFile.getPath(), KAnalyzeConstants.ERR_USAGE);
					return false;
				}
								
			} catch (FileNotFoundException ex) {
				error("Error reading interval file (" + option + "): " + ex.getMessage(), KAnalyzeConstants.ERR_USAGE);
				return false;
				
			} catch (IOException ex) {
				error("Error reading interval file (" + option + "): " + ex.getMessage(), KAnalyzeConstants.ERR_USAGE);
				return false;
				
			} catch (IntervalReaderInitException ex) {
				error("Error reading interval file (" + option + "): " + ex.getMessage(), KAnalyzeConstants.ERR_USAGE);
				return false;
			}
			
			return true;
		}
		
		/**
		 * Initialize this option.
		 */
		@Override
		public void init() {
			runnerBase.clearIntervals();
		}
	}
	
	/**
	 * Option: Set interval flank length.
	 */
	protected class OptSetFlankLength extends OptionSpecElement {
		
		/**
		 * Create option.
		 */
		public OptSetFlankLength() {
			super('\0', "flank",
					OptionArgumentType.REQUIRED,
					"LENGTH", String.format("k * %.2f", KestrelRunnerBase.DEFAULT_FLANK_LENGTH_MULTIPLIER),
					"When extracting intervals from reference sequences, this many bases are extracted " +
					"on both sides of the interval whenever possible. This gives Kestrel more bases for " +
					"active region detection, but it does not otherwise affect variant calls. Set to 0 " +
					"to disable flanks. If this option is not set, Kestrel will determine the appropriate " +
					"length of flank by multiplying " +
					String.format("%.2f", KestrelRunnerBase.DEFAULT_FLANK_LENGTH_MULTIPLIER) +
					" with the k-mer size."
			);
		}
		
		/**
		 * Invoke this option
		 * 
		 * @param option Option.
		 * @param argument Option argument.
		 */
		@Override
		public boolean invoke(String option, String argument) {
			
			int flankLength;
			
			try {
				flankLength = Integer.parseInt(argument);
				
				if (flankLength < 0) {
					error("Error setting flank length (" + option + "): Length is negative: " + flankLength, KAnalyzeConstants.ERR_USAGE);
					return false;
				}
				
				runnerBase.setFlankLength(flankLength);
				
			} catch (NumberFormatException ex) {
				error("Error setting flank length (" + option + "): Length is not a number: " + argument, KAnalyzeConstants.ERR_USAGE);
				return false;
				
			} catch (IllegalArgumentException ex) {
				error("Error setting flank length (" + option + "): " + ex.getMessage(), KAnalyzeConstants.ERR_USAGE);
				return false;
			}
			
			return true;
		}
		
		/**
		 * Initialize this option.
		 */
		@Override
		public void init() {
			runnerBase.setDefaultFlankLength();
		}
	}
	
	/**
	 * Option: Set interval flank length.
	 */
	protected class OptAutoFlankLength extends OptionSpecElement {
		
		/**
		 * Create option.
		 */
		public OptAutoFlankLength() {
			super('\0', "autoflank",
					OptionArgumentType.NONE,
					null, "",
					"When extracting intervals from reference sequences, some bases are extracted " +
					"on both sides of the interval whenever possible. This gives Kestrel more bases for " +
					"active region detection, but it does not otherwise affect variant calls. " +
					"This option tells Kestrel to pick the flank by multiplying " +
					String.format("%.2f", KestrelRunnerBase.DEFAULT_FLANK_LENGTH_MULTIPLIER) +
					" with the k-mer size."
			);
		}
		
		/**
		 * Invoke this option
		 * 
		 * @param option Option.
		 * @param argument Option argument.
		 */
		@Override
		public boolean invoke(String option, String argument) {
			
			runnerBase.setDefaultFlankLength();
			
			return true;
		}
		
		// Init by OptSetFLankLength
	}
	
	/**
	 * Option: Variant calls are relative to the reference sequence
	 */
	protected class OptVarCallRelativeReference extends OptionSpecElement {
		
		/**
		 * Create option.
		 */
		public OptVarCallRelativeReference() {
			super('\0', "byreference",
					OptionArgumentType.NONE,
					null, "",
					"If variant call regions were defined, variant call locations are still relative to the " +
					"reference sequence and not the region."
					);
		}
		
		/**
		 * Invoke this option
		 * 
		 * @param option Option.
		 * @param argument Option argument.
		 */
		@Override
		public boolean invoke(String option, String argument) {
			runnerBase.setVariantCallByReference();
			
			return true;
		}
		
		/**
		 * Initialize this option.
		 */
		@Override
		public void init() {
			runnerBase.setVariantCallByReference();
		}
	}
	
	/**
	 * Option: Variant calls are relative to the reference region
	 */
	protected class OptVarCallRelativeRegion extends OptionSpecElement {
		
		/**
		 * Create option.
		 */
		public OptVarCallRelativeRegion() {
			super('\0', "byregion",
					OptionArgumentType.NONE,
					null, null,
					"If variant call regions were defined, variant call locations are relative to the " +
					"region and not the reference sequence."
					);
		}
		
		/**
		 * Invoke this option
		 * 
		 * @param option Option.
		 * @param argument Option argument.
		 */
		@Override
		public boolean invoke(String option, String argument) {
			runnerBase.setVariantCallByRegion();
			
			return true;
		}
		
		// Init by OptVarCallRelativeReference
	}
	
	/**
	 * Option: Set maximum alignment states
	 */
	protected class OptMaxAlignStates extends OptionSpecElement {
		
		/**
		 * Create option.
		 */
		public OptMaxAlignStates() {
			super('\0', "maxalignstates",
					OptionArgumentType.REQUIRED,
					"STATES", "" + KmerAligner.DEFAULT_MAX_STATE,
					"Set the maximum number of alignment states. When haplotype assembly " +
					"branches into more than one possible sequence, the state of one is saved " +
					"while another is built. When the maximum number of saved states reaches " +
					"this value, the least likely one is discarded."
			);
		}
		
		/**
		 * Invoke this option
		 * 
		 * @param option Option.
		 * @param argument Option argument.
		 */
		@Override
		public boolean invoke(String option, String argument) {
			
			int maxState;
			
			try {
				maxState = Integer.parseInt(argument);
				
				if (maxState < 0) {
					error("Error setting maximum alignment states (" + option + "): Number of states is negative: " + maxState, KAnalyzeConstants.ERR_USAGE);
					return false;
				}
				
				runnerBase.setMaxAlignerState(maxState);
				
			} catch (NumberFormatException ex) {
				error("Error setting maximum alignment states (" + option + "): Argument is not a number: " + argument, KAnalyzeConstants.ERR_USAGE);
				return false;
				
			} catch (IllegalArgumentException ex) {
				error("Error setting maximum alignment states (" + option + "): " + ex.getMessage(), KAnalyzeConstants.ERR_USAGE);
				return false;
			}
			
			return true;
		}
		
		/**
		 * Initialize this option.
		 */
		@Override
		public void init() {
			runnerBase.setMaxAlignerState(KmerAligner.DEFAULT_MAX_STATE);
		}
	}
	
	/**
	 * Option: Set maximum alignment states
	 */
	protected class OptMaxHaplotypeStates extends OptionSpecElement {
		
		/**
		 * Create option.
		 */
		public OptMaxHaplotypeStates() {
			super('\0', "maxhapstates",
					OptionArgumentType.REQUIRED,
					"STATES", "" + KmerAlignmentBuilder.DEFAULT_MAX_HAPLOTYPES,
					"Set the maximum number of haplotypes for an active region. Alignments can generate " +
					"more than one haplotype, and with noisy sequence data or paralogues, many haplotypes " +
					"may be found. This options limits the amount of memory that can be consumed in these cases."
			);
		}
		
		/**
		 * Invoke this option
		 * 
		 * @param option Option.
		 * @param argument Option argument.
		 */
		@Override
		public boolean invoke(String option, String argument) {
			
			int maxHaplotypes;
			
			try {
				maxHaplotypes = Integer.parseInt(argument);
				
				if (maxHaplotypes < 0) {
					error("Error setting maximum number of haplotypes (" + option + "): Number is negative: " + maxHaplotypes, KAnalyzeConstants.ERR_USAGE);
					return false;
				}
				
				runnerBase.setMaxHaplotypes(maxHaplotypes);
				
			} catch (NumberFormatException ex) {
				error("Error setting maximum haplotypes (" + option + "): Argument is not a number: " + argument, KAnalyzeConstants.ERR_USAGE);
				return false;
				
			} catch (IllegalArgumentException ex) {
				error("Error setting maximum haplotypes (" + option + "): " + ex.getMessage(), KAnalyzeConstants.ERR_USAGE);
				return false;
			}
			
			return true;
		}
		
		/**
		 * Initialize this option.
		 */
		@Override
		public void init() {
			runnerBase.setMaxHaplotypes(KmerAlignmentBuilder.DEFAULT_MAX_HAPLOTYPES);
		}
	}
	
	/**
	 * Option: Set maximum alignment states
	 */
	protected class OptMaxRepeatCount extends OptionSpecElement {
		
		/**
		 * Create option.
		 */
		public OptMaxRepeatCount() {
			super('\0', "maxrepeat",
					OptionArgumentType.REQUIRED,
					"COUNT", "" + ActiveRegionDetector.DEFAULT_MAX_REPEAT_COUNT,
					"Cycles in the k-mer graph produce unreliable local assemblies. The default value for " +
					"this option (0) will terminate any local assembly that contains the same k-mer more " +
					"than once. For most applications, 0 is the recommended value. To attempt assemblies in " +
					"reptitive regions, this value can be increased, but the results may be variant calls on " +
					"haplotypes that do not exist in the sequence data."
			);
		}
		
		/**
		 * Invoke this option
		 * 
		 * @param option Option.
		 * @param argument Option argument.
		 */
		@Override
		public boolean invoke(String option, String argument) {
			
			int maxRepeatCount;
			
			try {
				maxRepeatCount = StringUtil.toInt(argument, false);
				
				if (maxRepeatCount < 0) {
					error("Error setting maximum repeat count (" + option + "): Number is negative: " + maxRepeatCount, KAnalyzeConstants.ERR_USAGE);
					return false;
				}
				
				runnerBase.setMaxRepeatCount(maxRepeatCount);
				
			} catch (NumberFormatException ex) {
				error("Error setting maximum repeat count (" + option + "): " + ex.getMessage(), KAnalyzeConstants.ERR_USAGE);
				return false;
			}
			
			return true;
		}
		
		/**
		 * Initialize this option.
		 */
		@Override
		public void init() {
			runnerBase.setMaxRepeatCount(ActiveRegionDetector.DEFAULT_MAX_REPEAT_COUNT);
		}
	}
	
	/**
	 * Option: Remove reference sequence descriptions
	 */
	protected class OptRemoveRefDescription extends OptionSpecElement {
		
		/**
		 * Create option.
		 */
		public OptRemoveRefDescription() {
			super('\0', "rmrefdesc",
					OptionArgumentType.NONE,
					null, (ReferenceReader.DEFAULT_REMOVE_SEQUENCE_DESCRIPTION ? "" : null),
					"When set, remove the description from reference sequence names. The descirption " +
					"is everything that occurs after the first whitespace character. FASTA files often " +
					"have a sequence name and a long description separated by whitespace. This " +
					"option ensures that the sequence name matches in the FASTA and an interval file, if used."
					);
		}
		
		/**
		 * Invoke this option
		 * 
		 * @param option Option.
		 * @param argument Option argument.
		 */
		@Override
		public boolean invoke(String option, String argument) {
			runnerBase.setRemoveReferencDescription(true);
			
			return true;
		}
		
		/**
		 * Initialize this option.
		 */
		@Override
		public void init() {
			runnerBase.setRemoveReferencDescription(ReferenceReader.DEFAULT_REMOVE_SEQUENCE_DESCRIPTION);
		}
	}
	
	/**
	 * Option: Do not remove reference sequence descriptions
	 */
	protected class OptNoRemoveRefDescription extends OptionSpecElement {
		
		/**
		 * Create option.
		 */
		public OptNoRemoveRefDescription() {
			super('\0', "normrefdesc",
					OptionArgumentType.NONE,
					null, (ReferenceReader.DEFAULT_REMOVE_SEQUENCE_DESCRIPTION ? null : ""),
					"Use the full sequence name as it appears in the reference sequence file. FASTA files " +
					"often include a description after the sequence name, and with this option, it becomes " +
					"part of the full sequence name. If using an interval file, the full sequence name and " +
					"description must match the sequence file."
					);
		}
		
		/**
		 * Invoke this option
		 * 
		 * @param option Option.
		 * @param argument Option argument.
		 */
		@Override
		public boolean invoke(String option, String argument) {
			runnerBase.setRemoveReferencDescription(false);
			
			return true;
		}
		
		// Init by OptRemoveRefDescription
	}
	
	
	
	
	/**
	 * Option: Reverse complement negative strand reference regions
	 */
	protected class OptRevComplNegRegStrand extends OptionSpecElement {
		
		/**
		 * Create option.
		 */
		public OptRevComplNegRegStrand() {
			super('\0', "revregion",
					OptionArgumentType.NONE,
					null, (KestrelRunnerBase.DEFAULT_REVERSE_COMPLEMENT_NEGATIVE_STRAND ? "" : null),
					"When set, reverse complement reference regions that occur on the negative strand. " +
					"Only itervals defined with on the negative strand are altered."
					);
		}
		
		/**
		 * Invoke this option
		 * 
		 * @param option Option.
		 * @param argument Option argument.
		 */
		@Override
		public boolean invoke(String option, String argument) {
			runnerBase.setRevComplementNegReferenceStrand(true);
			
			return true;
		}
		
		/**
		 * Initialize this option.
		 */
		@Override
		public void init() {
			runnerBase.setRevComplementNegReferenceStrand(KestrelRunnerBase.DEFAULT_REVERSE_COMPLEMENT_NEGATIVE_STRAND);
		}
	}
	
	/**
	 * Option: Do not reverse complement negative strand reference regions
	 */
	protected class OptNoRevComplNegRegStrand extends OptionSpecElement {
		
		/**
		 * Create option.
		 */
		public OptNoRevComplNegRegStrand() {
			super('\0', "norevregion",
					OptionArgumentType.NONE,
					null, (KestrelRunnerBase.DEFAULT_REVERSE_COMPLEMENT_NEGATIVE_STRAND ? "" : null),
					"When set, regions variants are called on are always in the same orientation as the reference " +
					"sequence. The stranded-ness of defined intervals is ignored."
					);
		}
		
		/**
		 * Invoke this option
		 * 
		 * @param option Option.
		 * @param argument Option argument.
		 */
		@Override
		public boolean invoke(String option, String argument) {
			runnerBase.setRemoveReferencDescription(false);
			
			return true;
		}
		
		// Init by OptRemoveRefDescription
	}
	
	/**
	 * Option: Haplotype output file name
	 */
	protected class OptHaplotypeOutputFileName extends OptionSpecElement {
		
		/**
		 * Create option.
		 */
		public OptHaplotypeOutputFileName() {
			super('p', "hapout",
					OptionArgumentType.REQUIRED,
					"OUT_FILE", null,
					"Set haplotype output file name."
					);
		}
		
		/**
		 * Invoke this option
		 * 
		 * @param option Option.
		 * @param argument Option argument.
		 */
		@Override
		public boolean invoke(String option, String argument) {
			
			argument = argument.trim();
			
			if (argument.isEmpty()) {
				error("Cannot set haplotype output file name (" + option + "): Argument is empty", KestrelConstants.ERR_USAGE);
				return false;
			}
			
			try {
				runnerBase.setHaplotypeOutputFile(argument);
				
			} catch (IllegalArgumentException ex) {
				error("Cannot set haplotype output file name (" + option + "): " + ex.getMessage(), KestrelConstants.ERR_SYSTEM, ex);
				return false;
			}
			
			return true;
		}
		
		/**
		 * Initialize this option.
		 */
		@Override
		public void init() {
			runnerBase.setHaplotypeOutputFile(null);
		}
	}

	/**
	 * Option: Haplotype output format
	 */
	protected class OptHaplotypeOutputFormat extends OptionSpecElement {
		
		/**
		 * Create option.
		 */
		public OptHaplotypeOutputFormat() {
			super('\0', "hapfmt",
					OptionArgumentType.REQUIRED,
					"OUT_FORMAT", KestrelRunnerBase.DEFAULT_HAPLOTYPE_OUTPUT_FORMAT,
					"Set haplotype output format. Ignored if a haplotype output file is not set."
					);
		}
		
		/**
		 * Invoke this option
		 * 
		 * @param option Option.
		 * @param argument Option argument.
		 */
		@Override
		public boolean invoke(String option, String argument) {
			
			if (argument.isEmpty()) {
				error("Cannot set haplotype output format (" + option + "): Argument is empty", KestrelConstants.ERR_USAGE);
				return false;
			}
			
			runnerBase.setHaplotypeOutputFormat(argument);
			
			return true;
		}
		
		/**
		 * Initialize this option.
		 */
		@Override
		public void init() {
			runnerBase.setHaplotypeOutputFormat(KestrelRunnerBase.DEFAULT_HAPLOTYPE_OUTPUT_FORMAT);
		}
	}
	
	//
	// Non-option argument (input files)
	//
	
	/**
	 * Option: Input files.
	 */
	protected class InputSourceFile extends NonOptionElement {
		
		/**
		 * Create option.
		 */
		public InputSourceFile() {
			super("SEQ_FILE", true, true,
					"Input sequence file.");
		}

		/**
		 * Invoke this option
		 * 
		 * @param argument Option argument.
		 */
		@Override
		public boolean invoke(String argument) {
			
			argument = argument.trim();
			
			if (argument.isEmpty()) {
				error("Cannot process empty file name argument from the command line", KestrelConstants.ERR_USAGE);
				return false;
			}
			
			// Automatically group files if set
			if (filesPerSample > 0 && inputSourceQueue.size() >= filesPerSample)
				flushSample();
			
			inputSourceQueue.add(new FileSequenceSource(new File(argument), formatType, charset, ++seqSourceIndex, filterSpec));
			
			return true;
		}
		
		/**
		 * Initialize this option.
		 */
		@Override
		public void init() {
			runnerBase.clearSamples();
			inputSourceQueue.clear();
			seqSourceIndex = 0;
		}
	}
	
	//////////////////////////////////
	//////////////////////////////////
	////  Help Topic Definitions  ////
	//////////////////////////////////
	//////////////////////////////////
	
	/**
	 * Help topic: Input file formats
	 */
	protected class HelpFormat extends HelpTopicElement {
		
		/**
		 * Create help topic.
		 */
		public HelpFormat() {
			super("options", "Command-line option conventions.",
					"The command-line option system is very robust and flexible. This topic outlines " +
					"the command-line conventions and how the help output (-h) can be interpreted.\n" +
					"\n" +
					"Options may have a short form, a long form, or both. Short options start with a " +
					"single dash (-) and are a single character (e.g. -v). Long options start with " +
					"two dashes (--) and are more than one character (e.g. --verbose). Where an " +
					"option has a short and long form, both are listed in the full help, but only the " +
					"short option is listed in the summary at the top of the help output.\n" +
					"\n" +
					"Options may have a required argument, an optional argument, and some do not " +
					"accept an argument. For a required argument, a short description is " +
					"listed after the option (e.g. -f <INPUT_FORMAT>). For an optional argument, the " +
					"short description is surrounded with square brackets (e.g. -r[<RMODE>]).\n" +
					"\n" +
					"On the command-line, arguments may be supplied in two ways. For short options, " +
					"the argument can appear after a space or appended to the option. For example, " +
					"\"-f fastq\" and \"-ffastq\" have the same meaning. Note that \"-ffastq\" cannot be a " +
					"long option because there is a single dash, so it must be a short option with an " +
					"argument. For long options, the argument may appear after a space or after an " +
					"equal sign (=). For example, \"--format fastq\" and \"--format=fastq\" have the same " +
					"meaning.\n" +
					"\n" +
					"For optional arguments, the argument cannot appear after a space because it " +
					"becomes ambiguous. For example, \"-rduplicate\" and \"--reverse=duplicate\" are the " +
					"only valid forms.\n" +
					"\n" +
					"Any arguments on the command-line that are not attached to an option are input " +
					"files containing sequences to be processed.\n" +
					"\n" +
					"Arguments are processed in the order they are found. The file format option,  " +
					"\"-f\", can set different formats for input files. For example, \"-f fastq f1.fq -f  " +
					"fastqgz f2.fq.gz\" inputs two FASTQ files where one is GZIP compressed and the  " +
					"other is not. Other options override each other, and the option listed last  " +
					"determines the behavior of the program. Where the KAnalyze command is setup as  " +
					"an alias with custom defaults, the defaults can be changed when the prgram is  " +
					"run. For example, if an alias runs KAnalyze in GUI mode with \"-g\", it can be  " +
					"disabled at runtime with \"-G\".\n" +
					"\n" +
					"For features with an inverse, the option that enables a feature is in " +
					"lower-case, and the option that disables the feature is in upper-case. The long " +
					"argument form of the option that disables the feature is prefixed with \"no\". For " +
					"example, \"-g\" and \"--gui\" turns on GUI mode while \"-G\" and \"--nogui\" disables " +
					"it.\n"
			);
		}
	}
	
	/**
	 * Help topic: Input file formats
	 */
	protected class HelpReader extends HelpTopicElement {
		
		/**
		 * Create help topic.
		 */
		public HelpReader() {
			super("reader", "Sequence reader options.",
					"Readers:"
			);
		}
		
		/**
		 * List readers that can be loaded at runtime.
		 */
		@Override
		public String dynamicContent() {
			StringBuilder builder = new StringBuilder();
			String description;
			URLClassLoader loader = runnerBase.getLoader();
			
			for (String name : SequenceReader.listReaders(loader)) {
				
				description = SequenceReader.getFormatFileDescription(name, loader);
				
				if (description == null)
					description = "";
				
				if (description.length() > 50)
					description = description.substring(0, 47) + "...";
				
				builder.append(String.format("  * %-15s %-50s%s", name, description, StringUtil.LINE_SEP));
			}
			
			return builder.toString();
		}
	}
	
	/**
	 * Help topic: Input file formats
	 */
	protected class HelpRefReader extends HelpTopicElement {
		
		/**
		 * Create help topic.
		 */
		public HelpRefReader() {
			super("refreader", "Reference sequence reader options.",
					"Readers:"
			);
		}
		
		/**
		 * List readers that can be loaded at runtime.
		 */
		@Override
		public String dynamicContent() {
			StringBuilder builder = new StringBuilder();
			String description;
			URLClassLoader loader = runnerBase.getLoader();
			
			for (String name : SequenceReader.listReaders(loader)) {
				
				description = SequenceReader.getFormatFileDescription(name, loader);
				
				if (description == null)
					description = "";
				
				if (description.length() > 50)
					description = description.substring(0, 47) + "...";
				
				builder.append(String.format("  * %-15s %-50s%s", name, description, StringUtil.LINE_SEP));
			}
			
			return builder.toString();
		}
	}
	
	/**
	 * Help topic: Cite Kestrel
	 */
	protected class HelpCite extends HelpTopicElement {
		
		/**
		 * Create help topic.
		 */
		public HelpCite() {
			super("cite", "Cite Kestrel",
					"Audano, P. A., Ravishankar, S., & Vannberg, F. O. (2017). Mapping-free variant calling using haplotype reconstruction from k-mer frequencies. Bioinformatics, (April), 1–7. https://doi.org/10.1093/bioinformatics/btx753"
			);
			
			return;
		}
	}
	
	/**
	 * Help topic: Cite Kestrel
	 */
	protected class HelpBib extends HelpTopicElement {
		
		/**
		 * Create help topic.
		 */
		public HelpBib() {
			super("citebib", "Cite Kestrel (BibTeX)",
					"@article{Audano2018,\n" +
					"abstract = {{\\textcopyright} The Author(s) 2017. Motivation The standard protocol for detecting variation in DNA is to map millions of short sequence reads to a known reference and find loci that differ. While this approach works well, it cannot be applied where the sample contains dense variants or is too distant from known references. De novo assembly or hybrid methods can recover genomic variation, but the cost of computation is often much higher. We developed a novel k-mer algorithm and software implementation, Kestrel, capable of characterizing densely packed SNPs and large indels without mapping, assembly or de Bruijn graphs. Results When applied to mosaic penicillin binding protein (PBP) genes in Streptococcus pneumoniae, we found near perfect concordance with assembled contigs at a fraction of the CPU time. Multilocus sequence typing (MLST) with this approach was able to bypass de novo assemblies. Kestrel has a very low false-positive rate when applied to the whole genome, and while Kestrel identified many variants missed by other methods, limitations of a purely k-mer based approach affect overall sensitivity. Availability and implementation Source code and documentation for a Java implementation of Kestrel can be found at https://github.com/paudano/kestrel. All test code for this publication is located at https://github.com/paudano/kescases.},\n" +
					"author = {Audano, P.A. and Ravishankar, S. and Vannberg, F.O.},\n" +
					"doi = {10.1093/bioinformatics/btx753},\n" +
					"issn = {14602059},\n" +
					"journal = {Bioinformatics},\n" +
					"number = {10},\n" +
					"title = {{Mapping-free variant calling using haplotype reconstruction from k-mer frequencies}},\n" +
					"volume = {34},\n" +
					"year = {2018}\n" +
					"}"
			);
			
			return;
		}
	}
	
	
	//
	// Condition listener and methods for handling errors
	//
	
	/**
	 * Intercepts configuration errors and stores them in a list to be handled after configuration
	 * returns.
	 */
	private class ConfigurationConditionListener implements ConditionListener {
		
		/**
		 * Intercept a condition and store it.
		 */
		@Override
		public void conditionOccurred(ConditionEvent condition) {
			
			if (condition == null)
				return;
			
			if (condition.type == ConditionEvent.TYPE_ERROR)
				errorList.add(condition);
			
			else if (condition.type == ConditionEvent.TYPE_WARNING)
				warningList.add(condition);
			
			return;
		}
	}
	
	/**
	 * If a configuration error has occurred, throw and exception here. This is how
	 * <code>KestrelRunnerBase.configure()</code> knows to throw an exception after processing
	 * command-line arguments. Calling this method does not reset the error condition. 
	 *  
	 * @throws ConfigurationErrorException If a configuration error has occurred in this parser, then
	 *   this exception is thrown, and it contains the details of the error. If multiple errors have
	 *   occurred, only the first one is thrown.
	 */
	public void throwConfigurationError()
		throws ConfigurationErrorException {
			
		if (errorList.size() > 0) {
			ConditionEvent condition = errorList.get(0);
			
			throw new ConfigurationErrorException(condition.message, condition.type, condition.cause);
		}
	}
	
	/**
	 * Get a list of all errors that have occurred. Calling this method does not reset the error condition.
	 * 
	 * @return A list of all errors that have occurred, or an empty list if there are no errors to report.
	 */
	public ConditionEvent[] getErrorConditions() {
		return errorList.toArray(new ConditionEvent[0]);
	}
	
	/**
	 * Get a list of all warnings that have occurred. Calling this method does not reset the warning condition.
	 * 
	 * @return A list of all warnings that have occurred, or an empty list if there are no warnings to report.
	 */
	public ConditionEvent[] getWarningConditions() {
		return warningList.toArray(new ConditionEvent[0]);
	}
	
	
	//
	// Other supporting methods
	//
	
	/**
	 * Flush the current sample and add it to <code>KestrelRunner</code>. This method clears
	 * <code>inputSourceQueue</code> and sets <code>sampleName</code> to <code>null</code>.
	 */
	private void flushSample() {
		
		if (! inputSourceQueue.isEmpty())
			runnerBase.addSample(new InputSample(sampleName, inputSourceQueue.toArray(new SequenceSource[0])));
		
		inputSourceQueue.clear();
		sampleName = null;
	}
}
