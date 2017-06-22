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

package edu.gatech.kestrel.interval.bed;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.regex.Pattern;

import edu.gatech.kanalyze.util.StringUtil;
import edu.gatech.kestrel.interval.IntervalReader;
import edu.gatech.kestrel.interval.RegionInterval;

/**
 * Read intervals from BED files.
 */
public class BedIntervalReader extends IntervalReader {
	
	/** Matches lines to be ignored. */
	private static final Pattern IGNORE_PATTERN = Pattern.compile("(#|(browser|track)\\s).*");
	
	/** Delimiter for the colmns of the BED files processed by this reader. */
	private String fieldDelim;
	
	/** Delimiter for files where columns are strictly tab delimited. This is the default delimiter. */
	private static final String TAB_DELIM_STRING = "\t";
	
	/** Delimiter for files where any whitespace delimits fields. */
	private static final String WS_DELIM_STRING = "\\s+";
	
	/** Name of this clas of reader. */
	public static final String NAME = "BED";
	
	/**
	 * Create a new BED interval reader.
	 */
	public BedIntervalReader() {
		
		super(NAME);
		
		fieldDelim = TAB_DELIM_STRING;
		
		return;
	}
	
	/**
	 * Set an argument for this reader.
	 * 
	 * @param attribute Attribute. Never <code>null</code> or empty, and never contains whitespace
	 *   or non-printable characters. Whitespace is trimmed, and the name is converted to lower case.
	 * @param value Value. May be <code>null</code> if no value was found or was empty. Whitespace
	 *   is trimmed.
	 * @param avpCount The number of this attribute-value pair in a list of pairs. The first pair is
	 *   <code>1</code>, the second is <code>2</code>, etc. 
	 * 
	 * @throws IllegalArgumentException If the attribute is not supported or there is any error setting
	 *   the value.
	 */
	@Override
	protected void setReaderArg(String attribute, String value, int avpCount)
			throws IllegalArgumentException {
		
		if (value == null)
			throw new IllegalArgumentException(String.format("%s interval reader: Attribute \"%s\" is missing a value (attribute/value pair %d): %d", name, attribute, avpCount));
		
		if (attribute.equals("tabonly")) {
			
			try {
				if (StringUtil.toBoolean(value))
					fieldDelim = TAB_DELIM_STRING;
				else
					fieldDelim = WS_DELIM_STRING;
				
			} catch (IllegalArgumentException ex) {
				throw new IllegalArgumentException(String.format("%s interval reader: Value for attribute \"tabonly\" is not boolean (attribute/value pair %d): %d", name, avpCount, value));
			}
			
		} else {
			throw new IllegalArgumentException(String.format("%s interval reader: Unknown attribute (attribute/value pair %d): %d", name, avpCount, attribute));
		}
		
		return;
	}
	
	/**
	 * Get a description of this reader.
	 * 
	 * @return A description of this reader.
	 */
	@Override
	public String getDescription() {
		
		return "BED File.";
	}
	
	/**
	 * Get a filename pattern.
	 * 
	 * @return A filename pattern.
	 */
	@Override
	public Pattern getFilenamePattern() {
		return Pattern.compile(".*\\.bed$", Pattern.CASE_INSENSITIVE);
	}
	
	/**
	 * Get intervals from a BED file.
	 * 
	 * @param intervalFile File to read.
	 * 
	 * @return An array of intervals found in the file.
	 * 
	 * @throws NullPointerException If <code>intervalFile</code> is <code>null</code>.
	 * @throws FileNotFoundException If <code>intervalFile</code> is not found.
	 * @throws IOException If there is any error reading the file or with the format
	 *   of the file contents.
	 */
	@Override
	public RegionInterval[] read(File intervalFile)
			throws NullPointerException, FileNotFoundException, IOException {
		
		// Declarations
		ArrayList<RegionInterval> intervals;  // Intervals to be returned
		
		String line;   // Line buffer
		String[] tok;  // Tokenized line (split on tabs)
		
		int lineNum;  // Line number
		
		String chrom;    // BED field: Chromosome name
		int chromStart;  // BED field: Start position (inclusive, first base is 0)
		int chromEnd;    // BED field: End position (exclusive, first base is 0)
		String name;     // BED field: Name (optional)
		boolean isFwd;   // BED field: Strand (true if +, false if -)
		
		// Check arguments
		if (intervalFile == null)
			throw new NullPointerException("Cannot open BED file for reading intervals: null");
		
		// Init
		intervals = new ArrayList<>();
		lineNum = 0;
		
		try (BufferedReader reader = new BufferedReader(new FileReader(intervalFile))) {
			
			// Read bed file
			READ_BED:
			while ((line = reader.readLine()) != null) {
				++lineNum;
				line = line.trim();
				
				// Check for lines to be ignored (blank lines, comments, track, and browser lines)
				if (line.isEmpty() || IGNORE_PATTERN.matcher(line).matches())
					continue READ_BED;
				
				// Split fields
				tok = line.split(fieldDelim);
				
				if (tok.length < 3)
					throw new IOException(String.format("Bad BED record on line %d in file %s: Must contain at least 3 fields, but only found %d", lineNum, intervalFile.getPath(), tok.length));
				
				// Get chromosome
				chrom = tok[0].trim();
				
				if (chrom.isEmpty())
					throw new IOException(String.format("Bad BED record on line %d in file %s: Chromosome name (first field) is empty", lineNum, intervalFile.getPath()));
				
				// Get start position
				try {
					chromStart = Integer.parseInt(tok[1]);
					
					if (chromStart < 0)
						throw new IOException(String.format("Bad BED record on line %d in file %s: Start location (second field) is negative: %d", lineNum, intervalFile.getPath(), chromStart));
					
				} catch (NumberFormatException ex) {
					throw new IOException(String.format("Bad BED record on line %d in file %s: Start location (second field) is not an integer: %s", lineNum, intervalFile.getPath(), tok[1]));
				}
				
				// Get end position
				try {
					chromEnd = Integer.parseInt(tok[2]);
					
					if (chromEnd < 0)
						throw new IOException(String.format("Bad BED record on line %d in file %s: End location (third field) is negative: %d", lineNum, intervalFile.getPath(), chromEnd));
					
					if (chromEnd <= chromStart)
						throw new IOException(String.format("Bad BED record on line %d in file %s: End location (%d) must be greater than the start location (%d)", lineNum, intervalFile.getPath(), chromEnd, chromStart));
					
				} catch (NumberFormatException ex) {
					throw new IOException(String.format("Bad BED record on line %d in file %s: End location (third field) is not an integer: %s", lineNum, intervalFile.getPath(), tok[2]));
				}
				
				// Get name
				if (tok.length < 4)
					name = null;
				else
					name = tok[3].trim();
				
				// Get strand
				if (tok.length >= 6) {
					
					tok[5] = tok[5].trim();
					
					if (tok[5].equals("+"))
						isFwd = true;
					else if (tok[5].equals("-"))
						isFwd = false;
					else
						throw new IOException(String.format("Bad BED record on line %d in file %s: Strand (sixth field) must be \"+\" or \"-\": %s", lineNum, intervalFile.getPath(), tok[5]));
					
				} else {
					isFwd = true;
				}
				
				// Add interval
				intervals.add(new RegionInterval(name, chrom, chromStart + 1, chromEnd, isFwd));
			}
			
			
		} catch (FileNotFoundException ex) {
			throw new FileNotFoundException("Interval BED file was not found: " + intervalFile.getPath());
			
		} catch (IOException ex) {
			throw new IOException("IO error reading interval BED file: " + ex.getMessage(), ex);
		}
		
		// Return intervals
		return intervals.toArray(new RegionInterval[0]);
	}
}
