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

package edu.gatech.kestrel.io;

import edu.gatech.kanalyze.comp.reader.SequenceSource;

/**
 * Represents one sample and all the files that contain it.
 * <p/>
 * Kestrel can support multiple samples, and samples may be in more than one input file
 * (e.g. paired-end reads). This class organizes all the input files for the sample.
 */
public class InputSample {
	
	/** Name of this sample. Will never be <code>null</code> or empty. */
	public final String name;
	
	/**
	 * An array of input sources for this sample. Will never be null or empty, and at
	 * instantiation, will not contain <code>null</code> references.
	 */
	public final SequenceSource[] sources;
	
	/**
	 * Create a new input sample.
	 * 
	 * @param name Name of this sample. If <code>null</code> or empty, the name is assigned
	 *   using the name of the first element of <code>sources</code>.
	 * @param sources An array of input sources.
	 * 
	 * @throws NullPointerException If <code>sources</code> is <code>null<code>.
	 * @throws IllegalArgumentException If <code>sources</code> is empty or contains
	 *   <code>null</code> references.
	 */
	public InputSample(String name, SequenceSource[] sources)
			throws NullPointerException, IllegalArgumentException {
		
		// Check sources
		if (sources == null)
			throw new NullPointerException("Cannot create sample with sources: null");
		
		if (sources.length == 0)
			throw new IllegalArgumentException("Cannot create sample with no input sources");
		
		for (SequenceSource source : sources)
			if (source == null)
				throw new IllegalArgumentException("Cannot create sample: Sequence source list contains null references");
		
		// Null name if empty
		if (name != null) {
			name = name.trim();
			
			if (name.isEmpty())
				name = null;
		}
		
		// Assign name if null
		if (name == null)
			name = sources[0].name;
		
		// Assign fields
		this.name = name;
		this.sources = sources;
		
		return;
	}
}
