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

package edu.gatech.kestrel.hapwriter;

import edu.gatech.kestrel.activeregion.Haplotype;

/**
 * Discards haplotypes. Used as a placeholder for a haplotype writer when haplotypes
 * are not output.
 */
public class NullHaplotypeWriter extends HaplotypeWriter {
	
	/**
	 * Create a new writer.
	 */
	public NullHaplotypeWriter() {
		super("EmptyHaplotypeWriter");
		
		return;
	}
	
	/**
	 * Ignore haplotype.
	 * 
	 * @param haplotype Haplotype to ignore.
	 */
	@Override
	public void add(Haplotype haplotype) {
		
		return;
	}
	
	/**
	 * Initialize.
	 */
	@Override
	protected void implInit(String writerSpec) {
		
		return;
	}
	
	/**
	 * Get a description of this writer.
	 * 
	 * @return Description of this writer.
	 */
	@Override
	public String getDescription() {
		
		return "Discards haplotypes without writing them.";
	}
}
