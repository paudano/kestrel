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

package edu.gatech.kestrel.counter;

import edu.gatech.kanalyze.module.count.CountModule;
import edu.gatech.kanalyze.util.KmerCounter;
import edu.gatech.kanalyze.util.kmer.KmerUtil;

/**
 * An in-memory k-mer counter.
 */
public class MemoryCountMap extends CountMap {
	
	/** In-memory K-mer count hash set. */
	private final KmerCounter counter;
	
	/**
	 * Create an in-memory count mapper.
	 * 
	 * @param kUtil K-mer utility.
	 * @param countModule A configured count module.
	 * 
	 * @throws NullPointerException If <code>kUtil</code> or <code>countModule</code>
	 *   is <code>null</code>.
	 */
	public MemoryCountMap(KmerUtil kUtil, CountModule countModule)
			throws NullPointerException {
		
		super(kUtil, countModule);  // throws NullPointerException
		
		counter = new KmerCounter(kUtil);
		countModule.setOutput(counter);
		
		return;
	}
	
	/**
	 * Called after the sample is set on the count module, but before it is run.
	 * 
	 * @see edu.gatech.kestrel.counter.CountMap#preModuleRun()
	 */
	@Override
	protected boolean preModuleRun() {
		
		counter.clear();
		
		return true;
	}
	
	/**
	 * Get a k-mer from this map.
	 * 
	 * @see edu.gatech.kestrel.counter.CountMap#get(int[])
	 */
	@Override
	public int get(int[] kmer) throws NullPointerException, IllegalArgumentException {
		
		return counter.get(kmer);  // throws NullPointerException, IllegalArgumentException
	}
}
