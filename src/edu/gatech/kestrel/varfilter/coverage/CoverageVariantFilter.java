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

package edu.gatech.kestrel.varfilter.coverage;

import edu.gatech.kestrel.varfilter.VariantFilter;
import edu.gatech.kestrel.variant.VariantCall;

/**
 * Filter variants based on their location in the reference.
 */
public class CoverageVariantFilter extends VariantFilter {
	
	/** Variants must end this many bases before the right end of the reference. */
	private double minCoverage;
	
	/** Variants must be supported by at least this many k-mers. */
	private int minDepth;
	
	
	
	/**
	 * Create a new variant filter.
	 */
	public CoverageVariantFilter() {
		super("CoverageVariantFilter");
		
		// Set defaults
		minCoverage = 0.0;
		minDepth = Integer.MAX_VALUE;
		
		return;
	}
	
	/**
	 * Initialize this filter. The argument string is a comma-separated list of key-value pairs. Some
	 * keywords may be used without a value.
	 * 
	 * @param argStr A comma-separated list of key-value pairs.
	 */
	@Override
	protected void init(String argStr)
			throws IllegalArgumentException {
		
		String[] argTok;  // argStr split on comma
		String[] attrTok;  // Attribute ([0]) and value ([1]) from one argument
		
		double minCoverage = 0.0;  // Minimum coverage
		int minDepth = 0;          // Minimum depth
		
		boolean arg0IsBare;  // Set to true if the first argument does not have an attribute name (bare number is coverage)
		
		// Check arguments
		argStr = argStr.trim();
		
		if (argStr.isEmpty())
			throw new IllegalArgumentException("Cannot configure coverage variant filter with empty arguments");
		
		// Init
		arg0IsBare = false;
		
		// Tokenize
		argTok = argStr.split("\\s*,\\s*");
		
		for (int argTokIndex = 0; argTokIndex < argTok.length; ++argTokIndex) {
			attrTok = argTok[argTokIndex].split("=", 2);
			
			if (attrTok[0].equals("coverage") || attrTok[0].equals("cov") || attrTok[0].equals("c")) {
				
				try {
					minCoverage = Double.parseDouble(attrTok[0]);
					
				} catch (NumberFormatException ex) {
					throw new IllegalArgumentException(String.format("Cannot configure coverage variant filter: Argument at positon %d is not a number: %s", argTokIndex + 1, argTok[argTokIndex]));
				}
				
				if (minCoverage < 0.0 || minCoverage > 1.0)
					throw new IllegalArgumentException(String.format("Cannot configure coverage variant filter: Argument at positon %d: Minimum coverage must be between 0.0 and 1.0 (inclusive): %f (%s)", argTokIndex + 1, minCoverage, argTok[argTokIndex]));
				
			} else if (attrTok[0].equals("depth") || attrTok[0].equals("dep") || attrTok[0].equals("d")) {
				
				try {
					minDepth = Integer.parseInt(attrTok[0]);
					
				} catch (NumberFormatException ex) {
					throw new IllegalArgumentException(String.format("Cannot configure depth variant filter: Argument at position %d: Minimum depth is not an integer: %s", argTokIndex + 1, argTok[argTokIndex]));
				}
				
				if (minDepth < 0)
					throw new IllegalArgumentException(String.format("Cannot configure depth variant filter: Argument at positon %d: Minimum depth must be greater than 0: %d (%s)", argTokIndex + 1, minDepth, argTok[argTokIndex]));
				
			} else {
				
				// The first two arguments may be a bare number
				if (attrTok.length == 1) {
					
					if (argTokIndex == 0) {
						// 1st position bare number: Coverage
						
						try {
							minCoverage = Double.parseDouble(attrTok[0]);
							
						} catch (NumberFormatException ex) {
							throw new IllegalArgumentException(String.format("Cannot configure coverage variant filter: Argument at positon 1 is not a number: %s", argTok[argTokIndex]));
						}
						
						if (minCoverage < 0.0 && minCoverage > 1.0)
							throw new IllegalArgumentException(String.format("Cannot configure coverage variant filter: Argument at positon 1: Minimum coverage must be between 0.0 and 1.0 (inclusive): %d", minCoverage));
						
						arg0IsBare = true;
						
					} else if (argTokIndex == 1) {
						// 1st position bare number: Depth
						
						if (! arg0IsBare)
							throw new IllegalArgumentException(String.format("Cannot configure coverage variant filter: Argument at positon 2 must have an attribute name if the first parameter does (the first two may be bare numbers): %s", argTok[argTokIndex]));
						
						// Get depth
						try {
							minDepth = Integer.parseInt(attrTok[0]);
							
						} catch (NumberFormatException ex) {
							throw new IllegalArgumentException(String.format("Cannot configure coverage variant filter: Argument at positon 2 is not an integer: %s", argTok[argTokIndex]));
						}
						
						if (minDepth < 0)
							throw new IllegalArgumentException(String.format("Cannot configure coverage variant filter: Argument at positon 2: Minimum depth must be at least 0: %d", minDepth));
						
					} else {
						throw new IllegalArgumentException(String.format("Cannot configure coverage variant filter: Unknown parameter at position %d: %s", argTokIndex + 1, argTok[argTokIndex]));
					}
					
				} else {
					throw new IllegalArgumentException(String.format("Cannot configure coverage variant filter: Unknown parameter at position %d: %s (%s)", argTokIndex + 1, attrTok[0], argTok[argTokIndex]));
				}
			}
		}
		
		// Set attributes
		this.minCoverage = minCoverage;
		this.minDepth = minDepth;
		
		return;
	}
	
	/**
	 * Get a description of this filter.
	 * 
	 * @return A description of this filter.
	 */
	@Override
	public String getDescription() {
		return "Filter variants by coverage and depth.";
	}
	
	/**
	 * Filter variants.
	 * 
	 * @return <code>variantCall</code> if it was not filtered, and <code>null</code> if it was.
	 */
	@Override
	public VariantCall filter(VariantCall variantCall) {
		
		if (variantCall == null ||
			variantCall.variantDepth < minDepth ||
			(double) variantCall.variantDepth / variantCall.locusDepth < minCoverage) {
			
			return null;
		}
		
		return variantCall;
	}
}
