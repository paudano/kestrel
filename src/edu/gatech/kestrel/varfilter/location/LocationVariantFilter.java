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

package edu.gatech.kestrel.varfilter.location;

import edu.gatech.kestrel.varfilter.VariantFilter;
import edu.gatech.kestrel.variant.VariantCall;

/**
 * Filter variants based on their location in the reference.
 */
public class LocationVariantFilter extends VariantFilter {
	
	/** Variants must start this many bases from the left end of the reference. */
	private int leftEndLength;
	
	/** Variants must end this many bases before the right end of the reference. */
	private int rightEndLength;
	
	/**
	 * Create a new variant filter.
	 */
	public LocationVariantFilter() {
		super("LocationVariantFilter");
		
		// Set defaults
		leftEndLength = 0;
		rightEndLength = 0;
		
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
		
		int leftEndLength = 0;           // Left end length attribute
		int rightEndLength = 0;          // Right end length attribute
		
		int intVal;  // Integer value
		
		// Check arguments
		argStr = argStr.trim();
		
		if (argStr.isEmpty())
			throw new IllegalArgumentException("Cannot configure location variant filter with empty arguments");
		
		// Tokenize
		argTok = argStr.split("\\s*,\\s*");
		
		for (int argTokIndex = 0; argTokIndex < argTok.length; ++argTokIndex) {
			attrTok = argTok[argTokIndex].split("=", 2);
			
			if (attrTok[0].equals("length") || attrTok[0].equals("len") || attrTok[0].equals("l")) {
				// Set both ends
				
				if (attrTok.length == 1)
					throw new IllegalArgumentException(String.format("Cannot configure location variant filter: End length value at position %d is missing the length value: %s", argTokIndex + 1, argTok[argTokIndex]));
				
				try {
					intVal = Integer.parseInt(attrTok[1]);
					
				} catch (NumberFormatException ex) {
					throw new IllegalArgumentException(String.format("Cannot configure location variant filter: End length value at position %d is not a number: %s", argTokIndex + 1, argTok[argTokIndex]));
				}
				
				if (intVal < 0)
					throw new IllegalArgumentException(String.format("Cannot configure location variant filter: End length value at position %d is negative: %s (value = %d)", argTokIndex + 1, argTok[argTokIndex], intVal));
				
				leftEndLength = intVal;
				rightEndLength = intVal;
				
			} else if (attrTok[0].equals("leftlength") || attrTok[0].equals("leftlen") || attrTok[0].equals("leftl") ||
					attrTok[0].equals("llength") || attrTok[0].equals("llen") || attrTok[0].equals("ll")) {
				// Set left end
				
				if (attrTok.length == 1)
					throw new IllegalArgumentException(String.format("Cannot configure location variant filter: Left end length value at position %d is missing the length value: %s", argTokIndex + 1, argTok[argTokIndex]));
				
				try {
					intVal = Integer.parseInt(attrTok[1]);
					
				} catch (NumberFormatException ex) {
					throw new IllegalArgumentException(String.format("Cannot configure location variant filter: Left end length value at position %d is not a number: %s", argTokIndex + 1, argTok[argTokIndex]));
				}
				
				if (intVal < 0)
					throw new IllegalArgumentException(String.format("Cannot configure location variant filter: Left end length value at position %d is negative: %s (value = %d)", argTokIndex + 1, argTok[argTokIndex], intVal));
				
				leftEndLength = intVal;
				
			} else if (attrTok[0].equals("rightlength") || attrTok[0].equals("rightlen") || attrTok[0].equals("rightl") ||
					attrTok[0].equals("rlength") || attrTok[0].equals("rlen") || attrTok[0].equals("rl")) {
				// Set right end
				
				if (attrTok.length == 1)
					throw new IllegalArgumentException(String.format("Cannot configure location variant filter: Right end length value at position %d is missing the length value: %s", argTokIndex + 1, argTok[argTokIndex]));
				
				try {
					intVal = Integer.parseInt(attrTok[1]);
					
				} catch (NumberFormatException ex) {
					throw new IllegalArgumentException(String.format("Cannot configure location variant filter: Right end length value at position %d is not a number: %s", argTokIndex + 1, argTok[argTokIndex]));
				}
				
				if (intVal < 0)
					throw new IllegalArgumentException(String.format("Cannot configure location variant filter: Right end length value at position %d is negative: %s (value = %d)", argTokIndex + 1, argTok[argTokIndex], intVal));
				
				rightEndLength = intVal;
				
			} else {
				throw new IllegalArgumentException(String.format("Cannot configure location variant filter: Unrecognized option at position %d: %s", argTokIndex + 1, argTok[argTokIndex]));
			}
		}
		
		// Check attributes
		if (leftEndLength == 0 && rightEndLength == 0)
			throw new IllegalArgumentException("Cannot configure location variant filter: Either the left end length, the right end length, or both must be set to a value greater than 0.");
		
		// Set attributes
		this.leftEndLength = leftEndLength;
		this.rightEndLength = rightEndLength;
		
		return;
	}
	
	/**
	 * Get a description of this filter.
	 * 
	 * @return A description of this filter.
	 */
	@Override
	public String getDescription() {
		return "Filter variants by distance from reference ends.";
	}
	
	/**
	 * Filter variants.
	 * 
	 * @return <code>variantCall</code> if it was not filtered, and <code>null</code> if it was.
	 */
	@Override
	public VariantCall filter(VariantCall variantCall) {
		
		// Check arguments
		if (variantCall == null)
			return null;
		
		if (variantCall.start <= leftEndLength)
			return null;
		
		if (variantCall.referenceEnd() >= variantCall.activeRegion.refRegion.size - rightEndLength)
			return null;
		
		return variantCall;
	}
}
