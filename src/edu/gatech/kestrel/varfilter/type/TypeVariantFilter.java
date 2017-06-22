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

package edu.gatech.kestrel.varfilter.type;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import edu.gatech.kestrel.varfilter.VariantFilter;
import edu.gatech.kestrel.variant.VariantCall;

/**
 * Filter variants based on their location in the reference.
 */
public class TypeVariantFilter extends VariantFilter {
	
	/** Logger. */
	private static final Logger logger = LoggerFactory.getLogger(TypeVariantFilter.class);
	
	/** Include SNP if <code>true</code>. */
	private boolean typeSnp;
	
	/** Include insertion if <code>true</code>. */
	private boolean typeIns;
	
	/** Include deletion if <code>true</code>. */
	private boolean typeDel;
	
	/**
	 * Create a new variant filter.
	 */
	public TypeVariantFilter() {
		super("TypeVariantFilter");
		
		// Set defaults
		typeSnp = true;
		typeIns = true;
		typeDel = true;
		
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
		
		boolean typeSnp = false;
		boolean typeIns = false;
		boolean typeDel = false;
		
		// Check arguments
		argStr = argStr.trim();
		
		if (argStr.isEmpty())
			throw new IllegalArgumentException("Cannot configure type variant filter with empty arguments");
		
		// Tokenize
		argTok = argStr.split("\\s*,\\s*");
		
		for (int argTokIndex = 0; argTokIndex < argTok.length; ++argTokIndex) {
			String tokStr = argTok[argTokIndex];
			
			if (tokStr.matches("snp|snv")) {
				typeSnp = true;
				
			} else if (tokStr.matches("ins|insertion")) {
				typeIns = true;
				
			} else if (tokStr.matches("del|deletion")) {
				typeDel = true;
				
			} else if (tokStr.matches("indel|insdel|insertiondeletion")) {
				typeIns = true;
				typeDel = true;
				
			} else {
				throw new IllegalArgumentException("Unrecognized type: " + argTok[argTokIndex]);
			}
		}
		
		// Set attributes
		this.typeSnp = typeSnp;
		this.typeIns = typeIns;
		this.typeDel = typeDel;
		
		return;
	}
	
	/**
	 * Get a description of this filter.
	 * 
	 * @return A description of this filter.
	 */
	@Override
	public String getDescription() {
		return "Filter variants by type.";
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
		
		// Check type
		switch (variantCall.type) {
		case SNP:
			if (typeSnp)
				return variantCall;
			
			break;
			
		case INSERTION:
			if (typeIns)
				return variantCall;
			
			break;
			
		case DELETION:
			if (typeDel)
				return variantCall;
			
			break;
		
		default:
			logger.error("Unrecognized variant type: {} for variant {}", variantCall.type, variantCall);
		}
		
		// No match
		return null;
	}
}
