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

package edu.gatech.kestrel.variant;

import java.util.Arrays;
import java.util.regex.Pattern;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import edu.gatech.kanalyze.util.RBTree;
import edu.gatech.kestrel.KestrelConstants;
import edu.gatech.kestrel.activeregion.ActiveRegion;
import edu.gatech.kestrel.activeregion.Haplotype;
import edu.gatech.kestrel.align.AlignNode;

/**
 * Scans haplotypes for variants and characterizes them.
 */
public class VariantCaller {
	
	/** Logger object. */
	private static final Logger logger = LoggerFactory.getLogger(VariantCaller.class);;
	
	/** Variant tree. */
	private RBTree<CallerVarNode,CallerVarNode> varTree;
	
	/** Region variants are called on. */
	private ActiveRegion activeRegion;
	
	/** Total depth of the active region (wild-type) and all haplotypes. */ 
	private int totalDepth;
	
	/** Call variants relative to the reference sequence. */
	private boolean variantCallByReference;
	
	/** Builder for inserted and deleted strings. */
	private StringBuilder strBuilder;
	
	/** If <code>true</code>, allow variants that span ambiguous regions of the reference. */
	private boolean callAmbiguousVariant;
	
	/** Default allow-ambiguous reference variant. */
	public static final boolean DEFAULT_CALL_AMBIGUOUS_VARIANT = true;
	
	/** Pattern for ambiguous bases. */
	private static final Pattern AMBIG_BASE_PATTERN = Pattern.compile(".*[^ACGTUacgtu].*");
	
	/**
	 * Create a new variant caller.
	 */
	public VariantCaller() {
		
		varTree = new RBTree<CallerVarNode,CallerVarNode>();
		activeRegion = null;
		totalDepth = 0;
		
		callAmbiguousVariant = DEFAULT_CALL_AMBIGUOUS_VARIANT;
		
		variantCallByReference = true;
		
		return;
	}
	
	/**
	 * Initialize this variant caller.
	 * 
	 * @param activeRegion Active region.
	 * 
	 * @throws NullPointerException If <code>activeRegion</code> is <code>null</code>.
	 */
	public void init(ActiveRegion activeRegion)
			throws NullPointerException {
		
		if (activeRegion == null)
			throw new NullPointerException("Cannot initialize a variant caller over active region: null");
		
		this.activeRegion = activeRegion;
		
		varTree.clear();
		totalDepth = 0;
		
		strBuilder = new StringBuilder();
		
		return;
	}
	
	/**
	 * When set, variant calls are relative to the reference sequence. This
	 * will undo a call to <code>setVariantCallByRegion()</code>.
	 */
	public void setVariantCallByReference() {
		
		variantCallByReference = true;
		
		return;
	}
	
	/**
	 * Variant calls are relative to the reference region. This
	 * will undo a call to <code>setVariantCallByReference()</code>.
	 */
	public void setVariantCallByRegion() {
		
		this.variantCallByReference = false;
		
		return;
	}
	
	/**
	 * Determine if variant calls are relative to the reference sequence.
	 * 
	 * @return <code>true</code> if variant calls are relative to the reference
	 *   sequence.
	 */
	public boolean isVariantCallByReference() {
		return variantCallByReference;
	}
	
	/**
	 * Add a haplotype to this variant caller.
	 * 
	 * @param haplotype Haplotype to scan for variants.
	 * 
	 * @throws NullPointerException If <code>haplotype</code> is <code>null</code>.
	 * @throws IllegalArgumentException If the active region in this haplotype does not match
	 *   the current active region in this variant caller (set with <code>init()</code>).
	 * @throws IllegalStateException If <code>init()</code> was not called before this method.
	 */
	public void add(Haplotype haplotype)
			throws NullPointerException, IllegalArgumentException, IllegalStateException {
		
		// Declarations
		AlignNode alignNode;  // Node from alignment
		int refPosition;      // Index in the reference sequence
		int altPosition;      // Index in the consensus sequence
		int startPosition;    // Position where a variant starts (used for deletions to save refPosition)
		
		CallerVarNode newNode;
		CallerVarNode existingNode;
		
		// Check arguments
		if (haplotype == null)
			throw new NullPointerException("Cannot add haplotype: null");
		
		if (activeRegion == null)
			throw new IllegalStateException("Cannot add haplotype to variant caller before init() is called on the variant caller");
		
		if (haplotype.activeRegion != activeRegion)
			throw new IllegalArgumentException(String.format("Haplotype active region does not match the active region active in this caller: current = %d, haplotype region = %d", activeRegion.toString(), haplotype.activeRegion.toString()));
		
		// Add depth
		totalDepth += haplotype.stats.min;
		
		// Find variants
		refPosition = activeRegion.startIndex;
		altPosition = 0;
		
		alignNode = haplotype.alignment;
		
//		System.out.println("Processing alignment: " + alignNode.getCigarString());  // DBGTMP
//		System.out.printf("\t* Ref pos=%d, length=%d\n", refPosition, activeRegion.refRegion.sequence.length);  // DBGTMP
//		System.out.printf("\t* Alt pos=%d, length=%d\n", altPosition, haplotype.sequence.length);  // DBGTMP
		
		while (alignNode != null) {
			
			if (alignNode.type == AlignNode.MATCH) {
				
//				System.out.printf("\tMATCH: (n=%d, refPos=%d, altPos=%d)\n", alignNode.n, refPosition, altPosition);  // DBGTMP
				
				// Move to the next position
				refPosition += alignNode.n;
				altPosition += alignNode.n;
				
			} else if (alignNode.type == AlignNode.MISMATCH) {
				
//				System.out.printf("\tMISMATCH: (n=%d, refPos=%d, altPos=%d)\n", alignNode.n, refPosition, altPosition);  // DBGTMP
				
				for (int index = 0; index < alignNode.n; ++index) {
					
					// Create variant node
					newNode = new CallerVarNode(VariantType.SNP,
							refPosition + 1,
							"" + (char) activeRegion.refRegion.sequence[refPosition],
							"" + (char) haplotype.sequence[altPosition],
							haplotype
					);
					
					++refPosition;
					++altPosition;
				
					// Add/search
					if (newNode.isFlank) {
						logger.trace("Dropping variant over a region flank: {}", newNode.toString());
						
					} else if (newNode.isAmbiguous && ! callAmbiguousVariant) {
						logger.trace("Dropping variant over an ambiguous region: {}", newNode.toString());
						
					} else {
						existingNode = varTree.put(newNode, newNode);  // Add new variant (returns a variant if an equivalent variant exists)
						
						if (existingNode != null)
							newNode.addVariant(existingNode);  // Add evidence from existing variant
					}
				}
				
			} else if (alignNode.type == AlignNode.INS) {
				
//				System.out.printf("\tINSERTION: (n=%d, refPos=%d, altPos=%d)\n", alignNode.n, refPosition, altPosition);  // DBGTMP
				
				// Find end of variant and build inserted sequence
				strBuilder.setLength(0);
				
				for (int index = 0; index < alignNode.n; ++index)
					strBuilder.append((char) haplotype.sequence[altPosition++]);
				
				// Create node
				newNode = new CallerVarNode(VariantType.INSERTION,
						refPosition,
						"",
						strBuilder.toString(),
						haplotype
				);
				
				// Add/search
				if (newNode.isFlank) {
					logger.trace("Dropping variant over a region flank: {}", newNode.toString());
					
				} else if (newNode.isAmbiguous && ! callAmbiguousVariant) {
					logger.trace("Dropping variant over an ambiguous region: {}", newNode.toString());
					
				} else {
					existingNode = varTree.put(newNode, newNode);  // Add new variant (returns a variant if an equivalent variant exists)
					
					if (existingNode != null)
						newNode.addVariant(existingNode);  // Add evidence from existing variant
				}
				
			} else if (alignNode.type == AlignNode.DEL) {
				
//				System.out.printf("\tDELETION: (n=%d, refPos=%d, altPos=%d)\n", alignNode.n, refPosition, altPosition);  // DBGTMP
				
				// Find end of variant and build deleted sequence
				strBuilder.setLength(0);
				
				startPosition = refPosition;
				
				for (int index = 0; index < alignNode.n; ++index)
					strBuilder.append((char) activeRegion.refRegion.sequence[refPosition++]);
				
				// Create node
				newNode = new CallerVarNode(VariantType.DELETION,
						startPosition + 1,
						strBuilder.toString(),
						"",
						haplotype
				);
				
				// Add/search
				if (newNode.isFlank) {
					logger.trace("Dropping variant over a region flank: {}", newNode.toString());
					
				} else if (newNode.isAmbiguous && ! callAmbiguousVariant) {
					logger.trace("Dropping variant over an ambiguous region: {}", newNode.toString());
					
				} else {
					existingNode = varTree.put(newNode, newNode);  // Add new variant (returns a variant if an equivalent variant exists)
					
					if (existingNode != null)
						newNode.addVariant(existingNode);  // Add evidence from existing variant
				}
				
			} else {
				throw new IllegalStateException("Unrecognized alignment node type: " + alignNode.type);
			}
			
			// Next alignment node
			alignNode = alignNode.next;
		}
		
		return;
	}
	
	/**
	 * Get variants from this caller.
	 * 
	 * @return An array of variants in this caller.
	 */
	public VariantCall[] getVariants() {
		
		VariantCall[] variants = new VariantCall[varTree.size()];
		int index = 0;
		
		for (CallerVarNode varNode : varTree)
			variants[index++] = varNode.toVariant();
		
		return variants;
	}
	
	/**
	 * Set the property to allow variants over ambiguous bases in the reference.
	 * 
	 * @param allowAmbiguousReference <code>true</code> to allow variants over ambiguous
	 *   bases in the reference.
	 * 
	 * @see #DEFAULT_CALL_AMBIGUOUS_VARIANT
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
	 * A variant node used for building evidence for variants.
	 */
	private class CallerVarNode implements Comparable<CallerVarNode>{
		
		/** Type of this variant. */ 
		public final VariantType type;
		
		/** Start position of this variant. */
		public final int start;
		
		/** Reference string. */
		public final String refString;
		
		/** Alternate string. */
		public final String altString;
		
		/** <code>true</code> if this variant spans an ambiguous base in the reference. */
		public final boolean isAmbiguous;
		
		/** <code>true</code> if this variant spans into the either flank of a reference region. */
		public final boolean isFlank;
		
		/** Approximate depth of reads supporting this variant. */
		public int varDepth;
		
		/** An array of haplotypes. */
		public Haplotype haplotype[];
		
		/** Number of elements in <code>haplotype</code>. */
		public int haploSize;
		
		/** Capacity of <code>haplotype</code>. */
		public int haploCapacity;
		
		/**
		 * Create a new variant node.
		 * 
		 * @param type Type of this variant.
		 * @param start Start position of this variant by coordinates (first base is 1).
		 * @param refString Reference string.
		 * @param altString Alternate string.
		 * @param haplotype Haplotype this variant was found in.
		 */
		public CallerVarNode(VariantType type, int start, String refString, String altString, Haplotype haplotype) {
			
			// Check arguments
			assert (type != null) :
				"Variant type is null";
			
			assert (start > 0) :
				"Start is less than 1: " + start;
			
			assert (refString != null) :
				"Reference string is null";
			
			assert (altString != null) :
				"Alternate string is null";
			
			assert (haplotype != null) :
				"Haplotype is null";
			
			// Set fields
			this.type = type;
			this.start = start;
			this.refString = refString;
			this.altString = altString;
			
			// Initialize variant depth
			varDepth = haplotype.stats.min;
			
			// Set haplotype capacity
			haploCapacity = 3;
			this.haplotype = new Haplotype[haploCapacity];
			this.haplotype[0] = haplotype;
			
			haploSize = 1;
			
			// Set ambiguous flag
			isAmbiguous = AMBIG_BASE_PATTERN.matcher(refString).matches();
			
			// Set flank flag
			if (altString.length() == 0)
				isFlank = haplotype.activeRegion.refRegion.isFlankByCoordinate(start, start);
			else
				isFlank = haplotype.activeRegion.refRegion.isFlankByCoordinate(start, start + altString.length() - 1);
			
			return;
		}
		
		/**
		 * Get the immutable variant call from this node.
		 * 
		 * @return Immutable variant call.
		 * 
		 * @throws IllegalStateException If <code>type</code> is an unrecognized variant type. This
		 *   indicates a software bug (invalid types, or types this caller was not updated to handle).
		 */
		public VariantCall toVariant()
				throws IllegalStateException {
			
			// Set start location
			int start = this.start;
			
			start -= haplotype[0].activeRegion.refRegion.leftFlankLength;
			
			if (variantCallByReference)
				start += haplotype[0].activeRegion.refRegion.interval.start - 1;
			
			// Create variant
			if (type == VariantType.SNP) {
				
				return new VariantSNP(start, varDepth, totalDepth, refString, altString, haplotype, haploSize, variantCallByReference);
				
			} else if (type == VariantType.INSERTION) {
				
				return new VariantInsertion(start, varDepth, totalDepth, altString, haplotype, haploSize, variantCallByReference);
				
			} else if (type == VariantType.DELETION) {
				
				return new VariantDeletion(start, varDepth, totalDepth, refString, haplotype, haploSize, variantCallByReference);
			}
			
			// type unknown
			throw new IllegalStateException("Variant caller cannot handle unknown variant type: " + type);
		}
		
		/**
		 * Add evidence from an equivalent variant to this variant.
		 * 
		 * @param o Other variant.
		 */
		public void addVariant(CallerVarNode o) {
			
			// Check arguments
			assert (o != null) :
				"Cannot add variant: null";
			
			assert (compareTo(o) == 0) :
				"Cannot add variant: Variants are not equal";
			
			// Expand haplotype array
			while (haploSize + o.haploSize >= haploCapacity) {
				int newCapacity = (int) (haploCapacity * KestrelConstants.ARRAY_EXPAND_FACTOR);
				
				if (newCapacity < 0 || newCapacity > KestrelConstants.MAX_ARRAY_SIZE) {
					
					if (haploCapacity == KestrelConstants.MAX_ARRAY_SIZE)
						throw new IllegalStateException("Exceeded maximum array size for haplotypes: " + KestrelConstants.MAX_ARRAY_SIZE);
					
					newCapacity = KestrelConstants.MAX_ARRAY_SIZE;
				}
				
				haplotype = Arrays.copyOf(haplotype, newCapacity);
				haploCapacity = newCapacity;
			}
			
			// Copy evidence
			for (int index = 0; index < o.haploSize; ++index)
				haplotype[haploSize++] = o.haplotype[index];
			
			varDepth += o.varDepth;
			
			return;
		}
		
		/**
		 * Get a hash code for this variant.
		 * 
		 * @return Hash code.
		 */
		@Override
		public int hashCode() {
			return type.hashCode ^ start ^ refString.hashCode() ^ altString.hashCode();
		}
		
		/**
		 * Compare this variant to another and return a negative number, <code>0</code>, or a positive
		 * number if this object is less than, equal to, or greater than <code>o</code> (respectively).
		 * 
		 * @return Compare value as described above, or a positive number of <code>o</code> is
		 *   <code>null</code>.
		 */
		@Override
		public int compareTo(CallerVarNode o) {
			
			int cmpVal;
			
			// Null first
			if (o == null)
				return 1;
			
			// Start position
			if (start != o.start)
				return start - o.start;
			
			// Type
			if (type != o.type)
				return type.compareTo(o.type);
			
			// Reference string
			cmpVal = refString.compareTo(o.refString);
			
			if (cmpVal != 0)
				return cmpVal;
			
			// Alternate string
			return altString.compareTo(o.altString);
		}
		
		/**
		 * Get a string describing thi call node.
		 * 
		 * @return A string describing this call node.
		 */
		@Override
		public String toString() {
			return String.format("CallerVarNode[type=%s, start=%d, depth=%d, ref=%s, alt=%s]", type, start, varDepth, refString, altString);
		}
	}
}
