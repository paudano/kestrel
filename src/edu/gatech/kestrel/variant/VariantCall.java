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

import edu.gatech.kestrel.activeregion.ActiveRegion;
import edu.gatech.kestrel.activeregion.Haplotype;
import edu.gatech.kestrel.util.SystemUtil;

/**
 * Base class of all variant calls.
 */
public abstract class VariantCall implements Comparable<VariantCall> {
	
	/** Active region this variant was called on. */
	public final ActiveRegion activeRegion;
	
	/** Haplotypes this variant was called in. */
	private final Haplotype[] haplotype;
	
	/** Type of this variant call. */
	public final VariantType type;
	
	/**
	 * Location in the reference sequence where this variant begins. Note that the first
	 * base in the reference is position <code>1</code> (this is not zero-based
	 * array index).
	 */
	public final int start;
	
	/** Estimated read depth of this variant. */
	public final int variantDepth;
	
	/** Estimated read depth of the locus this variant is in. */
	public final int locusDepth;
	
	/** Sequence of the reference at the point of this variant. */
	public final String ref;
	
	/** Sequence of the variant over <code>ref</code>. */
	public final String alt;
	
	/** <code>true</code> if this variant covers an ambiguous base in the reference sequence. */
	public final boolean isAmbiguous;
	
	/**
	 * <code>true</code> if the location of this variant is relative to the reference sequence
	 * and <code>false</code> if it is relative to a region within the reference sequence.
	 */
	public final boolean isReferenceAligned;
	
	/** Pattern for ambiguous bases. */
	private static final Pattern AMBIG_BASE_PATTERN = Pattern.compile(".*[^ACGTUacgtu].*");
	
	/**
	 * Create a variant call.
	 * 
	 * @param type Type of this variant.
	 * @param start Start position of this variant.
	 * @param variantDepth Depth of the haplotypes including the variant.
	 * @param locusDepth Depth of the wild-type active region and haplotypes not including
	 *   this variant.
	 * @param ref Reference sequence or <code>null</code> to use an empty string.
	 * @param alt Alternate sequence or <code>null</code> to use an empty string.
	 * @param haplotype A list of haplotypes the variant was found in.
	 * @param haploSize Number of non-null elements of <code>haplotype</code>.
	 * @param isReferenceAligned <code>true</code> if the location of this variant is relative
	 *   to the reference sequence and <code>false</code> if it is relative to a region within
	 *   the reference sequence.
	 * 
	 * @throws NullPointerException If <code>type</code> or <code>haplotype</code> is <code>null</code>. 
	 * @throws IllegalArgumentException If <code>haplotype</code> is empty or contains <code>null</code>
	 *   references within <code>haploSize</code> elements, <code>haploSize</code> is less than <code>0</code>
	 *   or greater than <code>haplotype.length</code>, <code>start</code> is less than <code>1</code> or does
	 *   not fit in the active region in the <code>haplotype</code> elements, <code>variantDepth</code> is
	 *   negative, <code>locusDepth</code> is less than <code>variantDepth</code>, or all <code>haplotypes</code>
	 *   do not reference the same active region.
	 */
	protected VariantCall(VariantType type, int start, int variantDepth, int locusDepth, String ref, String alt, Haplotype[] haplotype, int haploSize, boolean isReferenceAligned)
			throws NullPointerException, IllegalArgumentException {
		
		// Check arguments
		if (type == null)
			throw new NullPointerException("Variant type is null");
		
		if (haplotype == null)
			throw new NullPointerException("Cannot create variant with haplotypes: null");
		
		if (haplotype.length == 0)
			throw new IllegalArgumentException("Haplotype array is empty");
		
		if (start < 1)
			throw new IllegalArgumentException("Start position is less than 1: " + start);
		
		if (variantDepth < 0)
			throw new IllegalArgumentException("Variant depth is negative: " + variantDepth);
		
		if (locusDepth < variantDepth)
			throw new IllegalArgumentException(String.format("Locus depth (%d) is less than variant depth (%d)", locusDepth, variantDepth));
		
		if (ref == null)
			ref = "";
		
		if (alt == null)
			alt = "";
		
		if (haploSize < 1 || haploSize > haplotype.length)
			throw new IllegalArgumentException("Haplotype size (haploSize) is less than 1 or is greater than the haplotype length: " + haploSize + " (haplotype length = " + haplotype.length + ")");
		
		// Set fields
		this.type = type;
		this.start = start;
		this.variantDepth = variantDepth;
		this.locusDepth = locusDepth;
		this.ref = ref.trim();
		this.alt = alt.trim();
		this.isReferenceAligned = isReferenceAligned;
		
		// Determine if ambiguous bases are present
		isAmbiguous = AMBIG_BASE_PATTERN.matcher(this.ref).matches();
		
		// Copy haplotype
		haplotype = Arrays.copyOf(haplotype, haploSize);
		this.haplotype = haplotype;
		
		// Save and check active region
		activeRegion = haplotype[0].activeRegion;

		// TODO: May need to remove these commented lines
//		if (start <= activeRegion.startIndex)  // (start <= activeRegion.start): start counts from 1, activeRegion.start from 0
//			throw new IllegalArgumentException(String.format("Variant starts at a position before the active region: Variant start = %d, Active region start = %d", start, activeRegion.startIndex + 1));
		
		// Check for null in haplotype and that all active regions agree
		for (int index = 0; index < haplotype.length; ++index) {
			if (haplotype[index] ==  null)
				throw new IllegalArgumentException("Haplotype array contains a null reference at index " + index);
			
			if (haplotype[index].activeRegion != activeRegion)
				throw new IllegalArgumentException(String.format(
						"Haplotype array contains haplotypes from different active regions: index 0 = %s (%s) index %d = %s (%s)",
						activeRegion.toString(),
						SystemUtil.objectToString(activeRegion),
						index,
						haplotype[index].activeRegion.toString(),
						SystemUtil.objectToString(haplotype[index].activeRegion)
				));
		}
		
		return;
	}
	
	/**
	 * Get the location on the reference where this variant ends.
	 * 
	 * @return Location on the reference where this variant ends.
	 */
	public abstract int referenceEnd();
	
	/**
	 * Get an array of haplotypes supporting this variant.
	 * 
	 * @return An array of haplotypes supporting this variant.
	 */
	public Haplotype[] getHaplotypes() {
		return Arrays.copyOf(haplotype, haplotype.length);
	}
	
	/**
	 * Get the location of this variant as it would appear in a VCF file.
	 * 
	 * @return VCF location of this variant.
	 */
	public abstract int getVcfPos();
	
	/**
	 * Get the reference sequence of this variant.
	 * 
	 * @return Reference sequence of this variant.
	 */
	public abstract String getVcfRef();
	
	/**
	 * Get the alternate sequence of this variant.
	 * 
	 * @return Alternate sequence of this variant.
	 */
	public abstract String getVcfAlt();
	
	/**
	 * Get an HGVS-formatted string description of this variant.
	 * 
	 * @return String description of this variant.
	 * 
	 * @see <a href="http://www.hgvs.org/mutnomen/recs-DNA.html">http://www.hgvs.org/mutnomen/recs-DNA.html</a>
	 */
	@Override
	public abstract String toString();
	
	/**
	 * Get the variant string format with reference sequence name.
	 * 
	 * @return Variant with other useful debugging information.
	 */
	public String toStringFull() {
		return String.format("%s:%s", activeRegion.refRegion.referenceSequence.name, this.toString());
	}
	
	/**
	 * Compare this variant call to another.
	 * 
	 * @param other Other variant call.
	 * 
	 * @throws NullPointerException If <code>other</code> is <code>null</code>.
	 */
	@Override
	public int compareTo(VariantCall other)
		throws NullPointerException {
		
		int compareVal;  // Result of comparing subfields
		
		// Check arguments
		if (other == null)
			throw new NullPointerException("Cannot compare to variant: null");
		
		// Reference chromosome or contig
		compareVal = activeRegion.refRegion.referenceSequence.name.compareTo(other.activeRegion.refRegion.referenceSequence.name);
		
		if (compareVal != 0)
			return compareVal;
		
		// Position
		if (start != other.start)
			return start - other.start;
		
		// Variant type
		if (type != other.type)
			return type.compareTo(other.type);
		
		// Reference sequence
		compareVal = ref.compareTo(other.ref);
		
		if (compareVal != 0)
			return compareVal;
		
		// Alternate sequence (return 0 if equal)
		return alt.compareTo(other.alt);
	}
}
