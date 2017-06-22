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

package edu.gatech.kestrel.align;

import edu.gatech.kestrel.KestrelConstants;

/**
 * A trace-back matrix used for debugging. Kestrel uses a linked-graph structure for traces,
 * and the full matrix is not easily retrieved from it. For example, no reference is stored
 * to nodes that are not traversible, and they may be garbage collected. This matrix will
 * use more memory and should only be used for debugging.
 */
public class TraceMatrix {
	
	/** Matrix of trace elements. */
	public short[][] matrix;
	
	/** Number of rows. */
	public final int nRow;
	
	/** Current column. */
	private int column;
	
	/** Column capacity. The matrix must be expanded before <code>nCol</code> reaches this value. */
	private int colCapacity;
	
	/** The number of columns this matrix is initialized with. */
	public static final int DEFAULT_COLUMN_CAPACITY = 250;
	
	/**
	 * Create a new trace matrix.
	 * 
	 * @param nRow Number of rows in the matrix, which is the length of the
	 *   reference sequence active region.
	 * 
	 * @throws IllegalArgumentException If <code>nRow</code> is less than <code>1</code>.
	 */
	public TraceMatrix(int nRow)
			throws IllegalArgumentException {
		
		if (nRow < 1)
			throw new IllegalArgumentException("Cannot create trace matrix with less than 1 row: " + nRow);
		
		this.nRow = nRow;
		
		column = -1;
		colCapacity = DEFAULT_COLUMN_CAPACITY;
		
		matrix = new short[nRow][colCapacity];
		
		return;
	}
	
	/**
	 * Set a field in the trace-back matrix.
	 * 
	 * @param row Row of the field to set where the first row is <code>0</code>.
	 * @param typeFrom One of the <code>TAB_</code> constants indicating the table
	 *   the transition takes place from. If not a <code>TAB_</code> constant, the results
	 *   are undefined.
	 * @param typeTo One of the <code>TAB_</code> constants indicating the table
	 *   the transition moves to. If not a <code>TAB_</code> constant, the results
	 *   are undefined.
	 * 
	 * @throws IllegalArgumentException If <code>row</code> or <code>col</code> is less
	 *   than <code>0</code> or if <code>row</code> is not less than <code>this.nRow</code>.
	 * @throws IllegalStateException If <code>nextCol()</code> was not called at least once
	 *   before this method.
	 */
	public void set(int row, int typeFrom, int typeTo)
			throws IllegalArgumentException, IllegalStateException {
		
		// Check arguments
		if (row >= nRow || row < 0)
			throw new IllegalArgumentException(String.format("Row is out of range [0, %d]: %d", nRow - 1, row));
		
		if (column < 0)
			throw new IllegalStateException("No columns in matrix: must call nextCol() before set()");
		
		if (typeFrom > 1)
			typeFrom -= 1;  // Adjust so MATCH and MISMATCH are both 1, GAP_REF is 2, and GAP_CON is 3
		
		// Add to matrix
		matrix[row][column] |= (short) (0x1 << ((3 - typeFrom) + 3 * (3 - typeTo)));
		
		return;
	}
	
	/**
	 * Add a column to this matrix.
	 * 
	 * @return The current column index.
	 * 
	 * @throws IllegalStateException If the columns are at maximum capacity,
	 *   <code>KestrelConstants.MAX_ARRAY_SIZE</code>.
	 * 
	 * @see KestrelConstants#MAX_ARRAY_SIZE
	 */
	public int nextCol()
			throws IllegalStateException {
		
		++column;
		
		if (column >= colCapacity) {
			
			try {
				expandMatrix();
				
			} catch (IllegalStateException ex) {
				--column;  // Reset column
				throw ex;
			}
		}
		
		return column;
	}
	
	/**
	 * Expand the number of columns in this matrix.
	 * 
	 * @throws IllegalStateException If the matrix already has
	 *   <code>KestrelConstants.MAX_ARRAY_SIZE</code> columns.
	 */
	private void expandMatrix()
			throws IllegalStateException {
		
		int newCapacity = colCapacity;
		
		newCapacity = (int) (colCapacity * KestrelConstants.ARRAY_EXPAND_FACTOR);
		
		if (newCapacity < 0) {
			
			if (colCapacity == KestrelConstants.MAX_ARRAY_SIZE)
				throw new IllegalStateException("Cannot expand matrix: Column capacity is at its maximum size: " + KestrelConstants.MAX_ARRAY_SIZE);
			
			colCapacity = KestrelConstants.MAX_ARRAY_SIZE;
		}
		
		short[][] newMatrix = new short[nRow][newCapacity];
		
		for (int row = 0; row < nRow; ++row)
			for (int col = 0; col < column; ++col)  // Note: column is already incremented
				newMatrix[row][col] = matrix[row][col];
		
		matrix = newMatrix;
		colCapacity = newCapacity;
		
		return;
	}
	
	/**
	 * Get a string representation of this matrix. This string will be large with many lines.
	 * <p/>
	 * If <code>refSequence</code> and <code>conSequence</code> are not <code>null</code> and
	 * are long enough to contain this alignment, taking <code>refStart</code> into account,
	 * then the reference and consensus sequences are written with the matrix string. 
	 * 
	 * @param refSequence Reference sequence or <code>null</code>.
	 * @param conSequence Consensus sequence or <code>null</code>.
	 * @param refStart Start position in the reference sequence.
	 * 
	 * @return A string representation of this matrix.
	 */
	public String toString(byte[] refSequence, byte[] conSequence, int refStart) {
		
		int nWidth;
		int width;
		
		int rowSpace;  // Number of spaces for row numbers
		
		boolean writeRefSeq;
		boolean writeConSeq;
		
		StringBuilder builder;
		
		// Check arguments
		if (refStart < 0)
			refStart = 0;
		
		writeRefSeq = true;
		writeConSeq = true;
		
		if (refSequence == null || refSequence.length < nRow + refStart)
			writeRefSeq = false;
		
		if (conSequence == null || conSequence.length <= column || column < 0)
			writeConSeq = false;
		
		// Init
		builder = new StringBuilder();
		
		rowSpace = (int) Math.log10(nRow - 1) + 1;
		String rowIndexFormat = "%" + rowSpace + "d";
		
		// Write consensus sequence
		if (writeConSeq) {
			
			if (writeRefSeq)
				builder.append("  ");
			
			for (int index = 0; index < rowSpace; ++index)
				builder.append(" ");  // One space for each row digit
			
			for (int index = 0; index <= column; ++index) {
				builder.append("      ");
				builder.append((char) conSequence[index]);
				builder.append("     ");
			}
			
			builder.append("\n");
		}
		
		// Column indices
		if (writeRefSeq)
			builder.append("  ");
		
		for (int index = 0; index < rowSpace; ++index)
			builder.append(" ");  // One space for each row digit
		
		for (int col = 0; col <= column; ++col) {
			width = 11;
			
			if (col > 9)
				nWidth = (int) Math.log10(col) + 1;
			else
				nWidth = 1;
			
			width = (11 - nWidth) / 2;
			
			builder.append(" ");
			
			for (int index = 0; index < width; ++index)
				builder.append(" ");
			
			builder.append(col);
			
			width = 11 - width - nWidth;
			
			for (int index = 0; index < width; ++index)
				builder.append(" ");
		}
		
		builder.append("\n\n");
		
		// Add rows
		for (int row = 0; row < nRow; ++row) {
			
			if (writeRefSeq) {
				builder.append((char) refSequence[refStart + row]);
				builder.append(" ");
			}
			
			builder.append(String.format(rowIndexFormat, row));
			
			for (int col = 0; col <= column; ++col) {
				builder.append(" ");
				
				builder.append(((matrix[row][col] >> 8) & 0x1) != 0 ? "a" : "-");
				builder.append(((matrix[row][col] >> 7) & 0x1) != 0 ? "r" : "-");
				builder.append(((matrix[row][col] >> 6) & 0x1) != 0 ? "c" : "-");
				builder.append(",");
				builder.append(((matrix[row][col] >> 5) & 0x1) != 0 ? "a" : "-");
				builder.append(((matrix[row][col] >> 4) & 0x1) != 0 ? "r" : "-");
				builder.append(((matrix[row][col] >> 3) & 0x1) != 0 ? "c" : "-");
				builder.append(",");
				builder.append(((matrix[row][col] >> 2) & 0x1) != 0 ? "a" : "-");
				builder.append(((matrix[row][col] >> 1) & 0x1) != 0 ? "r" : "-");
				builder.append(((matrix[row][col] >> 0) & 0x1) != 0 ? "c" : "-");
			}
			
			builder.append("\n");
		}
		
		return builder.toString();
	}
	
	/**
	 * Get a string representation of this matrix. This string will be large with many lines.
	 * 
	 * @return A string representation of this matrix.
	 */
	@Override
	public String toString() {
		return toString(null, null, 0);
	}
}
