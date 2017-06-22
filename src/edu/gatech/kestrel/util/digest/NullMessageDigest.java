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

package edu.gatech.kestrel.util.digest;

import java.security.MessageDigest;

/**
 * Takes the place of a message digest when one is not used.
 */
public class NullMessageDigest extends MessageDigest {
	
	/** A name for the null algorithm. */
	public static final String ALGORITHM = "NO_DIGEST";
	
	/**
	 * Create a new null digest.
	 */
	public NullMessageDigest() {
		super(ALGORITHM);
		
		return;
	}

	@Override
	protected byte[] engineDigest() {
		return new byte[] {0};
	}

	@Override
	protected void engineReset() {
		return;
	}

	@Override
	protected void engineUpdate(byte arg0) {
		return;
	}

	@Override
	protected void engineUpdate(byte[] arg0, int arg1, int arg2) {
		return;
	}
}