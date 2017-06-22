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

package edu.gatech.kestrel.runner;

/**
 * Thrown when a reference reader cannot be found or initialized.
 */
public class ConfigurationErrorException extends Exception {
	
	/**
	 * Serial version ID.
	 */
	private static final long serialVersionUID = 1L;
	
	/** Error code associated with this exception. Codes are defined in <code>Constants</code>. */
	public final int errCode;
	
	/**
	 * Create a new reference reader initialization exception.
	 * 
	 * @param msg Message about why this exception was thrown.
	 * @param errCode Error code associated with this exception.
	 */
	public ConfigurationErrorException(String msg, int errCode) {
		super(msg);
		
		this.errCode = errCode;

		return;
	}

	/**
	 * Create a new reference reader initialization exception.
	 * 
	 * @param msg Message about why this exception was thrown.
	 * @param errCode Error code associated with this exception.
	 * @param cause Cause of this exception.
	 */
	public ConfigurationErrorException(String msg, int errCode, Throwable cause) {
		super(msg, cause);
		
		this.errCode = errCode;

		return;
	}
}
