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

package edu.gatech.kestrel.util;

import edu.gatech.kestrel.KestrelConstants;

/**
 * Output information about this software. This class was originally designed for
 * build utilities to get version numbers.
 */
public class InfoUtil {
	
	/**
	 * Valid information modes. 
	 */
	public enum InfoMode {
		
		/** Major version number. */
		VERSION_MAJOR("version.major"),
		
		/** Minor version number. */
		VERSION_MINOR("version.minor"),
		
		/** Development version number. If not <code>0<code> then this is not a release version. */
		VERSION_DEV("version.dev"),
		
		/** Full version string. */
		VERSION("version");
		
		/** User-friendly string version of this mode. */
		public final String modeName;
		
		/**
		 * Create mode.
		 * 
		 * @param modeName User-friendly mode name.
		 */
		private InfoMode(String modeName) {
			this.modeName = modeName;
			
			return;
		}
		
		/**
		 * Get the mode constant from a mode name.
		 * 
		 * @param modeName Mode name.
		 * 
		 * @return <code>InfoMode</code> constant or <code>null</code>
		 *   if <code>modeName</code> is <code>null</code> or not a valid
		 *   mode name.
		 */
		public static InfoMode getMode(String modeName) {
			if (modeName == null)
				return null;
			
			modeName = modeName.trim().toLowerCase();
			
			for (InfoMode mode : InfoMode.values()) {
				if (mode.modeName.equals(modeName))
					return mode;
			}
			
			return null;
		}
	}
	
	/**
	 * Hide the default constructor.
	 */
	private InfoUtil() {
		
		return;
	}
	
	/**
	 * Print information mode.
	 * 
	 * @param args A list of modes.
	 */
	public static void main(String[] args) {
		
		InfoMode mode;
		
		for (String arg : args) {
			
			arg = arg.trim().toLowerCase();
			
			mode = InfoMode.getMode(arg);
			
			if (mode == null) {
				System.err.println("Unrecognized mode: " + arg);
				System.exit(KestrelConstants.ERR_USAGE);
			}
			
			switch(mode) {
			case VERSION_MAJOR:
				System.out.println(KestrelConstants.VERSION_MAJOR);
				break;
				
			case VERSION_MINOR:
				System.out.println(KestrelConstants.VERSION_MINOR);
				break;
				
			case VERSION_DEV:
				System.out.println(KestrelConstants.VERSION_DEV);
				break;
				
			case VERSION:
				System.out.println(KestrelConstants.VERSION);
				break;
				
			default:
				System.err.println("Unknown output for mode: " + arg);
				System.exit(KestrelConstants.ERR_USAGE);
			}
		}
		
		return;
	}
}
