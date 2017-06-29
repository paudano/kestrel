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

package edu.gatech.kestrel.io;

import java.io.File;
import java.io.FileDescriptor;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;

/**
 * Represents an object that can be opened as an <code>OutputStream</code> object. This
 * is used for output files and log files.
 */
public class StreamableOutput {
	
	/**
	 * File this output object wraps, or <code>null</code> if the output object is
	 * a file descriptor. Exactly one of <code>file</code> or <code>fd</code> is
	 * <code>null</code>.
	 */
	public final File file;
	
	/**
	 * File descriptor this output object wraps, or <code>null</code> if the output object is
	 * a file. Exactly one of <code>file</code> or <code>fd</code> is <code>null</code>.
	 */
	public final FileDescriptor fd;
	
	/**
	 * Name of the output object. This is the file name or a short string describing the file
	 * descriptor, if it is known. This name can be used in error messages and logs.
	 */
	public final String name;
	
	/** Standard output file descriptor. */
	public static final StreamableOutput STDOUT = new StreamableOutput(null, FileDescriptor.out, "<STDOUT>");
	
	/** Standard output file descriptor. */
	public static final StreamableOutput STDERR = new StreamableOutput(null, FileDescriptor.err, "<STDERR>");
	
	/**
	 * Get a streamable output object.
	 * 
	 * @param file File.
	 * @param fd File descriptor.
	 * @param name Name for error, warnings, and documentation.
	 */
	private StreamableOutput(File file, FileDescriptor fd, String name) {
		
		// Exactly one of file or fd must be set
		assert (file != null) || (fd != null) :
			"File and file descriptor are both null (exactly one must be null)";
		
		assert (file == null) || (fd == null) :
			"Neither file nor file descriptor null (exactly one must be null)";
		
		// Set fields
		this.file = file;
		this.fd = fd;
		this.name = name;
		
		return;
	}
	
	/**
	 * Get object from file name.
	 * 
	 * @param fileName File name. Must not be <code>null</code>.
	 * @param name Name of this output object. If <code>null</code> (recommended), then the name is
	 *   automatically assigned from <code>outputObject</code>.
	 * 
	 * @return <code>StreamableOutput<code> object.
	 * 
	 * @throws IllegalArgumentException If <code>fileName</code> is empty.
	 */
	protected static StreamableOutput getFromFileName(String fileName, String name)
			throws IllegalArgumentException {
		
		assert (fileName != null) :
			"File name is null";
		
		fileName = fileName.trim();
		
		if (fileName.isEmpty())
			throw new IllegalArgumentException("Cannot set output for an empty file name");
		
		File file = new File(fileName);
		
		if (name == null)
			name = file.getName();
		
		return new StreamableOutput(file, null, name);
	}
	
	/**
	 * Get an object from a file.
	 * 
	 * @param file File. Must not be <code>null</code>.
	 * @param name Name of this output object. If <code>null</code> (recommended), then the name is
	 *   automatically assigned from <code>outputObject</code>.
	 * 
	 * @return <code>StreamableOutput<code> object.
	 */
	protected static StreamableOutput getFromFile(File file, String name) {
		
		assert (file != null) :
			"File is null";
		
		if (name == null)
			name = file.getName();
		
		return new StreamableOutput(file, null, file.getName());
	}
	
	/**
	 * Get an object from a file descriptor.
	 * 
	 * @param fd File descriptor. Must not be <code>null</code>.
	 * @param name Name of this output object. If <code>null</code> (recommended), then the name is
	 *   automatically assigned from <code>outputObject</code>.
	 * 
	 * @return <code>StreamableOutput<code> object.
	 */
	protected static StreamableOutput getFromFileDescriptor(FileDescriptor fd, String name) {
		
		assert (fd != null) :
			"File descriptor is null";
		
		// Return standard descriptors
		if (fd == FileDescriptor.out)
			return STDOUT;
		
		if (fd == FileDescriptor.err)
			return STDERR;
		
		// Set default name
		if (name == null)
			name = "<UNKNOWN_FILE_DESCRIPTOR>";
		
		// Return object
		return new StreamableOutput(null, fd, name);
	}
	
	/**
	 * Get a streamable output object from a <code>String</code> (file name), <code>File</code>,
	 * <code>FileDescriptor</code>, or <code>StreamableOutput</code> object.
	 * 
	 * @param outputObject A <code>String</code> (file name), <code>File</code>,
	 *   or <code>FileDescriptor</code> object.
	 * @param name Name of this output object. If <code>null</code> (recommended), then the name is
	 *   automatically assigned from <code>outputObject</code>.
	 * 
	 * @return <code>StreamableOutput<code> object. If <code>outputObject</code> was already
	 *   a <code>StreamableOutput</code>, then it is returned without modification.
	 * 
	 * @throws NullPointerException If <code>outputObject</code> is <code>null</code>.
	 * @throws IllegalArgumentException If <code>outputObject</code> is an empty string or is not
	 *   a valid object type.
	 */
	public static StreamableOutput get(Object outputObject, String name)
			throws NullPointerException, IllegalArgumentException {
		
		if (outputObject == null)
			throw new NullPointerException("Cannot create output object from null");
		
		if (name != null) {
			name = name.trim();
			
			if (name.isEmpty())
				name = null;
		}
		
		// Assign if object is an acceptable class
		if (outputObject instanceof String)
			return getFromFileName((String) outputObject, name);  // throws IllegalArgumentException
			
		else if (outputObject instanceof File)
			return getFromFile((File) outputObject, name);
			
		else if (outputObject instanceof FileDescriptor)
			return getFromFileDescriptor((FileDescriptor) outputObject, name);
		
		else if (outputObject instanceof StreamableOutput)
			return (StreamableOutput) outputObject;

		else
			throw new IllegalArgumentException("Output object must be of class String, File, FileDescriptor, or StreamableOutput");
	}
	
	/**
	 * Get a streamable output object from a <code>String</code> (file name), <code>File</code>,
	 * <code>FileDescriptor</code>, or <code>StreamableOutput</code> object.
	 * 
	 * @param outputObject A <code>String</code> (file name), <code>File</code>,
	 *   or <code>FileDescriptor</code> object.
	 * 
	 * @return <code>StreamableOutput<code> object. If <code>outputObject</code> was already
	 *   a <code>StreamableOutput</code>, then it is returned without modification.
	 * 
	 * @throws NullPointerException If <code>outputObject</code> is <code>null</code>.
	 * @throws IllegalArgumentException If <code>outputObject</code> is an empty string.
	 */
	public static StreamableOutput get(Object outputObject)
			throws NullPointerException, IllegalArgumentException {
		
		return get(outputObject, null);
	}
	
	/**
	 * Get a file output stream for this object.
	 * 
	 * @return A new file output stream.
	 * 
	 * @throws FileNotFoundException If a file or directory cannot be found.
	 */
	public FileOutputStream getStream()
			throws FileNotFoundException {
		
		// Exactly one of file or fd must be set
		assert (file != null) || (fd != null) :
			"File and file descriptor are both null (exactly one must be null)";
		
		assert (file == null) || (fd == null) :
			"Neither file nor file descriptor null (exactly one must be null)";
		
		// Return stream
		if (file != null)
			return new FileOutputStream(file);
		
		return new FileOutputStream(fd);
	}
	
	/**
	 * Determine if output goes to the screen.
	 * 
	 * @return <code>true</code> if output is directed to STOUT or STDERR.
	 */
	public boolean isScreen() {
		return fd == FileDescriptor.out || fd == FileDescriptor.err;
	}
}
