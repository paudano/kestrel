# Kestrel
Mapping-free variant caller for short-read Illumina data

Kestrel is a variant-caller for short-read Illumina data. It does not align sequence reads or perform a de novo assembly. Instead, it breaks reads into k-mers, which are short overlapping fragments of uniform length, and it counts the number of occurences of each unique k-mer in the sequence data. It also finds the ordered k-mers of a reference. With these two k-mer sets, Kestrel searches for patterns of variation. By using a novel local assembly method guided by the reference, it builds one or more sequences over the altered region and call variants from its alignment.

This work has been submitted to bioRxiv (BIORXIV/2017/153619) and Oxford Bioinformatics (BIOINF-2017-1097) and has not yet been peer-reviewed.

Check out releases for the latest version of Kestrel. The software can be run with the `kestrel` script or with `java -Xmx4G -jar kestrel.jar`. Use the `--help` (or `-h`) option for a list of available parameters.

Do not separate `kestrel.jar` from the other JAR files in its directory. You may create symbolic links in another direcotry to `kestrel.jar`. If the JAR file itself is moved, the other JAR files must be moved along with it. If the JAR files are separated, Kestrel will not be able to find libraries it depends on.

## Download a Release

The easiest way to get Kestrel is to download the [Latest Release](https://github.com/paudano/kestrel/releases/latest). Then, untar the downloaded file. There is nothing to compile or install.

`tar -zxvf kestrel-X.Y.Z-linux.tar.gz` where X, Y, and Z are the major, minor, and release version numbers.

All of Kestrel's dependencies are packaged in the directory, and there is nothing to build. Do not move `kestrel.jar` from the directory with the other JAR files. If it is separated form the other JAR files, Kestrel will not be able to find its dependenices. It can be symbolically linked into another directory without issues.

## Basic Usage

Java 1.7 or later is needed to run Kestrel.

From an un-tarred release, Kestrel can be run using the `kestrel` script or by running the JAR file directly with Java.

For example:
`./kestrel -r REFERENCE.fasta -o VARIANTS.vcf READS.fastq

See `./kestrel -h` for all of Kestrel's options. More documentation for these options will appear in the Kestrel manual, which is currently incomplete.

## Clone and Build

### Clone

Clone via SSH or HTTPS:

git clone git@github.com:paudano/kestrel.git

or

git clone https://github.com/paudano/kestrel.git

### Build

Kestrel's build system is coordinated by Apache Ant.

Build the current release:
`ant package`

A TAR/GZ file will appear in `dist`. This is the file that is uploaded as a release. To run a non-release version, the best way is to build the TAR and then un-tar it. This will place all of Kestrel's dependencies in their proper location within the directory.

For a complete list of build targets:
`ant -projecthelp`
