# kestrel
Mapping-free variant caller for short-read Illumina data

Kestrel is a variant-caller for short-read Illumina data. It does not align sequence reads or perform a de novo assembly. Instead, it breaks reads into k-mers, which are short overlapping fragments of uniform length, and it counts the number of occurences of each unique k-mer in the sequence data. It also finds the ordered k-mers of a reference. With these two k-mer sets, Kestrel searches for patterns of variation. By using a novel local assembly method guided by the reference, it builds one or more sequences over the altered region and call variants from its alignment.

This work has been submitted to bioRxiv (BIORXIV/2017/153619) and Oxford Bioinformatics (BIOINF-2017-1097) and has not yet been peer-reviewed.

Check out releases for the latest version of Kestrel. The software can be run with the `kestrel` script or with `java -Xmx4G -jar kestrel.jar`. Use the `--help` (or `-h`) option for a list of available parameters.

Do not separate `kestrel.jar` from the other JAR files in its directory. You may create symbolic links in another direcotry to `kestrel.jar`. If the JAR file itself is moved, the other JAR files must be moved along with it. If the JAR files are separated, Kestrel will not be able to find libraries it depends on.
