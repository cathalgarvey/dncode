# DNcode
## A 4x compression encoding for simple nucleotide sequences.
by Cathal Garvey, copyright 2014, released under the GNU AGPL: see COPYING.txt

### What
DNA sequence data can be quite large yet, in terms of computer-storage,
information-sparse. DNA has only four characters, for example, yet
compression algorithms are optimised for information-dense data like
human text or binary data and assume that text-only DNA sequences ought
to be treated as such.

DNcode converts plain DNA or RNA sequences into a binary format by mapping
four-base strings to single-byte values, adding a single byte to the end
to carry information such as sequence padding length and sequence type.

In benchmarks against Python's gzip library (which uses the GNU gzip
C library), DNcode is consistently faster, with a process time scalar with
the length of the input sequence rather than growing logarithmically.
It also compares well on compression, beating gzip on random DNA (which,
it must be said, is not representative or real-world use) in all cases
but producing significantly smaller compressions for smaller inputs.

DNcode started as a Gist consisting of increasingly incomprehensible
"terse python", which I did for amusement, but given its favourable
performance versus gzip I have rewritten and modularised it for further
use.

### How
For usage information, call `dncode_readable.py -h`.

### For What
Given the speed and excellent compression of small sequences, this would
be ideal for reducing the size of databases consisting of smaller elements,
such as barcodes.

I will shortly benchmark this on actual sequences; it is possible that
even in real-world use the compression remains competitive with gzip et al,
though the speed will certainly remain superior. This means if compression
ratio is not critical, but a 75% compression at low cost is desirable
(for example, with large databases where each entry is in regular use,
disincentivising the use of more advanced compression schemes), dncode may
be very useful.

More interesting would be the use of dncode as a pre-compression encoding
format for biological sequence data; that is, to encode DNA/RNA using dncode
prior to using Gzip or lzma. After increasing the information density of
the DNA through binary encoding, the ultimate compression ratio after gzip
may be superior to gzip alone, with (again) low additional overhead.

### Where From Here
* Clean up dncode.py a bit, maybe hack out tests to another file.
* Re-optimise for module usage and rearrange for setuptools/pip installation.
* Test and implement post-dncoding gzip, adding a "gzipped" flag to metabyte
* Optimise for further speed; use iterators where possible to enhance use on large sequence files.
