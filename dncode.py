#!/usr/bin/env python3
import os
import itertools
import random    # For generating random test cases.

def load_fasta(f, keep_non_acgtu = False):
    title = ''
    sequence = []
    for line in f.splitlines():
        line = line.strip()
        if line[0] == ">":
            title = line.lstrip(">").strip()
        elif line[0] == ";":
            continue
        else:
            sequence.append(line)
    sequence = ''.join(
        filter(None if not keep_non_acgtu else(lambda s:s in 'ACGTU'),itertools.chain(*sequence))
        )
    return title, sequence

def iter_fasta():
    'Should yield sequence data for first block (assume monolithic) and ignore title.'
    # Intended for large, monolithic sequence files.

# Create tables
acgt_encoding_table = {}
acgt_decoding_table = {}
# TODO:
iupac_encoding_table = {}
iupac_decoding_table = {}

# For every four-letter combination of ACGT (256 possibilities):
for n, n4 in enumerate((''.join(x) for x in itertools.product(*('ACGT',)*4))):
    nb = n.to_bytes(1,'big')
    acgt_encoding_table[n4] = nb
    acgt_decoding_table[nb] = n4

def get_encbyte(mask_len, rna=False, gzipped=False):
    '''The final byte of an encoded sequence is of form (binary):
    |unused|unused|unused|unused|gzipped|rna|mask|mask|
    ..that is, the length of the final n-quad mask is obtained using only
    the final two bits of the mask, and the third-last bit is 1 if RNA, 0 else.
    The remaining 5 bits are unspecified at present.
    '''
    if not 0 <= mask_len <= 3:
        raise ValueError("Mask length can only be 0-3.")
    mask  = 0
    mask |= mask_len
    if rna:
        mask |= 1<<2
    if gzipped:
        mask |= 1<<3    
    return mask.to_bytes(1, 'big')

def decode_encbyte(encbyte):
    'Decodes a byte as encoded by get_encbyte'
    #encbyte = int.from_byte(encbyte, 'big') # Be explicit on endianness
    if isinstance(encbyte, bytes):
        if not len(encbyte) == 1:
            raise ValueError("encbyte MAY be bytes, but must be length 1!")
        else:
            encbyte = int.from_bytes(encbyte, "big")
    if (not isinstance(encbyte, int)) or encbyte > 255:
        raise TypeError("decode_encbyte only accepts a single byte or a byte-size int. Value given was of type {} and value {}.".format(type(encbyte), encbyte))
    return dict(
        mask_len = encbyte & 3,
        rna      = True if encbyte & 1 << 2 else False,
        gzipped  = True if encbyte & 1 << 3 else False )

def encode(code, rna=False):
    '''Outputs a bytestring consisting of the encoded DNA, plus a mask byte.
    
    The encoding is static 4:1 compression, with quadrets of DNA translating to
    single-byte values evenly.
    The last encoded byte may be padded to permit encoding. To account for this
    a "mask" byte is always appended to the sequence, which consists simply of
    a big-endian integer specifying which of the last 4n are padding to be
    discarded on decoding. There is extra "room" in this last byte for other
    formatting details, such as alphabet, as a future improvement.'''
    # Iterate over the Dut DNA in four-letter chunks, maybe plus trailing.
    return b''.join(iter_encode(code, rna=rna))

def iter_encode(code_iter, table = acgt_encoding_table, rna=False):
    '''Iterates over sequences, preventing large sequence files being loaded into RAM.
    Yields encoded sequence bytes plus final encoding byte.'''
    # Create four iterators; one for each of four frames to iterate in steps of four.
    # Then zip these frames and use map to combine them into quadlets. Using
    # filter in the arguments to 'join' prevents a TypeError on last block if any
    # NoneType values are yielded from the zip_longest iterator.
    fr1, fr2, fr3, fr4 = [itertools.islice(code_iter, i, None, 4) for i in range(4)]
    framezip = itertools.zip_longest(fr1,fr2,fr3,fr4)
    zip_quads = map(lambda t:''.join(filter(None,t)), framezip)
    seql = 0
    for quad in zip_quads:
        try:
            enc_dna = acgt_encoding_table[quad]
            seql += 4
        except KeyError as E:
            # Should only evaluate on last quad, making all preceding lookups faster.
            # TODO: Make at least a token effort to ensure it's not a real KeyError.
            if len(quad) < 4:
                enc_dna = acgt_encoding_table[ quad + ( "A" * (4 - len(quad)) ) ]
                seql += len(quad)
            else:
                raise E
        yield enc_dna
    else:
        # Generate and yield final byte.
        mask = (4 - len(quad)) % 4 # Can only be 0-3, leaving unused byte-space..
        yield get_encbyte(mask, rna)

def iter_decode(enc_dna, encoding, table=acgt_decoding_table):
    'Returns an iterator to decode the input sequence or iterator over a sequence.'
    table = table.copy()
    if encoding['rna']:
        for k in table:
            table[k] = table[k].replace("T","U")
    for b in enc_dna:
        yield table[b]

def decode_all(enccode):
    'Straight decode of a sequence. Does not support iterators.'
    return b''.join(iter_decode(enccode[:-1], decode_encbyte(enccode[-1])))

def iter_decode_file(file_handle):
    'Returns an iterator from a file, after peeking at last byte for encoding.'
    fstat = os.stat(file_handle.name)
    # Get last byte for encoding.
    file_handle.seek(fstat.st_size - 1)
    enc_byte = file_handle.read()
    file_handle.seek(0)
    filesize_iterator = range( os.stat(file_handle.name).st_size - 1)
    seq_iterator = ( file_handle.read(1) for i in filesize_iterator )
    return iter_decode(seq_iterator, decode_encbyte(enc_byte))

if __name__ == "__main__":
    import argparse
    import sys
    P = argparse.ArgumentParser(description="")
    subP = P.add_subparsers(help="Sub commands")
 
    compressP = subP.add_parser("compress", help="Compress a sequence file")
    compressP.set_defaults(mode='compress')    
    compressP.add_argument("input", type=argparse.FileType("r"), default=sys.stdin,
        help="Input file to read; should be plain sequence data. Default stdin.")
    compressP.add_argument("-o","--output", type=argparse.FileType("wb"), default=sys.stdout.buffer,
        help="Output file to write compressed data to. Default stdout.")

    decompressP = subP.add_parser("decompress", help="Compress a sequence file")
    decompressP.set_defaults(mode='decompress')
    decompressP.add_argument("input", type=argparse.FileType("rb"), default=sys.stdin.buffer,
        help="Input file to read; should be plain sequence data. Default stdin.")
    decompressP.add_argument("-o","--output", type=argparse.FileType("w"), default=sys.stdout,
        help="Output file to write decompressed sequence to. Default stdout.")

    A = P.parse_args()
    # Check and run test
    if not hasattr(A, 'mode'):
        # Tried setting a default mode to catch modeless operation but somehow
        # that always ended up overwriting more specific modes. :s
        print("Error, no mode specified (mode was {}). For help, try -h".format(A.mode))
    
    if A.mode == 'compress':
        title, seq = load_fasta(A.input.read().strip())
        print("Compressing:",title,file=sys.stderr)
        #A.output.write(encode(seq))
        for b in iter_encode(seq):
            A.output.write(b)
    elif A.mode == 'decompress':
        for quad in iter_decode_file(A.input):
            A.output.write(quad)
        #A.output.write(decode(A.read().strip()))
