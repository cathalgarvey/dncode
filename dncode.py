#!/usr/bin/env python3
import itertools
import unittest
import binascii  # Prefer to store specified testcases as hex in this file.
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

# Create tables
encoding_table = {}
decoding_table = {}

# For every four-letter combination of ACGT (256 possibilities):
for n, n4 in enumerate((''.join(x) for x in itertools.product(*('ACGT',)*4))):
    n = n.to_bytes(1,'big')
    encoding_table[n4] = n
    decoding_table[n] = n4

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
    return dict(
    mask_len = encbyte & 3,
    rna      = True if encbyte & 1<<2 else False,
    gzipped  = True if encbyte & 1<<3 else False )

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
    output = []
    code = code.upper().replace("U","T")
    quadl = 0
    for quad in [ code[i: i+4] for i in range(0, len(code), 4) ]:
        # Pads anything less than 4n with "A"; this should always be
        # the last chunk anyway, so handling masking occurs in the "else" block
        quadl = len(quad)
        enc_dna = encoding_table[ quad + ('A'*(4- quadl)) ]
        output.append(enc_dna)
    else:
        mask = (4 - quadl) % 4 # Can only be 0-3, leaving unused byte-space..
        output.append( get_encbyte(mask, rna) )
    return b''.join(output)

def decode(enccode):
    encbyte = enccode[-1]
    encoding = decode_encbyte(encbyte)
    enccode = enccode[:-1]
    output = []
    for b in enccode:
        output.append( decoding_table[b.to_bytes(1,'big')] )
    # Make code equal to the original sequence minus the mask.
    # Must do "mask_len or None" because str[:-0] is the same as str[:0]: ''!
    code = ''.join(output[ : -encoding['mask_len'] or None ])
    if encoding['rna']:
        code = code.replace("T","U")
    return code

class TestEncoding(unittest.TestCase):
    @staticmethod
    def gen_rand_dna(length):
        return ''.join((random.choice("ACGT") for x in range(length)))
    @staticmethod
    def gen_rand_rna(length):
        return gen_rand_dna(length).replace("T","U")

    cases = dict( [ (k, binascii.unhexlify(v)) for k,v in {
        'ACGTCGTAGTGGTCGTGTAGTAGGCTTACGTGAGCGTTTCGGAATTCGTATGCTAGCTGTCGGGGATTATGCGTAC': 
            b'1b6cbadbb2ca7c6e26fda0f6ce727b6a8f39b100'
        }.items() ])

    randcases = [ ''.join((random.choice("ACGT") for x in range(i))) for i in range(0,1000,100) ]
    rnarandcases = [ x.replace("T","U") for x in randcases]

    def testEncoding(self):
        for case, encoded in self.cases.items():
            self.assertEqual(encoded, encode(case), 'Failed to encode test-case.')
            
    def testDecoding(self):
        for case, encoded in self.cases.items():
            self.assertEqual(case, decode(encoded), 'Failed to decode test-case.')

    def testRNAEncoding(self):
        for item in self.rnarandcases:
            encoded = encode(item,rna=True)
            self.assertEqual(item, decode(encoded))
      
    def testRoundTrip(self):
        for item in self.randcases:
            encoded = encode(item)
            self.assertEqual(item, decode(encoded))

def test():
    print("Running test-cases")
    unittest.main(argv=[__name__])

def bench():
    import gzip
    import time
    def rand_dna(length):
        return ''.join(random.choice("ACGT") for x in range(length))
    def cmpcmp(seq):
        slen = len(seq)
        print("Comparing gzip to dncode for seqlen {}".format(slen))
        gtime = time.time()
        gseq = gzip.compress(seq.encode())
        gtime = time.time() - gtime
        dtime = time.time()
        dseq = encode(seq)
        dtime = time.time() - dtime
        glen = len(gseq)
        dlen = len(dseq)
        print("Compression ratios: gzip={gz}; dncode={dn}".format(
            gz=len(gseq)/slen,dn=len(dseq)/slen))
        print("Time taken: gzip={}s, dncode={}s".format(gtime,dtime))

    for i in range(100, 10000000, 100000):
        cmpcmp(rand_dna(i))

if __name__ == "__main__":
    import argparse
    import sys
    P = argparse.ArgumentParser(description="")
    subP = P.add_subparsers(help="Sub commands")
    testP = subP.add_parser("test", help="Run test cases")
    testP.set_defaults(mode='test')
 
    benchP = subP.add_parser("bench", help="Run test cases")
    benchP.set_defaults(mode='bench')

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
    if A.mode == 'tests':
        test()
    elif A.mode == "bench":
        bench()
    elif A.mode == 'compress':
        title, seq = load_fasta(A.input.read().strip())
        print("Compressing:",title,file=sys.stderr)
        A.output.write(encode(seq))
    elif A.mode == 'decompress':
        A.output.write(decode(A.read().strip()))
