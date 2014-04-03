import sys, itertools
mode,output,code,i2b=sys.argv[1][0],sys.stdout.buffer,b'ACGT',lambda i:i.to_bytes(1,'big')
encode_table=dict(zip((b''.join((chr(y).encode()for y in x))for x in itertools.product(*(code,)*4)),range(256)))
decode_table=dict(((v,k) for k,v in encode_table.items()))
with open(sys.argv[2],'rb')as I:inp=b''.join([L.strip()for L in I]).strip()
if mode == "e":
  for q in [inp[i: i+4]for i in range(0,len(inp),4)]:output.write(i2b(encode_table[q+(b'A'*(4-len(q)))]))
  else:output.write(i2b(len(q)%4))
else:print(b''.join((decode_table[x]for x in inp[:-1]))[:-inp[-1]if inp[-1]else None].decode())