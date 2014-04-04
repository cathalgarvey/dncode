import sys, itertools
i2b,outs=lambda i:i.to_bytes(1,'big'),[sys.stdout.buffer.write,lambda s:0]
enct=dict(zip((b''.join((i2b(y)for y in x))for x in itertools.product(*(b'ACGT',)*4)),range(256)))
dect,_=dict(((v,k) for k,v in enct.items())),outs.reverse()if sys.argv[1][0]=="d"else 0
with open(sys.argv[2],'rb')as I:inp=b''.join([L.strip()for L in I]).strip()
for q in [inp[i: i+4]for i in range(0,len(inp),4)]:outs[0](i2b(enct.get(q+(b'A'*(4-len(q))),0)))
else:outs[0](i2b(len(q)%4))
outs[1](b''.join((dect.get(x,b'')for x in inp[:-1]))[:-inp[-1]if inp[-1]else None]+b'\n')