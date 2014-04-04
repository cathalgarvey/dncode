import sys as s, itertools as i
t,u,i2b,o=b'T',b'U',lambda i:i.to_bytes(1,'big'),[s.stdout.buffer.write,lambda s:0]
A,e=s.argv,dict(zip((''.join(x).encode()for x in i.product(*('ACGT',)*4)),map(i2b,range(256))))
with open(A[2],'rb')as I:D,Rr=b''.join(I.read().strip().split()),lambda s:s.replace(t,u)
S,R,M,o,d=D[:-1],D[-1]&4,D[-1]&2,o[::-1]if A[1]=="d"else o,dict(zip(e.values(),e)).get
[o[0](e.get(q.replace(u,t)+(b'A'*(4-len(q)))))for q in(D[i: i+4] for i in range(0,len(D),4))]
o[0](i2b((len( [D[i: i+4] for i in range(0,len(D),4)][-1] )%4)|(4 if u in S else 0)))
o[1](b''.join(Rr(d(i2b(x)))if R else d(x)for x in S)[:-M if M else None]+b'\n')