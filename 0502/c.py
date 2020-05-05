import itertools
import numpy as np

n,m,q = map(int,input().split())

a = np.array([],dtype = 'int64')
b = np.array([],dtype = 'int64')
c = np.array([],dtype = 'int64')
d = np.array([],dtype = 'int64')

for i in range(q):
    ai,bi,ci,di = map(int,input().split())
    a = np.append(a,ai)
    b = np.append(b,bi)
    c = np.append(c,ci)
    d = np.append(d,di)

nums = np.arange(1,m+1)
ans = 0
for balls in itertools.combinations_with_replacement(nums,n):
    tempans = 0
    for j in range(q):
        if(balls[b[j]-1] - balls[a[j]-1] == c[j]):
            #print(balls[b[j]-1],balls[a[j]-1],c[j])
            tempans = tempans + d[j]
    ans = max(ans,tempans)

print(ans)


