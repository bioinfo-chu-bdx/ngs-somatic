import random
import sys

def test(z):
  n=0
  for i in range(z):
    if problemathiques() == True :
      n+=1
  return n


def problemathiques(n=100, k=6):
    l = [random.randint(0, 1) for i in xrange(n)]
    return not any(all(l[((k-1)*i+j) % n] for j in xrange(k)) for i in xrange(n/(k-1)))


z = int(sys.argv[1])
count = test(z)
result = float(count)/float(z)

print result
