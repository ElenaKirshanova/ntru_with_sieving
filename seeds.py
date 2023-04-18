from keygen import NTRUKeyGenerator
import sys
import json
from multiprocessing import Pool

def getSeeds(useHRSS,n,q):
  generator = NTRUKeyGenerator(useHRSS,n,q)

  success = False
  while not success:
    seed = generator.newSeed()
    success = True
    try:
      f,g,h = generator.getKey(seed)
    except ZeroDivisionError:
      success = False
  
  return seed
  

n = int(sys.argv[1])
q = int(sys.argv[2])
processes = int(sys.argv[3])
fileName = sys.argv[4]

HRSS = (q==0)

pool = Pool( processes = processes )
results = pool.starmap(getSeeds, processes*[(HRSS,n,q,)] )

with open("seeds/"+fileName+".txt", "w+") as f:
  result = {"n":n,"q":q,"HRSS":useHRSS,"seeds":results}
  print( json.dumps( result ), file=f )