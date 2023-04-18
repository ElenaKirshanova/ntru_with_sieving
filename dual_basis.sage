import json

key = sys.argv[1]

with open(key, "r+") as f:
  fileContent = f.read()
  basis = json.loads(fileContent)
  
  B = Matrix(basis)
  B = B.pseudoinverse().transpose()
  B = B.denominator()*B
  
  basis = []
  for v in B:
    basis.append( v.list() )
  
  f.seek(0)
  f.write( str(basis) )
  f.truncate()