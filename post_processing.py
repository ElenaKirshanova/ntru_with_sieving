import numpy as np


"""
  Inverts the projection map.
"""
def projected(v):
  v = np.array( list(v) )
  d = int( len(v)/2 )
  
  ones_left = np.array( [1]*d + [0]*d )
  ones_right = np.array( [0]*d + [1]*d )
  all_zero = np.array( 2*d*[0] )
  
  foundPreImage = False
  i = 0
  
  while not foundPreImage and i < d:
    j = 0
    while not foundPreImage and j < d:
      candidate = v + i*ones_left + j*ones_right
      
      if np.array_equal(candidate % d, all_zero):
        foundPreImage = True
        v = (candidate/d).astype(int).tolist()
      else:
        j += 1
    
    i += 1
  
  return v
	
def challenge(v, h, sqnormH):
	h = np.array( list(h) )
	
	v = np.array( list(v) )
	d = int( len(v)/2 )
	
	f = v[d:]
	f = (f / sqnormH).astype(int).tolist()
	
	g = v[:d]
	
	for i in range( d ):
		g_i = g[i]
		h_i = h[i]
		if g_i % h_i == 0:
			candidate_g = g - int(g_i/h_i)*h
			
			valid = True
			i = 0
			while valid and i < d:
				valid = candidate_g[i] % sqnormH == 0
				i += 1
			
			if valid:
				g = (candidate_g / sqnormH).astype(int).tolist()
				break
	
	v = g + f
	
	return v