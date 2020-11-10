#------------------------------------------------------------------------------
#                             BEZIER CURVE
#------------------------------------------------------------------------------

# Load modules 
import numpy, math
				
def bezier(P, Np):
	
	# function that given P control points, degree is P-1, and Np is the curve
	# resolution, return the Bezier curve
	
	deg    = P.shape[0]-1	
	t      = numpy.linspace(0,1, Np)
	xb_sum = numpy.zeros((deg+1, Np)) 
	yb_sum = numpy.zeros((deg+1, Np))
	
	for i in range(0, deg+1):
		
		# Newton binomial coefficient (n k) = n!/(k!*(n-k)!)
		binom = math.factorial(deg)/(math.factorial(i)*math.factorial(deg-i)) 
		
		xb_sum[i,:] = binom*(1-t)**(deg-i)*t**i*P[i,0]
		yb_sum[i,:] = binom*(1-t)**(deg-i)*t**i*P[i,1]
		
	# Bezier curve
	xb = sum(xb_sum)
	yb = sum(yb_sum)
	
	return (xb, yb)
