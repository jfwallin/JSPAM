import math

class MathUtil:
  X_AXIS = 1
  Y_AXIS = 2
  Z_AXIS = 3
    
  # Builds a rotation matrix for rotating about the
  #specified axis by the specified angle in radians.
  def buildRot(angle, axis):
    mat = [[0] * 3 for i in range(3)]
    cosa = math.cos(angle)
    sina = math.sin(angle)
        
    if axis == X_AXIS:
      mat[0][0]=1.0
      mat[1][1]=cosa
      mat[1][2]=sina
      mat[2][1]=-sina
      mat[2][2]=cosa
    elif axis == Y_AXIS:
      mat[0][0]=cosa
      mat[0][2]=-sina
      mat[1][1]=1.0
      mat[2][0]=sina
      mat[2][2]=cosa
    elif axis == Z_AXIS:
      mat[0][0]=cosa
      mat[0][1]=sina
      mat[1][0]=-sina
      mat[1][1]=cosa
      mat[2][2]=1.0
        
    return mat

    
  # Calculates cross product of two vectors.
  def cross(v1,v2):
    v=[0,0,0]
        
    v[0] = v1[1]*v2[2]-v1[2]*v2[1]
    v[1] = v1[2]*v2[0]-v1[0]*v2[2]
    v[2] = v1[0]*v2[1]-v1[1]*v2[0]
        
    return v
   
 
  # Calculates the dot product of two vectors.
  def dot(v1,v2):
    size = len(v1)
    cp = 0.0
        
    for i in range(size):
      cp += v1[i]*v2[i]
        
    return cp


    
  # Calculates the magnitude of the vector.
  def mag(v):
    return math.sqrt(dot(v,v))
   
 
  # Add two vectors.
  def add(v1, v2):
    size = len(v1)
    v = [0]*size
   
    for i in range(size):
      v[i] = v1[i]+v2[i]
        
    return v


  # Sub two vectors.
  def sub(v1, v2):
    size = len(v1)
    v = [0]*size
   
    for i in range(size):
      v[i] = v1[i]-v2[i]
        
    return v

    
  # Scale the vector by the specified multiplier.
  def scale(sc, v1):
    size = len(v1)
    v = [0]*size
   
    for i in range(size):
      v[i] = sc*v1[i]
        
    return v


  # Determine the angle of rotation between two vectors.
  def angleBetween(v1,v2):
    m1 = mag(v1)
    m2 = mag(v2)
    ang = dot(v1,v2)/(m1*m2)
    ang = math.acos(ang)

    return ang
   
 
  # Multiplies a 3x3 matrix by a 3 vector.
  # It is mildly optimized by unrolling loops.
  def mult3(m,v):
    b = [0,0,0]

    if len(v) < 3 :
      return v

    b[0] = m[0][0]*v[0]+m[0][1]*v[1]+m[0][2]*v[2]
    b[1] = m[1][0]*v[0]+m[1][1]*v[1]+m[1][2]*v[2]
    b[2] = m[2][0]*v[0]+m[2][1]*v[1]+m[2][2]*v[2]
        
    return b
   
 
  # Multiplies an N vector by an MxN matrix.
  def mult(m,v):
    col = len( m[0])
    row = len(m)
    
    b = [0]*row
    
    for i in range(row):
      for j in range(col):
        b[i]+=m[i][j]*v[j]
        
    return b
