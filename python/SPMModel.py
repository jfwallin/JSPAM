import math
import numpy as np

from ForceModel import ForceModel

# This ForceModel represents a softened point-mass.

class SPMModel(ForceModel):

  def __init__(self):

    # temporary storage arrays dependent on n
    self.r22 = []
    self.r21 = []
    self.r1 = 0
    self.r2 = 0
    self.rn = 0
    self.a1 = 0
    self.a2 = 0
    self.a3 = 0

    # temporary storage independent of n
    self.xn = [0,0,0,0,0,0]
    self.xp = 0.0
    self.yp = 0.0
    self.zp = 0.0

    # force parameters
    self.params = None
    self.m1 = 0
    self.m2 = 0
    self.m3 = 0
    self.eps1 = 0
    self.eps2 = 0


  # Sets the simulation parameters to this model.
  def setParameters(self, p):
    self.params = p

    if p is not None:
      self.m1 = self.params.mass1
      self.m2 = self.params.mass2
      self.m3 = self.params.mass2
      self.eps1 = self.params.epsilon1 * self.params.epsilon1
      self.eps2 = self.params.epsilon2 * self.params.epsilon2


  # Initialize temporary storage.
  def initVars(self, n):
    self.r22 = np.zeros(n) 
    self.r21 = np.zeros(n)
    self.r1 = np.zeros(n)
    self.r2 = np.zeros(n)
    self.rn = np.zeros(n)
    self.a1 = np.zeros(n)
    self.a2 = np.zeros(n)
    self.a3 = np.zeros(n)


  # Cleanup temporary storage.
  def deallocateVars(self):
    self.r22 = None
    self.r21 = None
    self.r1 = None
    self.r2 = None
    self.rn = None
    self.a1 = None
    self.a2 = None
    self.a3 = None



  # Compute the circular velocity for a particle at a distance r from the specified mass.
  # The rout scale of the disk and softening length, eps, are provided.
  def circularVelocity(self,mass,r,rout,eps):

    ftotal = mass / ( r*r + eps )
    v = math.sqrt(ftotal*r)

    return v



  # For the given particle positions and velocities, calculate the accelerations.
  def diffeq(self,x,f):
  
    n = len(x)

    self.xn = x[n-1][:]
    #for i in range(6):
      #self.xn[i] = x[n-1][i]

    # make temps to handle perturber galaxy
    xp = self.xn[0]
    yp = self.xn[1]
    zp = self.xn[2]

    tmp4 = xp*xp+yp*yp+zp*zp   # r2n
    tmp5 = math.sqrt(tmp4)     # rn
    tmp4 = -self.m3/(tmp4+self.eps2)

    tmp1 = (x[:,0]-xp)
    tmp2 = (x[:,1]-yp)
    tmp3 = (x[:,2]-zp)
    self.r22 = tmp1*tmp1+tmp2*tmp2+tmp3*tmp3

    tmp1 = x[:,0]
    tmp2 = x[:,1]
    tmp3 = x[:,2]

    self.r21 = tmp1*tmp1+tmp2*tmp2 + tmp3*tmp3
    self.r2  = np.sqrt(self.r22)
    self.r1  = np.sqrt(self.r21)
    self.rn.fill(tmp5)

    self.a1 = -self.m1 / (self.r21 + self.eps1)
    self.a2 = -self.m2 / (self.r22 + self.eps2)
    self.a3.fill(tmp4)


    # this is a correction to prevent NaN errors in the vectorized
    # function evalution at the location of the second mass
    self.r2[n-1] = 1.0

    # calculate the RHS of the diffeq

        #// f(:,1) = x(:,4)
        #// f(:,2) = x(:,5)
        #// f(:,3) = x(:,6)

        #// f(:,4) = a1*x(:,1)/r1 + a2*(x(:,1)-xn(1))/r2 + a3*xn(1)/rn
        #// f(:,5) = a1*x(:,2)/r1 + a2*(x(:,2)-xn(2))/r2 + a3*xn(2)/rn
        #// f(:,6) = a1*x(:,3)/r1 + a2*(x(:,3)-xn(3))/r2 + a3*xn(3)/rn

    tmp1 = self.a1/self.r1
    tmp2 = self.a2/self.r2
    tmp3 = self.a3/self.rn

    f[:,0] = x[:,3]
    f[:,1] = x[:,4]
    f[:,2] = x[:,5]

    f[:,3] = tmp1*x[:,0] + tmp2*(x[:,0]-self.xn[0]) + tmp3*self.xn[0]
    f[:,4] = tmp1*x[:,1] + tmp2*(x[:,1]-self.xn[1]) + tmp3*self.xn[1]
    f[:,5] = tmp1*x[:,2] + tmp2*(x[:,2]-self.xn[2]) + tmp3*self.xn[2]


  #Calculate the acceleration felt by the secondary galaxy in orbit around the primary.
  def diffq1(self,x,f):
    r21 = 0
    r1 = 0
    a1 = 0
    a2 = 0

    r21 = x[0]*x[0]+x[1]*x[1]+x[2]*x[2]
    r1  = math.sqrt(r21)

    a1 = -self.m1 / (r21 + self.eps1)
    a2 = -self.m2 / (r21 + self.eps2)

    a1 = a1 + a2

    f[0] = x[3]
    f[1] = x[4]
    f[2] = x[5]
    f[3] = a1 * x[0] / r1
    f[4] = a1 * x[1] / r1
    f[5] = a1 * x[2] / r1
