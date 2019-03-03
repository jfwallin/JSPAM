import numpy as np

od6 = 1.0/6.0
od3 = 1.0/3.0

class Integrator:
 
  def __init__(self, force):
    self.force = force

    # temporary array used in integration 
    self.x = []

    # temporary array used in integration 
    self.xe = []

    # internal array to hold force,acceleration calculations 
    self.f = []

  # Initialize temporary storage arrays
  def initRKVar(self,n):
    self.x = np.zeros((n+1,6))
    self.xe = np.zeros((n+1,6))
    self.f = np.zeros((n+1,6))

    #self.x = [[0] * 6 for i in range(n+1)]
    #self.xe = [[0] * 6 for i in range(n+1)]
    #self.f = [[0] * 6 for i in range(n+1)]

    if self.force is not None:
      self.force.initVars(n+1)


  # Cleanup temporary storage.
  def deallocateRKVar(self):
    self.x = []
    self.xe = []
    self.f = []

    if self.force is not None:
      self.force.deallocateVars()

  # Using the provided ForceModel, compute the updated positions and velocities
  # by taking a step of size h.  These are computed for all particles.
  def rk4(self,x0, xout, h):
    n = len(x0)
    i=0

    # pre-compute products
    hod6 = h*od6
    h0p5 = h*0.5
    hod3 = h*od3
     
    # we unroll the second dimension of 6 elements to speed computation       
    np.copyto(self.x,x0)
          
    self.force.diffeq(self.x,self.f)

    self.xe = x0 + self.f*hod6
    self.x  = x0 + self.f*h0p5

    self.force.diffeq(self.x,self.f)

    self.xe = self.xe + self.f * hod3
    self.x  = x0 + self.f * h0p5

    self.force.diffeq(self.x,self.f)

    self.xe = self.xe + self.f * hod3
    self.x  = x0 + h * self.f

    self.force.diffeq(self.x,self.f)

    self.xe = self.xe + self.f * hod6
    np.copyto(xout,self.xe)

    # end rk4


  # Using the provided ForceModel, compute the updated position and velocity of
  # the secondary galaxy by taking a step of size h.
  def rk41(self,xx0,xxe,h):
    x = [0,0,0,0,0,0,0]
    f = [0,0,0,0,0,0,0]

    i=0
    n=6

    # pre-compute products
    hod6 = h*od6
    h0p5 = h*0.5
    hod3 = h*od3

    x=xx0[:]

    self.force.diffq1(x, f)

    xxe[0] = xx0[0] + f[0] * hod6
    x[0]   = xx0[0] + f[0] * h0p5
    xxe[1] = xx0[1] + f[1] * hod6
    x[1]   = xx0[1] + f[1] * h0p5
    xxe[2] = xx0[2] + f[2] * hod6
    x[2]   = xx0[2] + f[2] * h0p5
    xxe[3] = xx0[3] + f[3] * hod6
    x[3]   = xx0[3] + f[3] * h0p5
    xxe[4] = xx0[4] + f[4] * hod6
    x[4]   = xx0[4] + f[4] * h0p5
    xxe[5] = xx0[5] + f[5] * hod6
    x[5]   = xx0[5] + f[5] * h0p5

    self.force.diffq1(x, f)

    xxe[0] = xxe[0] + f[0] * hod3
    x[0]   = xx0[0] + f[0] * h0p5
    xxe[1] = xxe[1] + f[1] * hod3
    x[1]   = xx0[1] + f[1] * h0p5
    xxe[2] = xxe[2] + f[2] * hod3
    x[2]   = xx0[2] + f[2] * h0p5
    xxe[3] = xxe[3] + f[3] * hod3
    x[3]   = xx0[3] + f[3] * h0p5
    xxe[4] = xxe[4] + f[4] * hod3
    x[4]   = xx0[4] + f[4] * h0p5
    xxe[5] = xxe[5] + f[5] * hod3
    x[5]   = xx0[5] + f[5] * h0p5

    self.force.diffq1(x, f)

    xxe[0] = xxe[0] + f[0] * hod3
    x[0]   = xx0[0] + h*f[0]
    xxe[1] = xxe[1] + f[1] * hod3
    x[1]   = xx0[1] + h*f[1]
    xxe[2] = xxe[2] + f[2] * hod3
    x[2]   = xx0[2] + h*f[2]
    xxe[3] = xxe[3] + f[3] * hod3
    x[3]   = xx0[3] + h*f[3]
    xxe[4] = xxe[4] + f[4] * hod3
    x[4]   = xx0[4] + h*f[4]
    xxe[5] = xxe[5] + f[5] * hod3
    x[5]   = xx0[5] + h*f[5]

    self.force.diffq1(x, f)

    xxe[0] = xxe[0] + f[0] * hod6
    xxe[1] = xxe[1] + f[1] * hod6
    xxe[2] = xxe[2] + f[2] * hod6
    xxe[3] = xxe[3] + f[3] * hod6
    xxe[4] = xxe[4] + f[4] * hod6
    xxe[5] = xxe[5] + f[5] * hod6

    # // end rk41
