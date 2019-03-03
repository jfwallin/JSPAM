class Parameters:
  hbase = 0.001
  def __init__(self):

    self.potential_type = 0
    self.ndim = 3

    self.mass1 = 0.0
    self.epsilon1 = 0.0
    self.rin1 = 0.0
    self.rout1 = 0.0
    self.rscale1 = [0,0,0]
    self.theta1 = 0.0
    self.phi1 = 0.0
    self.heat1 = 0.0
    self.opt1 = 0.0

    self.mass2 = 0.0
    self.epsilon2 = 0.0
    self.rin2 = 0.0
    self.rout2 = 0.0
    self.rscale2 = [0,0,0]
    self.theta2 = 0.0
    self.phi2 = 0.0
    self.heat2 = 0.0
    self.opt2 = 0.0

    self.eps1 = 0.0
    self.eps2 = 0.0

    self.tcurrent = 0.0


    self.n = 0
    self.n1 = 0
    self.n2 = 0

    self.time = 0.0
    self.tstart = 0.0
    self.tend = 0.0
    self.tIsSet = False
    self.inclination_degree = 0.0
    self.omega_degree = 0.0
    self.rmin = 0.0
    self.velocity_factor = 0.0
    self.mratio = 0.0
    self.secondary_size = 0.0

    self.sec_vec = [0,0,0,0,0,0]
    self.use_sec_vec = False

    self.h = 0.0

    self.nstep = 0
    self.nout = 0
    self.iout = 0
    self.unit = 0
    self.istep = 0

    self.fname= ""
    self.iostat = 0
    self.header_on = False

    # info about time of closest approach
    self.tmin = 0.0
    self.vrmin = 0.0
    self.beta = 0.0

    self.defaultParameters()



  def defaultParameters(self):
    self.standardGalaxy1()
    self.standardGalaxy2()
    self.testCollision()
    self.customCollision()


  def standardGalaxy1(self):
    # disk profile #1
    self.mass1 = 1.0
    self.epsilon1  = 0.3
    self.rin1 = 0.05
    self.rout1 = 1.0
    self.rscale1 = [3.0,3.0,3.0]
    self.theta1 = 0.0
    self.phi1 = 0.0
    self.opt1 = 1
    self.heat1 = 0.0
    self.eps1 = self.epsilon1*self.epsilon1


  def standardGalaxy2(self):
    # disk profile #2
    self.mass2 = 1.0
    self.epsilon2  = 0.3
    self.rin2 = 0.05
    self.rout2 = 1.0
    self.rscale2 = [3.0,3.0,3.0]
    self.theta2 = 0.0
    self.phi2 = 0.0
    self.opt2 = 1
    self.heat2 = 0.0
    self.eps2 = self.epsilon2*self.epsilon2


  def testCollision(self):
    # collision parameters
    self.inclination_degree = 90.0
    self.omega_degree = 0.0
    self.rmin = 1.0
    self.velocity_factor = 1.0
    self.time = -3.0

    # time step
    self.h = self.hbase
    self.nout = 5
    self.nstep = 500

    # particle numbers
    self.n1 = 1000
    self.n2 = 1000
    self.n = self.n1 + self.n2


  def customCollision(self):
    self.phi1 = 5.0
    self.theta1 = 5.0
    self.rscale1 = [1.0,1.0,1.0]
    self.rout1 = 1.0
    self.mass1 = 1.0
    self.epsilon1 = 0.3
    self.eps1 = self.epsilon1*self.epsilon1
    self.n1 = 1000
    self.heat1 = 0.0
    self.opt1 = 1

    self.phi2 = 0.0
    self.theta2 = 0.0
    self.rscale2 = [0.3,0.3,0.3]
    self.rout2 = 0.5
    self.mass2 = 0.5
    self.epsilon2 = 0.3
    self.eps2 = self.epsilon2*self.epsilon2
    self.n2 = 500
    self.heat2 = 0.0
    self.opt2 = 1

    self.inclination_degree = 20.0
    self.omega_degree = 0.0
    self.rmin = 0.9
    self.velocity_factor = 0.9

    self.h = self.hbase
    self.time = -5

    self.n = self.n1 + self.n2
