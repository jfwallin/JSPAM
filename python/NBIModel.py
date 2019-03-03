import math
import numpy as np

from ForceModel import ForceModel

# This ForceModel represents an n-body halo/disk/bulge potential that is sampled
# and then interpolated to calculate the force.
class NBIModel(ForceModel):

  # static values for the mass distribution
  isInit = False
  nnn = 10000
  nnm1 = nnn-1
  rad=[]
  rho_halo=[]
  mass_halo=[]
  rho_disk=[]
  mass_disk=[]
  rho_bulge=[]
  mass_bulge=[]
  rho_total=[]
  mass_total=[]
  masses=[]
  radius = []
  density = []
  vr2 = []
  vr = []
  new_vr2 =[]
  new_vr = []
  acceleration=[]
  acceleration_particle=[]
  new_mass=[]
  new_rho=[]
  phi=[]


  rs_internal = 10.0
  rs_internal2 = rs_internal*rs_internal
  rs_internal3 = rs_internal*rs_internal*rs_internal
  rmax_scale = 100.0
  sqrtpi = math.sqrt(math.pi)
  pscale = 0
  tmpdfind = 0
  lnl = 0

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
    self.vxp = 0.0
    self.vyp = 0.0
    self.vzp = 0.0

    # force parameters
    self.params = None
    self.m1 = 0
    self.m2 = 0
    self.m3 = 0
    self.eps1 = 0
    self.eps2 = 0
    self.mm1rs2 = 0
    self.mm2rs2 = 0
    self.mm3rs2 = 0
    self.pn = 0
    self.pn1 = 0
    self.pn2 = 0
  
    # indices and df  
    self.ival11=[]
    self.ival22=[]
    self.ivaln=[]
    self.df_force11=[]
    self.df_force22=[]
    self.df_forcen=[]
    self.c3n=[] 
    self.rrout1 = 0
    self.rrout2 = 0



  # Sets the simulation parameters to this model.
  def setParameters(self, p):
    self.params = p

    if p is not None:
      self.m1 = self.params.mass1
      self.m2 = self.params.mass2
      self.m3 = self.params.mass2
      self.eps1 = self.params.epsilon1 * self.params.epsilon1
      self.eps2 = self.params.epsilon2 * self.params.epsilon2

      self.pn = self.params.n
      self.pn1 = self.params.n1
      self.pn2 = self.params.n2
      self.mm1rs2 = -self.m1*NBIModel.rs_internal2;
      self.mm2rs2 = -self.m2*NBIModel.rs_internal2;
      self.mm3rs2 = -self.m3*NBIModel.rs_internal2;
      self.rrout1 = self.params.rout1;
      self.rrout2 = self.params.rout2;


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

    self.ival11 = np.zeros(n,dtype=np.int8) 
    self.ival22 = np.zeros(n,dtype=np.int8) 
    self.ivaln = np.zeros(n,dtype=np.int8) 

    self.df_force11 = np.zeros(n)
    self.df_force22 = np.zeros(n)
    self.df_forcen = np.zeros(n)

    self.c3n = np.zeros(n)

    self.initDistribution()

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

    self.ival11 = None
    self.ival22 = None
    self.ivaln = None

    self.df_force11 = None
    self.df_force22 = None
    self.df_forcen = None

    self.c3n = None


  # Build an n-body halo/disk/bulge
  def initDistribution(self):
    if NBIModel.isInit:
      return

    NBIModel.isInit = True

    rmax=0
    mold=0
    dmold=0
    mtotal=0
    mtot=0
    rscale=0
    dx=0
    x=0
    alphahalo=0
    qhalo=0
    gammahalo=0
    mhalo=0
    rchalo=0
    rhalo=0
    epsilon_halo=0
    zdisk=0
    zdiskmass=0
    hdisk=0
    zdiskmax=0
    hbulge=0
    mbulge=0
    rho_tmp=0
    G=0
    factor=0
    r=0
    m=0
    pi=0
    p1=0 
    rd=0
    rho_local=0
    p=0
    rr=0
    dr=0
    rh=0
    dp=0
    mnew=0
    dm=0
    acc_merge=0
    rad_merge=0
    acc=0
    xmax=0

    j=0
    nmax=0
    k=0
    nmerge=0
    ntotal=0
    jj=0

    NBIModel.rad = np.zeros(NBIModel.nnn)
    NBIModel.rho_halo = np.zeros(NBIModel.nnn)
    NBIModel.mass_halo = np.zeros(NBIModel.nnn)
    NBIModel.rho_disk = np.zeros(NBIModel.nnn)
    NBIModel.mass_disk = np.zeros(NBIModel.nnn)
    NBIModel.rho_bulge = np.zeros(NBIModel.nnn)
    NBIModel.mass_bulge = np.zeros(NBIModel.nnn)
    NBIModel.rho_total = np.zeros(NBIModel.nnn)
    NBIModel.mass_total = np.zeros(NBIModel.nnn)
    NBIModel.masses = np.zeros(NBIModel.nnn)
    NBIModel.radius = np.zeros(NBIModel.nnn)
    NBIModel.density = np.zeros(NBIModel.nnn)
    NBIModel.vr2 = np.zeros(NBIModel.nnn)
    NBIModel.vr = np.zeros(NBIModel.nnn)
    NBIModel.new_vr2 = np.zeros(NBIModel.nnn)
    NBIModel.new_vr = np.zeros(NBIModel.nnn)
    NBIModel.acceleration = np.zeros(NBIModel.nnn)
    NBIModel.acceleration_particle = np.zeros(NBIModel.nnn)
    NBIModel.new_mass = np.zeros(NBIModel.nnn)
    NBIModel.new_rho = np.zeros(NBIModel.nnn)
    NBIModel.phi = np.zeros(NBIModel.nnn)

    # set the constant for dynamical friction
    NBIModel.lnl = 0.00

    # set up the parameters for the halo
    mhalo = 5.8
    rhalo = 10.0
    rchalo = 10.0
    gammahalo = 1.0
    epsilon_halo = 0.4
    pi = math.pi

    # derive additional constants
    qhalo = gammahalo / rchalo
    alphahalo = 1.0 / ( 1.0 - NBIModel.sqrtpi * qhalo * math.exp(qhalo*qhalo) * (1.0 - NBIModel.erf(qhalo)) )

    # set the integration limits and zero integration constants
    rmax = 20
    nmax = 2000
    dr = rmax / nmax
    mold = 0

    rscale = 5
    ntotal = NBIModel.nnn

    # set the limits for integration, and zero integration constants
    k = nmax / 2
    dx = 1.0  / k
    x = 0.0
    dmold = 0.0
    mtot = 0.0
    m = 0.0
    G = 1

    # set the fundamental disk parameters
    zdisk = 0.2
    zdiskmax = 3.5
    hdisk = 1.0


    # set the fundamental bulge parameters
    hbulge = 0.2
    mbulge = 0.3333


    #set up the radius array
    for j in range(nmax):
      x = x + dx
      NBIModel.rad[j]= x * rchalo


    dr = NBIModel.rad[1] - NBIModel.rad[0]
    dx = dr / rchalo

    for j in range(nmax):
      # set the position
      r = NBIModel.rad[j]
      x = r / rchalo

      # calculate the local rho based
      rho_tmp = alphahalo / (2*NBIModel.sqrtpi*NBIModel.sqrtpi*NBIModel.sqrtpi ) * (math.exp(-x*x) / (x*x + qhalo*qhalo))

      # renormalize for the new halo size
      rho_tmp = rho_tmp / ( rchalo * rchalo * rchalo)

      # calculate mass in local shell, and update total mass
      dm = rho_tmp * 4 * pi * r * r *dr
      mtot = mtot + dm

      # store values in an array
      NBIModel.rho_halo[j] = rho_tmp * mhalo
      NBIModel.mass_halo[j] = mtot * mhalo
    #}

    # now calculate the potential
    NBIModel.phi = -G*NBIModel.mass_halo/NBIModel.rad

    #for j in range(nmax):
    #  r = NBIModel.rad[j];
    #  m = mass_halo[j];
    #  p1 = -G * m / r;
    #  phi[j] = p1;


    # disk model

    # loop over the distribution
    for j in range(nmax):

      # set the radius
      rd = NBIModel.rad[j]

      # find the local density in the disk
      rho_local  = math.exp(-rd/hdisk)/ (8*pi*hdisk*hdisk)
      NBIModel.rho_disk[j] = rho_local

      # find the mass in the spherical shell
      mnew = 4 * pi * rho_local * rd *rd * dr

      NBIModel.mass_disk[j] = mnew + mold
      mold = NBIModel.mass_disk[j]
    #}


    # bulge model

    # loop over the distribution
    mold = 0.0
    for j in range(nmax):
      # set the radius
      rd = NBIModel.rad[j]

      # find the local density in the disk
      rho_local  = math.exp(-rd*rd/(hbulge*hbulge))
      NBIModel.rho_bulge[j] = rho_local

      # find the mass in the spherical shell
      mnew = 4 * pi * rho_local * rd *rd * dr

      NBIModel.mass_bulge[j] = mnew + mold
      mold = NBIModel.mass_bulge[j]



    # renormalize distribution
    factor = mbulge / NBIModel.mass_bulge[nmax-1]
    NBIModel.mass_bulge = np.multiply(NBIModel.mass_bulge,factor)
    NBIModel.rho_bulge = np.multiply(NBIModel.rho_bulge,factor)

    #for j in range(nmax):
    #  mass_bulge[j] = mass_bulge[j] * factor;
    #  rho_bulge[j]  = rho_bulge[j]  * factor;


    dr = NBIModel.rad[1] - NBIModel.rad[0];



    j = 0
    NBIModel.mass_total[j]=  (NBIModel.mass_halo[j] + NBIModel.mass_disk[j] + NBIModel.mass_bulge[j])
    r = NBIModel.rad[j]
    NBIModel.rho_total[j] = NBIModel.mass_total[j] /  (4.0/3.0 * pi * r * r * r)
    dr = NBIModel.rad[1] - NBIModel.rad[0]

    for j in range(1,nmax):
      r = NBIModel.rad[j]
      NBIModel.mass_total[j]=  (NBIModel.mass_halo[j] + NBIModel.mass_disk[j] + NBIModel.mass_bulge[j])

      dm = NBIModel.mass_total[j] - NBIModel.mass_total[j-1]
      NBIModel.rho_total[j] = dm / (4 * pi * r * r * dr)



    # find the velocity dispersion v_r**2

    NBIModel.masses = NBIModel.mass_total;
    NBIModel.radius = NBIModel.rad;
    NBIModel.density = NBIModel.rho_total;


    for j in range(nmax):
      p = 0.0
      rr = NBIModel.radius[j]
      dr = NBIModel.radius[nmax-1] / nmax

      for jj in range(j,nmax):
        m  = NBIModel.masses[jj]
        rh = NBIModel.density[jj]
        rr = rr + dr

        dp = rh * G * m / (rr*rr) * dr
        p = p + dp

      NBIModel.vr2[j] = 1/NBIModel.density[j] * p
      NBIModel.vr[j] = math.sqrt(NBIModel.vr2[j])




    # find the accelerations felt by the particles and center of mass
    NBIModel.masses = NBIModel.mass_total;
    NBIModel.radius = NBIModel.rad;
    NBIModel.density = NBIModel.rho_total;

    for j in range(nmax):
      rr = NBIModel.radius[j]
      m  = NBIModel.masses[j]
      NBIModel.acceleration[j] = G * m / (rr*rr)
      NBIModel.acceleration_particle[j] = NBIModel.acceleration[j]

    nmerge = 50
    acc_merge = NBIModel.acceleration[nmerge]
    rad_merge = NBIModel.rad[nmerge]

    for j in range(nmax):
      rr = NBIModel.radius[j]
      m  = NBIModel.masses[j]

      # smoothed acceleration
      acc = G * m / (rr*rr + 0.1* (rad_merge -rr))
      NBIModel.acceleration_particle[j] = acc


    # rederive the masses from the new particle acceleration
    NBIModel.radius = NBIModel.rad
    dr = NBIModel.rad[1] - NBIModel.rad[0]

    # find the accelerations felt by the particles and center of mass
    NBIModel.radius = NBIModel.rad

    for j in range(nmax):
      rr = NBIModel.radius[j]
      NBIModel.new_mass[j] = rr*rr * NBIModel.acceleration_particle[j] / G
      NBIModel.new_rho[j]  = NBIModel.new_mass[j] / (4 * pi * rr * rr * dr)


    # find the velocity dispersion v_r**2 using the new density and masses

    NBIModel.masses = NBIModel.new_mass;
    NBIModel.radius = NBIModel.rad;
    NBIModel.density = NBIModel.new_rho;


    for j in range(nmax):
      p = 0.0
      rr = NBIModel.radius[j]
      dr = NBIModel.radius[nmax-1] / nmax

      for jj in range(j,nmax):
        m  = NBIModel.masses[jj]
        rh = NBIModel.density[jj]
        rr = rr + dr
        dp = rh * G * m / (rr*rr) * dr
        p = p + dp

      NBIModel.new_vr2[j] = 1/NBIModel.density[j] * p;
      NBIModel.new_vr[j] = math.sqrt(NBIModel.new_vr2[j])


    # extend the values to large rmax
    for j in range(nmax,ntotal):
      NBIModel.mass_total[j] = NBIModel.mass_total[nmax-1]
      NBIModel.mass_halo[j] = NBIModel.mass_halo[nmax-1]
      NBIModel.mass_disk[j] = NBIModel.mass_disk[nmax-1]
      NBIModel.mass_bulge[j] = NBIModel.mass_bulge[nmax-1]
      NBIModel.new_mass[j] = NBIModel.new_mass[nmax-1]

      NBIModel.rho_total[j] = 0.0
      NBIModel.new_rho[j]   = 0.0

      NBIModel.vr[j]      = 1e-6
      NBIModel.vr2[j]     = 1e-6
      NBIModel.new_vr[j]  = 1e-6
      NBIModel.new_vr2[j] = 1e-6

      m = NBIModel.mass_total[nmax-1]
      rr = NBIModel.rad[nmax-1] + dr*(j - nmax)
      NBIModel.rad[j] = rr
      acc = G * m / (rr*rr)
      NBIModel.acceleration_particle[j] = acc
      NBIModel.acceleration[j] = acc



    # normalize to the unit mass
    for j in range(ntotal):
      NBIModel.mass_total[j]  = NBIModel.mass_total[j] / 7.13333
      NBIModel.mass_halo[j]   = NBIModel.mass_halo[j]  / 7.13333
      NBIModel.mass_disk[j]   = NBIModel.mass_disk[j]  / 7.13333
      NBIModel.mass_bulge[j]  = NBIModel.mass_bulge[j] / 7.13333
      NBIModel.new_mass[j]    = NBIModel.new_mass[j]   / 7.13333

      NBIModel.rho_total[j]   = NBIModel.rho_total[j]  / 7.13333
      NBIModel.new_rho[j]     = NBIModel.new_rho[j]    / 7.13333

      NBIModel.vr[j]          = NBIModel.vr[j]      / 7.13333
      NBIModel.vr2[j]         = NBIModel.vr2[j]     / 7.13333
      NBIModel.new_vr[j]      = NBIModel.new_vr[j]  / 7.13333
      NBIModel.new_vr2[j]     = NBIModel.new_vr2[j] / 7.13333

      NBIModel.rad[j]         = NBIModel.rad[j];

      NBIModel.acceleration_particle[j] = NBIModel.acceleration_particle[j] / 7.13333
      NBIModel.acceleration[j]          = NBIModel.acceleration[j]  / 7.13333

    NBIModel.pscale = 1.0

    NBIModel.tmpdfind = NBIModel.pscale*NBIModel.rs_internal/NBIModel.rmax_scale*NBIModel.nnn;
  # end initDistribution





  # Determine the index to use for interpolating the force.
  def dfIndex(self,rin,rs):
    ival = 0
    tmp = rin*NBIModel.tmpdfind
    ival = int(tmp)
    if ival > NBIModel.nnm1:
      ival = NBIModel.nnm1

    return ival


  # Compute the circular velocity for a particle at a distance r from the specified mass.
  # The rout scale of the disk and softening length, eps, are provided.
  def circularVelocity(self, mass, r, rout, eps):
    ftotal = 0

    ival = self.dfIndex(r, rout)

    ftotal = mass * NBIModel.acceleration_particle[ival] * NBIModel.rs_internal * NBIModel.rs_internal

    v = math.sqrt(ftotal*r)

    return v



  # For the given particle positions and velocities, calculate the accelerations.
  def diffeq(self,x,f):

    tmp1 = 0
    tmp2 = 0
    tmp3 = 0
    tmp4 = 0
    tmp5 = 0
    n = len(x)
    i = 0;

    df_sigma = 0;
    df_rho = 0;
    c1 = 0;
    c2 = 0;
    xvalue = 0;
    v1 = 0;
    v21 = 0;

    self.xn = x[n-1][:]
    #for i in range(6):
      #self.xn[i] = x[n-1][i]

    # make temps to handle perturber galaxy
    xp = self.xn[0]
    yp = self.xn[1]
    zp = self.xn[2]
    vxp = self.xn[3]
    vyp = self.xn[4]
    vzp = self.xn[5]

    tmp4 = xp*xp+yp*yp+zp*zp   # r2n
    tmp5 = math.sqrt(tmp4)     # rn


    iv1 = 0
    iv2 = 0
    ivn = 0
    df1 = 0
    df2 = 0
    dfn = 0
    tx = None
    r22d = 0
    r21d = 0
    aa1 = 0
    aa2 = 0 
    aa3 = 0


    # get the forces, sigma and rho, and rescale them
    ivn = self.dfIndex(tmp5,self.rrout2);

    df_sigma  =  NBIModel.new_vr2[ivn] * NBIModel.rs_internal2;
    df_rho    =  NBIModel.new_rho[ivn] * NBIModel.rs_internal3;

    # df
    v21 = vxp*vxp + vyp*vyp + vzp*vzp;
    v1  = math.sqrt(v21);

    xvalue = v1 / df_sigma;
    c1 = NBIModel.erf(xvalue) - 2.0 * xvalue / NBIModel.sqrtpi * math.exp(-xvalue*xvalue);

    # df formula with G=1
    c2  = 4.0 * math.pi * self.m2 * NBIModel.lnl / v21;
    self.c3n.fill(0)

    for i in range(self.pn1,n):
      self.c3n[i] = c1*c2*df_rho


    tf = None
    tv1 = 1.0/v1
    c3tv = 0
    dv = 0

    ivn  = self.dfIndex(tmp5, self.rrout2)

    for i in range(n):
      # distance between the companion and the particle
      tx = x[i];

      tmp1 = (tx[0]-xp);
      tmp2 = (tx[1]-yp);
      tmp3 = (tx[2]-zp);
      r22d = tmp1*tmp1 + tmp2*tmp2 + tmp3*tmp3;

      # distance between the main galaxy and the particle
      tmp1 = tx[0];
      tmp2 = tx[1];
      tmp3 = tx[2];
      r21d = tmp1*tmp1 + tmp2*tmp2 + tmp3*tmp3;

      r22d = math.sqrt(r22d);
      r21d = math.sqrt(r21d);

      dv = r21d*NBIModel.tmpdfind;
      iv1 = int(dv)
      if(iv1>NBIModel.nnm1):
        iv1=NBIModel.nnm1;

      dv = r22d*NBIModel.tmpdfind;
      iv2 = int(dv)
      if(iv2>NBIModel.nnm1):
        iv2=NBIModel.nnm1;

      aa1 = NBIModel.acceleration_particle[iv1]*self.mm1rs2;
      aa2 = NBIModel.acceleration_particle[iv2]*self.mm2rs2;
      aa3 = NBIModel.acceleration_particle[ivn]*self.mm3rs2;

      tmp1 = aa1/r21d;
      tmp2 = aa2/r22d;
      tmp3 = aa3/tmp5;
      c3tv = self.c3n[i]*tv1;


      tf = f[i];
      tf[0] = tx[3];
      tf[1] = tx[4];
      tf[2] = tx[5];

      tf[3] = tmp1*tx[0] + tmp2*(tx[0]-self.xn[0]) + tmp3*self.xn[0] - self.xn[3] * c3tv;
      tf[4] = tmp1*tx[1] + tmp2*(tx[1]-self.xn[1]) + tmp3*self.xn[1] - self.xn[4] * c3tv;
      tf[5] = tmp1*tx[2] + tmp2*(tx[2]-self.xn[2]) + tmp3*self.xn[2] - self.xn[5] * c3tv;

    # end loop

    #redo acceleration for last particle
    r22d=1.0;
    tmp2 = aa2/r22d;
    tf[3] = tmp1*tx[0] + tmp2*(tx[0]-self.xn[0]) + tmp3*self.xn[0] - self.xn[3] * c3tv;
    tf[4] = tmp1*tx[1] + tmp2*(tx[1]-self.xn[1]) + tmp3*self.xn[1] - self.xn[4] * c3tv;
    tf[5] = tmp1*tx[2] + tmp2*(tx[2]-self.xn[2]) + tmp3*self.xn[2] - self.xn[5] * c3tv;


  # this is an attempt at a pythonic implementation of arrays
  def diffeqpy(self,x,f):
    n = len(x) 
    xn = x[n-1,:]

    # distance between the main galaxy and the particle
    r21 = x[:,0]*x[:,0]  + x[:,1]*x[:,1] + x[:,2]*x[:,2]
    r1 = np.sqrt(r21)

    # distance between the companion and the particle
    r22 = (x[:,0]-xn[0])*(x[:,0]-xn[0])+ (x[:,1]-xn[1])*(x[:,1]-xn[1])+ (x[:,2]-xn[2])*(x[:,2]-xn[2])
    r2 = np.sqrt(r22)


    # distance between the two galaxies - the tidal force
    r2n = xn[0]*xn[0]+xn[1]*xn[1]+xn[2]*xn[2]
    rn = math.sqrt(r2n)
   
    for i in range(n):
      self.ival11[i] =  self.dfIndex(r1[i], self.rrout1)
      self.ival22[i] =  self.dfIndex(r2[i], self.rrout2)
      self.ivaln[i]  =  self.dfIndex(rn, self.rrout2)


    df_force11 = NBIModel.acceleration_particle[self.ival11] * NBIModel.rs_internal2
    df_force22 = NBIModel.acceleration_particle[self.ival22] * NBIModel.rs_internal2
    df_forcen  = NBIModel.acceleration_particle[self.ivaln]  * NBIModel.rs_internal2



    # get the forces, sigma and rho, and rescale them
    df_sigma  =  NBIModel.new_vr2[self.ivaln[0]] * NBIModel.rs_internal2
    df_rho    =  NBIModel.new_rho[self.ivaln[0]] * NBIModel.rs_internal3


    # interpolated forces
    a1 = -self.m1 * df_force11
    a2 = -self.m2 * df_force22
    a3 = -self.m3 * df_forcen

    # df
    v21 = xn[3]*xn[3]+xn[4]*xn[4]+xn[5]*xn[5]
    v1  = math.sqrt(v21)


    xvalue = v1 / df_sigma
    c1 = self.erf(xvalue) - 2.0 * xvalue / NBIModel.sqrtpi * math.exp(-xvalue*xvalue)

    # df formula with G=1
    c2  = 4.0 * math.pi * self.m2 * NBIModel.lnl / v21
    self.c3n[0:n-1] = 0.0
    self.c3n[self.pn1+1:n]  = c1 * c2 * df_rho

    # this is a correction to prevent NaN errors in the vectorized
    # function evalution at the location of the second mass
    r2[n-1] = 1.0


    # calculate the RHS of the diffeq

    f[:,0] = x[:,3]
    f[:,1] = x[:,4]
    f[:,2] = x[:,5]

    f[:,3] = a1*x[:,0]/r1 + a2*(x[:,0]-xn[0])/r2 + a3*xn[0]/rn - self.c3n * xn[3]/ v1
    f[:,4] = a1*x[:,1]/r1 + a2*(x[:,1]-xn[1])/r2 + a3*xn[1]/rn - self.c3n * xn[4]/ v1
    f[:,5] = a1*x[:,2]/r1 + a2*(x[:,2]-xn[2])/r2 + a3*xn[2]/rn - self.c3n * xn[5]/ v1







  # Calculate the acceleration felt by the secondary galaxy in orbit around the primary.
  def diffq1(self,x,f):

    r21 = 0
    r1 = 0
    a1 = 0
    a2 = 0
    at = 0
    c1 = 0
    c2 = 0
    c3 = 0
    v21 = 0
    v1 = 0
    xvalue = 0

    ival2 = 0

    df_force1 = 0
    df_force2 = 0
    df_sigma = 0
    df_rho = 0

    i = 0
    rr = 0
    rlocal = 0
    ivalue = 0
    dr = 0
    mmm = 0
    dm = 0

    r21 = x[0]*x[0] + x[1]*x[1] + x[2]*x[2];
    r1  = math.sqrt(r21);

    # get the index for the calculations
    ival  =  self.dfIndex(r1, self.rrout1   );
    ival2 =  self.dfIndex(r1, self.rrout2   );

    # get the forces, sigma and rho, and rescale them
    df_force1 = NBIModel.acceleration_particle[ival] * NBIModel.rs_internal * NBIModel.rs_internal;
    df_force2 = NBIModel.acceleration_particle[ival2]* NBIModel.rs_internal * NBIModel.rs_internal;

    df_sigma  = NBIModel.new_vr2[ival] * NBIModel.rs_internal * NBIModel.rs_internal;
    df_rho    = NBIModel.new_rho[ival] * ( NBIModel.rs_internal * NBIModel.rs_internal * NBIModel.rs_internal );

    # interpolated forces
    a1 = -self.m1 * df_force1;
    a2 = -self.m2 * df_force2;
    at = a1 + a2;

    # df
    v21 = x[3]*x[3] + x[4]*x[4] + x[5]*x[5];
    v1  = math.sqrt(v21);

    xvalue = v1 / df_sigma;
    c1 = self.erf(xvalue) - 2.0 * xvalue / NBIModel.sqrtpi * math.exp(-xvalue*xvalue);

    # df formula with G=1
    c2 = -4.0 * math.pi * self.m2 * NBIModel.lnl / v21;
    c3 = c1 * c2 * df_rho;

    f[0] = x[3];
    f[1] = x[4];
    f[2] = x[5];

    f[3] = at * x[0]/r1 - c3 * x[3]/v1;
    f[4] = at * x[1]/r1 - c3 * x[4]/v1;
    f[5] = at * x[2]/r1 - c3 * x[5]/v1;


   
  # Computes the error function.  Based on code from Numerical Recipies and the following URL: 
  # http://www.cs.princeton.edu/introcs/21function/ErrorFunction.java.html
  # fractional error in math formula less than 1.2 * 10 ^ -7.
  # although subject to catastrophic cancellation when z in very close to 0
  # from Chebyshev fitting formula for erf(z) from Numerical Recipes, 6.2
  @staticmethod
  def erf(z) :
    t = 1.0 / (1.0 + 0.5 * math.fabs(z))

    # use Horner's method
    ans = 1.0 - t * math.exp( -z*z   -   1.26551223 + \
                                            t * ( 1.00002368 + \
                                            t * ( 0.37409196 +  \
                                            t * ( 0.09678418 +  \
                                            t * (-0.18628806 +  \
                                            t * ( 0.27886807 +  \
                                            t * (-1.13520398 +  \
                                            t * ( 1.48851587 +  \
                                            t * (-0.82215223 +  \
                                            t * ( 0.17087277))))))))))
    if (z >= 0):
      return ans
    else:
      return -ans;

    return ans
