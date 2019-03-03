import random
import math

from MathUtil import MathUtil
from Parameters import Parameters
from SPMModel import SPMModel
from NBIModel import NBIModel
from MONDModel import MONDModel

class SetupUtil:
  def __init__(self):
    self.forceModel = None
    self.integrator = None
    self.params = None


  def createForceModel(self,potential_type, apply):
    self.force = None

    if potential_type == 0:
      force = SPMModel()
    elif potential_type == 1:
      force = NBIModel()
    elif potential_type == 2:
      force = MONDModel()
    else:
      force = SPMModel()

    if apply :
      self.forceModel = force
      if self.integrator is not None:
        self.integrator.force=force

    return force



  def setHelpers(self, forceModel, integrator, params):
    self.forceModel = forceModel
    self.params = params
    self.integrator = integrator



  #Determine if the caller provides a parameter file or a parameter string
  def customCollision(self,args):

    self.params.tIsSet = False

    if len(args) > 0:
      if len(args) > 1 and args[0] == "-f":
        #parse the file
        IOUtil.readParameterFile(self.params,args[1])

      else:
        # parse the state string
        IOUtil.parseStateInfoString(self.params,args[0])
        self.params.potential_type=1
        self.params.h = Parameters.hbase
        self.params.tstart = -5
        self.params.tend = 0
        self.params.time = -5

        if len(args) > 1:
          t = IOUtil.parseDouble(args[1])
          if(t != 0):
            self.params.tstart = t
            self.params.time = t 
            self.params.tIsSet = true

    else:
      self.params.phi1   = 5.0
      self.params.theta1 = 5.0
      self.params.rscale1 = [1.0,1.0,1.0]
      self.params.rout1   = 1.0
      self.params.mass1   = 1.0
      self.params.epsilon1 = 0.3
      self.params.eps1 = self.params.epsilon1*self.params.epsilon1
      self.params.n1 = 1000
      self.params.heat1 = 0.0
      self.params.opt1 = 1

      self.params.phi2   = 0.0
      self.params.theta2 = 0.0
      self.params.rscale2 = [0.3,0.3,0.3]
      self.params.rout2   = 0.50
      self.params.mass2   = 0.50
      self.params.epsilon2 = 0.3
      self.params.eps2 = self.params.epsilon2*self.params.epsilon2
      self.params.n2 = 500
      self.params.heat2 = 0.0
      self.params.opt2 = 1

      self.params.inclination_degree = 20.0
      self.params.omega_degree = 0.0
      self.params.rmin = 0.90
      self.params.velocity_factor = 0.90

      self.params.h = Parameters.hbase
      self.params.time = -5
      self.params.tstart = -5
      self.params.tIsSet = True


    self.params.n = self.params.n1 + self.params.n2
    self.params.eps1 = self.params.epsilon1*self.params.epsilon1
    self.params.eps2 = self.params.epsilon2*self.params.epsilon2

  # end customCollision







  #Initialize the disks and set the initial positions.
  def createCollision(self):
    self.profile(self.params.rin1, self.params.rout1, self.params.rscale1, 0, self.params.n1, self.params.mass1,
                self.params.eps1, self.params.theta1, self.params.phi1, self.params.opt1, self.params.heat1,
                self.integrator.x)        

    self.profile(self.params.rin2, self.params.rout2, self.params.rscale2, self.params.n1, self.params.n, self.params.mass2,
                self.params.eps2, self.params.theta2, self.params.phi2, self.params.opt2, self.params.heat2,
                self.integrator.x)        

    # determine if we need to calculate tStart
    if( not self.params.tIsSet ) :
      sec_vec = self.params.sec_vec
      rv4min = [sec_vec[0],sec_vec[1],sec_vec[2],-sec_vec[3],-sec_vec[4],-sec_vec[5],0.0]
      tminVals = self.getTStart(rv4min, self.params.mass1, self.params.mass2, self.params.eps1, self.params.eps2 , self.params.h,-30.0 ,10.0*self.params.rout1, self.params.rout1, self.params.rout2)

      tmpT = tminVals[0]
      if ( tmpT < -12.0 ) :
        tmpT = -5

      if ( math.fabs(tmpT) < self.params.h):
        tmpT = -5

      self.params.tstart = tmpT
      self.params.time = tmpT
      self.params.tIsSet  = true


    # end !tIsSet

    mins = None
    if(self.params.use_sec_vec):
      mins = self.perturberPositionVec(self.params.sec_vec, self.params.mass1, self.params.mass2, 
                                        self.params.eps1, self.params.eps2,
                                        self.params.h, self.params.n, self.params.n1, self.params.time, self.integrator.x)
    else:
      mins = self.perturberPosition(self.params.inclination_degree, self.params.omega_degree, 
                                     self.params.rmin, self.params.velocity_factor, self.params.rout1,
                                     self.params.mass1, self.params.mass2, self.params.eps1, self.params.eps2,
                                     self.params.h, self.params.n, self.params.n1, self.params.time, 
                                     self.integrator.x, self.params.sec_vec)

    self.params.tmin = mins[0]
    self.params.rmin = mins[1]
    self.params.vrmin = mins[2]
    self.params.beta = mins[3]

  # end createCollision

  def randm(self):
    return random.random() 


  def distrb(self,r1, opt, rscale):
    distrb = 0
        
    if (opt == 1):
      distrb = 1.0/r1
    elif (opt == 2):
      distrb = math.exp(-r1/rscale[0])
    elif (opt == 3):
      distrb = math.exp(-r1*r1*rscale[0] - rscale[1]*r1 - rscale[2] )
        
    return distrb
    #  end function distrb


  def profile(self,rin, rout, rscale, nstart, ntot, mass, eps, theta, phi, opt, heat, x0):

  #     double stheta,ctheta,sphi,cphi,pi
  #     double x3,y3,z3,xv3,yv3,zv3,x2,y2,z2,xv2,yv2,zv2
  #     double x,y,z,xv,yv,zv

  #     int i, j, n, nprof

  #     double rnorm
  #     double[] rp, r, angle, v, p_ring, cp_ring
  #     double st, ct, dr, ran, r1, r2, ptot
  #     int[] n_ring
  #     int nring, dnring, is, ie, iring, tct

  #     double xp, yp, zp
  #     double fx, fy, fz, ftotal
  #     double tmp

  #     int ival

    n = len(x0)
    if(ntot > n):
      ntot = n
      print("ntot exceeds n")

    r = [0]*n
    angle = [0]*n
    v = [0]*n

    pi = math.pi

    stheta = math.sin(theta*pi/180.0)    
    ctheta = math.cos(theta*pi/180.0)    
    sphi   = math.sin(phi*pi/180.0)    
    cphi   = math.cos(phi*pi/180.0)    

    # set up the probability distribution for the disk
    nprof = 1000
    nring = int(nprof / 10)

    dnring = nprof/nring
    rp = [0]*nprof
    n_ring = [0]*nprof
    p_ring = [0]*nprof
    cp_ring = [0]*nprof

    # set the differential sum of the probability funtion into a vector
    rnorm = 0.0
    dr = (rout - rin)/(nprof)

    for i in range(nprof):
      r1 = i*dr + rin
      rp[i] = self.distrb(r1, opt, rscale) * r1 * dr * 2.0 * pi
      rnorm = rnorm + rp[i]


    # normalize the vector
    for i in range(nprof):
      rp[i] /= rnorm

    # take the fine bins and put them into the selection bins
    tct = 0
    for iring in range(nring):
      isi = int((iring) * dnring + 1)
      ie = int((iring+1) * dnring )

      ptot = 0.0

      for i in range(isi,ie):
        ptot = ptot + rp[i]

      p_ring[iring] = ptot


    # formulative cumulative distribution function
    cp_ring[0] = p_ring[0]
    for iring in range(1,nring):
      cp_ring[iring] = cp_ring[iring -1] + p_ring[iring]


    # find the number of particles in each bin
    n_ring = [0]*nprof

    for i in range(nstart-1,ntot):

      # find the radial position bin
      ran = self.randm()
      j = 0
      while(j<nring and ran > cp_ring[j]):
        j = j + 1

      j = j-1
      if(j<0):
          j=0

      n_ring[j] = n_ring[j] + 1
    # end for


    tct = 0
    i = nstart
    for iring in range(nring):
      isi = (iring) * dnring + 1
      ie = (iring+1) * dnring

      r1 = (isi)*dr + rin
      r2 = (ie)*dr + rin

      for j in range(n_ring[iring]):
        if i>=n:
          break
        ran = self.randm()
        r[i] = r1 + ran * (r2 - r1)
        i = i + 1

    #enddo

    # set the angular positions and orbital velocities
    for i in range(nstart,ntot):
      angle[i] = 2.0 * pi * self.randm()
      v[i] = self.forceModel.circularVelocity(mass,r[i],rout,eps)

    # set position and velocity based on the distribution parameters
    for i in range(nstart,ntot):
      st  =  math.sin(angle[i])
      ct  =  math.cos(angle[i])

      x   =  ct*r[i]
      y   =  st*r[i]
      z   =  0.0

      xv  = -v[i]*st
      yv  =  v[i]*ct
      zv  =  0.0

      x2  =   x * ctheta +  z * stheta
      y2  =   y
      z2  =  -x * stheta +  z * ctheta
      xv2 =  xv * ctheta + zv * stheta
      yv2 =  yv
      zv2 = -xv * stheta + zv * ctheta

      x3  =  x2  * cphi -  y2 * sphi
      y3  =  x2  * sphi +  y2 * cphi
      z3  =  z2
      xv3 =  xv2 * cphi - yv2 * sphi
      yv3 =  xv2 * sphi + yv2 * cphi
      zv3 =  zv2
      

      x0[i][0] = x3
      x0[i][1] = y3
      x0[i][2] = z3
      x0[i][3] = xv3  + self.randm()*heat
      x0[i][4] = yv3  + self.randm()*heat
      x0[i][5] = zv3  + self.randm()*heat


  # end profile




  def perturberPosition(self,inclinationDegree, omegaDegree, rMin, velocityFactor, rout1, mass1, mass2, eps1, eps2, h, n, n1, t0, x0, sec_vec):
    xx0 = [0,0,0,0,0,0]

    omega = math.pi*(omegaDegree)/180.0
    inc = math.pi*(inclinationDegree)/180.0

    v = math.sqrt(2.0)*self.forceModel.circularVelocity(mass1+mass2,rMin,rout1,eps1)

    v = -v * velocityFactor

    #      setup the transformaton based on the matrix in
    #      fundamentals of astrodynamics p. 82 by
    #      bate, mueller, and white (1971)

    xx0[0] = math.cos(omega) * rMin
    xx0[1] = math.sin(omega) * math.cos(inc) * rMin
    xx0[2] = math.sin(omega) * math.sin(inc) * rMin

    xx0[3] = -math.sin(omega) * v
    xx0[4] =  math.cos(omega) * math.cos(inc) * v
    xx0[5] =  math.cos(omega) * math.sin(inc) * v

    sec_vec[0] = xx0[0]
    sec_vec[1] = xx0[1]
    sec_vec[2] = xx0[2]
    sec_vec[3] = -xx0[3]
    sec_vec[4] = -xx0[4]
    sec_vec[5] = -xx0[5]

    return self.perturberPositionVec(sec_vec, mass1, mass2, eps1, eps2, h, n, n1, t0, x0)

  # end perturberPosition



  def perturberPositionVec(self,xx0, mass1, mass2, eps1, eps2, h, n, n1, t0, x0):
    xxe = [0,0,0,0,0,0]
    tcurrent = 0
    i = 0


    # reverse the velocity for backward integration
    xx0[3] = -xx0[3]
    xx0[4] = -xx0[4]
    xx0[5] = -xx0[5]

    tmin = tcurrent
    tmpr = 0
    tmpv = 0

    # avoid multiple calls to sqrt, save it until the end
    rmin = xx0[0]*xx0[0]+xx0[1]*xx0[1]+xx0[2]*xx0[2]
    vrmin = xx0[3]*xx0[3]+xx0[4]*xx0[4]+xx0[5]*xx0[5]

    # now move position back to t0 from t=0.0
    while( t0 < tcurrent ):
      self.integrator.rk41(xx0, xxe, h)
      xx0[0] = xxe[0]
      xx0[1] = xxe[1] 
      xx0[2] = xxe[2]
      xx0[3] = xxe[3]  
      xx0[4] = xxe[4] 
      xx0[5] = xxe[5]

      tcurrent = tcurrent - h

      tmpr = xx0[0]*xx0[0]+xx0[1]*xx0[1]+xx0[2]*xx0[2]

      if(tmpr < rmin):
        rmin = tmpr
        vrmin = xx0[3]*xx0[3]+xx0[4]*xx0[4]+xx0[5]*xx0[5]
        tmin = tcurrent

    # end while


    # reverse the velocity for forward integration
    xx0[3] = -xx0[3]
    xx0[4] = -xx0[4]
    xx0[5] = -xx0[5]

    # now adjust the test particles from the second disk
    # to the proper velocity and positions
    if(self.params.n > self.params.n1):
      for i in range(self.params.n1,self.params.n):
        x0[i][0] += xx0[0]
        x0[i][1] += xx0[1]
        x0[i][2] += xx0[2]
        x0[i][3] += xx0[3]
        x0[i][4] += xx0[4]
        x0[i][5] += xx0[5]



    # include the perturbing galaxy
    i = self.params.n
    x0[i][0] += xx0[0]
    x0[i][1] += xx0[1]
    x0[i][2] += xx0[2]
    x0[i][3] += xx0[3]
    x0[i][4] += xx0[4]
    x0[i][5] += xx0[5]

    vrmin = math.sqrt(vrmin)
    beta = (mass1+mass2)/(rmin*vrmin)
    rmin = math.sqrt(rmin)

    return [tmin,rmin,vrmin,beta]

  #end perturberPositionVec


  # Convert the position and velocity to classical orbital
  # elements.
   
  # coe[0] = p
  # coe[1] = ecc
  # coe[2] = inc
  # coe[3] = LAN
  # coe[4] = w
  # coe[5] = v
  # coe[6] = u        

  def rvToCoe(r,v,mu):
    muInv = 1.0/mu

    rmag = MathUtil.mag(r)
    vmag = MathUtil.mag(v)
    
    h = MathUtil.cross(r,v)
    hmag = MathUtil.mag(h)
   
    K = [0.0,0.0,1.0] 
    n = MathUtil.cross(K,h)
    nmag = MathUtil.mag(n)
    
    tmp1 = vmag*vmag - mu/rmag
    tmp2 = MathUtil.dot(r,v)
    
    v1 = MathUtil.scale(tmp1,r)
    v2 = MathUtil.scale(tmp2,v)
    
    ev = MathUtil.sub(v1,v2)
    ev = MathUtil.scale(muInv,ev)
    
    p = hmag*hmag*muInv
    ecc = MathUtil.mag(ev)
    cosi = h[2]/hmag
    cosO = n[0]/nmag
    cosw = MathUtil.dot(n,ev)/(nmag*ecc)
    cosv = MathUtil.dot(ev,r)/(ecc*rmag)
    cosu = MathUtil.dot(n,r)/(nmag*rmag)

    outv = [0]*7
    outv[0] = p
    outv[1] = ecc
    outv[2] = math.acos(cosi)
    
    tmp1 = math.acos(cosO)
            
    if (n[0] < 0):
      tmp1 = 2.0*math.PI-tmp1

    outv[3] = tmp1 
        
    tmp1 = math.acos(cosw)
   
    if (ev[1] < 0 ) :
      tmp1 = 2.0*math.PI-tmp1

    outv[4] = tmp1 
    
    tmp1 = math.acos(cosv)
    
    if (MathUtil.dot(r,v) < 0):
      tmp1 = 2.0*math.PI-tmp1

    outv[5] = tmp1 
        
    if(cosu > 1.0 or cosu < -1.0) :
      outv[6] = Double.NaN
    else:
      tmp1 = math.acos(cosu)
      if(r[1]>0) :
        tmp1 = 2.0*math.PI-tmp1
      outv[6] = tmp1

    return outv
    # end rvToCoe    

    
  # * Find the time of rmin, assuming earlier than now, given
  # * the r and v values.  Returns r and v at time of rmin
  # * by replacing r and v.  r and v are given as
  # * {rx,ry,rz,vx,vy,vz}.
  def getTStart(self,rv,mass1,mass2,eps1,eps2,h,tmin,mind,rout1,rout2):

    i = 0

    mu = mass1+mass2
    t = 0.0
   
    r =[0,0,0]
    r[0] = rv[0]
    r[1] = rv[1]
    r[2] = rv[2]
    v = [-rv[3],-rv[4],-rv[5]] 


    coe = self.rvToCoe(r,v,mu)

    ecc = coe[1]

    a = coe[0]/(1.0-ecc*ecc)
    period = 0.0
    apocenter = a*(1.0+ecc)
    a2 = apocenter*apocenter
    tApp = 0.0
      
    isEllipse = False
            
    if (ecc < 1.0):
      isEllipse = True
      period = 2.0*math.pi/math.sqrt(mu)*math.pow(a,1.5)
      period = period * 1.0
    
      
    xxe = [0]*7
    rvmin = [0]*7
     
    for i in range(7):
      rvmin[i] = rv[i] 

      
    distNew = rv[0]*rv[0] + rv[1]*rv[1] + rv[2]*rv[2]
    distOld = 2.0*distNew
      
    distNearApp = -1e30
   
    # keep looping as long as distance is decreasing
    while(tmin < t) :
      r[0] = rv[0] 
      r[1] = rv[1] 
      r[2] = rv[2]
      v[0] = rv[3] 
      v[1] = rv[4] 
      v[2] = rv[5]

      coe = self.rvToCoe(r,v,mu)
      xxe[6] =t +h

      self.integrator.rk41(rv, xxe, h)
         
      distNew = xxe[0]*xxe[0] + xxe[1]*xxe[1] + xxe[2]*xxe[2]
          
      # if it's ellipse and it's near apocenter, take this time
      if ( isEllipse and (math.fabs(distNew-a2)/a2 < 0.05) ):
        if(distNew > distNearApp):
          distNearApp = distNew
          tApp = t


          
      if (distNew < distOld):
        distOld = distNew
        for i in range(7):
          rvmin[i] = xxe[i]

        rvmin[6] = rvmin[6] - h

      for i in range(7):
        rv[i] = xxe[i]

      rv[6] = xxe[6] - h * 2.0
      t = t - h
    # end loop

    for i in range(7):
      rv[i] = rvmin[i]
      
    minDist = math.sqrt(rv[0]*rv[0] + rv[1]*rv[1] + rv[2]*rv[2])
    minVel = math.sqrt(rv[3]*rv[3] + rv[4]*rv[4] + rv[5]*rv[5])
      
    t = rv[6]
      
    if(isEllipse and tApp < 0.0):
      t = tApp
    else:
      t = t - mind/minVel
      
    outStuff = [t,minDist,minVel,rv[6]]
      
      
    return outStuff

  # end getTStart
