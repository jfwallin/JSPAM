import numpy as np
import sys
import datetime

from Parameters import Parameters
from SetupUtil import SetupUtil
from SPMModel import SPMModel
from MONDModel import MONDModel
from NBIModel import NBIModel
from Integrator import Integrator
from IOUtil import IOUtil

class Run:

  def __init__(self):
    self.params = Parameters()
    #self.forceModel = SPMModel()
    #self.forceModel = MONDModel()
    #self.params.potential_type=2
    self.forceModel = NBIModel()
    self.params.potential_type=1
    self.integrator = Integrator(self.forceModel)
    x0 = None
    xout = None

  def initRun(self,args):
    su = SetupUtil()

    su.setHelpers(self.forceModel, self.integrator, self.params)
    su.customCollision(args)

    # update the forceModel based upon passed in args
    self.forceModel = su.createForceModel(self.params.potential_type,True)

    self.forceModel.setParameters(self.params)
    self.integrator.initRKVar(self.params.n)

    n = self.params.n
    self.x0 = np.zeros((n+1,6))
    self.xout = np.zeros((n+1,6))

    su.createCollision()

    self.copyParticles(self.integrator.x,self.x0)

    self.forceModel.setParameters(self.params)

  # end initRun


  def initRunP(self,params):
    su = SetupUtil()

    self.params = params
    su.setHelpers(self.forceModel, self.integrator, self.params)

    # update the forceModel based upon passed in args
    self.forceModel = su.createForceModel(self.params.potential_type,true)

    self.forceModel.setParameters(self.params)
    self.integrator.initRKVar(self.params.n)

    n = self.params.n
    #self.x0 = [[0] * 6 for i in range(n+1)]
    #self.xout = [[0] * 6 for i in range(n+1)]
    self.x0 = np.zeros((n+1,6))
    self.xout = np.zeros((n+1,6))

    su.createCollision()

    self.copyParticles(self.integrator.x,self.x0)

    self.forceModel.setParameters(self.params)
  # end initRunP



  def getMinValues(self,params):
    su = SetupUtil()

    su.setHelpers(self.forceModel, self.integrator, self.params)

    # update the forceModel based upon passed in args
    self.forceModel = su.createForceModel(self.params.potential_type,true)

    self.forceModel.setParameters(self.params)
    self.integrator.initRKVar(self.params.n)

    n = self.params.n
    #self.x0 = [[0] * 6 for i in range(n+1)]
    #self.xout = [[0] * 6 for i in range(n+1)]
    self.x0 = np.zeros((n+1,6))
    self.xout = np.zeros((n+1,6))

    mins = su.perturberPositionVec(self.params.sec_vec, self.params.mass1, self.params.mass2, 
                                   self.params.eps1, self.params.eps2,
                                   self.params.h, self.params.n, self.params.n1, self.params.time, self.integrator.x)

    return mins;
  # end getMinValues

  #Copies the particles from x1 to x2.
  def copyParticles(self,x1,x2):
    n = len(x1)
    np.copyto(x2,x1)
    # assuming 6 members
    #for i in range(n):
      #x2[i][0] = x1[i][0]
      #x2[i][1] = x1[i][1]
      #x2[i][2] = x1[i][2]
      #x2[i][3] = x1[i][3]
      #x2[i][4] = x1[i][4]
      #x2[i][5] = x1[i][5]

  # end copyParticles

  def takeAStep(self):
    self.integrator.rk4(self.x0,self.xout,self.params.h)
    self.copyParticles(self.xout,self.x0)
    self.params.time = self.params.time + self.params.h;



  def getFilename(self,i):

    st = str(i)
    while(len(st)<3):
      st = "0" + st

    return "a."+st



  # Run the simulation with current parameters from tstart to tend
  def calculate(self,tstart,tend):
    t0 = 0
    time_interval = 0

    nstep_local = 7500
    self.params.tstart = tstart
    self.params.tend = tend
    t0 = self.params.tstart
    self.params.nstep = ((self.params.tend-t0)/self.params.h)+2
    nstep_local = self.params.nstep
    time_interval = (self.params.tend-t0)*2

    for i in range(nstep_local):
      self.takeAStep()


def main():
  args = sys.argv[1:]
  run = Run()
  run.initRun(args)
  params = run.params;

  t0 = 0
  time_interval = 0

  nstep_local = 7500;

  t0 = params.tstart;
  params.nstep = ((params.tend-t0)/params.h)+2;
  nstep_local = params.nstep;
  time_interval = (params.tend-t0)*2;
  #IOUtil.writeParameterFile(params,"tmp.p")

  d1 = datetime.datetime.now()
  for i in range(1,int(nstep_local+1)):
    run.takeAStep()
    if(i % 10 == 5):
      run.params.iout = run.params.iout+1
      print(run.params.iout)
      #IOUtil.outputParticles(run.getFilename(run.params.iout), run.integrator.x)

  d2 = datetime.datetime.now()
  IOUtil.outputParticles(run.getFilename(run.params.iout), run.integrator.x)

  print(d1)
  print(d2)



if __name__ == "__main__": main()

