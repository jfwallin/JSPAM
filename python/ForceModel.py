
# This interface defines the methods required for defining a force model.
class ForceModel:
    def __init__(self):
      pass

    # Sets the simulation parameters to this model.
    def setParameters(self,p):
      pass

    # Initialize temporary storage.
    def initVars(self,n):
      pass

    # Cleanup temporary storage.
    def deallocateVars(self,):
      pass

    # Compute the circular velocity for a particle at a distance r from the specified mass.
    # The rout scale of the disk and softening length, eps, are provided.
    def circularVelocity(self,mass, r, rout, eps):
      pass

    # For the given particle positions and velocities, calculate the accelerations.
    def diffeq(self,x,f):
      pass

    # Calculate the acceleration felt by the secondary galaxy in orbit around the primary.
    def diffq1(self,x,f):
      pass
