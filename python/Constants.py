import math

G = 6.673e-11
solarMassToKg = 1.98892e30
metersPerMpc = 3.08568025e22
kmPerMpc = 3.08568025e19 # or m per kpc
degToRad = math.pi/180.0
    
MU_kg = 1.0e11*solarMassToKg
DU_m = 15.0*kmPerMpc # meters in a kpc
TU_s = math.sqrt(DU_m*DU_m*DU_m/(G*MU_kg))

MU_sm = 1.0e11
DU_mpc = 15.0/1000.0 

VEL_KMS = DU_m/1000.0/TU_s
A_MSS = DU_m/TU_s/TU_s

A0_MKS = 1e-10
A0 = A0_MKS/A_MSS

DEFAULT_EPS = math.sqrt(0.1)

# dump constants
def main():
  print(G)
  print(MU_kg)
  print(DU_m)
  print(TU_s)
  print(MU_sm)
  print(DU_mpc)
  print(VEL_KMS)
  print(A_MSS)
  print(A0)

if __name__ == "__main__": main()
