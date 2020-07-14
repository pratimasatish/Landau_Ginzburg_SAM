import numpy as np

randopt=0

tmax=10000
Nmax=11

periodx=3.0
periody=1.0

AmpMean=5.0
AmpStdev=2.0
AmpCos=1.0
AmpSin=0.0

alpha = 2.0*np.pi*periodx/Nmax
beta  = 2.0*np.pi*periody/Nmax
fp=open("fftTest.xyz", "w")
for t in range(tmax):
  AmpCos = np.random.normal(AmpMean, AmpStdev)
#   AmpSin = np.random.normal(AmpMean, AmpStdev)
  for i in range(Nmax):
    for j in range(Nmax):
      fp.write("{} {} {}\n".format(i, j, AmpCos*np.cos(alpha*i + beta*j) + AmpSin*np.sin(alpha*i + beta*j)))
fp.close()  
