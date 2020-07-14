import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rc('font', size=28)
plt.rc('font', weight='bold')
plt.rc('xtick', labelsize=28)
plt.rc('ytick', labelsize=28)

data = np.genfromtxt('stripped-300-big.xyz')
# data = np.genfromtxt('fftTest.xyz')
N = 21
eta_txz = data[:,2].reshape((-1, N, N))
T = eta_txz.shape[0]
X = eta_txz.shape[1]
Z = eta_txz.shape[2]
kB = 1.3806503 * 6.0221415 / 4184.0
temp = 300
beta = 1/(kB * temp)
a = 0.305
kappa = 0.321

step = T / 10000
eta_tk_real = np.zeros((T/2, X, Z))
eta_tk_imag = np.zeros((T/2, X, Z))
# eta_tk_real = np.zeros((X, Z))
# eta_tk_imag = np.zeros((X, Z))
kvecs = np.zeros((2, X, Z))

for kx in range(-(N-1)/2, (N+1)/2):
#     print kx
    xind = kx + (N-1)/2
    for kz in range(-(N-1)/2, (N+1)/2):
        zind = kz + (N-1)/2
        kvecs[0][xind][zind] = 2 * np.pi * kx / N
        kvecs[1][xind][zind] = 2 * np.pi * kz / N

for tind in range(T/2,T,step):
    t = (tind-T/2) / step
    for kx in range(-(N-1)/2, (N+1)/2):
        xind = kx + (N-1)/2
        for kz in range(-(N-1)/2, (N+1)/2):
            zind = kz + (N-1)/2
            for x in range(N):
                for z in range(N):
                    eta_tk_real[t][xind][zind] += eta_txz[tind][x][z] * ( np.cos(2*np.pi*x*kx/N + 2*np.pi*z*kz/N) )
                    eta_tk_imag[t][xind][zind] += eta_txz[tind][x][z] * ( np.sin(2*np.pi*x*kx/N + 2*np.pi*z*kz/N) )

# tind = 0
# t = 0
# for kx in range(-(N-1)/2, (N+1)/2):
#     xind = kx + (N-1)/2
#     for kz in range(-(N-1)/2, (N+1)/2):
#         zind = kz + (N-1)/2
#         for x in range(N):
#             for z in range(N):
#                 eta_tk_real[t][xind][zind] += eta_txz[tind][x][z] * ( np.cos(2*np.pi*x*kx/N + 2*np.pi*z*kz/N) )
#                 eta_tk_imag[t][xind][zind] += eta_txz[tind][x][z] * ( np.sin(2*np.pi*x*kx/N + 2*np.pi*z*kz/N) )
#         print kx, kz, eta_tk_real[t][xind][zind]/(N*N), np.mean(eta_txz[t, :, :]) 
       

# eta_tk_real /= (N*N/2)
# eta_tk_imag /= (N*N/2)

eta_var_real = []
eta_var_imag = []
eta_var = []
eta_mean = []
norm_k = []
predicted_eta_var = []

for kx in range(N):
    for kz in range(N):
        eta_mean.append( 0.5 * np.mean(eta_tk_real[:, kx, kz]) + 0.5 * np.mean(eta_tk_imag[:, kx, kz]) )
        eta_var_real.append( np.var(eta_tk_real[:, kx, kz]) )
        eta_var_imag.append( np.var(eta_tk_imag[:, kx, kz]) )
        kx_term = kx - (N-1)/2
        kz_term = kz - (N-1)/2
        eta_var.append( 0.5 * ( np.var(eta_tk_real[:, kx, kz]) + np.var(eta_tk_imag[:, kx, kz]) ) )
        knorm = 2 * np.pi / N * np.sqrt( kx_term**2 + kz_term**2 )
        predicted_eta_var.append( (N*N) / (2 * beta * ( a + kappa * ( np.sin(kx_term)**2 + np.sin(kz_term)**2 ) )) )
        norm_k.append( knorm )
#         norm_k.append( np.sqrt( kvecs[0][kx][kz]**2 + kvecs[1][kx][kz]**2 ) )

eta_mean = np.array(eta_mean)
eta_var = np.array(eta_var)
eta_var_real = np.array(eta_var_real)
eta_var_imag = np.array(eta_var_imag)
norm_k = np.array(norm_k)
predicted_eta_var = np.array(predicted_eta_var)

print norm_k
# print eta_mean
print eta_var

plt.plot(norm_k, predicted_eta_var, 'k', linestyle='dashed', linewidth=3)
plt.plot(norm_k, 2*eta_var, 'ro')
plt.plot(norm_k, 2*eta_var, 'r', linewidth=4, alpha=0.7)
plt.xlabel(r'$|k|$')
plt.ylabel(r'$\langle\hat\eta_k^2\rangle$')
plt.show()

exit(1)

norm_k_unique = np.unique(norm_k)
eta_var_unique = np.zeros(len(norm_k_unique))
predicted_eta_var_unique = np.zeros(len(norm_k_unique))

for i in range(len(norm_k_unique)):
    eta_var_unique[i] = np.mean( eta_var[ np.where( np.abs(norm_k-norm_k_unique[i]) < 1.0e-2) ] )
    predicted_eta_var_unique[i] = np.mean( predicted_eta_var[ np.where( np.abs(norm_k-norm_k_unique[i]) < 1.0e-2) ] )

# print eta_var_unique
print norm_k_unique

plt.plot(norm_k_unique, eta_var_unique, 'ro')
plt.plot(norm_k_unique, eta_var_unique, 'r', linewidth=4, alpha=0.7)
plt.xlabel(r'$|k|$')
plt.ylabel(r'$\langle\hat\eta_k^2\rangle$')
plt.show()



# # lattice spacing in each direction
# dx = 1
# dz = 1
# 
# def calc_kvecs(N, dx, dz):
#     kvec_lat = np.meshgrid( np.fft.fftfreq(N, d=dx), np.fft.fftfreq(N, d=dz) )
#     kvec_lat[0] *= 2*np.pi
#     kvec_lat[1] *= 2*np.pi
#     return np.array(kvec_lat)
# 
# kvecs_py = calc_kvecs(N, dx, dz)
# 
# print (kvecs - kvecs_py)[ np.where(np.abs(kvecs - kvecs_py) > 1.0e-2) ]

exit(1)

# print data.shape
# print eta_txz.shape

# real_space_lat = np.zeros((X, Z, 2))
# t = 0
# for i in range(X):
#     for j in range(Z):
#         real_space_lat[i][j][0] = data[t][0]
#         real_space_lat[i][j][1] = data[t][1]
#         t += 1
# 
# # print real_space_lat[1][9]

kx = kvecs[1]
ky = kvecs[0]

# print kx.shape
# print ky.shape
# print kvecs[0][5][5]
# print kvecs[1][5][5]

# start fourier transform of order parameter field
eta_tk = []
for t in range(1000, T):
    eta_tk.append(np.fft.fft2(eta_txz[t,:,:]))

eta_tk = np.array(eta_tk)
# eta_tk = np.array(eta_tk)/(N*N/2)
eta_tk0_real = np.copy(eta_tk[0,:,:].real)
eta_tk0_imag = np.copy(eta_tk[0,:,:].imag)
# print np.where(np.abs(eta_tk0_real)>1.0e-2)
# print kx
# print ky
# print eta_tk[0,:,:].shape
eta_tk0_real[np.abs(eta_tk0_real) < 1.0e-2] = np.nan
eta_tk0_imag[np.abs(eta_tk0_imag) < 1.0e-2] = np.nan

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(kx, ky, eta_tk0_real, s=30, c='r', marker='o', alpha=1.0)
ax.scatter(kx, ky, eta_tk0_imag, s=30, c='b', marker='v', alpha=1.0)
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
plt.xlim(-N,N)
plt.ylim(-N,N)
plt.show()

kx_vec = kx[:,0]
ky_vec = ky[0,:]

# print kx_vec, ky_vec
eta_k_realav = np.mean(eta_tk.real, axis=0)
eta_k_realvar = np.var(eta_tk.real, axis=0)
# print eta_k_realav
# print eta_k_realav[ np.where(np.abs(eta_k_realav) > 1.0e-4) ]
# print eta_k_realvar[ np.where(np.abs(eta_k_realav) > 1.0e-4) ]

eta_k_imagav = np.mean(eta_tk.imag, axis=0)
eta_k_imagvar = np.var(eta_tk.imag, axis=0)
# print eta_k_imagav
# print eta_k_imagav[ np.where(np.abs(eta_k_imagav) > 1.0e-4) ]
# print eta_k_imagvar[ np.where(np.abs(eta_k_imagav) > 1.0e-4) ]

# calculate variance as a function of norm

var_eta_k = []
norm_k = []
for i in range(X):
    for j in range(Z):
#         var_eta_k.append( ( np.var(eta_tk[:,i,j].real) + np.var(eta_tk[:,i,j].imag) ) / 2 )
        var_eta_k.append( ( np.var(eta_tk[:,i,j].real) + np.var(eta_tk[:,i,j].imag) )  )
        norm_k.append( np.sqrt( kx[i,j]*kx[i,j] + ky[i,j]*ky[i,j] ) )

var_eta_k = np.array(var_eta_k)
norm_k = np.array(norm_k)
# print var_eta_k.shape, norm_k.shape
# print norm_k[np.where(var_eta_k > 1.0e-2)]
# print var_eta_k[np.where(var_eta_k > 1.0e-2)]
plt.clf()
# sort_indices = norm_k.argsort()
# var_eta_k_sorted = var_eta_k[sort_indices]
# norm_k_sorted = norm_k[sort_indices]

# plt.plot(np.sort(var_eta_k), 'ro')
plt.plot(norm_k, var_eta_k, 'ro')
plt.hlines(np.mean(var_eta_k), 0, 5, 'k', linestyle='dashed', linewidth=3)
# plt.plot(norm_k_sorted, var_eta_k_sorted, 'ro')
# plt.loglog(norm_k_sorted, var_eta_k_sorted, 'r', linewidth=4, alpha=0.5)
plt.xlabel(r'$|k|$')
plt.ylabel(r'$\langle\hat\eta_k^2\rangle$')
plt.show()

# plot probability distributions of eta_k at specified (kx, ky)
i = 4
j = 6
kx1 = kx[i,j]
ky1 = ky[i,j]
norm_k1 = np.sqrt(kx1*kx1 + ky1*ky1)
bin_centres = np.linspace(eta_tk[:,i,j].real.min(), eta_tk[:,i,j].real.max(), 50)
# print bin_centres

plt.hist(eta_tk[:,i,j].real, bins=bin_centres, histtype='stepfilled', alpha=0.7, color='r')
plt.xlabel(r'$\boldsymbol{\eta_k}$', fontsize=28)
plt.ylabel(r'$\boldsymbol{P(\eta_k)}$', fontsize=28)
plt.title(r'$\textbf{{Probability distribution at site {:.0f}, {:.0f} with }} \boldsymbol{{\|k\|}} \textbf{{ = {:1.3f} }}$'.format(i, j, norm_k1))
plt.xlim(-60, 60)
plt.show()

i = 0
j = 0
kx1 = kx[i,j]
ky1 = ky[i,j]
norm_k1 = np.sqrt(kx1*kx1 + ky1*ky1)
bin_centres = np.linspace(eta_tk[:,i,j].real.min(), eta_tk[:,i,j].real.max(), 50)
# print bin_centres

plt.hist(eta_tk[:,i,j].real, bins=bin_centres, histtype='stepfilled', alpha=0.7, color='r')
plt.xlabel(r'$\boldsymbol{\eta_k}$', fontsize=28)
plt.ylabel(r'$\boldsymbol{P(\eta_k)}$', fontsize=28)
plt.title(r'$\textbf{{Probability distribution at site {:.0f}, {:.0f} with }} \boldsymbol{{\|k\|}} \textbf{{ = {:1.3f} }}$'.format(i, j, norm_k1))
plt.xlim(-60, 60)
plt.show()

i = 5
j = 5
kx1 = kx[i,j]
ky1 = ky[i,j]
norm_k1 = np.sqrt(kx1*kx1 + ky1*ky1)
bin_centres = np.linspace(eta_tk[:,i,j].real.min(), eta_tk[:,i,j].real.max(), 50)
# print bin_centres

plt.hist(eta_tk[:,i,j].real, bins=bin_centres, histtype='stepfilled', alpha=0.7, color='r')
plt.xlabel(r'$\boldsymbol{\eta_k}$', fontsize=28)
plt.ylabel(r'$\boldsymbol{P(\eta_k)}$', fontsize=28)
plt.title(r'$\textbf{{Probability distribution at site {:.0f}, {:.0f} with }} \boldsymbol{{\|k\|}} \textbf{{ = {:1.3f} }}$'.format(i, j, norm_k1))
plt.xlim(-60, 60)
plt.show()

i = 8
j = 9
kx1 = kx[i,j]
ky1 = ky[i,j]
norm_k1 = np.sqrt(kx1*kx1 + ky1*ky1)
bin_centres = np.linspace(eta_tk[:,i,j].real.min(), eta_tk[:,i,j].real.max(), 50)
# print bin_centres

plt.hist(eta_tk[:,i,j].real, bins=bin_centres, histtype='stepfilled', alpha=0.7, color='r')
plt.xlabel(r'$\boldsymbol{\eta_k}$', fontsize=28)
plt.ylabel(r'$\boldsymbol{P(\eta_k)}$', fontsize=28)
plt.title(r'$\textbf{{Probability distribution at site {:.0f}, {:.0f} with }} \boldsymbol{{\|k\|}} \textbf{{ = {:1.3f} }}$'.format(i, j, norm_k1))
plt.xlim(-60, 60)
plt.show()

i = 0
j = 9
kx1 = kx[i,j]
ky1 = ky[i,j]
norm_k1 = np.sqrt(kx1*kx1 + ky1*ky1)
bin_centres = np.linspace(eta_tk[:,i,j].real.min(), eta_tk[:,i,j].real.max(), 50)
# print bin_centres

plt.hist(eta_tk[:,i,j].real, bins=bin_centres, histtype='stepfilled', alpha=0.7, color='r')
plt.xlabel(r'$\boldsymbol{\eta_k}$', fontsize=28)
plt.ylabel(r'$\boldsymbol{P(\eta_k)}$', fontsize=28)
plt.title(r'$\textbf{{Probability distribution at site {:.0f}, {:.0f} with }} \boldsymbol{{\|k\|}} \textbf{{ = {:1.3f} }}$'.format(i, j, norm_k1))
plt.xlim(-60, 60)
plt.show()



