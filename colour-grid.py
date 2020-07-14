import numpy as np
import argparse
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from mpl_toolkits.mplot3d import Axes3D
import math

parser = argparse.ArgumentParser(description="")
parser.add_argument("-temp", default='360', type=str, help="temp value to analyse")
parser.add_argument("-fil", default='360', type=str, help="file keyword to analyse")
parser.add_argument("-boxsize", default=11, type=int, help="length of each dimension")
parser.add_argument("-savefigs", action='store_true', help="whether to save figures of CG lattice to make movie or not")
args = parser.parse_args()

Lx = args.boxsize
Lz = args.boxsize
data_tot = np.genfromtxt('stripped-' + args.fil + '.xyz', delimiter=' ')
data_txz = data_tot[:,2].reshape((-1,Lx,Lz))

plt.rc('text', usetex=True)
plt.rc('font', family='serif', weight='bold')
plt.rcParams['text.latex.preamble']=[r"\usepackage{amsmath}"]
plt.rc('font', size=28)
plt.rc('xtick', labelsize=28)
plt.rc('ytick', labelsize=28)

T = len(data_txz)

# get eta values at each lattice site as a function of time
eta_lat = []
for i in range(T):
    eta_lat.append( data_txz[i, :, :].flatten() )
eta_lat = np.array(eta_lat)
eta_lat = eta_lat.reshape((-1, Lx, Lz))
# print np.mean(eta_lat)
max_eta = eta_lat.max()
min_eta = eta_lat.min()
# print eta_lat.shape
bins = np.linspace(-0.6, 0.6, 101)
bin_centres = 0.5 * (bins[1:] + bins[:-1])

a = -31.08
c = 701.3
kappa = 32.71
T0 = 375.6
h = -0.0488

temp = float(args.temp)
# print 'temp = {}'.format(temp)

def hamiltonian(eta, T):
    return h * (T - T0) * eta + a * eta**2 / 2.0 + c * eta**2 * eta**2 / 4.0


# hist, bins = np.histogram(eta_lat[4000:,3,4], bins=bins, density=True)
# # dummy = hist[hist != 0]
# # min_log = (-np.log(dummy)).min()
# plt.plot(bin_centres, -np.log(hist), 'ro')
# plt.plot(bin_centres, -np.log(hist), 'r', linewidth=4, alpha=0.6)
# plt.plot(bin_centres, hamiltonian(bin_centres, temp), 'b', linewidth=4, alpha=0.8)
# plt.title('histogram for site (4,5) at T = {} K'.format(temp))
# plt.show()
# 
# hist, bins = np.histogram(eta_lat[4000:,2,9], bins=bins, density=True)
# # dummy = hist[hist != 0]
# # min_log = (-np.log(dummy)).min()
# plt.plot(bin_centres, -np.log(hist), 'ro')
# plt.plot(bin_centres, -np.log(hist), 'r', linewidth=4, alpha=0.7)
# plt.plot(bin_centres, hamiltonian(bin_centres, temp), 'b', linewidth=4, alpha=0.8)
# plt.title('histogram for site (3,10) at T = {} K'.format(temp))
# plt.show()
# 
# hist, bins = np.histogram(eta_lat[4000:,7,10], bins=bins, density=True)
# # dummy = hist[hist != 0]
# # min_log = (-np.log(dummy)).min()
# plt.plot(bin_centres, -np.log(hist), 'ro')
# plt.plot(bin_centres, -np.log(hist), 'r', linewidth=4, alpha=0.7)
# plt.plot(bin_centres, hamiltonian(bin_centres, temp), 'b', linewidth=4, alpha=0.8)
# plt.title('histogram for site (8,11) at T = {} K'.format(temp))
# plt.show()
# 
# eta_t = np.mean(eta_lat, axis=(1,2))
# hist_full, bins = np.histogram(eta_t[4000:], bins=bins, density=True)
# plt.plot(bin_centres, -np.log(hist_full), 'ro')
# plt.plot(bin_centres, -np.log(hist_full), 'r', linewidth=4, alpha=0.7)
# plt.plot(bin_centres, hamiltonian(bin_centres, temp), 'b', linewidth=4, alpha=0.8)
# plt.title('histogram for spatially averaged eta at T = {} K'.format(temp))
# plt.show()

eta_mean = np.mean(eta_lat, axis=0)

plt.imshow(eta_mean, aspect=1.0, cmap="seismic_r", origin="lower", interpolation="none", vmin=min_eta, vmax=max_eta)
plt.xticks(np.arange(0, Lx, 1))
plt.yticks(np.arange(0, Lz, 1))
plt.xlim(-0.5,Lx-0.5)
plt.ylim(-0.5,Lz-0.5)
for i in np.arange(-0.5,Lx-1,1.0):
    plt.vlines(i, -0.5, Lz-0.5, linestyle='solid', linewidth=2)
for i in np.arange(-0.5,Lz-1,1.0):
    plt.hlines(i, -0.5, Lx-0.5, linestyle='solid', linewidth=2)
plt.colorbar()
plt.show()

plt.imshow(eta_lat[9,:,:], aspect=1.0, cmap="seismic_r", origin="lower", interpolation="none", vmin=-0.4, vmax=0.4)
plt.xticks(np.arange(0, Lx, 1))
plt.yticks(np.arange(0, Lz, 1))
plt.xlim(-0.5,Lx-0.5)
plt.ylim(-0.5,Lz-0.5)
for i in np.arange(-0.5,Lx-1,1.0):
    plt.vlines(i, -0.5, Lz-0.5, linestyle='solid', linewidth=2)
for i in np.arange(-0.5,Lz-1,1.0):
    plt.hlines(i, -0.5, Lx-0.5, linestyle='solid', linewidth=2)
plt.title('Snapshot at t=10')
plt.colorbar()
plt.show()

eta_x = np.mean(eta_mean, axis=0)
x_axis = np.arange(-Lx/2+1, Lx/2+1, 1)
def surf_profile(x, eta0, l):
    return eta0 * np.tanh(x / l)

popt, pcov = curve_fit(surf_profile, x_axis, eta_x, p0=(0.3, 1.5))
print 'best fit is eta0 = {} and l = {}'.format(popt[0], popt[1])
eta0 = popt[0]
l = popt[1]
surf_fit = surf_profile(x_axis, eta0, l)
# surf_fit = surf_profile(x_axis, 0.3, 1.45)

plt.plot(x_axis, eta_x, 'ro')
plt.plot(x_axis, eta_x, 'r', linewidth=4, alpha=0.4, label='simulation data')
plt.plot(x_axis, surf_fit, 'bv')
plt.plot(x_axis, surf_fit, 'b', linewidth=4, alpha=0.4, label='fit')
plt.legend(loc='best')
plt.show()


exit(1)


if args.savefigs:
    name_arr = range(0, eta_lat.shape[0], 100)
    name_arr = np.array(name_arr)
    for j in range(len(name_arr)):
    # for j in range(0, 100, 20):
        matr = eta_lat[name_arr[j]].transpose()
#         plt.imshow(matr, aspect=1.0, cmap="seismic_r", origin="lower", interpolation="none", vmin=min_eta, vmax=max_eta)
        plt.imshow(matr, aspect=1.7, cmap="PRGn_r", origin="lower", interpolation="none", vmin=-0.4, vmax=0.4)
        plt.yticks(np.arange(0, Lz, 1))
        plt.xticks(np.arange(0, Lx, 1))
        plt.ylim(-0.5,10.5)
        plt.xlim(-0.5,10.5)
        for i in np.arange(-0.5,10,1.0):
            plt.hlines(i, -0.5, 10.5, linestyle='solid', linewidth=2)
        for i in np.arange(-0.5,10,1.0):
            plt.vlines(i, -0.5, 10.5, linestyle='solid', linewidth=2)
        plt.colorbar()
#         plt.title('timestep = {}'.format(j*100))
        plt.title('T = {} K'.format(temp))
        plt.savefig('lat-' + args.fil + '-{:05d}.png'.format(j))
#         plt.savefig('lat-0.7250-{:05d}.png'.format(j))
        plt.clf()

eta_lat -= eta_lat.mean()
all_corr_xz = []
samplez = 10
# DT = t_steps / samplez
X = eta_lat.shape[1]
Z = eta_lat.shape[2]
DT = eta_lat.shape[0] / samplez
fig = plt.figure()
ax  = fig.add_subplot(111, projection='3d')
var_list = []
for sample in xrange(samplez):
    To = DT*sample
    Tf = DT*(sample+1)
    sub_txz = eta_lat[To:Tf, :, :]
    
    cov_xz = np.zeros((X/2+1, Z/2+1))
    for t in xrange(DT):
       for x in xrange(X/2 + 1):
            for z in xrange(Z/2 + 1):
#                 print cov_xz.shape
#                 print sub_txz[t,x,z].shape
#                 print sub_txz[t,x : x + X/2,z : z + Z/2].shape
                cov_xz += sub_txz[t, x, z] * sub_txz[t, x : x + X/2 + 1, z : z + Z/2 + 1]
    
    cov_xz /= (DT * (X/2+1) * (Z/2+1) )
    var_list.append(cov_xz[0,0])
    corr_xz = cov_xz / cov_xz[0,0]
#     corr_xz = cov_xz
    
    x = range(X/2 + 1)
    z = range(Z/2 + 1)
    xv, zv = np.meshgrid(x, z)
    
    all_corr_xz.append(corr_xz)

all_corr_xz = np.array(all_corr_xz)
m_corr_xz = np.mean(all_corr_xz, axis=0)
d_corr_xz = np.std (all_corr_xz, axis=0) / np.sqrt(samplez)
ax.plot_surface(xv, zv, m_corr_xz.T)
plt.show()

# print np.mean(np.array(var_list))
# print m_corr_xz[1:,0].shape, range(1,X/2)
 
# correlation plots with first data point removed
plt.clf()
# plt.errorbar(range(1, X/2), m_corr_xz[1:,0], d_corr_xz[1:,0], c='b', label="X", linewidth=2)
plt.errorbar(range(0, X/2+1), m_corr_xz[:,0], d_corr_xz[:,0], c='b', label="X", linewidth=2)
plt.hlines(0, 0, X/2+1, linestyles="dashed")
plt.xlim([0.0,X/2+1])
plt.ylim(-0.1,1.1)
# plt.errorbar(range(1, Z/2), m_corr_xz[0,1:], d_corr_xz[0,1:], c='g', label="Z", linewidth=2)
plt.errorbar(range(0, Z/2+1), m_corr_xz[0,:], d_corr_xz[0,:], c='g', label="Z", linewidth=2)
plt.hlines(0, 0, Z/2+1, linestyles="dashed")
plt.xlabel('x or z', fontsize=30)
plt.ylabel('G(x, z)', fontsize=30)
plt.legend(loc='upper right', fontsize=30)
plt.show()






