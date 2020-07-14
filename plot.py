import numpy as np
import matplotlib.pyplot as plt

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rcParams['text.latex.preamble']=[r"\usepackage{amsmath}"]
plt.rc('font', size=28)
plt.rc('xtick', labelsize=28)
plt.rc('ytick', labelsize=28)

temp_list = []
# for i in range(365, 387, 2):
# for i in range(200, 400, 50):
for i in range(320, 425, 5):
    temp_list.append(i)
temp_list = np.array(temp_list)
# print temp_list

data_list = []
for i in range(len(temp_list)):
    data_list.append( np.genfromtxt('log-{}.txt'.format(temp_list[i]), skip_header=1, delimiter='  ') )
data_list = np.array(data_list)
eta_list = np.copy(data_list[:,:,6])

mean_eta_vs_temp = np.zeros(len(temp_list))
for i in range(len(temp_list)):
    mean_eta_vs_temp[i] = np.mean(eta_list[i][5000:])
# print mean_eta_vs_temp

plt.figure(0)
plt.plot(temp_list, mean_eta_vs_temp, 'bv')
plt.plot(temp_list, mean_eta_vs_temp, 'b', linewidth=4, alpha=0.6)
plt.xlabel('T(K)', fontsize=28)
plt.ylabel(r'$\langle\eta\rangle$', fontsize=28)
plt.title(r'$\langle\eta\rangle(T) \textrm{ for } \mathcal{H} = h(T-T_0)\eta + \dfrac{a}{2}\eta^2 + \dfrac{c}{4}\eta^4 + \dfrac{\kappa}{2}|\nabla\eta|^2')
plt.show()

temps = np.array( [365, 367, 369, 371, 373, 375, 377, 379, 381, 383, 385] )
temps = np.array( [200, 250, 300, 350] )
temps = np.array( [320, 330, 340, 350, 360, 365, 370, 375, 380, 390, 400, 410, 420] )

a = -31.08
c = 701.3
kappa = 32.71
kappa = 327.1
T0 = 375.6
h = -0.0488
def hamiltonian(eta, T):
    return h * (T - T0) * eta + a * eta**2 / 2.0 + c * eta**2 * eta**2 / 4.0

bins = np.linspace(-0.5, 0.5, 251)
eta_bins = 0.5 * (bins[1:] + bins[:-1])
# plt.plot(eta_bins, hamiltonian(eta_bins, 400), 'ro')
# plt.show()

for T in temps:
    index = np.where(temp_list == T)[0]
    current_eta = np.copy(eta_list[index].flatten()[5000:])
    eta_hist, bins = np.histogram(current_eta, bins=bins, density=True)
    bins_hist = eta_bins[eta_hist != 0]
    eta_hist = eta_hist[eta_hist != 0]

    plt.plot(bins_hist, -np.log(eta_hist), 'go', label='MC data at T = {}K'.format(T))
    plt.plot(bins_hist, -np.log(eta_hist), 'g', linewidth=4, alpha=0.8)
    plt.plot(eta_bins, hamiltonian(eta_bins, T), 'r', linewidth=4, alpha=0.5, label='hamiltonian at T = {}K'.format(T))
    plt.legend(loc='best')
    plt.show()



