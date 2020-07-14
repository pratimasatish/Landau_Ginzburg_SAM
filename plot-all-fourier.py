import numpy as np
import matplotlib.pyplot as plt

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rc('font', size=28)
plt.rc('font', weight='bold')
plt.rc('xtick', labelsize=28)
plt.rc('ytick', labelsize=28)

data_short = np.genfromtxt('eta_var_short.txt')
data_long = np.genfromtxt('eta_var_long.txt')
data_noflucs = np.genfromtxt('eta_var_noflucs.txt')

N = 11 * 11
temp = 300
kB = 1.3806503 * 6.0221415 / 4184.0
beta = 1/(kB * temp)
a = 0.305
kappa = 0.321

norm_k = data_short[:,0]
eta_var_short = data_short[:,1]
eta_var_long = data_long[:,1]
eta_var_noflucs = data_noflucs[:,1]
norm_k_unique = np.unique(data_short[:,0])
eta_var_short_unique = np.zeros(len(norm_k_unique))
eta_var_long_unique = np.zeros(len(norm_k_unique))
eta_var_noflucs_unique = np.zeros(len(norm_k_unique))

for i in range(len(norm_k_unique)):
    eta_var_short_unique[i] = np.mean( eta_var_short[ np.where( np.abs(norm_k-norm_k_unique[i]) < 1.0e-2) ] )
    eta_var_long_unique[i] = np.mean( eta_var_long[ np.where( np.abs(norm_k-norm_k_unique[i]) < 1.0e-2) ] )
    eta_var_noflucs_unique[i] = np.mean( eta_var_noflucs[ np.where( np.abs(norm_k-norm_k_unique[i]) < 1.0e-2) ] )

predicted_var_flucs = 0.5 * N / (beta * (a + kappa * norm_k_unique**2))
predicted_var_noflucs = 0.5 * N / (beta * a) * np.ones(len(norm_k_unique))

plt.plot(norm_k_unique, predicted_var_flucs, 'k', linestyle='dashed', linewidth=3, label=r'$\kappa \neq 0$')
plt.plot(norm_k_unique, predicted_var_noflucs, 'b', linestyle='dashed', linewidth=3, label=r'$\kappa = 0$')

plt.plot(norm_k_unique, eta_var_short_unique, 'ro', label='short run')
plt.plot(norm_k_unique, eta_var_short_unique, 'r', linewidth=4, alpha=0.7)

plt.plot(norm_k_unique, eta_var_long_unique, 'gv', label='long run')
plt.plot(norm_k_unique, eta_var_long_unique, 'g', linewidth=4, alpha=0.7)

plt.plot(norm_k_unique, eta_var_noflucs_unique, 'm^', label='no fluctuations')
plt.plot(norm_k_unique, eta_var_noflucs_unique, 'm', linewidth=4, alpha=0.7)

plt.legend(loc='best', fontsize=24)
plt.xlabel(r'$|k|$')
plt.ylabel(r'$\langle\hat\eta_k^2\rangle$')
plt.show()




