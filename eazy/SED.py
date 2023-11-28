import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
plt.style.use('classic')

c=299792458
k_B=1.380649e-23
h=6.62607015e-34

M_sun=1.98847e30
pc=3.085677581e16

def B_nu(nu,T):
    return 2*h*nu**3/c**2/(np.exp(h*nu/k_B/T)-1)

def tau_nu(nu):
    
    kappa_0=0.015
    nu_0=250e9
    beta=2
    
    M_d=1e10*M_sun
    
    return kappa_0*(nu/nu_0)**beta*M_d

def D_A(z):
    def E(z):
        return (0.315*(1+z)**3+0.685)**(-0.5)
    chi=quad(E,0,z)[0]
    return 2.998e8/(67.3*1e3/1e6/3.086e16)*chi/(1+z)

nu=10**np.linspace(np.log10(100*1e9),np.log10(5000*1e9),300)

#Model with perfect blackbody dust grains single temperature
T_d1=30
z1=1
S_nu_gal1=tau_nu(nu*(1+z1))*B_nu(nu*(1+z1),T_d1)/D_A(z1)**2/(1+z1)**3

T_d2=30
z2=2
S_nu_gal2=tau_nu(nu*(1+z2))*B_nu(nu*(1+z2),T_d1)/D_A(z2)**2/(1+z2)**3

T_d3=30
z3=4
S_nu_gal3=tau_nu(nu*(1+z3))*B_nu(nu*(1+z3),T_d1)/D_A(z3)**2/(1+z3)**3

fig=plt.figure()
ax=[]
ax.append(fig.add_subplot(111))
ax.append(ax[0].twiny())

ax[0].plot(nu,S_nu_gal1,linewidth=5,marker=None,label='$z = 1$',color='cornflowerblue')

ax[0].plot(nu,S_nu_gal2,linewidth=5,marker=None,label='$z = 2$',color='green')

ax[0].plot(nu,S_nu_gal3,linewidth=5,marker=None,label='$z = 4$',color='salmon')

ax[0].set_xscale('log')
ax[0].set_yscale('log')

ax[0].axis([100e9,5000e9,5e-30,8e-27])

ax[0].set_xlabel('$\\nu_{\mathrm{obs}}\,[\mathrm{GHz}]$',fontsize=25)
ax[0].set_ylabel('$S_{\\nu_{\mathrm{obs}}}\,[\mathrm{mJy}]$',fontsize=25)

ax[0].set_xticks([100e9,300e9,1000e9,3000e9])
ax[0].set_xticklabels(['$100$','$300$','$1000$','$3000$'],fontsize=20)
ax[0].set_yticks([1*1e-3*1e-26,10*1e-3*1e-26,100*1e-3*1e-26])
ax[0].set_yticklabels(['$1$','$10$','$100$'],fontsize=20)

ax[1].set_xscale('log')
ax[1].axis([100e9,5000e9,5e-30,8e-27])

ax[1].set_xlabel('$\lambda_{\mathrm{obs}}\,[\mu\mathrm{m}]$',fontsize=25, labelpad=15)

ax[1].set_xticks([3e8/3000e-6,3e8/2000e-6,3e8/1100e-6,3e8/850e-6,3e8/500e-6,3e8/350e-6,3e8/250e-6,3e8/100e-6])
ax[1].set_xticklabels(['$3000$','$2000$','$1100$','$850$','$500$','$350$','$250$','$100$'],fontsize=15)
ax[1].tick_params(axis='x', which='minor', top=False)

handles,labels=ax[0].get_legend_handles_labels()
lg=plt.legend(handles,labels,loc=2,fontsize=20,frameon=False)
ax[1].add_artist(lg)

plt.savefig('/System/Volumes/Data/Users/ptanouri/Desktop/Flux_example.png',dpi=300,bbox_inches='tight',pad_inches=0.1)



