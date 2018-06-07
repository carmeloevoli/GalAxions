import matplotlib.pyplot as plt
import numpy as np
import plot_lib as pl

def plot_ne(filename, title):
    d, ne = np.loadtxt(filename, skiprows=1, usecols=(0, 4), unpack=True)
    
    plt.plot(d, ne)

    plt.xlabel(r'd [kpc]')
    plt.ylabel(r'n$_e$ [cm$^{-3}$]')
    plt.title(title, fontsize = 25)
    
    plt.legend(['YMW'])
    
    return 'ne_los.pdf'

def plot_B(filename, title):
    d, B_perp, B_tot = np.loadtxt(filename, skiprows=1, usecols=(0, 1, 2), unpack=True)

    plt.plot(d, B_perp, 'r')
    plt.plot(d, B_tot, 'r:')

    plt.xlabel(r'd [pc]')
    plt.ylabel(r'B [$\mu$G]')
    plt.title(title, fontsize = 25)

    return 'B_los.pdf'

pl.set_plot_style()

filename = 'output/HESSJ1640-465.los'

#plotname = plot_ne(filename, 'HESSJ1640-465')
plotname = plot_B(filename, 'HESSJ1640-465')

#plt.show()
plt.savefig(plotname)
