import matplotlib.pyplot as plt
import numpy as np
import plot_lib as pl

def plot_Iav(filename, title):
    E, Iav = np.loadtxt(filename, skiprows=1, usecols=(0, 1), unpack=True)
    
    plt.plot(E, Iav)

    plt.xlabel(r'E [TeV]')
    plt.xscale('log')
    plt.ylabel(r'I$_e$ [cm$^{-3}$]')

    plt.title(title, fontsize = 25)
    
#    plt.legend(['YMW'])
    
    return 'Iav.pdf'

def plot_Pag(filename, title):
    E, P = np.loadtxt(filename, skiprows=1, usecols=(0, 2), unpack=True)

    plt.plot(E, P, 'r')

    plt.xlabel(r'E [TeV]')
    plt.xscale('log')
    
    plt.ylabel(r'P []')

    plt.title(title, fontsize = 25)
    
    return 'Pag.pdf'

pl.set_plot_style()

filename = 'output/HESSJ1640-465.txt'

#plotname = plot_Iav(filename, 'm$_a$ = 10$^{-9}$ eV - g$_{ag}$ = 10$^{-10}$ GeV$^{-1}$')  
plotname = plot_Pag(filename, 'm$_a$ = 10$^{-9}$ eV - g$_{ag}$ = 10$^{-10}$ GeV$^{-1}$')  

#plt.show()
plt.savefig(plotname)
