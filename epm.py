from __future__ import division
import itertools as it
import numpy as np
import scipy.linalg as sl
import matplotlib.pyplot as plt
import csv
import sys
from zdict import zdict

# hbar^2/2m in eV-Ang^2
hb2m0 = 3.807

# hcR_infty (eV per Rydberg)
Ry = 13.60569253

def mag2(V):
    """ Return the magnitude squared of a tuple, list, or array """
    return sum([v**2 for v in V])

class EPM(object):
    def __init__(self,m,a0,VS,VA,bands):

        self.a0 = a0
        self.bands = bands        

        # Range of Fourier modes
        d = range(-m,1+m)

        # Construct all G-vectors in reciprocal lattice basis
        G = [np.array((-i+j+k,i-j+k,i+j-k)) for i in d for j in d for k in d]

        # Restrict to vectors of squared magnitude 12 or less
        self.G = [g for g in G if mag2(g) < 13]

        # Number of G vectors 
        ng = len(self.G)

        # Assemble potential part of Hamiltonian
        self.H = np.zeros((ng,ng),dtype=complex)       

        # Loop over all pairs of G vectors
        for ii in range(ng):
            for jj in range(ii,ng):
 
                # Difference between two G vectors
                dg = self.G[jj]-self.G[ii]
                dgmag2 = mag2(dg)

                # Dot product of dg and tau
                theta = np.pi*sum(dg)/4 
                c = np.cos(theta)
                s = np.sin(theta)

                self.H[ii,jj] = (VS[dgmag2]*c-1j*VA[dgmag2]*s)*Ry
                self.H[jj,ii] = self.H[ii,jj].conj()


    def solve(self,k):

        # Incorporate kinetic (main diagonal) part of Hamiltonian
        kpg2 = np.array([mag2(k-g) for g in self.G])
        kinetic = hb2m0*kpg2*(2*np.pi/self.a0)**2 

        # Insert diagonal elements
        np.fill_diagonal(self.H,kinetic)  

        # Calculate eigenvalues of a Hermitian matrix
        E = sl.eigvalsh(self.H)[:self.bands]
         
        return E


if __name__ == '__main__':

    # Name of semiconductor, e.g. Si, GaAs, InP, ZnS...
    material = sys.argv[1]   

    reader = csv.reader(open('form_factors.csv','r'),delimiter=',')
    param = [[entry.split()[0] for entry in row] for row in reader]
 
    # Store form factors in dictionaries
    VS = zdict()
    VA = zdict()
    a0 = None

    # Read form factors and lattice constant from file
    row = [p[0] for p in param].index(material)
    for i in range(1,len(param[0])):
        exec(param[0][i] + '=' + param[row][i])

    # Symmetry points in the FCC/Zincblende Brillouin zone
    bz = {r'$\Gamma$':np.array((0,0,0)),
          r'X':np.array((0,1,0)),
          r'L':np.array((1/2,1/2,1/2)),
          r'W':np.array((1/2,1,0)),
          r'U':np.array((1/4,1,1/2)),
          r'K':np.array((3/4,3/4,0))}

    # Follow this path through the Brillouin zone to construct
    # the band diagram
    path = [r'L',r'$\Gamma$',r'X',r'W',r'K',r'$\Gamma$']

    path_dex = range(len(path)-1)

    # Highest Fourier mode to use
    fmodes = 3

    # Number of energy bands to compute
    bands = 8 
 
    # Number of k-point along each path to evaluate
    kpts = 40

    # k-space path parametric variable
    t = np.linspace(0,1,kpts)

    # Construct solver object
    epm = EPM(fmodes,a0,VS,VA,bands)

    # Sequence of path directions in k-space
    kdir = np.diff(np.vstack([bz[p] for p in path]),n=1,axis=0)
    
    # Lengths of k-space path segments
    path_length = [np.sqrt(mag2(k)) for k in kdir]

    # Relative positions of k-space symmetry points along x axis
    xticks = np.cumsum([0]+path_length)
    x=np.hstack([xticks[i]*(1-t)+xticks[i+1]*t for i in path_dex])
    
    # Parameterize path between two Brilluoin zone symmetry points
    K = lambda d: (np.outer((1-t),bz[path[d]])+np.outer(t,bz[path[d+1]]))

    # Compute eigenvalues along k-space path
    E = np.vstack([np.vstack([epm.solve(k) for k in K(j)]) for j in path_dex])
    
    Emin = np.min(E)-1
    Emax = np.max(E)+1
  
    # Display E-k diagram
    fig = plt.figure(1,(10,6))   
    plt.plot(x,E,'r-',lw=2)
    plt.xlim(x[0],x[-1])
    plt.ylim(Emin,Emax)
    plt.xticks(xticks,path,fontsize=20)
    plt.ylabel('Energy (eV)',fontsize=20)
    plt.title(material + ' bandstructure by EPM without S-O')
    plt.vlines(xticks,Emin,Emax)
    plt.show()
    




