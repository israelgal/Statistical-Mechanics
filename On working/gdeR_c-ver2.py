import scipy as sp
import numpy as np
import time
import ctypes

from matplotlib import pyplot as plt
from numpy.ctypeslib import ndpointer 



 
class RDF_obj:
    def __init__(self, filename, X, Y, Z, resolution):
        
         
        with open(filename, 'r') as f:
            data = f.readlines()
        
        self.x = X
        self.y = Y
        self.z = Z

        self.n_atoms = int(data[0].split()[0])
        self.coordinates = np.zeros((self.n_atoms, 3),dtype=np.float32)
        self.resolution = resolution

        for j, line in enumerate(data[2 : self.n_atoms + 2]):
            self.coordinates[j, :] = [float(value) for value in line.split()]
         
        self.density_number()

    def volume(self,r):
        """ volumen de una esfera de radio r """
        volume = 4.0 / 3.0 * sp.pi * r**3
        return volume
     
    def distance(self,a, b):
        """ distancia minima entre dos particulas, considerando las dimensiones de la celda primaria """
        dx = abs(a[0] - b[0])
        x = min(dx, abs(self.x - dx))
         
        dy = abs(a[1] - b[1])
        y = min(dy, abs(self.y - dy))
         
        dz = abs(a[2] - b[2])
        z = min(dz, abs(self.z - dz))
         
        return np.sqrt(x**2 + y**2 + z**2)
     
    def density_number(self):
        """ calcula la densidad numérica"""
        self.dn = self.n_atoms /(self.x * self.y * self.z)  
     
    def compute_rdf(self):
        """ el radio de corte es la mitad de la longitud minima de las dimensiones de la celda """
        r_cutoff = min(min(self.x, self.y),self.z) / 2.0
        dr = r_cutoff / self.resolution
        volumes = np.zeros(self.resolution)
         
        self.radii = np.linspace(0.0, self.resolution * dr, self.resolution)
        self.rdf = np.zeros(self.resolution,dtype=np.float32)
         
        
        print('Calculando g(r) para {:4d} particulas...'.format(self.n_atoms))
        start = time.time()
        
        """ corremos sobre cada par de particulas, calculamos su distancia, construimos un histograma
        cada par de particulas contribuye dos veces al valor del histograma """
        
        clibs = ctypes.cdll.LoadLibrary("./gdeR_libs.so")

        c_rdf = clibs.RDF
        c_rdf.restype =  None
        c_rdf.argtypes = [ndpointer(ctypes.c_float, flags="C_CONTIGUOUS"),
                ndpointer(ctypes.c_float, flags="C_CONTIGUOUS"),
                ctypes.c_int,
                ctypes.c_int,
                ctypes.c_int,
                ctypes.c_int,
                ctypes.c_int,
                ctypes.c_float]

        c_rdf(self.rdf, self.coordinates, self.n_atoms, self.x, self.y, self.z, self.resolution,dr)

         
        for j in range(self.resolution):
                r1 = j * dr
                r2 = r1 + dr
                v1 = self.volume(r1)
                v2 = self.volume(r2)
                volumes[j] += v2 - v1

        self.rdf = self.rdf/self.n_atoms
        """ normalizamos con respecto al volumen del cascaron esferico que pertene a cada radio """
        for i, value in enumerate(self.rdf):
            self.rdf[i] = value/ (volumes[i] * self.dn) 
        
        end = time.time()
        print("Tiempo total de computo: {:.3f} segundos".format(end - start))

    def plot(self, rdf_filename):
         
        plt.xlabel('r (Å)')
        plt.ylabel('g(r)')
        plt.plot(self.radii, self.rdf)
         
        if rdf_filename:
            plt.savefig(rdf_filename, dpi=300, bbox_inches='tight', format='pdf')
         
        plt.show()
 
 

 
particles_rdf = RDF_obj('coor100e.dat', 25, 25, 25,200)
particles_rdf.compute_rdf()
particles_rdf.plot("rdf.pdf")