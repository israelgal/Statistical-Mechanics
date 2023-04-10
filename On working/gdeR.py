import scipy as sp
import numpy as np
from matplotlib import pyplot as plt
 
def volume(r):
    """ volumen de una esfera de radio r """
    volume = 4.0 / 3.0 * sp.pi * r**3
    return volume
 
def distance(a, b):
    """ distancia minima entre dos particulas, considerando las dimensiones de la celda primaria """
    dx = abs(a[0] - b[0])
    x = min(dx, abs(A - dx))
     
    dy = abs(a[1] - b[1])
    y = min(dy, abs(B - dy))
     
    dz = abs(a[2] - b[2])
    z = min(dz, abs(C - dz))
     
    return np.sqrt(x**2 + y**2 + z**2)
 
class RDF_obj:
    def __init__(self, filename, X, Y, Z, resolution):
        
         
        with open(filename, 'r') as f:
            data = f.readlines()
        
        self.x = X
        self.y = Y
        self.z = Z

        self.n_atoms = int(data[0].split()[0])
        self.coordinates = np.zeros((self.n_atoms, 3))
        self.resolution = resolution

        for j, line in enumerate(data[2 : self.n_atoms + 2]):
            self.coordinates[j, :] = [float(value) for value in line.split()]
         
        self.density_number()
     
    def density_number(self):
        """ calcula la densidad numérica"""
        self.dn = self.x * self.y * self.z  / self.n_atoms
     
    def compute_rdf(self):
        """ el radio de corte es la mitad de la longitud minima de las dimensiones de la celda """
        r_cutoff = min(min(A, B),C) / 2.0
        dr = r_cutoff / self.resolution
        volumes = np.zeros(self.resolution)
         
        self.radii = np.linspace(0.0, self.resolution * dr, self.resolution)
        self.rdf = np.zeros(self.resolution)
         
        
        print('Calculando g(r) para {:4d} particulas...'.format(self.n_atoms))
        
        """ corremos sobre cada par de particulas, calculamos su distancia, construimos un histograma
        cada par de particulas contribuye dos veces al valor del histograma """
        for i, part_1 in enumerate(self.coordinates):
            for j in range(self.resolution):
                r1 = j * dr
                r2 = r1 + dr
                v1 = volume(r1)
                v2 = volume(r2)
                volumes[j] += v2 - v1
             
            for part_2 in self.coordinates[i:]:
                dist = distance(part_1, part_2)
                index = int(dist / dr)
                if 0 < index < self.resolution:
                    self.rdf[index] += 2.0
         
        """ normalizamos con respecto al volumen del cascaron esferico que pertene a cada radio """
        for i, value in enumerate(self.rdf):
            self.rdf[i] = value * self.dn / volumes[i]
         
    def plot(self, rdf_filename):
         
        plt.xlabel('r (Å)')
        plt.ylabel('g(r)')
        plt.plot(self.radii, self.rdf)
         
        if rdf_filename:
            plt.savefig(rdf_filename, dpi=300, bbox_inches='tight', format='pdf')
         
        plt.show()
 
 
A = 25
B = 25
C = 25

 
particles_rdf = RDF_obj('coor100e.dat', A, B, C,200)
particles_rdf.compute_rdf()
particles_rdf.plot("rdf.pdf")
plt.savefig('try.pdf')