import scipy as sp
import numpy as np
import time
from rdfpy import rdf

from matplotlib import pyplot as plt





 
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
         
        

    
     
    
     
    def compute_rdf(self):
        """ el radio de corte es la mitad de la longitud minima de las dimensiones de la celda """
        r_cutoff = min(min(self.x, self.y),self.z) / 2.0
        delta_r = r_cutoff / self.resolution
        
         
         
        
        print('Calculando g(r) para {:4d} particulas...'.format(self.n_atoms))
        start = time.time()


        """ corremos sobre cada par de particulas, calculamos su distancia, construimos un histograma
        cada par de particulas contribuye dos veces al valor del histograma """
        
        self.rdf, self.radii = rdf(self.coordinates, dr=delta_r)
        
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