import scipy as sp
import numpy as np
import time

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
         
        self.density_number()

    def volume(self,r):
        """ volumen de una esfera de radio r """
        volume = 4.0 / 3.0 * sp.pi * r**3
        return volume
     
    def distance(self,cx,cy,cz,px,py,pz):
        """ distancia minima entre dos particulas, considerando las dimensiones de la celda primaria """
        dx = np.absolute(px - cx)
        x = np.minimum(dx, abs(self.x - dx))
         
        dy = np.absolute(py - cy)
        y = np.minimum(dy, abs(self.y - dy))
         
        dz = np.absolute(pz - cz)
        z = np.minimum(dz, abs(self.z - dz))
         
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
        self.rdf = np.zeros(self.resolution)
         
        
        print('Calculando g(r) para {:4d} particulas...'.format(self.n_atoms))
        start = time.time()


        """ corremos sobre cada par de particulas, calculamos su distancia, construimos un histograma
        cada par de particulas contribuye dos veces al valor del histograma """

        vec_distance = np.vectorize(self.distance)
        
        
        for i, part_1 in enumerate(self.coordinates):

            indexes = (vec_distance(self.coordinates[i:,0],self.coordinates[i:,1],self.coordinates[i:,2],part_1[0],part_1[1],part_1[2])/dr).astype(int)
            indexes = indexes[(indexes > 0) & (indexes < self.resolution)]
            np.add.at(self.rdf,indexes,2)
            
                    

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