import scipy as sp
import numpy as np
from matplotlib import pyplot as plt
 
def volume(r):
    """ volumen de una esfera de radio r """
    volume = 4.0 / 3.0 * sp.pi * r**3
    return volume

def distance(a,b);
	
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

		self.n_atoms = int(data[0].split()[0]) # cadena a la cual la dividimos y obtenemos una lista, el valor 0 del indice
		self.coordinates = np.zeros((self.n_atoms, 3))
		self.resolution = resolution

		for j, line in enumerate(data[2 : self.n_atoms + 2]): #seleccionamos la linea dos de los datos 
            self.coordinates[j, :] = [float(value) for value in line.split()] #divide cada línea de los archivos como si fueran listas nrng
         
        self.density = density_number()

    def density_number(self)

    	self.dn = self.x * self.y * self.z / self.n_atoms



