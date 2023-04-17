import numpy as np

""" Abre documento .cif y lee sus lineas """

with open('1esh.cif') as f:
    data = f.readlines()


""" Crea una lista y añade la información del átomo
     en cada entrada   """

Ats_inf = []

for line in range(len(data)):
    if 'ATOM' in data[line]:
        if data[line].split()[5] == 'A' or 'C' or 'G' or 'U':
            Ats_inf.append(data[line])


with open('1esh_raw.pdb','w') as g:
    for line in Ats_inf:
        g.write(line)

N_atoms = len(Ats_inf)

Hidrogen = 1.00784
Oxigen = 15.999
Nitrogen = 14.0067
Carbon = 12.011
Fosforo = 30.973762


Matrix = np.zeros( ( N_atoms , 6 ) )

""" Función que calcula el centro de masa de los nucleotidos """

def CoM(stack):
    M = np.sum(stack[:,0])
    X = 0
    Y = 0
    Z = 0
    for i in range(len(stack)):
        X += stack[i,0] * stack[i,1]
        Y += stack[i,0] * stack[i,2]
        Z += stack[i,0] * stack[i,3]
    return [X, Y, Z]/M


""" Se sustituyen los pesos atómicos de H, O, N, C, P
        y además se sustituye ACGU por 1234 """

for j, line in enumerate(Ats_inf):
    if  line.split()[2] == 'H':
        Matrix[j, 0] = Hidrogen
    if  line.split()[2] == 'O':
        Matrix[j, 0] = Oxigen
    if  line.split()[2] == 'N':
        Matrix[j, 0] = Nitrogen
    if  line.split()[2] == 'C':
        Matrix[j, 0] = Carbon
    if  line.split()[2] == 'P':
        Matrix[j, 0] = Fosforo
    if  line.split()[5] == 'A':
        Matrix[j, 1] = 1
    if  line.split()[5] == 'C':
        Matrix[j, 1] = 2
    if  line.split()[5] == 'G':
        Matrix[j, 1] = 3
    if  line.split()[5] == 'U':
        Matrix[j, 1] = 4
    Matrix[j,2:] = [line.split()[8],
            line.split()[10], line.split()[11],line.split()[12]]


""" Encontrar las dimensiones de los nucleotidos """

x = 0
y = 0
z = 0

for i in range(N_atoms):

    if x < abs(Matrix[i,3]):
        x = abs(Matrix[i,3])
    if y < abs(Matrix[i,4]):
        y = abs(Matrix[i,4])
    if z < abs(Matrix[i,5]):
        z = abs(Matrix[i,5])
print(x,y,z)

A = []
C = []
G = []
U = []
stack = []

""" Calcula de los centros de masas y de guardar los en listas para cada ACGU """

for i in range(1, int(Matrix[N_atoms-1,2]) + 1):
    for j in range(len(Matrix)):
        if int(Matrix[j,2]) == i:
            stack.append([Matrix[j,0], Matrix[j,3], Matrix[j,4], Matrix[j,5]])
            ACGU = int(Matrix[j,1])

    stack = np.array(stack)
    CM = CoM(stack)

    if ACGU == 1:
        A.append(CM)
    if ACGU == 2:
        C.append(CM)
    if ACGU == 3:
        G.append(CM)
    if ACGU == 4:
        U.append(CM)
    stack = []


def distance(a, b):
        """ distancia minima entre dos particulas, considerando las dimensiones de la celda primaria """
        dx = abs(a[0] - b[0])
        x = min(dx, abs(x - dx))

        dy = abs(a[1] - b[1])
        y = min(dy, abs(y - dy))

        dz = abs(a[2] - b[2])
        z = min(dz, abs(z - dz))

def density_number():
        """ calcula la densidad numérica"""
        dn = N_atoms /(x * y * z)

def volume(r):

        volume = ( 4.0 * sp.pi * r**3) / 3.0
        return volume

def compute_rdf(SA, SB):
        """ el radio de corte es la mitad de la longitud minima de las dimensiones de la celda """
        r_cutoff = min( min(x, y), z ) / 2.0
        dr = r_cutoff / resolution
        volumes = np.zeros(resolution)

        radii = np.linspace(0.0, resolution * dr, resolution)
        rdf = np.zeros(resolution)


        print('Calculando g(r) para {:4d} particulas...'.format(N_atoms))
        start = time.time()


        """ corremos sobre cada par de particulas, calculamos su distancia, construimos un histograma
        cada par de particulas contribuye dos veces al valor del histograma """


        for i, part_1 in enumerate(SA):


            for j, part_2 in enumerate(SB):

                dist = distance(part_1, part_2)

                index = int(dist / dr)
                if 0 < index < self.resolution:
                    self.rdf[index] += 2.0


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

resolution = 1
