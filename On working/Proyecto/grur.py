import scipy as sp
import numpy as np
from matplotlib import pyplot as plt


def Read_RNA_file( filename ):

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
        return (X, Y, Z)/M

    """ Encontrar las dimensiones de los nucleotidos """

    def Walls(Matrix_inf):

        X_max = 0
        Y_max = 0
        Z_max = 0

        for i in range(N_atoms):

            if X_max < abs(Matrix[i,3]):
                X_max = abs(Matrix[i,3])
            if Y_max < abs(Matrix[i,4]):
                Y_max = abs(Matrix[i,4])
            if Z_max < abs(Matrix[i,5]):
                Z_max = abs(Matrix[i,5])

        return X_max, Y_max, Z_max

    def distance(a, b):
    """ Calculo de la distancia considerando las dimensiones de X_max, Y_max y Z_max  """

        dx = abs(a[0] - b[0])
        x = min(dx, abs(X_max - dx))

        dy = abs(a[1] - b[1])
        y = min(dy, abs(Y_max - dy))

        dz = abs(a[2] - b[2])
        z = min(dz, abs(Z_max - dz))

    return np.sqrt(x**2 + y**2 + z**2)

    def density_number(N_A, N_B):

    """ calculo de densidad númerica """
        dn = N_A * N_B /(X_max * Y_max * Z_max)
        return dn

    def volume(r):

        volume = ( 4.0 * sp.pi * r**3) / 3.0
        return volume

#################################################################################

    with open( filename ) as f:
    data = f.readlines()

    Ats_inf = []

    for line in range(len(data)):
        if 'ATOM' in data[line]:
            l_inf = data[line].split()[5]
        if l_inf == 'A' or l_inf == 'C' or l_inf == 'G' or l_inf == 'U':
            Ats_inf.append(data[line])

    with open( filename + 'raw.cif','w') as g:
        for line in Ats_inf:
            g.write(line)

    """ Datos de las masas y el número de átomos totales """

    N_atoms = len(Ats_inf)

    Hidrogen = 1.00784
    Oxigen = 15.999
    Nitrogen = 14.0067
    Carbon = 12.011
    Fosforo = 30.973762

    """ Matriz con los datos de las masas, nucleotido, residuo, X, Y, Z """

    Matrix = np.zeros( ( N_atoms , 6 ) )

    """ Se sustituyen los pesos atómicos de H, O, N, C, P
            y además se sustituye ACGU por 1234 """

    for j, line in enumerate(Ats_inf):
        if  line.split()[2] == 'H':
            Matrix[j, 0] = Hidrogen
        elif  line.split()[2] == 'O':
            Matrix[j, 0] = Oxigen
        elif  line.split()[2] == 'N':
            Matrix[j, 0] = Nitrogen
        elif  line.split()[2] == 'C':
            Matrix[j, 0] = Carbon
        elif  line.split()[2] == 'P':
            Matrix[j, 0] = Fosforo
        else:
            print(0)
        if  line.split()[5] == 'A':
            Matrix[j, 1] = 1
        elif  line.split()[5] == 'C':
            Matrix[j, 1] = 2
        elif  line.split()[5] == 'G':
            Matrix[j, 1] = 3
        elif  line.split()[5] == 'U':
            Matrix[j, 1] = 4
        else:
            print(0)
        Matrix[j,2:] = [line.split()[8],
                line.split()[10], line.split()[11],line.split()[12]]


    A = []
    C = []
    G = []
    U = []
    stack = []

    """ Calcula de los centros de masas y de guardar los en listas para cada ACGU """
    l = 0
    for i in range(1, int(Matrix[N_atoms-1,2]) + 1):
        for j in range(len(Matrix)):
            if int(Matrix[j,2]) == i:
                stack.append([Matrix[j,0], Matrix[j,3], Matrix[j,4], Matrix[j,5]])
                ACGU = int(Matrix[j,1])

        if int(len(stack)) == 0:
            l = l + 1
        else:
            stack = np.array(stack)
            CM = CoM(stack)

            if ACGU == 1:
                A.append([CM[0], CM[1], CM[2]])
            elif ACGU == 2:
                C.append([CM[0], CM[1], CM[2]])
            elif ACGU == 3:
                G.append([CM[0], CM[1], CM[2]])
            elif ACGU == 4:
                U.append([CM[0], CM[1], CM[2]])


    stack = []

    A = np.array(A)
    C = np.array(C)
    G = np.array(G)
    U = np.array(U)
