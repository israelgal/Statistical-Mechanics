import numpy as np

""" Abre documento .cif y lee sus lineas """

with open('1esh.cif') as f:
    data = f.readlines()


""" Crea una lista y añade la información del átomo
     en cada entrada   """

Atoms_information = []

for line in range(len(data)):
    if 'ATOM' in data[line]:
        Atoms_information.append(data[line])

#list0 = Atoms_information[0].split()[:]
#print(list0)

#print(type(list0))
#print(Atoms_information[0], Atoms_information[len(Atoms_information)-1])

""" Crea un archivo de datos, con los datos de sus entradas línea por línea  """

with open('1esh.dat','w') as g:

    g.write(str(len(Atoms_information)))
    for line in Atoms_information:
        #print(type(line))
        #print(line)
        #data_splited = []
        #data_splited[0] = Atoms_information[line].split()[:]
        g.write('\n')
        g.write( line.split()[5] + ' '
         + line.split()[6] + ' ' + line.split()[10] + ' ' +
         line.split()[11] + ' ' + line.split()[12] )
        #g.write(line)



""" Masas de los aminoacidos y Ribonucleotidos en g/mol"""

# Adenine
A = 135.13

# Citosina / Cytosine
C = 111.1

#Guanina / Guanine
G = 151.13

#Uracil
U = 112.09

#Alanine
ALA = 89.09

#Arginine
ARG = 174.2

#Asparagine
ASN = 132.12

#Aspartic
ASP = 133.1

#Cysteine
CYS = 121.16

#Glutamine
GLN = 146.14

#Glutamic acid
GLU = 147.13

#Glycine
GLY = 75.07

#Histidine
HIS = 155.1546

#Isoleucine
ILE = 131.07

#Leucine
LEU = 131.17

#Lysine
LYS = 146.19

#Methionine
MET = 149.21

#Phenylalanine
PHE = 165.19

#Proline
PRO = 115.13

#Serine
SER = 105.09

#Threonine
THR = 119.1192

#Tryptophan
TRP = 204.23

#Tyrosine
TYR = 181.19

#Valine
VAL = 117.151

""" Leemos el archivo de datos creado anteriormente"""

with open('1esh.dat') as h:
    data_Mol = h.readlines()

#print(len(data_Mol))


# Definimos el número de átomos

n_atoms = len(data_Mol) - 1

""" Creamos una matriz 418x4 de ceros """

coordenadas = np.zeros((len(data_Mol) -1 , 4))

""" Creamos una lista de 418 elementos con los valores de las masas """

chain = []

for line in range(1, n_atoms + 1  ):
    if data_Mol[line].split()[0] == 'A':
        chain.append(A)
    elif data_Mol[line].split()[0] == 'C':
        chain.append(C)
    elif data_Mol[line].split()[0] == 'G':
        chain.append(G)
    elif data_Mol[line].split()[0] == 'U':
        chain.append(U)
    elif data_Mol[line].split()[0] == 'ALA':
        chain.append(ALA)
    elif data_Mol[line].split()[0] == 'ARG':
        chain.append(ARG)
    elif data_Mol[line].split()[0] == 'ASN':
        chain.append(ASN)
    elif data_Mol[line].split()[0] == 'ASP':
        chain.append(ASP)
    elif data_Mol[line].split()[0] == 'CYS':
        chain.append(CYS)
    elif data_Mol[line].split()[0] == 'GLN':
        chain.append(GLN)
    elif data_Mol[line].split()[0] == 'GLU':
        chain.append(GLU)
    elif data_Mol[line].split()[0] == 'GLY':
        chain.append(GLY)
    elif data_Mol[line].split()[0] == 'HIS':
        chain.append(HIS)
    elif data_Mol[line].split()[0] == 'ILE':
        chain.append(ILE)
    elif data_Mol[line].split()[0] == 'LEU':
        chain.append(LEU)
    elif data_Mol[line].split()[0] == 'LYS':
        chain.append(LYS)
    elif data_Mol[line].split()[0] == 'MET':
        chain.append(MET)
    elif data_Mol[line].split()[0] == 'PHE':
        chain.append(PHE)
    elif data_Mol[line].split()[0] == 'PRO':
        chain.append(PRO)
    elif data_Mol[line].split()[0] == 'SER':
        chain.append(SER)
    elif data_Mol[line].split()[0] == 'THR':
        chain.append(THR)
    elif data_Mol[line].split()[0] == 'TRP':
        chain.append(TRP)
    elif data_Mol[line].split()[0] == 'TYR':
        chain.append(TYR)
    elif data_Mol[line].split()[0] == 'VAL':
        chain.append(VAL)
    else:
        chain.append(0)
        print(data_Mol[line].split()[0])

#print(chain[-1:])

""" Llenamos la matriz creada de 418x4, en la primera columna están las masas
    y en las siguientes las coordenadas XYZ de la molécula """

for j, line in enumerate(data_Mol[1 :  n_atoms + 1 ]):
            coordenadas[j, 0] = chain[j]
            coordenadas[j, 1:] = [float(value) for value in line.split()[2:]]

""" Creamos una función que encuentra el centro de masa """

def CM(XYZ, mass):
    mass = mass.reshape((n_atoms,1))
    #print(mass)
    M_T = mass.sum()
    #print(mass[0])
    #print(XYZ*mass[0])
    #print(M_T)
    Vec = XYZ*mass



    #X = np.sum(Vec,axis = 0) / M_T
    #Y = np.sum(Vec,axis = 1) / M_T
    #Z = np.sum(Vec,axis = 2) / M_T
    return np.sum(Vec, axis = 0)/M_T

Com = CM(coordenadas[:,1:4], coordenadas[:,0])
print(Com)


""" Creamos un método iterativo para encontrar el centro de masa """

X_ = 0
Y_ = 0
Z_ = 0
Masa = 0

for j in range(0,n_atoms):
    X_ = X_ + coordenadas[j,1]*coordenadas[j,0]
    Y_ = Y_ + coordenadas[j,2]*coordenadas[j,0]
    Z_ = Z_ + coordenadas[j,3]*coordenadas[j,0]
    Masa = Masa + coordenadas[j,0]

#Masa = np.sum(coordenadas, axis = 0)[0]

#print(Masa)

#CM2 = [X_/Masa, Y_/ Masa, Z_/ Masa]
#print(X_)

CM2 = [X_, Y_, Z_] / Masa
print(CM2)









