import numpy as np
import matplotlib.pyplot as plt
import scipy as sp

""" Abre documento .cif y lee sus lineas """

with open('6awc.cif') as f:
    data = f.readlines()

""" Crea una lista y añade la información del átomo
     en cada entrada   """

Ats_inf = []

for line in range(len(data)):
    if 'ATOM' in data[line]:
        l_inf = data[line].split()[5]
        if l_inf == 'A' or l_inf == 'C' or l_inf == 'G' or l_inf == 'U':
            Ats_inf.append(data[line])

""" Creo un archivo con los datos de los átomos """

with open('6awc_raw.pdb','w') as g:
    for line in Ats_inf:
        g.write(line)

""" Datos de las masas y el número de átomos totales """

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
    return (X, Y, Z)/M


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
        #print(stack)
        CM = CoM(stack)
        #print(stack)
        if ACGU == 1:
            A.append([CM[0], CM[1], CM[2]])
        elif ACGU == 2:
            C.append([CM[0], CM[1], CM[2]])
        elif ACGU == 3:
            G.append([CM[0], CM[1], CM[2]])
        elif ACGU == 4:
            U.append([CM[0], CM[1], CM[2]])

        if l == 0:
            #print(i, j, ACGU)
            l = l+1

    stack = []


A = np.array(A)
C = np.array(C)
G = np.array(G)
U = np.array(U)

N_A = len(A)
N_C = len(C)
N_G = len(G)
N_U = len(U)

print(N_A, N_U)

Concatenate = np.concatenate( ( A , U ) )

Max_pos = np.amax( Concatenate, axis = 0 )
Min_pos = np.amin( Concatenate, axis = 0 )

X_max = Max_pos[0] - Min_pos[0]
Y_max = Max_pos[1] - Min_pos[1]
Z_max = Max_pos[2] - Min_pos[2]

print(X_max, Y_max, Z_max)

##############################################################################################


def reduced_box( species_A, species_B , N_A , N_B , percent ):

    if int(percent) == 1:
        return species_A, species_B

    copyA = np.copy( species_A )
    copyB = np.copy( species_B )
    Iter_N = int( ( N_A + N_B ) * ( 1 - percent ) )
    #print(Iter_N)

    for i in range(Iter_N):
        SAmax = np.max( copyA )
        SBmax = np.max( copyB )
        SAmin = np.min( copyA )
        SBmin = np.min( copyB )

        Max_value =  max( abs( max(  SAmax , SBmax ) ), abs( min(  SAmin , SBmin ) ) )
        #print(Max_value)

        if Max_value == abs( SAmax ):
            index = np.where( copyA == SAmax )[0][0]
            copyA = np.delete( copyA , index , axis = 0)
        elif Max_value == abs( SBmax ):
            index =  np.where( copyB == SUmax )[0][0]
            copyB = np.delete( copyB , index , axis = 0)
        elif Max_value == abs( SAmin ):
            index =  np.where( copyA == SAmin )[0][0]
            copyA = np.delete( copyA , index , axis = 0)
        else:
            index = np.where( copyB == SBmin)[0][0]
            copyB = np.delete( copyB , index, axis = 0)

    return copyA, copyB


def distance(a, b):
    """ Calculo de la distancia  """

    dx = abs(a[0] - b[0])
    x = min(dx, abs(X_max - dx))

    dy = abs(a[1] - b[1])
    y = min(dy, abs(Y_max - dy))

    dz = abs(a[2] - b[2])
    z = min(dz, abs(Z_max - dz))

    return np.sqrt(x**2 + y**2 + z**2)

def density_number(N_A, N_B):
    """ calcula la densidad numérica"""
    dn = N_A * N_B /(X_max * Y_max * Z_max)

    return dn

def volume(r):

    volume = ( 4.0 * sp.pi * r**3) / 3.0
    return volume

def compute_rdf(Species_A, Species_B, resolution):
    """ el radio de corte es la mitad de la longitud minima de las dimensiones de la celda """

    #resolution = 300
    N_A = len(Species_A)
    N_B = len(Species_B)
    N_species = N_A + N_B

    r_cutoff =  15 #min( min(X_max, Y_max ), Z_max ) / 2.0
    dr = r_cutoff / resolution
    volumes = np.zeros(resolution)

    radii = np.linspace(0.0, resolution * dr, resolution)
    rdf = np.zeros((int(len(Species_A)), resolution))

    #print('Calculando g(r) para {:4d} particulas...'.format(N_atoms))
    #start = time.time()

    """ corremos sobre cada par de particulas, calculamos su distancia, construimos un histograma
    cada par de particulas contribuye dos veces al valor del histograma """

    for i, part_1 in enumerate(Species_A):

        for j, part_2 in enumerate(Species_B):

            dist = distance(part_1, part_2)
            index = int(dist / dr)
            if 0 < index < resolution:
                rdf[i,index] += 2.0

    for j in range(resolution):
        r1 = j*dr
        r2 = r1 + dr
        v1 = volume(r1)
        v2 = volume(r2)
        volumes[j] += v2 -v1

        #rdf = rdf / N_species

    rdf = np.mean(rdf, axis = 0)


    for i, value in enumerate(rdf):
        rdf[i] = value/ (volumes[i] * density_number(N_A,N_B))

    return radii, rdf

def plotrdf( species_A, species_B, resolution, species_name):

    plt.xlabel( 'r (Å)')
    plt.ylabel(  'g(r)' )

    radii, rdf = compute_rdf(species_A, species_B, resolution)

    plt.plot(radii, rdf, label = species_name)
    plt.legend()
    plt.savefig( species_name + '.pdf' , dpi=resolution, bbox_inches='tight', format='pdf')

def plot3d(species_A, species_B, species_name):

    X_1 = species_A[:,0]
    Y_1 = species_A[:,1]
    Z_1 = species_A[:,2]

    X_2 = species_B[:,0]
    Y_2 = species_B[:,1]
    Z_2 = species_B[:,2]

    fig3d = plt.figure( figsize = ( 10, 10) )
    ax = plt.axes( projection = '3d')
    ax.grid()

    ax.scatter( X_1, Y_1, Z_1, c = 'r', s = 50)
    ax.scatter( X_2, Y_2, Z_2, c = 'b', s = 50)
    ax.set_title( 'Nucleotidos ' + species_name )

    ax.set_xlabel('x', labelpad = 20)
    ax.set_ylabel('y', labelpad = 20)
    ax.set_zlabel('z', labelpad = 20)

    #plt.show()
    #plt.savefig(species_name + '_3d' + '.png')






compute_rdf(A,U, 100)
plotrdf(A, U, 100, 'AU')

#plotrdf(C, G, 300, 'CG')
#plot3d(A, U, 'AU')
