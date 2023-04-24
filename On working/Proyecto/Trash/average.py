import numpy
masses = numpy.array([[0,  0,  0,  0],
[0,  1,  0,  0],
[0,  2,  0,  0],
[1,  0,  0,  0],
[1,  1,  0,  1],
[1,  2,  0,  1],
[2,  0,  0,  0],
[2,  1,  0,  0],
[2,  2,  0,  0]])

nonZeroMasses = masses[numpy.nonzero(masses[:,3])] # Not really necessary, can just use masses because 0 mass used as weight will work just fine.

CM = numpy.average(nonZeroMasses[:,:3], axis=0, weights=nonZeroMasses[:,3])
print(CM)
