import numpy as np

with open('prueba.txt') as f:
	data = f.readlines()


	#print('linea 0', data[0].split()[0])
	print(len(data))

coordenadas = np.zeros((4, 4))

chain = []

for line in range(4):
	if data[line].split()[0] == 'a':
		chain.append(1)
	elif data[line].split()[0] == 'c':
		chain.append(2)
	elif data[line].split()[0] == 'g':
		chain.append(3)
	elif data[line].split()[0] == 'u':
		chain.append(4)
	else:
		chain.append(0)

#print(chain)

for j, line in enumerate(data[0 : 4]):
			coordenadas[j, 0] = chain[j] 
			coordenadas[j, 1:] = [float(value) for value in line.split()[1:]]
			
print(coordenadas)

#print(coordenadas[j,0])
#coordenadas[0,1] = [1: float(value) for value in line.split()[1:]]

"""for line in enumerate(data[0:4]):
	for j in enumerate(data[0:4]):
		if data[j:0] == 1:
			coordenadas[j, -1]"""
