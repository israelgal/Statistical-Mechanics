import tensorflow as tf 
import matplotlib.pyplot as plt 
from matplotlib.pyplot import figure 
import numpy as np

np.random.seed(123)
tf.random.set_seed(123)

#condicion inicial 
u0 = 2
#infinitesimal 
inf_s = np.sqrt(np.finfo(np.float32).eps)

#parametros 
learning_rate = 0.01
training_steps = 500

display_step = training_steps/10

#parametros de NN
n_input = 1 
n_hidden_1 = 32
n_hidden_2 = 32 
n_output = 1 

weights = {
	'h1' : tf.Variable(tf.random.normal([n_input, n_hidden_1])),
	'h2' : tf.Variable(tf.random.normal([n_hidden_1,n_hidden_2])),
	'out' : tf.Variable(tf.random.normal([n_hidden_2, n_output]))
}

biases = {
	'b1' : tf.Variable(tf.random.normal([n_hidden_1])),
	'b2' : tf.Variable(tf.random.normal([n_hidden_2])),
	'out' : tf.Variable(tf.random.normal([n_output]))
}

#optimizador 
optimizador = tf.optimizers.SGD(learning_rate)

#modelo de la red 

def NN(x):
	x = np.array([[[x]]], dtype = 'float32')

	layer_1 = tf.add(tf.matmul(x,weights['h1']), biases['b1'])
	layer_1 = tf.nn.sigmoid(layer_1)

	layer_2 = tf.add(tf.matmul(layer_1,weights['h2']), biases['b2'])
	layer_2 = tf.nn.sigmoid(layer_2)

	output = tf.matmul(layer_2,weights['out']) + biases['out']

	return output

#aproximador 

def g(x):
	return x*NN(x) + u0

def f(x):
	return 3*x*x

#funcion de perdida 

def perdida():
	summation = []
	for x in np.linspace(0,1,10):
		dNN = (g(x + inf_s) - g(x - inf_s))/(2*inf_s)
		summation.append((dNN-f(x)))
	return tf.reduce_sum(tf.abs(summation)**2)#aquiborar**2


def entrenamiento():
	with tf.GradientTape() as tape:
		loss = perdida()
	variables_entrenamiento = list(weights.values()) + list(biases.values())
	gradientes = tape.gradient(loss,variables_entrenamiento)
	optimizador.apply_gradients(zip(gradientes,variables_entrenamiento))

for i in range(training_steps):
	entrenamiento()
	if i % display_step == 0:
		print('funcion de perdida: %f' % (perdida()))

def analitica(x):
	return x**3 + 2

figure(figsize = (10,10))

X = np.linspace(0,1,100)
S = analitica(X)

result = []

for i in X:
	result.append(g(i).numpy()[0][0][0])

plt.plot(X,S,label = 'Solucion analitica')
plt.plot(X,result, label = 'Solucion aproximada')
plt.legend(loc=2, prop={'size' :20})
plt.savefig('figura.pdf')
plt.show()