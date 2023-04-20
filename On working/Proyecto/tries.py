import matplotlib.pyplot as plt
import numpy as np
#%matplotlib inline

t = np.arange(0.0, 2.0, 0.01)
s = 1 + np.sin(2*np.pi*t)
plt.plot(t, s)
print(5)

plt.title('About as simple as it gets, folks')
plt.show()
