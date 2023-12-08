import matplotlib.pyplot as plt
import numpy as np

dados = np.loadtxt("dados.txt", delimiter='\t')

load = dados[:,1]
depth = dados[:,2]

ymin, xmin, ymax, xmax = list(), list(), list(), list()
f0 = load[0]

for i in range(0, len(load)):

    f = load[i]
    a = f - f0
    if(a > 0  and load[i] < load[i-1]):
        ymax.append(load[i-1])
        xmax.append(depth[i-1])
        f0 = load[i-1]
    if(a < 0  and load[i] > load[i-1]):
        ymin.append(load[i-1])
        xmin.append(depth[i-1])
        f0 = load[i-1]
ymin.append(load[len(load)-1])
xmin.append(depth[len(load)-1])  

plt.plot(depth, load)
plt.title('Gráfico de Força por Profundidade')
plt.scatter(xmin, ymin, color="r", marker="D")
plt.scatter(xmax, ymax, color="g", marker="D")
plt.xlabel('h (µm)')
plt.ylabel('F (N)')
plt.grid(True)
plt.show()
