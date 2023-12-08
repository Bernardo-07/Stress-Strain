import matplotlib.pyplot as plt
import numpy as np

def pontos_criticos(x, y):

    ymin, xmin, ymax, xmax = list(), list(), list(), list()
    indices_min, indices_max = list(), list()
    y0 = y[0]

    for i in range(0, len(y)):

        a = y[i] - y0

        if(a > 0  and y[i] < y[i-1]):
            y0 = y[i-1]
            ymax.append(y[i-1])
            xmax.append(x[i-1])
            indices_max.append(i-1)

        if(a < 0  and y[i] > y[i-1]):
            y0 = y[i-1]
            ymin.append(y[i-1])
            xmin.append(x[i-1])
            indices_min.append(i-1)

    ymin.append(y[len(y)-1])
    xmin.append(x[len(y)-1])  
    indices_min.append(len(y)-1)

    return ymin, xmin, ymax, xmax, indices_min, indices_max

def reta_inicial(x, y, xmax, ymax):
    yf = ymax[0]
    xf = xmax[0]
    a = (yf - y[0])/(xf - x[0])
    b = yf - a*xf
    return a, b

dados = np.loadtxt("dados.txt", delimiter='\t')

load = dados[:,1]
depth = dados[:,2]

Fmin, hmin, Fmax, hmax, i_min, i_max = pontos_criticos(depth, load)
coef_ang, coef_lin = reta_inicial(depth, load, hmax, Fmax)

newdepth = depth - (-coef_lin/coef_ang)

plt.plot(depth, load, color='b', label='Dados originais')
plt.plot(newdepth, load, label='Curva Calibrada')
plt.title('Gráfico de Força por Profundidade')
plt.scatter(hmin, Fmin, color="r", marker="D")
plt.scatter(hmax, Fmax, color="g", marker="D")
plt.xlabel('h (µm)')
plt.ylabel('F (N)')
plt.grid(True)
plt.show()
