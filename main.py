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

def coeficientes(x, y):
    coef = np.polyfit(x, y, 1)
    a = coef[0]
    b = coef[1]
    return a, b

def regressao_linear(x, y):
    coef = np.polyfit(x, y, 1)
    polinomio = np.poly1d(coef)
    y_fit = polinomio(x)
    return y_fit

dados = np.loadtxt("dados.txt", delimiter='\t')

load = dados[:,1]
depth = dados[:,2]

coef_ang, coef_lin = coeficientes(depth, load)
newdepth = depth - (-coef_lin/coef_ang)
Fmin, hmin, Fmax, hmax, i_min, i_max = pontos_criticos(newdepth, load)

for i in range(0, len(Fmax)):
    j = i_max[i]
    k = i_min[i]

    h_descarga = newdepth[j:k]
    F_descarga = load[j:k]
    F_fit = regressao_linear(h_descarga, F_descarga)
    plt.plot(h_descarga, F_fit, color='y')

#plt.plot(depth, load, color='b', label='Dados originais')
plt.plot(newdepth, load, label='Curva Calibrada')
plt.title('Gráfico de Força por Profundidade')
plt.scatter(hmin, Fmin, color="r", marker="D")
plt.scatter(hmax, Fmax, color="g", marker="D")
plt.xlabel('h (µm)')
plt.ylabel('F (N)')
plt.legend()
plt.grid(True)
plt.show()
