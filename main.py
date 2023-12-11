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

def rigidez(h, F, Fmax):
    inicial = Fmax*1
    final = Fmax*0.75

    for i in range(0, len(F)):
        if (F[i] == inicial):
            i_inicial = i
        else:
            menor = F[0]
            for j in range(0, len(F)):
                dif = abs(F[j] - inicial)
                if (dif < menor):
                    menor = dif
                    i_inicial = j

    for i in range(0, len(F)):
        if (F[i] == final):
            newh = h[i_inicial:i]
            newF = F[i_inicial:i]
        else:
            menor = F[i_inicial]
            for j in range(0, len(F)):
                dif = abs(F[j] - final)
                if (dif < menor):
                    menor = dif
                    i_final = j
            newh = h[i_inicial:i_final]
            newF = F[i_inicial:i_final]

    a, b = coeficientes(newh, newF)
    newF_fit = regressao_linear(newh, newF)
    hp = -b/a

    return a, newh, newF_fit, hp


dados = np.loadtxt("dados.txt", delimiter='\t')

load = dados[:,1]
depth = dados[:,2]

compliance = 0
depth = depth - compliance*load

coef_ang, coef_lin = coeficientes(depth, load)
newdepth = depth - (-coef_lin/coef_ang)
Fmin, hmin, Fmax, hmax, i_min, i_max = pontos_criticos(newdepth, load)

S = [None] * (len(Fmin))
hp = [None] * (len(Fmin))
for i in range(0, len(Fmax)):
    j = i_max[i]
    k = i_min[i]

    h_descarga = newdepth[j:k]
    F_descarga = load[j:k]
    S[i], h_fit, F_fit, hp[i] = rigidez(h_descarga, F_descarga, Fmax[i])
    plt.plot(h_fit, F_fit, color='y')

print(S)
print(hp) #deformação plástica

plt.plot(depth, load, color='b', label='Dados originais')
plt.plot(newdepth, load, label='Curva Calibrada')
plt.title('Gráfico de Força por Profundidade')
plt.scatter(hmin, Fmin, color="r", marker="D")
plt.scatter(hmax, Fmax, color="g", marker="D")
plt.xlabel('h (µm)')
plt.ylabel('F (N)')
plt.legend()
plt.grid(True)
plt.show()
