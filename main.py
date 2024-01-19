import matplotlib.pyplot as plt
import numpy as np

def critical_point(x, y):

    ymin, xmin, ymax, xmax = list(), list(), list(), list()
    index_min, index_max = list(), list()
    y0 = y[0]

    for i in range(0, len(y)):

        a = y[i] - y0

        if(a > 0  and y[i] < y[i-1]):
            y0 = y[i-1]
            ymax.append(y[i-1])
            xmax.append(x[i-1])
            index_max.append(i-1)

        if(a < 0  and y[i] > y[i-1]):
            y0 = y[i-1]
            ymin.append(y[i-1])
            xmin.append(x[i-1])
            index_min.append(i-1)

    ymin.append(y[len(y)-1])
    xmin.append(x[len(y)-1])  
    index_min.append(len(y)-1)

    return ymin, xmin, ymax, xmax, index_min, index_max

def coefficients(x, y):
    coef = np.polyfit(x, y, 1)
    a = coef[0]
    b = coef[1]
    return a, b

def linear_regression(x, y):
    coef = np.polyfit(x, y, 1)
    polinomio = np.poly1d(coef)
    y_fit = polinomio(x)
    return y_fit

def stiffness(h, F, Fmax):
    initial = Fmax*1
    final = Fmax*0.75

    for i in range(0, len(F)):
        if (F[i] == initial):
            i_initial = i
        else:
            lowest = F[0]
            for j in range(0, len(F)):
                dif = abs(F[j] - initial)
                if (dif < lowest):
                    lowest = dif
                    i_initial = j

    for i in range(0, len(F)):
        if (F[i] == final):
            newh = h[i_initial:i]
            newF = F[i_initial:i]
        else:
            lowest = F[i_initial]
            for j in range(0, len(F)):
                dif = abs(F[j] - final)
                if (dif < lowest):
                    lowest = dif
                    i_final = j
            newh = h[i_initial:i_final]
            newF = F[i_initial:i_final]

    a, b = coefficients(newh, newF)
    newF_fit = linear_regression(newh, newF)
    hp = -b/a

    return a, newh, newF_fit, hp

def stress_strain(n, hmax, Fmax, S):
    R = 0.5
    hc = hmax - (0.75*(Fmax/S))
    a = np.sqrt((5*(2-n)/(2*(4+n)))*((2*R*hc)-(hc*hc)))
    stress = (1/3)*(Fmax/(np.pi*a*a))
    aux = np.sqrt(1-((a/R)*(a/R)))
    strain = (0.14/aux)*(a/R)
    return stress, strain


dados = np.loadtxt("dados.txt", delimiter='\t')

load = dados[:,1]
depth = dados[:,2]

compliance = 0
depth = depth - compliance*load

coef_ang, coef_lin = coefficients(depth, load)
newdepth = depth - (-coef_lin/coef_ang)
Fmin, hmin, Fmax, hmax, i_min, i_max = critical_point(newdepth, load)

S = [None] * (len(Fmin))
hp = [None] * (len(Fmin))
for i in range(0, len(Fmax)):
    j = i_max[i]
    k = i_min[i]

    h_unload = newdepth[j:k]
    F_unload = load[j:k]
    S[i], h_fit, F_fit, hp[i] = stiffness(h_unload, F_unload, Fmax[i])
    plt.plot(h_fit, F_fit, color='y')

print(S) #rigidez
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
