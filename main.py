from scipy.optimize import least_squares
import numpy as np
import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt

# Constantes físicas globais
R = 0.5              # raio do indentador
v = 0.3              # Poisson do material
vi = 0.21            # Poisson do indentador
Ei = 630000          # módulo do indentador (MPa)

# Parâmetros de método
UNLOAD_MIN = 0.75    # fração de Fmax para rigidez
CAL_MIN = 0.70       # calibração: 70% do 1º ciclo
CAL_MAX = 0.95       # calibração: 95% do 1º ciclo


# Funções matemáticas 
def function(strain, K, n):
    """
    Modelo constitutivo tensão-deformação do material, em que 
    'n' é o expoente de encruamento e 'K' é uma constante.
    """ 
    return K * strain**n 

def contact_radius(R, n, hc):
    """
    Calcula o raio de contato entre o indentador e o material, 
    em que 'hc' é a profundidade de contato.
    """
    return np.sqrt((5*(2-n)/(2*(4+n))) * (2*R*hc - hc**2))

def young_modulus(a, S, v, vi, Ei):
    """
    Calcula o módulo de elasticidade do material a partir da rigidez
    de descarregamento 'S' e do raio de contato 'a'.
    """
    return (1 - (v*v))/((2*a/S) - ((1-(vi*vi))/Ei))

def representative_stress_strain(R, n, hmax, Fmax, S):
    """
    Calcula a tensão e deformação representativas de um ciclo
    de indentação instrumentada, em que 'hmax' é a profundidade 
    máxima do ciclo e 'Fmax' é a força máxima do ciclo.
    """
    hc = hmax - 0.75 * (Fmax / S)
    a = contact_radius(R, n, hc)
    
    stress = (1/3) * Fmax / (np.pi * a**2)
    strain = (0.14 / np.sqrt(1 - (a/R)**2)) * (a/R)
    return stress, strain


# Função de otimização
def error(params, strain, stress):
    """
    Função erro para ajuste por mínimos quadrados do modelo constitutivo.
    """
    K, n = params
    calc_stress = function(strain, K, n)
    return stress - calc_stress


# Processamento de dados
def load_data(filename):
    """
    Carrega dados experimentais de um arquivo de texto.
    """
    data = np.loadtxt(filename, delimiter='\t')
    load = data[:, 1]
    depth = data[:, 2]
    return load, depth

def apply_compliance(depth, load, compliance):
    """
    Corrige a profundidade considerando a complacência do sistema.
    """
    return depth - compliance * load

def calibrate_depth(depth, load, i_max):
    """
    Calibra a profundidade deslocando a curva para que o prolongamento
    do primeiro ciclo de carregamento passe pela origem.
    """
    j = i_max[0]
    start = int(CAL_MIN * j)
    end = int(CAL_MAX * j)

    a, b = np.polyfit(depth[start:end], load[start:end], 1)
    return depth - (-b/a)

def unloading_stiffness(h, F, Fmax):
    """
    Calcula a rigidez de descarregamento (S) a partir de
    uma regressão linear no trecho superior do ciclo,
    conforme metodologia.
    """
    mask = (F >= UNLOAD_MIN*Fmax) & (F <= Fmax)
    a, b = np.polyfit(h[mask], F[mask], 1)
    return a, -b/a

def extract_cycles(x, y):
    """
    Identifica pontos críticos de carga e descarga
    a partir da curva força x profundidade.
    """
    ymin, xmin, ymax, xmax = list(), list(), list(), list()
    index_min, index_max = list(), list()
    y0 = y[0]

    for i in range(0, len(y)):

        a = y[i] - y0

        if (a > 0 and y[i] < y[i-1]):
            y0 = y[i-1]
            ymax.append(y[i-1])
            xmax.append(x[i-1])
            index_max.append(i-1)

        if (a < 0 and y[i] > y[i-1]):
            y0 = y[i-1]
            ymin.append(y[i-1])
            xmin.append(x[i-1])
            index_min.append(i-1)

    ymin.append(y[len(y)-1])
    xmin.append(x[len(y)-1])
    index_min.append(len(y)-1)

    return ymin, xmin, ymax, xmax, index_min, index_max


# Leitura e pré-processamento 
load, depth = load_data("dados.txt")
depth = apply_compliance(depth, load, compliance=0)

depth = depth * 1e-3

Fmin, hmin, Fmax, hmax, i_min, i_max = extract_cycles(depth, load)

depth_cal = calibrate_depth(depth, load, i_max)

Fmin, hmin, Fmax, hmax, i_min, i_max = extract_cycles(depth_cal, load)


# Rigidez e deformação plástica
n_cycles = len(i_max)
S = np.zeros(n_cycles)
hp = np.zeros(n_cycles)

for i in range(n_cycles):
    j = i_max[i]   
    k = i_min[i]   

    h_unload = depth_cal[j:k]
    F_unload = load[j:k]

    S[i], hp[i] = unloading_stiffness(
        h_unload,
        F_unload,
        Fmax[i]
    )


# Tensão-deformação
stress = []
strain = []
n_guess = 0.14    

for i in range(n_cycles):
    s, e = representative_stress_strain(
        R,
        n_guess,
        hmax[i],
        Fmax[i],
        S[i]
    )
    stress.append(s)
    strain.append(e)

stress = np.array(stress)
strain = np.array(strain)

initial_params = [1406.048, 0.14]
bounds = ([703.24, 0.07], [2109.072, 0.21]) 
result = least_squares(
    error,
    initial_params,
    bounds=bounds,
    args=(strain, stress)
)
K, n = result.x

x = [None] * len(hmin)
y = [None] * len(Fmin)
a = [None] * len(Fmin)

for i in range(0, len(Fmax)):
    hc = hmax[i] - 0.75 * (Fmax[i] / S[i])
    a[i] = contact_radius(R, n, hc)
    if (i == 0):
        E = young_modulus(a[i], S[i], v, vi, Ei)
    y[i] = Fmax[i]/(4*a[i]*a[i])
    x[i] = a[i]/R

newx = np.log(x)
newy = np.log(y)
A, B = np.polyfit(newx, newy, 1)
M = A + 2
material_param = np.exp(B)
yield_stress = 0.2285*material_param + 0

stress_opt = function(strain, K, n)
yield_strain = yield_stress/E
x_line = np.linspace(0, yield_strain, 100)
y_line = E*x_line

# Gráficos
plt.figure()
plt.plot(depth, load, color='b', label='Dados originais')
plt.plot(depth_cal, load, label='Curva Calibrada')
plt.title('Gráfico de Força por Profundidade')
plt.scatter(hmin, Fmin, color="r", marker="D")
plt.scatter(hmax, Fmax, color="g", marker="D")
plt.xlim(left=0)
plt.ylim(bottom=0)
plt.xlabel('h (mm)')
plt.ylabel('F (N)')
plt.grid(True)
plt.show()

plt.figure()
plt.scatter(strain, stress_opt)
plt.plot(x_line, y_line, color='orange')
plt.title('Gráfico de Tensão-Deformação')
plt.xlabel('Deformação')
plt.ylabel('Tensão (Mpa)')
plt.grid(True)
plt.show()