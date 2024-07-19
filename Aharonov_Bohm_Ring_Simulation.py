import math
import numpy as np
import matplotlib.pyplot as plt
from numpy import sqrt

# Initialization
def f(x):
    return (1/(1+math.exp(x/(kb*300))))

n = 52 # Number of atoms at each ring of AGNR
Hex_N = 5 # Number of hexagon columns in the channel
N_Ch = 2*n*Hex_N # Number of atoms in the device channel

# Basic Distances in Carbon Honey-Comb Structure
a0 = 1.42
d1 = a0
d2 = 2.45951214
d3 = 2.84
mrg = 0.01

# Hopping Parameters
eps = 0
t1 = -2.65
t2 = 0
t3 = 0

# Device Size
D_w = 14 * d1 # Device width
Ch_l = 25 * sqrt(3) * d1 / 10 # Channel (device) length

print(f"Device width = {D_w}, Device length = {Ch_l}")

# Make Site
print("Generating atom sites")

atom_x = np.zeros(N_Ch+1)
atom_y = np.zeros(N_Ch+1)

atom_x[0] = -10
atom_y[0] = -10

# Assigning geometric position to each atom
for j in range(Hex_N+1):
    if j == 0: continue
    for i in range(n // 2 + 1):
        if i == 0: continue
        atom_x[2 * (i - 1) + 1 + (j - 1) * 2 * n] = d1 / 2 + (j - 1) * d1 * 3
        atom_x[2 * i + (j - 1) * 2 * n] = 0 + (j - 1) * d1 * 3
        atom_y[2 * (i - 1) + 1 + (j - 1) * 2 * n] = (i - 1) * d1 * math.sqrt(3)
        atom_y[2 * i + (j - 1) * 2 * n] = (i - 0.5) * d1 * math.sqrt(3)

        atom_x[2 * (i - 1) + 1 + n + (j - 1) * 2 * n] = d1 * 1.5 + (j - 1) * d1 * 3
        atom_x[2 * i + n + (j - 1) * 2 * n] = d1 * 2 + (j - 1) * d1 * 3
        atom_y[2 * (i - 1) + 1 + n + (j - 1) * 2 * n] = (i - 1) * d1 * math.sqrt(3)
        atom_y[2 * i + n + (j - 1) * 2 * n] = (i - 0.5) * d1 * math.sqrt(3)

    if (n % 2 == 1):
        atom_x[n + (j - 1) * 2 * n] = 0 + (j - 1) * d1 * 3
        atom_y[n + (j - 1) * 2 * n] = n // 2 * d1 * math.sqrt(3)
        atom_x[2 * n + (j - 1) * 2 * n] = d1 * 2 + (j - 1) * d1 * 3
        atom_y[2 * n + (j - 1) * 2 * n] = n // 2 * d1 * math.sqrt(3)

# Flip the entire structure before making the cuts
temp = atom_x
atom_x = atom_y
atom_y = temp

# Initial Cut (to get rid of the extra row of atoms)
N_Cut = 0

for i in range(N_Ch + 1):
    if i == 0: continue
    if (atom_x[i] > 61.5):
        atom_x[i] = -10
        atom_y[i] = -10
        N_Cut += 1

# Make Cut
Square_N = 7
S_x = [0.0] * (Square_N + 1)
S_y = [0.0] * (Square_N + 1)
S_l = [0.0] * (Square_N + 1)
S_w = [0.0] * (Square_N + 1)

S_x[1] = 13
S_y[1] = 0
S_l[1] = 15
S_w[1] = 15

'''
# You may replace the code block above with this to make the exact structure
# proposed by Motohiko Ezawa, but keep in mind that the Band Gap will become 0.55.
S_x[1] = 0
S_y[1] = 0
S_l[1] = 40
S_w[1] = 15
'''

S_x[2] = 20.906
S_y[2] = 0
S_l[2] = 1
S_w[2] = 11

S_x[3] = 22.135
S_y[3] = 0
S_l[3] = 1
S_w[3] = 9

S_x[4] = 23.365
S_y[4] = 0.3
S_l[4] = 1
S_w[4] = 1

S_x[5] = 30.744
S_y[5] = 9.94
S_l[5] = 14
S_w[5] = 4

S_x[6] = 30.744
S_y[6] = 9.94
S_l[6] = 10
S_w[6] = 10

S_x[7] = 30.744
S_y[7] = 9.94
S_l[7] = 8
S_w[7] = 12

# Cutting out the antidotes using symmetry
for j in range (Square_N + 1):
    if j==0: continue
    for i in range (N_Ch + 1):
        if i==0: continue
        if (atom_x[i] < (S_x[j] + S_l[j]/2)) and (atom_x[i] > (S_x[j] - S_l[j]/2)) and (atom_y[i] < (S_y[j] + S_w[j]/2)) and (atom_y[i] > (S_y[j] - S_w[j]/2)):
            atom_x[i] = -10
            atom_y[i] = -10
            N_Cut += 1
        if (atom_x[i] < ((Ch_l*10 - S_x[j]) + S_l[j]/2)) and (atom_x[i] > ((Ch_l*10 - S_x[j]) - S_l[j]/2)) and (atom_y[i] < (S_y[j] + S_w[j]/2)) and (atom_y[i] > (S_y[j] - S_w[j]/2)):
            atom_x[i] = -10
            atom_y[i] = -10
            N_Cut += 1
        if (atom_x[i] < ((Ch_l*10 - S_x[j]) + S_l[j]/2)) and (atom_x[i] > ((Ch_l*10 - S_x[j]) - S_l[j]/2)) and (atom_y[i] < ((D_w - S_y[j]) + S_w[j]/2)) and (atom_y[i] > ((D_w - S_y[j]) - S_w[j]/2)):
            atom_x[i] = -10
            atom_y[i] = -10
            N_Cut += 1
        if (atom_x[i] < ((S_x[j]) + S_l[j]/2)) and (atom_x[i] > ((S_x[j]) - S_l[j]/2)) and (atom_y[i] < ((D_w - S_y[j]) + S_w[j]/2)) and (atom_y[i] > ((D_w - S_y[j]) - S_w[j]/2)):
            atom_x[i] = -10
            atom_y[i] = -10
            N_Cut += 1

# Flip the entire structure again after finishing the cuts
temp = atom_x
atom_x = atom_y
atom_y = temp

# Plotting the atomic structure
for i in range(len(atom_x)):
    if i == 0: continue
    if i >= len(atom_x): break # Bug fix
    if atom_x[i] > -0.5:
        plt.scatter(atom_x[i], atom_y[i], c='blue')
plt.title("Atomic Structure")
plt.xlabel("X Coordinate")
plt.ylabel("Y Coordinate")
plt.axis('equal')
plt.show()

print("Done!")

# Make Hamiltonians
print("Making Hamiltonians")

# Generate H of Channel
H_Ch = np.zeros((N_Ch, N_Ch))

for j in range(1, N_Ch+1):
    for i in range(1, N_Ch+1):
        dis = sqrt(((atom_x[j]-atom_x[i])**2) + ((atom_y[j]-atom_y[i])**2))

        if dis < mrg: H_Ch[j-1, i-1] = eps
        elif d1-mrg < dis < d1+mrg: H_Ch[j-1, i-1] = t1
        elif d2-mrg < dis < d2+mrg: H_Ch[j-1, i-1] = t2
        elif d3-mrg < dis < d3+mrg: H_Ch[j-1, i-1] = t3
        else: H_Ch[j-1, i-1] = 0

# Generate alpha and beta matrices to be used in NEGF calculation
alpha = np.zeros((2*n, 2*n))
beta = np.zeros((2*n, 2*n))

atom_x0 = np.zeros(2*n+1)
atom_y0 = np.zeros(2*n+1)

atom_xr = np.zeros(2*n+1)
atom_yr = np.zeros(2*n+1)

for i in range(1, 2*n+1):
    atom_x0[i] = atom_x[i]
    atom_y0[i] = atom_y[i]

    atom_xr[i] = atom_x[i] + (3*d1)
    atom_yr[i] = atom_y[i]

# Generate alpha
for j in range(1, 2*n+1):
    for i in range(1, 2*n+1):
        dis = sqrt(((atom_x0[j]-atom_x0[i])**2) + ((atom_y0[j]-atom_y0[i])**2))

        if dis < mrg: alpha[j-1, i-1] = eps
        elif d1-mrg < dis < d1+mrg: alpha[j-1, i-1] = t1
        elif d2-mrg < dis < d2+mrg: alpha[j-1, i-1] = t2
        elif d3-mrg < dis < d3+mrg: alpha[j-1, i-1] = t3
        else: alpha[j-1, i-1] = 0

# Generate beta

for j in range(1, 2*n+1):
    for i in range(1, 2*n+1):
        dis = sqrt(((atom_x0[j]-atom_xr[i])**2) + ((atom_y0[j]-atom_yr[i])**2))

        if dis < mrg: beta[j-1, i-1] = eps
        elif d1-mrg < dis < d1+mrg: beta[j-1, i-1] = t1
        elif d2-mrg < dis < d2+mrg: beta[j-1, i-1] = t2
        elif d3-mrg < dis < d3+mrg: beta[j-1, i-1] = t3
        else: beta[j-1, i-1] = 0

print("Done!")

# NEGF
# Define Constants
q = 1.6e-19
hbar = 1.06e-34
zplus = 1j*1e-3
kb = 8.3333e-05
IE = (q*q)/(2*np.pi*hbar)

# Temperature (in Kelvin)
Temp = 300

# Fermi level at equilibrium
Band_Gap = 0.346668 #Band gap value is obtained from a seperate band structure calculation
Ef = (Band_Gap/2)+0.02 #Fermi level is chosen slightly higher than half of band gap.

# NEGF error
errmax = 0.001

# Variables dependent on Channel Length (Ch_l)
V_min = 50
V_max = 2000
div_num = 60


# Back-up lists
E_list = []
T_list = []

V = np.linspace(V_min/1000,V_max/1000,div_num)
E = np.linspace((Ef - 4*kb*Temp),(Ef + 4*kb*Temp),100)

I = np.zeros(len(V))

for i in range(len(E)): E_list.append(E[i])

print("Running NEGF...please wait.")
print(f"Energy range = {min(E)} eV, {max(E)} eV")

N_Contact = 2*n

# Define matrices for NEGF
T = np.zeros(len(E))
g1 = np.linalg.inv(E[0]*np.eye(N_Contact)-alpha)
g2 = np.linalg.inv(E[0]*np.eye(N_Contact)-alpha)
U = np.zeros((N_Ch, N_Ch))
H_Buff = H_Ch
Ch_l = Ch_l * 10
b = np.zeros(N_Ch)

for p in range(len(V)):
    print("Loop Number (out of 60) : ", p+1)
    for i in range(1, N_Ch+1): b[i-1] = atom_x[i]/Ch_l

    U = -1 * np.diag(b) * V[p]
    H_Ch = H_Buff + U

    # Energy Loop
    for k in range(len(E)):
        err = 100
        while err > errmax:
            g1new = np.linalg.inv((E[k]+zplus)*np.eye(N_Contact)-alpha-np.conj(beta).T @ g1 @ beta)
            err = np.sum(np.sum(np.abs(g1new-g1)))/np.sum(np.sum(np.abs(g1new+g1)))
            g1 = g1new

        sigma1 = np.conj(beta).T @ g1 @ beta

        err = 100

        while err > errmax:
            g2new = np.linalg.inv((E[k] + zplus) * np.eye(N_Contact) - alpha - (beta @ g2 @ np.conj(beta).T))
            err = np.sum(np.sum(np.abs(g2new - g2))) / np.sum(np.sum(np.abs(g2new + g2)))
            g2 = g2new

        sigma2 = beta @ g2 @ np.conj(beta).T

        # Calculate self energy matrices

        diag1 = np.diag([1] + [0] * (Hex_N - 1))
        diag2 = np.diag([0] * (Hex_N - 1) + [1])
        E1 = np.kron(diag1, sigma1)
        E2 = np.kron(diag2, sigma2)

        # Calculate broadening and transmission

        G1 = 1j*(E1 - np.conj(E1).T)
        G2 = 1j*(E2 - np.conj(E2).T)

        G = np.linalg.inv((E[k] + zplus)*np.eye(N_Ch)-H_Ch-E1-E2)
        T[k] = np.real(np.trace(G1 @ G @ G2 @ np.conj(G).T))

    dE = E[1] - E[0]
    dI = [0] * len(E)
    term1 = []
    term2 = []
    for i in range(0, len(E)):
        term1.append(E_list[i] - Ef)
        term2.append(E_list[i] - Ef + V[p])
        a = f(term1[i])
        c = f(term2[i])
        dI[i] = IE * (a - c) * T[i] * dE

    I[p] = np.sum(dI)

# Plotting the relationships
plt.plot(E_list, T)
plt.xlabel("E(eV)")
plt.ylabel("T(E)")
plt.title("Transmission without Logarithm")
plt.show()

T_log10 = np.log10(T)

plt.plot(E_list, T_log10)
plt.xlabel("E(eV)")
plt.ylabel("T(E) - log 10")
plt.title("Transmission with Logarithm")
plt.show()

#Plotting the Current - Voltage Graph
plt.plot(V, I)
plt.xlabel("Voltage(V)")
plt.ylabel("Current(I)")
plt.title("V - I Plot")
plt.show()

print("Done!")
print("Simulation is complete.")