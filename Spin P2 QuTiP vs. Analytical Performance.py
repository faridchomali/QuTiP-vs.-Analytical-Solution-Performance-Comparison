number_of_points = 10000 # Number of points to compute. Adjust as needed


####################################################################
# PART 1: ANALYTIC SIMULATION

import numpy as np

# Define constants
Case = 2 #1 FOR CONSTANT MAGNETIC FIELD, 2 for OSCILLATING MAGNETIC FIELD
omega = 1  # Adjust as needed
gammaBx = 10  # Adjust as needed
t_0 = 0
t_f = 2 * np.pi / omega  # Total time for simulation

#Initial conditions for case 2:
a_0 = 0
b_0 = 1

# Function definitions for analytic solutions
def f(t):
    return gammaBx / (2 * omega) * np.sin(omega * t)

def alpha(t):
    if (Case == 1):
        return 1j * np.sin(1/2*t)
    if (Case == 2):
        #Esta es nuestra solución analítica original, que dejo comentada:
        #return a_0 * np.cos(np.abs(f(t))) + 1j * b_0 * np.sin(np.abs(f(t)))
        #La siguiente, por otro lado, es sin el valor absoluto:
        #Para el primer periodo, no genera ninguna diferencia, pero cuando
        #comienza el segundo período, ahí comienza a retornar el opuesto del <Sy>
        #que retorna QuTiP.
        return a_0 * np.cos((f(t))) + 1j * b_0 * np.sin((f(t)))

def beta(t):
    if (Case == 1):
        return np.cos(1/2*t)
    if (Case == 2):
        #return 1j * a_0 * np.sin(np.abs(f(t))) + b_0 * np.cos(np.abs(f(t)))
        return 1j * a_0 * np.sin((f(t))) + b_0 * np.cos((f(t)))
    
# Define the vector components
def Sx(t):
    return np.real(np.conjugate(alpha(t)) * beta(t)) / (1 / 2)  # 1/2 as in hbar/2

def Sy(t):
    return np.imag(np.conjugate(alpha(t)) * beta(t) - alpha(t) * np.conjugate(beta(t)))

def Sz(t):
    return ((np.abs(alpha(t)))**2 - (np.abs(beta(t)))**2)

# Generate time points
t_points = np.linspace(t_0, t_f, number_of_points)

# Compute <Sx>, <Sy>, <Sz> at each time point and store as list of tuples
S_x_list = [(t, Sx(t)) for t in t_points]
S_y_list = [(t, Sy(t)) for t in t_points]
S_z_list = [(t, Sz(t)) for t in t_points]

def print_analytic_results():
    print("S_x_list:")
    print(S_x_list)
    print("S_y_list:")
    print(S_y_list)
    print("S_z_list:")
    print(S_z_list)

#print_analytic_results()


####################################################################
# PART 2: QUTIP NUMERICAL SIMULATION

from qutip import *
from numpy import sqrt, pi, array, sin, cos, arange
import matplotlib.pyplot as plt
import numpy as np

pz = basis(2,0) # i.e. the column vector (1,0)
mz = basis(2,1) # i.e. the column vector (0,1)
px = 1/sqrt(2)*(pz + mz)
mx = 1/sqrt(2)*(pz - mz)
py = 1/sqrt(2)*(pz + 1j*mz)
my = 1/sqrt(2)*(pz - 1j*mz)

Sx = 1/2.0*sigmax() / (1/2)
Sy = 1/2.0*sigmay() / (1/2)
Sz = 1/2.0*sigmaz() / (1/2)

# Note that \hbar has been set = 1 (natural units).
# Also, we are reduntantly dividing by 1/2 to make it clear that we are normalizing.

# CONSTANT MAGNETIC FIELD SETUP:
if (Case == 1):
    Hx = -omega*Sx/2
    t = np.linspace(0, t_f, number_of_points)
    expect_ops = [Sx,Sy,Sz]
    psi0 = mz
    #print(psi0)
    result = sesolve(Hx, psi0, t, expect_ops)

# OSCILLATING MAGNETIC FIELD SETUP
if Case == 2:
    t = np.linspace(0, t_f, number_of_points)
    
    def H1_coeff(t, args):
        return np.cos(args['omega'] * t)
    
    H = [[-1/2 * gammaBx * Sx, H1_coeff]]
    args = {'omega': omega}
    expect_ops = [Sx, Sy, Sz]
    psi0 = mz
    #print(psi0)
    result = sesolve(H, psi0, t, expect_ops, args=args)

def print_qutip_results():
    print(result.expect[0])
    print(result.expect[1])
    print(result.expect[2])

#print_qutip_results


####################################################################
# PART 3: GRAPHING BOTH RESULTS

# Assuming S_x_list, S_y_list, S_z_list are already defined and filled with data
# result is assumed to be the QuTiP simulation result containing expect[0], expect[1], expect[2]

# Time array for analytic data (assuming S_x_list, S_y_list, S_z_list are tuples of (time, value))
time_analytic = np.array([t[0] for t in S_x_list])

# Analytic data
Sx_analytic = np.array([t[1] for t in S_x_list])
Sy_analytic = np.array([t[1] for t in S_y_list])
Sz_analytic = np.array([t[1] for t in S_z_list])

# QuTiP data
Sx_qutip = result.expect[0]
Sy_qutip = result.expect[1]
Sz_qutip = result.expect[2]

# Plot 1: Analytic <Sx>, <Sy>, <Sz> against time
plt.figure(figsize=(12, 6))
plt.subplot(131)
plt.plot(time_analytic, Sx_analytic, 'b-', label='<Sx> Analytic', linewidth=1)
plt.plot(time_analytic, Sy_analytic, 'g-', label='<Sy> Analytic', linewidth=1)
plt.plot(time_analytic, Sz_analytic, 'r-', label='<Sz> Analytic', linewidth=1)
plt.xlabel('Time')
plt.ylabel('Expectation Values')
plt.title('Analytic <Sx>, <Sy>, <Sz> vs Time')
plt.legend()

# Plot 2: QuTiP <Sx>, <Sy>, <Sz> against time
plt.subplot(132)
plt.plot(result.times, Sx_qutip, 'b-', label='<Sx> QuTiP')
plt.plot(result.times, Sy_qutip, 'g-', label='<Sy> QuTiP')
plt.plot(result.times, Sz_qutip, 'r-', label='<Sz> QuTiP')
plt.xlabel('Time')
plt.ylabel('Expectation Values')
plt.title('QuTiP <Sx>, <Sy>, <Sz> vs Time')
plt.legend()

# Plot 3: Combined plot of Analytic (dotted) and QuTiP (continuous)
plt.subplot(133)
plt.plot(time_analytic, Sx_analytic, 'bo-', label='<Sx> Analytic', markersize=4, linewidth=0.5)
plt.plot(time_analytic, Sy_analytic, 'go-', label='<Sy> Analytic', markersize=4, linewidth=0.5)
plt.plot(time_analytic, Sz_analytic, 'ro-', label='<Sz> Analytic', markersize=4, linewidth=0.5)
plt.plot(result.times, Sx_qutip, 'b-', label='<Sx> QuTiP')
plt.plot(result.times, Sy_qutip, 'g-', label='<Sy> QuTiP')
plt.plot(result.times, Sz_qutip, 'r-', label='<Sz> QuTiP')
plt.xlabel('Time')
plt.ylabel('Expectation Values')
plt.title('Analytic vs QuTiP Comparison')
plt.legend()

plt.tight_layout()
plt.show()


####################################################################
# PART 4: PERFORMANCE COMPARATION

# Step 1: Define the array
array_0 = result.expect[0]
array_1 = result.expect[1]
array_2 = result.expect[2]

# Step 2: Define the list of tuples
list_of_tuples_x = S_x_list
list_of_tuples_y = S_y_list
list_of_tuples_z = S_z_list

# Step 3: Extract the second value of each tuple to create a list of those values
extracted_values_x = [t[1] for t in list_of_tuples_x][:number_of_points]
extracted_values_y = [t[1] for t in list_of_tuples_y][:number_of_points]
extracted_values_z = [t[1] for t in list_of_tuples_z][:number_of_points]

array_list_0 = array_0.tolist()[:number_of_points]
array_list_1 = array_1.tolist()[:number_of_points]
array_list_2 = array_2.tolist()[:number_of_points]

print(array_list_1)
print(extracted_values_y)

# Step 5: Calculate the mean square error between the two lists
def mean_square_error(list1, list2):
    if len(list1) != len(list2):
        raise ValueError("Lists must have the same length.")
    return np.mean((np.array(list1) - np.array(list2))**2)

mse_x = mean_square_error(array_list_0, extracted_values_x)
mse_y = mean_square_error(array_list_1, extracted_values_y)
mse_z = mean_square_error(array_list_2, extracted_values_z)

#print(f"Array: {array_list}")
#print(f"Extracted Values: {extracted_values}")
print(f"Mean Square Error in X: {mse_x}")
print(f"Mean Square Error in Y: {mse_y}")
print(f"Mean Square Error in Z: {mse_z}")
print(f"Average Mean Square Error: {(mse_x+mse_y+mse_z)/3}")