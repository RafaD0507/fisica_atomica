# -*- coding: utf-8 -*-
"""
Created on Fri Apr  6 15:29:48 2018

@author: Rafael Sanabria
"""

import numpy as np
import matplotlib.pyplot as plt

def d2r_helium_atom_initial(u, r, dr, constants):
    energy, qn_l = constants
    return -(2/u)*dr-(energy+2/u)*r

def d2r_helium_atom_initial_2(u, r, dr, constants):
    energy, qn_l = constants
    return -(2/u)*dr-(energy+1/u)*r

def d2r_helium_atom(u, r, dr, constants, potential):
    energy, qn_l = constants
    return -(2/u)*dr-(energy-potential)*r

def solhde_numsolution(initconditions, x, d2f_function, constants=None, ans='f', div_val=None, potential = None):
    f, df = initconditions
    # list of lists - [f, df, d2f]
    solutions = [[], [], []]
    dx = x[1]-x[0]
    div_val = (len(x)-1) if div_val is None else div_val
    for i in range(len(x)):
        solutions[0].append(f)
        solutions[1].append(df)
        if potential is None:
            d2f = d2f_function(x[i], f, df, constants)
        else:
            d2f = d2f_function(x[i], f, df, constants, potential[i])
        solutions[2].append(d2f)
        df = df + d2f*dx
        f = f + df*dx
    if ans == 'f':
        return np.asarray(solutions[0])
    elif ans == 'div':
        return solutions[0][div_val]
    
def energies_location(initconditions, x, d2f_function, e_values, constants, potential = None):
    energy_locations = []
    e1 = e_values[0]
    constants1 = [e1]
    constants1.extend(constants[1:])
    div_1 = solhde_numsolution(initconditions, x, d2f_function, constants=tuple(constants1), ans='div', potential=potential)
    for i in range(1, len(e_values)-1): #Revisar este -1
        e2 = e_values[i]
        constants2 = [e2]
        constants2.extend(constants[1:])
        div_2 = solhde_numsolution(initconditions, x, d2f_function, constants=tuple(constants2), ans='div', potential=potential)
        sign = div_1*div_2
        if sign <= 0:
            energy_locations.append((e1, e2))
        e1 = e2
        div_1 = div_2
    return energy_locations

def energy_search(initconditions, x, d2f_function, constants, energy_location, epsilon, potential = None):
    down, up = energy_location
    mid = (down+up)/2.0
    count = 0
    f = 1000
    constants1 = [down]
    constants1.extend(constants[1:])
    sign = solhde_numsolution(initconditions, x, d2f_function, constants=tuple(constants1), ans='div', potential=potential)
    while abs(f) >= epsilon and count < 500:
        count += 1
        f = solhde_numsolution(initconditions, x, d2f_function, constants=tuple(constants1), ans='div', potential=potential)
        if sign > 0:
            if f > 0:
                down = mid
            elif f < 0:
                up = mid
            else:
                break
        else:
            if f > 0:
                up = mid
            elif f < 0:
                down = mid
            else:
                break
        mid = (up+down)/2.0
        constants1 = [mid]
        constants1.extend(constants[1:])
    if count == 500:
        pass
        # print('Maximum number of iterations reached at energy: '+ str(mid))
    return solhde_numsolution(initconditions, x, d2f_function, constants=tuple(constants1), potential=potential), mid

def wave_function_Schrodinger(initconditions, x, quantum_l, energy_range):
    l_names = {0: 's', 1: 'p', 2: 'd', 3: 'f'}
    plt.figure(figsize=(10, 10))
    for q_l in quantum_l:
        constants = [0, q_l]
        energyloc = energies_location(initconditions, x, d2r_hydrogen_atom, energy_range, constants)
        for energy_pair in energyloc:
            r_func, energy = energy_search(initconditions, x, d2r_hydrogen_atom, constants, energy_pair, 0.001)
            temp = str(int(round((1/abs(energy))**0.5))) + l_names[int(q_l)]
            plt.plot(x, r_func/abs(max(r_func)), label=temp)
    plt.title('Funcion de onda', fontsize=25)
    plt.xlabel('u', fontsize=20)
    plt.xlim(0, max(u))
    plt.ylabel('Funcion radial R(u)', fontsize=20)
    plt.axhline(0, color='k')
    plt.legend(fontsize=15)
    plt.savefig(str(quantum_l[len(quantum_l)-1])+'.png')
    plt.close()
    
initial_conditions = (1, -1)
u = np.linspace(0.01, 10, 2000)
constants_He = [0, 0]
du = u[1]-u[0]
r_A, e_A = None, None
r_B, e_B = None, None
V_total_B = None
V_iniciales = [[0, 0.13, 0.14],[0.51, 0.14, 0.145]]

for i in range(3):
    if i == 0:     
        plt.figure()
        #Para A1
        e_v = np.linspace(-1, -0.8, 100)
        energyloc_He_A1 = energies_location(initial_conditions, u, d2r_helium_atom_initial, e_v, constants_He)
        for energy_pair in energyloc_He_A1:
            r_A, e_A = energy_search(initial_conditions, u, d2r_helium_atom_initial, constants_He, energy_pair, 0.001)
            temp = 'A, '+'e = ' + str(round(e_A,2))
            plt.plot(u, r_A, label=temp)
        #Para B1
        e_v = np.linspace(-0.3,-0.2, 100)
        energyloc_He_B1 = energies_location(initial_conditions, u, d2r_helium_atom_initial_2, e_v, constants_He)
        r_B, e_B = None, None
        for energy_pair in energyloc_He_B1:
            r_B, e_B = energy_search(initial_conditions, u, d2r_helium_atom_initial_2, constants_He, energy_pair, 0.001)
            temp = 'B, '+'e = ' + str(round(e_B,2))
            plt.plot(u, r_B, label= temp)
        plt.legend()
        plt.title("Iteracion 0")
        plt.show()
    else:
        #Para A
        plt.figure()
        plt.title("Iteracion "+ str(i))
        e_v = np.linspace(-1.0, -0, 100)
        energyloc_He_A = energies_location(initial_conditions, u, d2r_helium_atom, e_v, constants_He, potential=V_total_B)
        for enery_pair in energyloc_He_A[:1]:
            print(i)
            r_A, e_A = energy_search(initial_conditions, u, d2r_helium_atom, constants_He, energy_pair, 0.001, potential=V_total_B)
            temp = 'A, '+'e = ' + str(round(e_A,4))
            plt.plot(u, r_A, label=temp)
        P_A = (r_A**2)*(u**2) #Verificar longitud
        S_A = [0]
        for probabilidad in P_A[:-1]:
            S_A.append(probabilidad*du + S_A[-1])
        Q = S_A/S_A[-1]
        V_A = [V_iniciales[0][i]] #Definir el potencial inicial
        for j in range(len(Q)-1):
            V_A.append(V_A[-1]-Q[j]*du/u[j]**2)
        V_total_A = V_A+(-2/u)
        
        e_v = np.linspace(-1, -0, 100)
        energyloc_He_B = energies_location(initial_conditions, u, d2r_helium_atom, e_v, constants_He, potential=V_total_A)
        for enery_pair in energyloc_He_B[:1]:
            print(i*10)
            r_B, e_B = energy_search(initial_conditions, u, d2r_helium_atom, constants_He, energy_pair, 0.001, potential=V_total_A)
            temp = 'B, '+'e = ' + str(round(e_B,4))
            plt.plot(u, r_B, label=temp)
        plt.legend()
        plt.figure()
        plt.title("V inicial A, iteracion="+str(i))
        plt.loglog(u, 1/u)
        plt.loglog(u, 2/u)
        plt.loglog(u, abs(V_total_A))
        
    
    P_B = (r_B**2)*(u**2) #Verificar longitud
    S_B = [0]
    for probabilidad in P_B[:-1]:
        S_B.append(probabilidad*du + S_B[-1])
    Q = S_B/S_B[-1]
    V_B = [V_iniciales[1][i]] #Definir el potencial inicial
    for j in range(len(Q)-1):
        V_B.append(V_B[-1]-Q[j]*du/u[j]**2)
    V_total_B = V_B+(-2/u)
    plt.figure()
    plt.title("V inicial B, iteracion=" + str(i))
    plt.loglog(u, 1/u)
    plt.loglog(u, 2/u)
    plt.loglog(u, abs(V_total_B))



