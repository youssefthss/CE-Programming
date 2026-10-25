#El script está incompleto

from CoolProp.CoolProp import PropsSI
from CoolProp.CoolProp import PhaseSI
import numpy as np
import pandas as pd

list_sustancias = []
Resultados = {}

R = 8.3145 #Constante gases ideales

def model_GI():
    Z = np.array([1])
    
    return Z

def model_VdW(P, T, Tc, Pc):
    a = (27/64)*(R*Tc)**2/Pc
    b = R*Tc/(8*Pc)
    A = a*P/(R*T)**2
    B = b*P/(R*T)
    Z = np.roots([1 , -(1 + B) , A , -A*B]) 
    #Z es un np.array. 
    #Si como mínimo una raíz es compleja, cada elemento de Z es un numpy.complex128, cosa diferente al data type 'complex' de Python
    #Si cada raíz es real, cada elemento de Z es numpy.float64, cosa diferente al data type 'float' de Python
    
    for i in range(3): #Si Z[0] tiene un valor sin sentido físico, Z[0] = 0

        if Z[i].imag != 0 or Z[i].real < 0 or Z[i].real < B or Z[i].real == 0:
            Z[i] = 0    

    i = 0
    while i < len(Z):
        if Z[i] == 0:
            Z = np.delete(Z,i)

        else:
            i = i + 1

    if len(Z) == 3: #Si Z aún contiene 3 raíces, sabemos que la intermedia no tiene significado físico
        Z = np.array([Z[np.argmax(Z)], Z[np.argmin(Z)]])

    Z = Z.real

    return Z

def model_RK(P, T, Tc, Pc):
    a = 0.42748*R**2*Tc**2.5/Pc
    b = 0.08664*R*Tc/Pc
    A = a*P/(R**2*T**2.5)
    B = b*P/(R*T)

    Z = np.roots([1 , -1 , (A - B - B**2) , -A*B])

    for i in range(3):

        if Z[i].imag != 0 or Z[i].real < 0 or Z[i].real < B or Z[i].real == 0:

            Z[i] = 0

    i = 0
    while i < len(Z):
        if Z[i] == 0:
            Z = np.delete(Z,i)

        else:
            i = i + 1

    if len(Z) == 3: 
        Z = np.array([Z[np.argmax(Z)], Z[np.argmin(Z)]])
            
    Z = Z.real

    return Z

def model_SRK(P, T, Tc, Pc, omega):
    a = 0.42748*R**2*Tc**2/Pc
    b = 0.08664*R*Tc/Pc
    m = 0.48 + 1.574*omega - 0.176*omega**2
    alfa = (1 + m*(1 - np.sqrt(T/Tc)))**2
    A = a*alfa*P/(R*T)**2
    B = b*P/(R*T)

    Z = np.roots([1 , -1 , (A - B - B**2) , -A*B])

    for i in range(3):

        if Z[i].imag != 0 or Z[i].real < 0 or Z[i].real < B or Z[i].real == 0:

            Z[i] = 0

    i = 0
    while i < len(Z):
        if Z[i] == 0:
            Z = np.delete(Z,i)

        else:
            i = i + 1

    if len(Z) == 3: 
        Z = np.array([Z[np.argmax(Z)], Z[np.argmin(Z)]])
    
    Z = Z.real

    return Z

def model_PR(P, T, Tc, Pc, omega):
    a = 0.45724*((R*Tc)**2/Pc)
    b = 0.0778*R*Tc/Pc
    k = 0.37464 + 1.54226*omega - 0.26992*omega**2
    alfa = (1 + k*(1 - np.sqrt(T/Tc)))**2
    A = a*alfa*P/(R*T)**2
    B = b*P/(R*T)

    Z = np.roots([1 , -(1 - B) , A - 3*B**2 - 2*B , -(A*B - B**2 - B**3)])

    for i in range(3):

        if Z[i].imag != 0 or Z[i].real < 0 or Z[i].real < B or Z[i].real == 0:

            Z[i] = 0

    i = 0
    while i < len(Z):
        if Z[i] == 0:
            Z = np.delete(Z,i)

        else:
            i = i + 1

    if len(Z) == 3: 
        Z = np.array([Z[np.argmax(Z)], Z[np.argmin(Z)]])
 
    Z = Z.real

    return Z

#Las funciones descartan los valores de Z que no tienen sentido físico. 
#Puede ocurrir que Z tenga uno, dos o tres valores con sentido físico
#Si dada una Z hay 3 valores con sentido físico, se descarta la intermedia
#Si dada una Z hay 2 valores con sentido físico, el mayor entre los dos se corresponde a un estado vapor, el menor al líquido. 


input_ok = False
while not input_ok:
    try:
        num_sustancias = int(input('Cuantas sustancias quieres analizar? '))

        while not num_sustancias > 0:
            print('Introduzca un número válido de sustancias (1,2,3...): ')
            num_sustancias = int(input('Cuantas sustancias quieres analizar? '))
            
        input_ok = True
        
    except ValueError:
        print('Has introducido un número de sustancias incorrecto')

contador = 1

for i in range(num_sustancias):  

    try:  
        
        V_real = []

        sustancia = str(input('Que sustancia quieres analizar? Usa un nombre válido: ')).strip()

        nombre = sustancia
        repetidos = 2

        while nombre in list_sustancias: #Si se cumple esta condición, la sustancia 'nombre' ya se encuentra en list_sustancias

            nombre = sustancia + ' ' f'{repetidos}' 
            #Si se introduce la sustancia 'Water' por tercera vez, se renombra a 'Water 3'.
            repetidos = repetidos + 1 

        list_sustancias.append(nombre)

        P = float(input('A que presión está el sistema? (Pa): '))
        while P <= 0:
            print('Introduce un valor correcto de la presión')
            P = float(input('A que presión está el sistema? (Pa): '))
            
        T = float(input('A que temperatura está el sistema? (K): '))
        while T <= 0:
            print('Introduce un valor correcto de la temperatura')
            T = float(input('A que temperatura está el sistema? (K): '))

        try: 

            fase = PhaseSI('P',P,'T',T,sustancia)

#Por como funciona Phase SI, cuando el par (P,T) escogidos sean (Psat,Tsat), PhaseSI devuelve un ValueError
#Es decir, cada vez que PhaseSI devuelve un ValueError surgirá debido a la existéncia de dos fases en el sistema

        except ValueError:

            fase = 'twophase'

        Tc = PropsSI('Tcrit',sustancia) #Temperatura crítica
        Pc = PropsSI('Pcrit',sustancia) #Presión crítica
        
        omega = PropsSI('acentric',sustancia)

        ZGI = model_GI()
        ZVDW = model_VdW(P,T,Tc,Pc)
        ZRK = model_RK(P,T,Tc,Pc)
        ZSRK = model_SRK(P,T,Tc,Pc,omega)
        ZPR = model_PR(P,T,Tc,Pc,omega)

        Z_modelos = [ZGI, ZVDW, ZRK, ZSRK, ZPR]
        V_modelos = []
        e_rel_modelos = []

        if fase == 'gas':

            V_real.append(1/PropsSI('Dmolar','P',P,'T',T,sustancia))

        elif fase == 'liquid':
            
            V_real.append(1/PropsSI('Dmolar','P',P,'T',T,sustancia))

        elif fase == 'twophase':

            Tit_vapor = str(input('Quiere introducir el título de vapor? Y/N: ')).upper()

            while Tit_vapor not in ['Y', 'N']:

                print('Entrada inválida. Introduzca Y/N o y/n.')
                Tit_vapor = str(input('Quiere introducir el título de vapor? Y/N: ')).upper()

            if Tit_vapor == 'Y':

                Q = float(input('Cuál es el título de vapor del sistema?\n'
                                '0 ---> 100% Líquido\n'
                                '1 ---> 100% Vapor/Gas\n' 
                                'Mayor que 0 y menor que 1 ---> Mezcla bifásica: '))
                
                while Q < 0 or Q > 1:

                    print('Valor de Q incorrecto. Q está contenido en el intervalo [0,1].')
                    Q = float(input('Introduzca el título de vapor del sistema: '))

                V_real.append(1/PropsSI('Dmolar','T',T,'Q',Q,sustancia))

            else:

                V_real.append(1/PropsSI('Dmolar','T',T,'Q',0,sustancia))
                V_real.append(1/PropsSI('Dmolar','T',T,'Q',1,sustancia))

        else: #Si la sustancia no es ni líquido, gas, ni mezcla bifásica ---> Si la sustancia está en condiciones supercríticas:

            V_real.append(1/PropsSI('Dmolar','P',P,'T',T,sustancia))

        for j , Z_modelo in enumerate(Z_modelos):

            if len(Z_modelo) > 1: #Si el modelo no es el de los Gases Ideales

                if fase == 'gas':

                    Z_modelo = np.delete(Z_modelo, np.argmin(Z_modelo)) #Elimina el Z del líquido (el de mínimo valor)

                elif fase == 'liquid':

                    Z_modelo = np.delete(Z_modelo, np.argmax(Z_modelo)) #Elimina el Z del gas (el de máximo valor)

                elif fase == 'twophase':

                    Z_modelo = np.sort(Z_modelo) 
                    #El primer elemento de Z_modelo, después de hacer el np.sort, es Z líquido, y el segundo es Z gas.

                    if Tit_vapor == 'Y':

                        Z_modelo = Z_modelo[0] + Q*(Z_modelo[1] - Z_modelo[0])
                        Z_modelo = np.array([Z_modelo])

            Z_modelos[j] = Z_modelo 

            if len(Z_modelo) == 0: #Evitamos que se quede vacío y rompa la lógica de las listas
                V_modelos.append('N/A (Sin sentido físico)')
                e_rel_modelos.append('N/A')
            else:
                V_modelo = (Z_modelo*R*T/P)
                V_modelos.append(V_modelo)

                e_rel_modelo = 100 * ((np.array(V_real) - np.array(V_modelo)) / np.array(V_real))
                e_rel_modelos.append(e_rel_modelo) 

    #En este punto, para sustancias con Q = 0 o Q = 1 tenemos una única raíz de Z. Si Q != 0 y Q != 1, tenemos 2 raíces de Z

    #V es el Volumen calculado según el modelo usado, mientras V_real es el volumen medido en el laboratorio
    #V_modelo contiene los datos calculados de 1 modelo para una sustancia
    #V_modelos es una lista que contiene cada V_modelo para una sustancia
                
        Resultados[list_sustancias[i]] = {
            'Z teórica:':Z_modelos,
            'Volumen teórico':V_modelos,
            'Error relativo:':e_rel_modelos
            }

    except ValueError:

        print('Has introducido algun input incorrectamente. La sustancia afectada es:',sustancia,)