
#TRATAMIENTO DE INFORMACION#
import numpy as np
from numpy.core.fromnumeric import transpose
import pandas as pd 
import matplotlib.pyplot as plt
from scipy import signal

data = pd.read_csv("t3.csv", delimiter=",")  

print(data)


data = data.replace(',','.', regex=True)

## SEÑALES DE CORRIENTE Y TENSIÓN#

# Vector de tiempos
time_vector = data.iloc[:,0].astype(float).to_numpy()

# Señales de tensión/corriente
current_deltaA = data.iloc[:, 1].astype(float).to_numpy()
voltage_deltaA =data.iloc[ :,4].astype(float).to_numpy()

current_deltaB = data.iloc[:, 2].astype(float).to_numpy()
voltage_deltaB =data.iloc[ :,5].astype(float).to_numpy()

current_deltaC = data.iloc[:, 3].astype(float).to_numpy()
voltage_deltaC =data.iloc[ :,6].astype(float).to_numpy()



# ADC
delta_t = 0.641e-06
k = 1/ (1.56e06/2)
theta = 0

# IMPEDANCIA CARACTERISTICA
Z = (20.631 + 89.93j)/110
Y = (676.75e-06j)/110


ZC = np.sqrt(Z/Y)

# print(ZC)

#  TRANSFORMACION MODAL

w = 2*np.pi*60
t = time_vector


# output_modal = output*T_dq0

# print(output_modal)



# print(transform_Matrix)
# # ONDAS VIAJERAS DE LLEGADA Y SALIDA

# current_delta_transform = (1/3)*current_delta.dot(transform_Matrix)
# voltage_delta_transform = (1/3)*voltage_delta.dot(transform_Matrix)

# S1 = voltage_delta_transform + ZC*current_delta_transform 
# S2 = voltage_delta_transform - ZC*current_delta_transform 

# print(voltage_delta_transform)

# plt.plot(current_deltaA[:,0],current_deltaA[:,1])
# plt.plot(outputA[:,0],outputA[:,1])

plt.plot(t, current_deltaB)
plt.plot(t, current_deltaA)
plt.plot(t, current_deltaC)
# plt.plot(t,outputB[:,1])
plt.xlabel('Time (s)')
plt.ylabel('Current (kA)')
plt.show()

# AGG 2
# print(time_vector)
# contador = np.size(current_deltaA)
contador = len(time_vector)
contador_arr = np.linspace(1,contador,contador)
print('Hola')
print(contador_arr)
Ad_k_vector = np.zeros((contador))
Aq_k_vector = np.zeros((contador))
A0_k_vector = np.zeros((contador))

Vd_k_vector = np.zeros((contador))
Vq_k_vector = np.zeros((contador))
V0_k_vector = np.zeros((contador))

# print('Loop')
for i in range(contador):
#   print('Iteración')
  # Vector ABC
#   print('Impresion de Aabc_k')
  Aabc_k = np.array([[current_deltaA[i]], [current_deltaB[i]], [current_deltaC[i]]])
  Vabc_k = np.array([[voltage_deltaA[i]], [voltage_deltaB[i]], [voltage_deltaC[i]]])
#   print(Aabc_k)

  # Construcción Matriz Tdq0
#   print('Impresion de phi')

  phi = i*delta_t*2*np.pi*60 + theta
#   print(phi)
#   print('Impresion de Matriz')
  T_dq0 =(2/3) * np.array([[np.sin(phi),np.sin(phi-(2/3)*np.pi), np.sin(phi+(2/3)*np.pi)],[np.cos(phi),np.cos(phi-(2/3)*np.pi), np.cos(phi+(2/3)*np.pi) ],[0.5, 0.5,0.5]])
#   print(T_dq0)
  # Cálculo de la trasnformada
#   print('Impresion de Park')
  Adq0_k = np.matmul(T_dq0, Aabc_k)
  Vdq0_k = np.dot(T_dq0, Vabc_k)
#   print(Adq0_k)

  # Almacenamiento de valores en vectores
#   print('Impresion de Almacenamiento vector')
  Ad_k = Adq0_k[0]
  Aq_k = Adq0_k[1]
  A0_k = Adq0_k[2]


  Vd_k = Vdq0_k[0]
  Vq_k = Vdq0_k[1]
  V0_k = Vdq0_k[2]
#   print(Ad_k)

  Ad_k_vector[i] = Ad_k
  Aq_k_vector[i] = Aq_k
  A0_k_vector[i] = A0_k

  Vd_k_vector[i] = Vd_k
  Vq_k_vector[i] = Vq_k
  V0_k_vector[i] = V0_k
#   print(Ad_k_vector)
  if i == (contador - 1):
        print("The end")
        break
#   print('Número de iteración')    
#   print(i)
#   print('Impresion contador')
#   print(contador)

# print('Impresion componente directa')
# print(Ad_k_vector)
# print('Impresion componente cuadratura')
# print(Aq_k_vector)
# print('Impresion componente eje 0')
# print(A0_k_vector)

print(Aq_k_vector)
print(Ad_k_vector)
print(A0_k_vector)

plt.clf()

plt.plot(t, Ad_k_vector)
plt.plot(t, Aq_k_vector)
plt.plot(t, A0_k_vector)
# plt.plot(t,outputB[:,1])
plt.xlabel('Time (s)')
plt.ylabel('Current (kA)')
plt.show()

S_forward = Vd_k_vector + ZC*Ad_k_vector
S_backward = Vd_k_vector - ZC*Ad_k_vector

#FILTER
  # Cut-off frequency of the filter

tau = 1560000
tau_arr = np.array(tau)
w = 20e03/tau
# Normalize the frequency
b1, a1 = signal.butter(3, w, 'high')
S_f = signal.filtfilt(b1, a1, S_forward,axis=0)

b2, a2 = signal.butter(3, w, 'high')
S_b = signal.filtfilt(b2, a2, S_backward,axis=0)


plt.plot(t, S_f)
plt.plot(t, S_b)
# plt.plot(t,outputB[:,1])
plt.xlabel('Time (s)')
plt.ylabel('Current (kA)')
plt.show()

index_S = np.argwhere(time_vector == 0.05).astype(int)
index_S = index_S[0,0]
print(index_S)

S_template = S_f[index_S:round(index_S+contador/2)]

# def corr_fn(tau):
#   phi = 0
#   delta_t = 1e-04
#   k = 1
#   M = round(contador/2)
#   while k<M:
#     phi = phi+S_backward[k+tau]*S_template[k]
#     k = k+1
#     print(k)
#   return(phi)


corr = np.correlate(S_template,S_b,mode="same")
corr = corr*np.abs(corr)
print(np.argmax(np.abs(corr),axis=0))
tau_f = np.argmax(np.abs(corr),axis=0)

d = ((298.863e03*tau_f*delta_t)/(2))*1e-03
print(d)

ventana = tkinter.Tk()
ventana.geometry("400x300")
ventana.mainloop()









