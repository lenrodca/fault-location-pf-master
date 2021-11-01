from datetime import time
from zipfile import ZIP_MAX_COMMENT
import numpy as np
from numpy.core.fromnumeric import transpose
from numpy.lib.npyio import zipfile_factory
import pandas as pd 
import matplotlib.pyplot as plt
from scipy import signal
import tkinter
from tkinter import StringVar, filedialog, Canvas, ttk
from matplotlib.backends.backend_tkagg import (
    FigureCanvasTkAgg, NavigationToolbar2Tk)
# Implement the default Matplotlib key bindings.
from matplotlib.backend_bases import key_press_handler
from matplotlib.figure import Figure
from scipy.signal import find_peaks
import heapq

import csv
global filename
filename = ""
global delta_t
global velprop
global Zc
global faultdistance
global combo


win = tkinter.Tk()
d= ""

def UploadAction():
    global filename
    filename = filedialog.askopenfilename()
    

def FaultLoc():
    global time_vector
    global current_deltaA, current_deltaB,current_deltaC
    global voltage_deltaA, voltage_deltaB, voltage_deltaC
    global Vd_k, Vq_k,V0_k
    global Ad_k, Aq_k, A0_k
    global S_b, S_f
    global corr
    global d 

    data = pd.read_csv(filename, delimiter=",") 

    print(data)

    data = data.replace(',','.', regex=True)

    ## SEÑALES DE CORRIENTE Y TENSIÓN#

    # Vector de tiempos
    time_vector = data.iloc[:,0].astype(float).to_numpy()

    # Señales de tensión/corriente
    current_deltaA = data.iloc[:, 4].astype(float).to_numpy()
    voltage_deltaA =data.iloc[ :,1].astype(float).to_numpy()

    current_deltaB = data.iloc[:, 5].astype(float).to_numpy()
    voltage_deltaB =data.iloc[ :,2].astype(float).to_numpy()

    current_deltaC = data.iloc[:, 6].astype(float).to_numpy()
    voltage_deltaC =data.iloc[ :,3].astype(float).to_numpy()

    # Transformada de PARK
    contador = len(time_vector)

    Ad_k_vector = np.zeros((contador))
    Aq_k_vector = np.zeros((contador))
    A0_k_vector = np.zeros((contador))

    Vd_k_vector = np.zeros((contador))
    Vq_k_vector = np.zeros((contador))
    V0_k_vector = np.zeros((contador))
    
    theta = 0

    for i in range(contador):

        Aabc_k = np.array([[current_deltaA[i]], [current_deltaB[i]], [current_deltaC[i]]])
        Vabc_k = np.array([[voltage_deltaA[i]], [voltage_deltaB[i]], [voltage_deltaC[i]]])
        
        phi = i*delta_t*2*np.pi*60 + theta

        T_dq0 = (2/3) * np.array([[np.cos(phi),np.cos(phi-(2/3)*np.pi), np.cos(phi+(2/3)*np.pi)],[-np.sin(phi),-np.sin(phi-(2/3)*np.pi), -np.sin(phi+(2/3)*np.pi)],[0.5, 0.5,0.5]])
  
        Adq0_k = np.dot(T_dq0, Aabc_k)
        Vdq0_k = np.dot(T_dq0, Vabc_k)

        Ad_k = Adq0_k[0]
        Aq_k = Adq0_k[1]
        A0_k = Adq0_k[2]


        Vd_k = Vdq0_k[0]
        Vq_k = Vdq0_k[1]
        V0_k = Vdq0_k[2]

        Ad_k_vector[i] = Ad_k
        Aq_k_vector[i] = Aq_k
        A0_k_vector[i] = A0_k

        Vd_k_vector[i] = Vd_k
        Vq_k_vector[i] = Vq_k
        V0_k_vector[i] = V0_k

    # Ondas Viajeras
    S_forward = -Vd_k_vector + Zc*Ad_k_vector
    S_backward = Vd_k_vector + Zc*Ad_k_vector

    #Filtro
    tau = 1.56e06
    w = 20e03/tau
    b1, a1 = signal.butter(3, w, 'high')
    S_f = signal.filtfilt(b1, a1, S_forward,axis=0)

    b2, a2 = signal.butter(3, w, 'high')
    S_b = signal.filtfilt(b2, a2, S_backward,axis=0)

    #Correlacion

    corr = np.correlate(S_b,S_f,mode="same")/contador
    peaks, _ = find_peaks(corr)

    taus = np.zeros(len(peaks))
    for i in range(0,len(peaks)):
        taus[i] = corr[peaks[i]]

    print(taus)

    max_taus = heapq.nlargest(2, taus)
   
    tau_1 = np.argwhere(taus == max_taus[0]).astype(int)
    tau_2 = np.argwhere(taus== max_taus[1]).astype(int)
    # tau_f = np.argmax(corr)

    d = ((velprop*(np.abs(time_vector[peaks[tau_1]]-time_vector[peaks[tau_2]])))/(2))
    d = round(d[0,0],4)
    print(d)
    labelText.set(d)

    error = (np.abs(faultdistance-d)/(faultdistance))*100
    error = round(error,4)
    labelText2.set(error)


win.geometry("1280x720")
win.title("Aplicación para localización de fallas - PF202130")


label1 = tkinter.Label(win,text = "Ingrese el valor de la frecuencia de muestreo:").place(x = 20, y = 50)  
label2 = tkinter.Label(win,text = "Ingrese el valor de la velocidad de propagación:").place(x = 20, y = 80) 
label3 = tkinter.Label(win,text = "Ingrese el valor de la impedancia característica:").place(x = 20, y = 110) 
label4 = tkinter.Label(win,text = "Ingrese la distancia de falla:").place(x = 20, y = 140) 
label5 = tkinter.Label(win,text = "Programa para la localización de fallas utilizando ondas viajeras:", font=('Helvetica', 18, 'bold')).place(x = 360, y = 10)  

labelText = tkinter.StringVar()
labelText2 = tkinter.StringVar()

label_d1 = tkinter.Label(win,text = "La distancia de la falla calculada es:", fg = "black", font =('Helvetica', 12)).place(x = 500, y = 50 )
label_d2 = tkinter.Label(win,textvariable= labelText, fg = "green", font =('Helvetica', 12, 'bold')).place(x = 755, y = 50 )

label_e1 = tkinter.Label(win,text = "El error en (%) calculado es:", fg = "black", font =('Helvetica', 12)).place(x = 500, y = 70 )
label_e2 = tkinter.Label(win,textvariable= labelText2, fg = "blue", font =('Helvetica', 12, 'bold')).place(x = 755, y = 70 )

entry = tkinter.Entry(win)
entry.place(x = 280, y = 50)
entry2 = tkinter.Entry(win)
entry2.place(x = 280, y = 80)
entry3 = tkinter.Entry(win)
entry3.place(x = 280, y = 110)
entry4 = tkinter.Entry(win)
entry4.place(x = 280, y = 140)

def VarReading():
    global delta_t
    global velprop
    global Zc
    global faultdistance
    
    # delta_t = float(entry.get())
    delta_t = 0.641E-6
    # velprop = float(entry2.get())
    velprop = 2.92051e05
    # Zc = float(entry3.get())
    Zc = 369.498
    faultdistance = float(entry4.get())
    
button = tkinter.Button(win, text='Abrir archivo .csv', command=UploadAction)
button.place(x = 1120, y = 50)
buttonrun = tkinter.Button(win, text='Ejecutar Script', command=lambda:[VarReading(), FaultLoc()])
buttonrun.place(x = 1130, y = 80)


widget_var = tkinter.StringVar()
combo = ttk.Combobox(win,textvariable=widget_var,values=["Señales de entrada", "Señales transformadas", "Ondas Viajeras", "Correlación"], state="readonly")
combo.set("Seleccione la gráfica...")
combo.place(x = 1100, y = 140)



def Graph():
    value = combo.get()

    if value == "Señales de entrada":
        print(value)
        
        fig1 = Figure(figsize=(6,8), dpi=100)
        canvas = FigureCanvasTkAgg(fig1, master=win)  # A tk.DrawingArea.
        
        canvas.get_tk_widget().pack(padx=50, pady=170)
        toolbar = NavigationToolbar2Tk(canvas, win)
        toolbar.update()
        canvas.get_tk_widget().pack(side=tkinter.TOP, fill=tkinter.BOTH, expand=1)
       
        plot1 = fig1.add_subplot(1,2,2)
        plot1.plot(time_vector, current_deltaA, 'r')
        plot1.plot(time_vector, current_deltaB, 'b')
        plot1.plot(time_vector, current_deltaC,'g')
         # Agregar punto de falla en inputs
        plot1.axvline(x=0.016, color='k', linestyle='--')
        plot1.legend(labels=['Ph-A', 'Ph-B', 'Ph-C','Falla'],loc='lower left')
        plot1.title.set_text("Señales de corriente - Entrada")
        plot1.set_xlabel('Tiempo (s)')
        plot1.set_ylabel('Corriente (A)')


        plot2 = fig1.add_subplot(1,2,1)
        plot2.plot(time_vector, voltage_deltaA, 'r')
        plot2.plot(time_vector, voltage_deltaB, 'b')
        plot2.plot(time_vector, voltage_deltaC,'g')
        # Agregar punto de falla en inputs
        plot2.axvline(x=0.016, color='k', linestyle='--')
        plot2.legend(labels=['Ph-A', 'Ph-B', 'Ph-C','Falla'],loc='lower left')
        plot2.title.set_text("Señales de tensión - Entrada")
        plot2.set_xlabel('Tiempo (s)')
        plot2.set_ylabel('Voltaje (V)')

        fig1.canvas.draw_idle()
        
        

        

    
    elif value == "Señales transformadas":
        print(value)
        
        fig2 = Figure(figsize=(6,8), dpi=100)
        canvas2 = FigureCanvasTkAgg(fig2, master=win)  # A tk.DrawingArea.
       
        canvas2.get_tk_widget().pack(padx=50, pady=170)
        toolbar = NavigationToolbar2Tk(canvas2, win)
        toolbar.update()
        canvas2.get_tk_widget().pack(side=tkinter.TOP, fill=tkinter.BOTH, expand=1)

        plot1 = fig2.add_subplot(1,2,2)
        plot1.plot(Ad_k, 'r')
        plot1.plot(Aq_k, 'b')
        plot1.plot(A0_k,'g')
         # Agregar punto de falla en inputs
        plot1.axvline(x=0.016, color='k', linestyle='--')
        plot1.legend(labels=['Ph-A', 'Ph-B', 'Ph-C','Falla'],loc='lower left')
        plot1.title.set_text("Señales de corriente - Transformadas")
        plot1.set_xlabel('Tiempo (s)')
        plot1.set_ylabel('Corriente (A)')


        plot2 = fig2.add_subplot(1,2,1)
        plot2.plot(Vd_k, 'r')
        plot2.plot(Vq_k, 'b')
        plot2.plot(V0_k,'g')
        # Agregar punto de falla en inputs
        plot2.axvline(x=0.016, color='k', linestyle='--')
        plot2.legend(labels=['Ph-A', 'Ph-B', 'Ph-C','Falla'],loc='lower left')
        plot2.title.set_text("Señales de tensión - Transformadas")
        plot2.set_xlabel('Tiempo (s)')
        plot2.set_ylabel('Voltaje (V)')

        fig2.canvas.draw_idle()

        

    


    elif value == "Ondas Viajeras":
        print(value)
        
        fig3 = Figure(figsize=(6,8), dpi=100)
        canvas3 = FigureCanvasTkAgg(fig3, master=win)  # A tk.DrawingArea.
        
        canvas3.get_tk_widget().pack(padx=50, pady=170)
        toolbar = NavigationToolbar2Tk(canvas3, win)
        toolbar.update()
        canvas3.get_tk_widget().pack(side=tkinter.TOP, fill=tkinter.BOTH, expand=1)

        plot1 = fig3.add_subplot(111)
        plot1.plot(time_vector,S_f)
        plot1.plot(time_vector,S_b)
        plot1.legend(labels=['Onda de ida', 'Onda de llegada'],loc='lower left')
        plot1.title.set_text("Ondas viajeras")
        plot1.set_xlabel('Tiempo (s)')
        plot1.set_ylabel('Voltaje (V)')

        fig3.canvas.draw_idle()

       

   
        
    elif value == "Correlación":
        print(value)

buttongraph = tkinter.Button(win, text='Realizar Gráfica', command=Graph)
buttongraph.place(x = 1125, y = 110)

win.mainloop()




