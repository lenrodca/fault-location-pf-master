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
import matplotlib.backends.backend_tkagg as tkagg
# Implement the default Matplotlib key bindings.
from matplotlib.backend_bases import key_press_handler
from matplotlib.figure import Figure
from scipy.signal import find_peaks
import heapq
from scipy.signal import argrelextrema
from PIL import ImageTk, Image
from modwt import modwt

import csv
global filename
filename = ""
global delta_t
global velprop
global Zc
global faultdistance
global combo,combo2,combo22,combo3
global canvas, canvas2


win = tkinter.Tk()
tabControl = ttk.Notebook(win)

tab1 = ttk.Frame(tabControl)
tab2 = ttk.Frame(tabControl)
tab3 = ttk.Frame(tabControl)

tabControl.add(tab1, text ='Localización de fallas')
tabControl.add(tab2, text ='Validación del método')
tabControl.add(tab3, text ='Manual de Usuario')
tabControl.pack(expand = 1, fill ="both")


def UploadAction():
    global filename
    filename = filedialog.askopenfilename()

def UploadAction2():
    global filename_arrays
    filename_arrays = filedialog.askopenfilenames()
    print(len(filename_arrays))
    

def FaultLoc():
    global d 
    global delta_t
    global velprop
    global Zc

    global line_distance, error 

    global time_vector
    global current_deltaA, current_deltaB,current_deltaC
    global voltage_deltaA, voltage_deltaB, voltage_deltaC
    global Vd_k_vector, Vq_k_vector,V0_k_vector
    global Ad_k_vector, Aq_k_vector, A0_k_vector
    global S_b, S_f,S_backward_modwt,S_forward_modwt
    global corr,peaks,tau_1,tau_2,max_taus,corr2,peaks2,tau_11,tau_22,max_taus2

    data = pd.read_csv(filename, delimiter=",") 

    data = data.replace(',','.', regex=True)

    print(data)

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

    #MODWT 
    wt = modwt(Ad_k_vector, 'db4', 4)
    wt_A = wt[3,:]
    wt = modwt(Vd_k_vector, 'db4', 4)
    wt_V = wt[3,:]

    # Ondas Viajeras
    S_forward = -Vd_k_vector + Zc*Ad_k_vector
    S_backward = Vd_k_vector + Zc*Ad_k_vector

    S_forward_modwt = -wt_V + Zc*wt_A
    S_backward_modwt = wt_V + Zc*wt_A

    #Filtro
    tau = 1.56e06
    w = 20e03/tau
    b1, a1 = signal.butter(3, w, 'high')
    S_f = signal.filtfilt(b1, a1, S_forward,axis=0)
    S_forward_modwt = signal.filtfilt(b1, a1, S_forward_modwt,axis=0)
    b2, a2 = signal.butter(3, w, 'high')
    S_b = signal.filtfilt(b2, a2, S_backward,axis=0)
    S_backward_modwt = signal.filtfilt(b2, a2, S_backward_modwt,axis=0)

    #Correlacion

    corr = np.correlate(S_b,S_f,mode="same")/contador
    corr2 = np.correlate(S_backward_modwt, S_forward_modwt,mode="same")/contador
    peaks, _ = find_peaks(corr)
    peaks2, _ = find_peaks(corr2,distance=50)

    taus = np.zeros(len(peaks))
    taus2 = np.zeros(len(peaks2))
    for i in range(0,len(peaks)):
        taus[i] = corr[peaks[i]]
        
    for i in range(0,len(peaks2)):
        taus2[i] = corr2[peaks2[i]]

    # MAXIMOS CORRELACION CONVENCIONAL
    max_taus = heapq.nlargest(2, taus)

    tau_1 = np.argwhere(taus == max_taus[0]).astype(int)
    tau_2 = np.argwhere(taus == max_taus[1]).astype(int)

    # MAXIMOS MODWT
    max_taus2 = heapq.nlargest(2, taus2)

    tau_11 = np.argwhere(taus2 == max_taus2[0]).astype(int)
    tau_22 = np.argwhere(taus2 == max_taus2[1]).astype(int)

    print(time_vector[peaks[tau_1]])
    print(time_vector[peaks[tau_2]])

    #DISTANCIA CORRELACION CONVENCIONAL
    d = ((velprop*(np.abs(time_vector[peaks[tau_1]]-time_vector[peaks[tau_2]])))/(2))
    d = d[0,0]

    #DISTANCIA MODWT
    d2 = ((velprop*(np.abs(time_vector[peaks2[tau_11]]-time_vector[peaks2[tau_22]])))/(2))
    d2 = d2[0,0]

    d2 = str(d2) + ' ' + 'km'
    d = str(d) + ' ' + 'km'
    print (d)
    print(d2)
    
    value_d = combo2.get()
    if (value_d == "Correlación convencional"):
        labelText.set(d)
    elif(value_d == "MODWT"):
        labelText.set(d2)
    global label_d1    
    label_d1 = tkinter.Label(tab1,text = "La distancia de la falla calculada es...", fg = "black", font =('Helvetica', 12)).place(x = 550, y = 65)    

    

# error = (np.abs(faultdistance-d)/(faultdistance))*100
# error = round(error,4)
# labelText2.set(error)

def FaultLoc2():
    global d
    global line_distance, error 
    global filename_arrays
    global delta_t, Zc,velprop
    global error_array,vector_distancia
    value_met = combo22.get()
    calculo = int(100/(len(filename_arrays)))
    print('calculo')
    print(calculo)
    count = 0
    vector_distancia = np.zeros(len(filename_arrays))
    for i in range(len(vector_distancia)):
        vector_distancia[i] = 100/(len(filename_arrays))*(i+1)
    
    error_array = np.zeros(len(filename_arrays))

    
    
    for x in filename_arrays:
        data_cyc = pd.read_csv(x, delimiter=",") 

        data_cyc = data_cyc.replace(',','.', regex=True)

        ## SEÑALES DE CORRIENTE Y TENSIÓN#

        # Vector de tiempos
        time_vector_cyc = data_cyc.iloc[:,0].astype(float).to_numpy()

        # Señales de tensión/corriente
        current_deltaA_cyc = data_cyc.iloc[:, 4].astype(float).to_numpy()
        voltage_deltaA_cyc = data_cyc.iloc[ :,1].astype(float).to_numpy()

        current_deltaB_cyc = data_cyc.iloc[:, 5].astype(float).to_numpy()
        voltage_deltaB_cyc =data_cyc.iloc[ :,2].astype(float).to_numpy()

        current_deltaC_cyc = data_cyc.iloc[:, 6].astype(float).to_numpy()
        voltage_deltaC_cyc =data_cyc.iloc[ :,3].astype(float).to_numpy()

        # Transformada de PARK
        contador_cyc = len(time_vector_cyc)

        Ad_k_vector_cyc = np.zeros((contador_cyc))
        Aq_k_vector_cyc = np.zeros((contador_cyc))
        A0_k_vector_cyc = np.zeros((contador_cyc))

        Vd_k_vector_cyc = np.zeros((contador_cyc))
        Vq_k_vector_cyc = np.zeros((contador_cyc))
        V0_k_vector_cyc = np.zeros((contador_cyc))
        
        theta_cyc = 0

        for i in range(contador_cyc):

            Aabc_k_cyc = np.array([[current_deltaA_cyc[i]], [current_deltaB_cyc[i]], [current_deltaC_cyc[i]]])
            Vabc_k_cyc = np.array([[voltage_deltaA_cyc[i]], [voltage_deltaB_cyc[i]], [voltage_deltaC_cyc[i]]])
            
            phi_cyc = i*delta_t*2*np.pi*60 + theta_cyc

            T_dq0_cyc = (2/3) * np.array([[np.cos(phi_cyc),np.cos(phi_cyc-(2/3)*np.pi), np.cos(phi_cyc+(2/3)*np.pi)],[-np.sin(phi_cyc),-np.sin(phi_cyc-(2/3)*np.pi), -np.sin(phi_cyc+(2/3)*np.pi)],[0.5, 0.5,0.5]])
    
            Adq0_k_cyc = np.dot(T_dq0_cyc, Aabc_k_cyc)
            Vdq0_k_cyc = np.dot(T_dq0_cyc, Vabc_k_cyc)

            Ad_k_cyc = Adq0_k_cyc[0]
            Aq_k_cyc = Adq0_k_cyc[1]
            A0_k_cyc = Adq0_k_cyc[2]


            Vd_k_cyc = Vdq0_k_cyc[0]
            Vq_k_cyc = Vdq0_k_cyc[1]
            V0_k_cyc = Vdq0_k_cyc[2]

            Ad_k_vector_cyc[i] = Ad_k_cyc
            Aq_k_vector_cyc[i] = Aq_k_cyc
            A0_k_vector_cyc[i] = A0_k_cyc

            Vd_k_vector_cyc[i] = Vd_k_cyc
            Vq_k_vector_cyc[i] = Vq_k_cyc
            V0_k_vector_cyc[i] = V0_k_cyc


        #MODWT 
        wt_cyc = modwt(Ad_k_vector_cyc, 'db4', 4)
        wt_A_cyc = wt_cyc[3,:]
        wt_cyc = modwt(Vd_k_vector_cyc, 'db4', 4)
        wt_V_cyc = wt_cyc[3,:]

        # Ondas Viajeras
        S_forward_cyc = -Vd_k_vector_cyc + Zc*Ad_k_vector_cyc
        S_backward_cyc = Vd_k_vector_cyc + Zc*Ad_k_vector_cyc

        S_f_cyc_2 = -wt_V_cyc + Zc*wt_A_cyc
        S_b_cyc_2 = wt_V_cyc + Zc*wt_A_cyc

        #Filtro
        tau = 1.56e06
        w = 20e03/tau
        b1, a1 = signal.butter(3, w, 'high')
        S_f_cyc = signal.filtfilt(b1, a1, S_forward_cyc,axis=0)
        S_f_cyc_2 = signal.filtfilt(b1, a1, S_f_cyc_2,axis=0)

        b2, a2 = signal.butter(3, w, 'high')
        S_b_cyc = signal.filtfilt(b2, a2, S_backward_cyc,axis=0)
        S_b_cyc_2 = signal.filtfilt(b2, a2, S_b_cyc_2,axis=0)

        #Correlacion

        corr_cyc = np.correlate(S_b_cyc,S_f_cyc,mode="same")/contador_cyc
        peaks_cyc, _ = find_peaks(corr_cyc)

        corr_cyc_2 = np.correlate(S_b_cyc_2,S_f_cyc_2,mode="same")/contador_cyc
        peaks_cyc_2, _ = find_peaks(corr_cyc_2,distance=50)

        taus_cyc = np.zeros(len(peaks_cyc))
        for i in range(0,len(peaks_cyc)):
            taus_cyc[i] = np.abs(corr_cyc[peaks_cyc[i]])

        taus_cyc_2 = np.zeros(len(peaks_cyc_2))
        for i in range(0,len(peaks_cyc_2)):
            taus_cyc_2[i] = np.abs(corr_cyc_2[peaks_cyc_2[i]])


        max_taus_cyc = heapq.nlargest(2, taus_cyc)
    
        tau_1_cyc = np.argwhere(taus_cyc == max_taus_cyc[0]).astype(int)
        tau_2_cyc = np.argwhere(taus_cyc== max_taus_cyc[1]).astype(int)

        max_taus_cyc_2 = heapq.nlargest(2, taus_cyc_2)
    
        tau_1_cyc_2 = np.argwhere(taus_cyc_2 == max_taus_cyc_2[0]).astype(int)
        tau_2_cyc_2 = np.argwhere(taus_cyc_2== max_taus_cyc_2[1]).astype(int)
        # tau_f = np.argmax(corr)

        d = ((velprop*(np.abs(time_vector_cyc[peaks_cyc[tau_1_cyc]]-time_vector_cyc[peaks_cyc[tau_2_cyc]])))/(2))
        d = round(d[0,0],4)
        print(d)

        d2 = ((velprop*(np.abs(time_vector_cyc[peaks_cyc_2[tau_1_cyc_2]]-time_vector_cyc[peaks_cyc_2[tau_2_cyc_2]])))/(2))
        d2 = round(d2[0,0],4)
        print(d2)

        if (value_met == "Correlación convencional"):
            d = d
        elif (value_met == "MODWT"):
            d = d2
        
        error = ((np.abs(((vector_distancia[count])/100)*line_distance - d))/(((vector_distancia[count])/100)*line_distance))*100
        error_array[count] = error
    
        count = count + 1

    print(vector_distancia)
    print(error_array)

win.geometry("1280x720")
win.title("Aplicación para localización de fallas - PF202130")


label1 = tkinter.Label(tab1,text = "Ingrese el valor de la frecuencia de muestreo:").place(x = 20, y = 65)  
label2 = tkinter.Label(tab1,text = "Ingrese el valor de la velocidad de propagación:").place(x = 20, y = 100) 
label3 = tkinter.Label(tab1,text = "Ingrese el valor de la impedancia característica:").place(x = 20, y = 140) 
label4 = tkinter.Label(tab1,text = "Seleccione método de cálculo:").place(x = 20, y = 180) 
label5 = tkinter.Label(tab1,text = "Herramienta para la localización de fallas utilizando teoría de ondas viajeras", font=('Helvetica', 18, 'bold')).place(x = 200, y = 10)  
labelalt = tkinter.Label(tab2,text = "Validación/Desempeño de la Herramienta Computacional", font=('Helvetica', 18, 'bold')).place(x = 320, y = 10)  

labelText = tkinter.StringVar()
labelText2 = tkinter.StringVar()

# label_d1 = tkinter.Label(tab1,text = "La distancia de la falla calculada es...", fg = "black", font =('Helvetica', 12)).place(x = 550, y = 65)
label_d2 = tkinter.Label(tab1,textvariable= labelText, fg = "green", font =('Helvetica', 30, 'bold')).place(x = 550, y = 100 )

# img = ImageTk.PhotoImage(Image.open("FaultLocApp.png"))
# panel = tkinter.Label(tab1, image = img).place(x = 500, y = 10) 
# label_e1 = tkinter.Label(win,text = "El error en (%) calculado es:", fg = "black", font =('Helvetica', 12)).place(x = 500, y = 70 )
# label_e2 = tkinter.Label(win,textvariable= labelText2, fg = "blue", font =('Helvetica', 12, 'bold')).place(x = 755, y = 70 )

entry = tkinter.Entry(tab1)
entry.place(x = 280, y = 65)
entry2 = tkinter.Entry(tab1)
entry2.place(x = 280, y = 100)
entry3 = tkinter.Entry(tab1)
entry3.place(x = 280, y = 140)

label1 = tkinter.Label(tab2,text = "Ingrese la distancia de la línea (en km):").place(x = 20, y = 65)  
entry4 = tkinter.Entry(tab2)
entry4.place(x = 280, y = 65)

label6 = tkinter.Label(tab2,text = "Seleccione método de cálculo:").place(x = 20, y = 105) 

def VarReading():
    global delta_t
    global velprop
    global Zc
    global line_distance, error 
    # global faultdistance
    
    # delta_t = float(entry.get())
    delta_t = 0.641E-6
    # velprop = float(entry2.get())
    velprop = 2.92051e05
    # Zc = float(entry3.get())
    Zc = 369.498

def VarReading2():
    global line_distance, error,delta_t,Zc,velprop

    line_distance = float(entry4.get())

    delta_t = 0.641E-6
    # velprop = float(entry2.get())
    velprop = 2.92051e05
    # Zc = float(entry3.get())
    Zc = 369.498

    print(delta_t)

    
    
button = tkinter.Button(tab1, text='Abrir archivo .csv', command=UploadAction)
button.place(x = 1120, y = 65)
button = tkinter.Button(tab2, text='Abrir archivos .csv', command=UploadAction2)
button.place(x = 1000, y = 65)
buttonrun = tkinter.Button(tab1, text='Ejecutar Script', command=lambda:[VarReading(), FaultLoc()])
buttonrun.place(x = 1130, y = 100)


widget_var = tkinter.StringVar()
combo = ttk.Combobox(tab1,values=["Señales de entrada", "Señales transformadas", "Ondas Viajeras", "Correlación"])
combo.set("Seleccione la gráfica...")
combo.place(x = 1100, y = 180)

widget_var2 = tkinter.StringVar()
combo2 = ttk.Combobox(tab1,values=["Correlación convencional", "MODWT"])
combo2.set("Seleccione método...")
combo2.place(x = 280, y = 180)

widget_var22 = tkinter.StringVar()
combo22 = ttk.Combobox(tab2,values=["Correlación convencional", "MODWT"])
combo22.set("Seleccione método...")
combo22.place(x = 280, y = 105)

fig1 = plt.figure(1,figsize=(6,4))
fig1.patch.set_facecolor('#F0F0F0')
canvas = FigureCanvasTkAgg(fig1, master=tab1)  # A tk.DrawingArea.
widg = canvas.get_tk_widget()
# widg.pack(padx=50, pady=170)
widg.place(x=80,y=215,height=440,width=1080)
    
fig2 = plt.figure(2,figsize=(6,4))
fig2.patch.set_facecolor('#F0F0F0')
canvas2 = FigureCanvasTkAgg(fig2, master=tab2)  # A tk.DrawingArea.
widg2 = canvas2.get_tk_widget()
# widg2.pack(padx=50, pady=170)
widg2.place(x=80,y=215,height=440,width=1080)

def Graph2():
    global vector_distancia, error_array, canvas2

    print(vector_distancia)
    print(error_array)

    fig2 = plt.figure(2,figsize=(6,4))
    fig2.patch.set_facecolor('#F0F0F0')
    canvas2 = FigureCanvasTkAgg(fig2, master=tab2)  # A tk.DrawingArea.
    widg2 = canvas2.get_tk_widget()
    # widg2.pack(padx=50, pady=170)
    widg2.place(x=80,y=160,height=440,width=1080)
    canvas2.draw()

    toolbar2 = NavigationToolbar2Tk(canvas2, tab2, pack_toolbar=False)
    toolbar2.update()

    plt.clf()
    plot3 = fig2.add_subplot(1,1,1)

    plot3.plot(vector_distancia, error_array, 'r')
    plot3.legend(labels=['Validación de fallas'],loc='lower left')
    plot3.title.set_text("Validación de fallas")
    plot3.set_xlabel('Distancia de la falla (%)')
    plot3.set_ylabel('Error de estimación(%)')

    
def clearCanv():
    for item in canvas.get_tk_widget().find_all():
       canvas.get_tk_widget().delete(item)

    for item in canvas2.get_tk_widget().find_all():
       canvas2.get_tk_widget().delete(item)

def Graph():
    value = combo.get()
    value2 = combo2.get()
   
    global canvas
    global Vd_k_vector, Vq_k_vector,V0_k_vector
    global Ad_k_vector, Aq_k_vector, A0_k_vector
    global corr, tau_1, tau_2,max_taus,peaks

    global time_vector
    global current_deltaA, current_deltaB,current_deltaC
    global voltage_deltaA, voltage_deltaB, voltage_deltaC
    global Vd_k_vector, Vq_k_vector,V0_k_vector
    global Ad_k_vector, Aq_k_vector, A0_k_vector
    global S_b, S_f
    global corr,peaks,tau_1,tau_2,max_taus
    
    

    if value == "Señales de entrada":
        print(value)
        # canvas.get_tk_widget().pack_forget()
       
        fig1 = plt.figure(1,figsize=(6,4))
        fig1.patch.set_facecolor('#F0F0F0')
        canvas = FigureCanvasTkAgg(fig1, master=tab1)  # A tk.DrawingArea.
        widg = canvas.get_tk_widget()
        # widg.pack(padx=50, pady=170)
        widg.place(x=80,y=215,height=440,width=1080)
        canvas.draw()

        toolbar = NavigationToolbar2Tk(canvas, tab1, pack_toolbar=False)
        toolbar.update()
        
        plt.clf()

        plot2 = fig1.add_subplot(1,2,1)
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


        plot2.plot(time_vector, voltage_deltaA, 'r')
        plot2.plot(time_vector, voltage_deltaB, 'b')
        plot2.plot(time_vector, voltage_deltaC,'g')
        # Agregar punto de falla en inputs
        plot2.axvline(x=0.016, color='k', linestyle='--')
        plot2.legend(labels=['Ph-A', 'Ph-B', 'Ph-C','Falla'],loc='lower left')
        plot2.title.set_text("Señales de tensión - Entrada")
        plot2.set_xlabel('Tiempo (s)')
        plot2.set_ylabel('Voltaje (V)')

        

    elif value == "Señales transformadas":
        print(value)
    
        # canvas.get_tk_widget().pack_forget()

        fig2 = plt.figure(2,figsize=(6,4))
        fig2.patch.set_facecolor('#F0F0F0')
        canvas = FigureCanvasTkAgg(fig2, master=tab1)  # A tk.DrawingArea.
        widg = canvas.get_tk_widget()
        # widg.pack(padx=50, pady=170)
        widg.place(x=80,y=215,height=440,width=1080)
        canvas.draw()

        toolbar = NavigationToolbar2Tk(canvas, tab1, pack_toolbar=False)
        toolbar.update()

        plt.clf()

        plot2 = fig2.add_subplot(1,2,1)
        plot1 = fig2.add_subplot(1,2,2)
        plot1.plot(time_vector,Ad_k_vector, 'r')
        plot1.plot(time_vector,Aq_k_vector, 'b')
        plot1.plot(time_vector,A0_k_vector,'g')
         # Agregar punto de falla en inputs
        plot1.axvline(x=0.016, color='k', linestyle='--')
        plot1.legend(labels=['Ph-A', 'Ph-B', 'Ph-C','Falla'],loc='lower left')
        plot1.title.set_text("Señales de corriente - Transformadas")
        plot1.set_xlabel('Tiempo (s)')
        plot1.set_ylabel('Corriente (A)')

        plot2.plot(time_vector,Vd_k_vector, 'r')
        plot2.plot(time_vector,Vq_k_vector, 'b')
        plot2.plot(time_vector,V0_k_vector,'g')
        # Agregar punto de falla en inputs
        plot2.axvline(x=0.016, color='k', linestyle='--')
        plot2.legend(labels=['Ph-A', 'Ph-B', 'Ph-C','Falla'],loc='lower left')
        plot2.title.set_text("Señales de tensión - Transformadas")
        plot2.set_xlabel('Tiempo (s)')
        plot2.set_ylabel('Voltaje (V)')

        
    elif value == "Ondas Viajeras":
        print(value)
        
        # canvas.get_tk_widget().pack_forget()

        fig3 = plt.figure(3,figsize=(6,4))
        fig3.patch.set_facecolor('#F0F0F0')
        canvas = FigureCanvasTkAgg(fig3, master=tab1)  # A tk.DrawingArea.
        widg = canvas.get_tk_widget()
        # widg.pack(padx=50, pady=170)
        widg.place(x=80,y=215,height=440,width=1080)
        canvas.draw()

        toolbar = NavigationToolbar2Tk(canvas, tab1, pack_toolbar=False)
        toolbar.update()

        plt.clf()

        if (value2 == 'Correlación convencional'):
            plot1 = fig3.add_subplot(111)
            plot1.plot(time_vector,S_f)
            plot1.plot(time_vector,S_b)
            plot1.legend(labels=['Onda de ida', 'Onda de llegada'],loc='lower left')
            plot1.title.set_text("Ondas viajeras - Correlación convencional")
            plot1.set_xlabel('Tiempo (s)')
            plot1.set_ylabel('Voltaje (V)')
        
        elif(value2 == 'MODWT'):
            plot1 = fig3.add_subplot(111)
            plot1.plot(time_vector,S_forward_modwt)
            plot1.plot(time_vector,S_backward_modwt)
            plot1.legend(labels=['Onda de ida', 'Onda de llegada'],loc='lower left')
            plot1.title.set_text("Ondas viajeras - MODWT")
            plot1.set_xlabel('Tiempo (s)')
            plot1.set_ylabel('Voltaje (V)')
        
        
    elif value == "Correlación":
        print(value)
        
        # canvas.get_tk_widget().pack_forget()

        fig4 = plt.figure(4,figsize=(6,4))
        fig4.patch.set_facecolor('#F0F0F0')
        canvas = FigureCanvasTkAgg(fig4, master=tab1)  # A tk.DrawingArea.
        widg = canvas.get_tk_widget()
        # widg.pack(padx=50, pady=170)
        widg.place(x=80,y=215,height=440,width=1080)
        canvas.draw()

        toolbar = NavigationToolbar2Tk(canvas, tab1, pack_toolbar=False)
        toolbar.update()

        plt.clf()

        if (value2 == 'Correlación convencional'):
            plot1 = fig4.add_subplot(111)
            plot1.plot(corr[0:round(len(corr))])
            plot1.scatter(peaks[tau_1],max_taus[0],c="red")
            plot1.scatter(peaks[tau_2],max_taus[1],c="black")
            # plot1.invert_xaxis()
            plot1.legend(labels=['Correlación cruzada - Convencional','Punto máximo #1', 'Punto máximo #2'],loc='lower right')
            plot1.title.set_text("Correlación cruzada - Convencional")
            plot1.set_xlabel('Número de muestras')
            plot1.set_ylabel('Puntos de correlación')

        elif(value2 == 'MODWT'):
            plot1 = fig4.add_subplot(111)
            plot1.plot(corr2[0:round(len(corr2))])
            plot1.scatter(peaks2[tau_11],max_taus2[0],c="red")
            plot1.scatter(peaks2[tau_22],max_taus2[1],c="black")
            # plot1.invert_xaxis()
            plot1.legend(labels=['Correlación cruzada - MODWT','Punto máximo #1', 'Punto máximo #2'],loc='lower right')
            plot1.title.set_text("Correlación cruzada - MODWT")
            plot1.set_xlabel('Número de muestras')
            plot1.set_ylabel('Puntos de correlación')



buttongraph = tkinter.Button(tab1, text='Realizar Gráfica', command=lambda:[clearCanv(),Graph()])
buttongraph.place(x = 1125, y = 140)

buttongraph2 = tkinter.Button(tab2, text='Realizar Validación', command=lambda:[VarReading2(),FaultLoc2(),Graph2()])
buttongraph2.place(x = 1000, y = 105)

win.resizable(0,0)
win.mainloop()