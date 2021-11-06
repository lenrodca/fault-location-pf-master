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
from PIL import ImageTk, Image

import csv
global filename
filename = ""
global delta_t
global velprop
global Zc
global faultdistance
global combo
global canvas


win = tkinter.Tk()
tabControl = ttk.Notebook(win)

tab1 = ttk.Frame(tabControl)
tab2 = ttk.Frame(tabControl)

tabControl.add(tab1, text ='Localización de fallas')
tabControl.add(tab2, text ='Validación del método')
tabControl.pack(expand = 1, fill ="both")


def UploadAction():
    global filename
    filename = filedialog.askopenfilename()

def UploadAction2():
    global filename_arrays
    filename_arrays = filedialog.askopenfilenames()
    

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
    global S_b, S_f
    global corr,peaks,tau_1,tau_2,max_taus

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

    # error = (np.abs(faultdistance-d)/(faultdistance))*100
    # error = round(error,4)
    # labelText2.set(error)
    
    def FaultLoc2():
        global d,error_array,vector_distancia
        global line_distance, error 

        count = 0
        vector_distancia = np.zeros(100/len(filename_arrays))
        for i in range(len(vector_distancia)):
            vector_distancia[i] = 100/len(filename_arrays)*(i+1)

        for x in filename_arrays:
            data_cyc = pd.read_csv(x, delimiter=",") 
            
            print(data)

            data = data.replace(',','.', regex=True)

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

            for i in range(contador):

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

            # Ondas Viajeras
            S_forward_cyc = -Vd_k_vector_cyc + Zc*Ad_k_vector_cyc
            S_backward_cyc = Vd_k_vector_cyc + Zc*Ad_k_vector_cyc

            #Filtro
            tau = 1.56e06
            w = 20e03/tau
            b1, a1 = signal.butter(3, w, 'high')
            S_f_cyc = signal.filtfilt(b1, a1, S_forward_cyc,axis=0)

            b2, a2 = signal.butter(3, w, 'high')
            S_b_cyc = signal.filtfilt(b2, a2, S_backward_cyc,axis=0)

            #Correlacion

            corr_cyc = np.correlate(S_b_cyc,S_f_cyc,mode="same")/contador
            peaks_cyc, _ = find_peaks(corr_cyc)

            taus_cyc = np.zeros(len(peaks_cyc))
            for i in range(0,len(peaks_cyc)):
                taus_cyc[i] = corr_cyc[peaks_cyc[i]]


            max_taus_cyc = heapq.nlargest(2, taus_cyc)
        
            tau_1_cyc = np.argwhere(taus_cyc == max_taus_cyc[0]).astype(int)
            tau_2_cyc = np.argwhere(taus_cyc== max_taus_cyc[1]).astype(int)
            # tau_f = np.argmax(corr)

            d = ((velprop*(np.abs(time_vector_cyc[peaks_cyc[tau_1_cyc]]-time_vector_cyc[peaks_cyc[tau_2_cyc]])))/(2))
            d = round(d[0,0],4)
            print(d)

            error = ((((vector_distancia[count])/100)*line_distance - d)/((vector_distancia[count])/100)*line_distance)*100

            error_array = np.zeros(len(filename_arrays))
            
            for i in range (len(error_array)):
                error_array[i] = error

            count = count + 1

            

            



            


win.geometry("1280x720")
win.title("Aplicación para localización de fallas - PF202130")


label1 = tkinter.Label(tab1,text = "Ingrese el valor de la frecuencia de muestreo:").place(x = 20, y = 50)  
label2 = tkinter.Label(tab1,text = "Ingrese el valor de la velocidad de propagación:").place(x = 20, y = 80) 
label3 = tkinter.Label(tab1,text = "Ingrese el valor de la impedancia característica:").place(x = 20, y = 110) 
label5 = tkinter.Label(tab1,text = "Programa para la localización de fallas utilizando ondas viajeras", font=('Helvetica', 18, 'bold')).place(x = 360, y = 10)  

labelText = tkinter.StringVar()
labelText2 = tkinter.StringVar()

label_d1 = tkinter.Label(tab1,text = "La distancia de la falla calculada es:", fg = "black", font =('Helvetica', 12)).place(x = 500, y = 50 )
label_d2 = tkinter.Label(tab1,textvariable= labelText, fg = "green", font =('Helvetica', 12, 'bold')).place(x = 755, y = 50 )

# img = ImageTk.PhotoImage(Image.open("FaultLocApp.png"))
# panel = tkinter.Label(tab1, image = img).place(x = 500, y = 10) 
# label_e1 = tkinter.Label(win,text = "El error en (%) calculado es:", fg = "black", font =('Helvetica', 12)).place(x = 500, y = 70 )
# label_e2 = tkinter.Label(win,textvariable= labelText2, fg = "blue", font =('Helvetica', 12, 'bold')).place(x = 755, y = 70 )

entry = tkinter.Entry(tab1)
entry.place(x = 280, y = 50)
entry2 = tkinter.Entry(tab1)
entry2.place(x = 280, y = 80)
entry3 = tkinter.Entry(tab1)
entry3.place(x = 280, y = 110)

label1 = tkinter.Label(tab2,text = "Ingrese el valor de la distancia de la linea:").place(x = 20, y = 50)  
entry4 = tkinter.Entry(tab2)
entry4.place(x = 280, y = 50)


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
    global line_distance, error 

    line_distance = float(entry4.get())

    
    
button = tkinter.Button(tab1, text='Abrir archivo .csv', command=UploadAction)
button.place(x = 1120, y = 50)
button = tkinter.Button(tab2, text='Abrir grupo de archivos .csv', command=UploadAction2)
button.place(x = 1120, y = 50)
buttonrun = tkinter.Button(tab1, text='Ejecutar Script', command=lambda:[VarReading(), FaultLoc()])
buttonrun.place(x = 1130, y = 80)


widget_var = tkinter.StringVar()
combo = ttk.Combobox(tab1,values=["Señales de entrada", "Señales transformadas", "Ondas Viajeras", "Correlación"])
combo.set("Seleccione la gráfica...")
combo.place(x = 1100, y = 140)




fig1 = plt.figure(1,figsize=(6,4))
canvas = FigureCanvasTkAgg(fig1, master=tab1)  # A tk.DrawingArea.
widg = canvas.get_tk_widget()
widg.pack(padx=50, pady=170,side=tkinter.TOP, fill=tkinter.BOTH, expand=1)
    
fig2 = plt.figure(2,figsize=(6,4))
canvas2 = FigureCanvasTkAgg(fig2, master=tab2)  # A tk.DrawingArea.
widg2 = canvas2.get_tk_widget()
widg2.pack(padx=50, pady=170,side=tkinter.TOP, fill=tkinter.BOTH, expand=1)

# fig3 = plt.figure(3,figsize=(6,4))
# canvas3 = FigureCanvasTkAgg(fig1, master=tab1)  # A tk.DrawingArea.
# widg3 = canvas3.get_tk_widget()
# widg3.pack(padx=50, pady=170,side=tkinter.TOP, fill=tkinter.BOTH, expand=1)

def Graph2():
        global vector_distancia, error_array
        canvas2.get_tk_widget().pack_forget()
    
        fig2 = plt.figure(1,figsize=(6,4))
        canvas = FigureCanvasTkAgg(fig1, master=tab1)  # A tk.DrawingArea.
        widg = canvas.get_tk_widget()
        widg.pack(padx=50, pady=170,side=tkinter.TOP, fill=tkinter.BOTH, expand=1)
        canvas.draw()

        plot1 = fig2.add_subplot(1,1,1)
    
        plot1.plot(vector_distancia, error_array, 'r')
        plot1.legend(labels=['Validación de fallas'],loc='lower left')
        plot1.title.set_text("Validación de fallas")
        plot1.set_xlabel('Distancia de la falla (%)')
        plot1.set_ylabel('Error de estimación(%)')

    

def Graph():
    value = combo.get()
   
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
        canvas.get_tk_widget().pack_forget()
       
        fig1 = plt.figure(1,figsize=(6,4))
        canvas = FigureCanvasTkAgg(fig1, master=tab1)  # A tk.DrawingArea.
        widg = canvas.get_tk_widget()
        widg.pack(padx=50, pady=170,side=tkinter.TOP, fill=tkinter.BOTH, expand=1)
        canvas.draw()

        
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
    
        canvas.get_tk_widget().pack_forget()

        fig2 = plt.figure(2,figsize=(6,4))
        canvas = FigureCanvasTkAgg(fig2, master=tab1)  # A tk.DrawingArea.
        widg = canvas.get_tk_widget()
        widg.pack(padx=50, pady=170,side=tkinter.TOP, fill=tkinter.BOTH, expand=1)
        canvas.draw()

        
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
        
        canvas.get_tk_widget().pack_forget()

        fig3 = plt.figure(3,figsize=(6,4))
        canvas = FigureCanvasTkAgg(fig3, master=tab1)  # A tk.DrawingArea.
        widg = canvas.get_tk_widget()
        widg.pack(padx=50, pady=170,side=tkinter.TOP, fill=tkinter.BOTH, expand=1)
        canvas.draw()

        plot1 = fig3.add_subplot(111)
        plot1.plot(time_vector,S_f)
        plot1.plot(time_vector,S_b)
        plot1.legend(labels=['Onda de ida', 'Onda de llegada'],loc='lower left')
        plot1.title.set_text("Ondas viajeras")
        plot1.set_xlabel('Tiempo (s)')
        plot1.set_ylabel('Voltaje (V)')
        
        
    elif value == "Correlación":
        print(value)
        
        canvas.get_tk_widget().pack_forget()

        fig4 = plt.figure(4,figsize=(6,4))
        canvas = FigureCanvasTkAgg(fig4, master=tab1)  # A tk.DrawingArea.
        widg = canvas.get_tk_widget()
        widg.pack(padx=50, pady=170,side=tkinter.TOP, fill=tkinter.BOTH, expand=1)
        canvas.draw()

        plot1 = fig4.add_subplot(111)
        plot1.plot(time_vector[0:round(len(time_vector)/2)],corr[0:round(len(corr)/2)])
        plot1.scatter(time_vector[peaks[tau_1]],max_taus[0])
        plot1.scatter(time_vector[peaks[tau_2]],max_taus[1])
        plot1.invert_xaxis()
        plot1.legend(labels=['Correlación cruzada','Punto máximo #1', 'Punto máximo #2'],loc='lower right')
        plot1.title.set_text("Correlación cruzada")
        plot1.set_xlabel('Tiempo (s)')
        plot1.set_ylabel('Puntos de correlación')

buttongraph = tkinter.Button(tab1, text='Realizar Gráfica', command=Graph)
buttongraph.place(x = 1125, y = 110)

buttongraph2 = tkinter.Button(tab2, text='Realizar Gráfica', command=Graph2)
buttongraph2.place(x = 1125, y = 110)


win.mainloop()




