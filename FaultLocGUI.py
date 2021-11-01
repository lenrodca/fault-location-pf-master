from tkinter import *
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
root = Tk()
root.geometry('600x500+10+10')

def do_plot(x, y):
    [ax[x].clear() for x in range(4)]
    ax[1].plot(x,y)
    canvas.draw()

frame1 = Frame(root); frame1.place(x=0, y=0, width=400, height=400)
figure = plt.Figure(figsize=(5,5), facecolor='yellow')
canvas = FigureCanvasTkAgg(figure, frame1)
canvas.get_tk_widget().place(x=0,y=0,width=500,height=500)
ax = [figure.add_subplot(2, 2, x+1) for x in range(4)]

frame2 = Frame(root); frame2.place(x=500, y=0, width=100, height=400)
btplot1 = Button(frame2, text='plot 1', command= lambda: do_plot([0,1,2],[5,3,7]))
btplot1.place(x=0, y=50, width=50, height=20)
btplot2 = Button(frame2, text='plot 2', command= lambda: do_plot([5,6,7],[3,8,2]))
btplot2.place(x=0, y=100, width=50, height=20)

root.mainloop()
