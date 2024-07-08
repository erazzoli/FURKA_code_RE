# First things, first. Import the wxPython package.
import wx
from bsread import source
from collections import deque
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.animation import FuncAnimation
from threading import Thread
from time import sleep
from scipy.optimize import curve_fit
import random





from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigureCanvas

class MyApp(wx.App):
    def OnInit(self):
        frame = MyFrame(None, -1, "wxPython Matplotlib Example")
        frame.Show(True)
        self.SetTopWindow(frame)
        return True

class MyFrame(wx.Frame):
    def __init__(self, parent, ID, title):
        wx.Frame.__init__(self, parent, ID, title, size=(600, 400))
        self.figure = plt.figure(figsize=(10, 5))
        self.canvas = FigureCanvas(self, -1, self.figure)
        self.draw_figure()

    def draw_figure(self):
        ax = self.figure.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
        ax.stock_img()

app = MyApp(0)
app.MainLoop()

