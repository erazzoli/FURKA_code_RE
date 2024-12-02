#!/usr/bin/env python

import os
import wx

from slic.gui.persist import Persistence

from bstrd import bsstream

from spc_wx import MainPanel


def run(*args, **kwargs):
    app = wx.App()
    wx.Yield = app.Yield
    frame = MainFrame(*args, **kwargs)
    frame.Show()
    app.MainLoop()
    app.Yield() #TODO: without this, wxPython segfaults locking a mutex


def get_wx_icon(fname="FurkaLogo.png"):
    iname = os.path.dirname(__file__)
    iname = os.path.join(iname, fname)
    return wx.Icon(iname)



class MainFrame(wx.Frame):

    def __init__(self, title="Cmos RIXS"):
        super().__init__(None, title=title)
        self.SetIcon(get_wx_icon())

        main = MainPanel(self)

        # make sure the window is large enough
        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(main, proportion=1, flag=wx.EXPAND)
        self.SetSizerAndFit(sizer)

        self.persist = persist = Persistence(".coffee365", self)
        persist.load()

        self.Bind(wx.EVT_CLOSE, self.on_close)


    def on_close(self, event):
        print("bye")
        bsstream.stop()
        self.persist.save()
        event.Skip() # forward the close event
        wx.GetApp().Yield()



if __name__ == "__main__":
    run()




