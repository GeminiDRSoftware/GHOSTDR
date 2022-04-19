#!/usr/bin/env python3

import wx
import numpy as np
 
class Mywin(wx.Frame): 
   def __init__(self, parent, title, params): 
      super(Mywin, self).__init__(parent, title = title)
      
      panel = wx.Panel(self) 
      vbox = wx.BoxSizer(wx.VERTICAL) 

      for aspect in params.keys():
          nm = wx.StaticBox(panel, -1, aspect) 
          nmSizer = wx.StaticBoxSizer(nm, wx.VERTICAL)

          parshape=params[aspect].shape
          
          for row in parshape[0]:
              nmbox = wx.BoxSizer(wx.HORIZONTAL) 
              fn = wx.StaticText(panel, -1, "Parameter") 

              nmbox.Add(fn, 0, wx.ALL|wx.CENTER, 5) 
              nm1 = wx.TextCtrl(panel, -1, style = wx.ALIGN_LEFT) 
              nm2 = wx.TextCtrl(panel, -1, style = wx.ALIGN_LEFT) 
              ln = wx.StaticText(panel, -1, "Last Name") 

              nmbox.Add(nm1, 0, wx.ALL|wx.CENTER, 5)
              nmbox.Add(ln, 0, wx.ALL|wx.CENTER, 5) 
              nmbox.Add(nm2, 0, wx.ALL|wx.CENTER, 5) 
              nmSizer.Add(nmbox, 0, wx.ALL|wx.CENTER, 10)  





          sbox = wx.StaticBox(panel, -1, 'buttons:') 
          sboxSizer = wx.StaticBoxSizer(sbox, wx.VERTICAL) 

          hbox = wx.BoxSizer(wx.HORIZONTAL) 
          okButton = wx.Button(panel, -1, 'ok') 

          hbox.Add(okButton, 0, wx.ALL|wx.LEFT, 10) 
          cancelButton = wx.Button(panel, -1, 'cancel') 

          hbox.Add(cancelButton, 0, wx.ALL|wx.LEFT, 10) 
          sboxSizer.Add(hbox, 0, wx.ALL|wx.LEFT, 10) 
          vbox.Add(nmSizer,0, wx.ALL|wx.CENTER, 5) 
          vbox.Add(sboxSizer,0, wx.ALL|wx.CENTER, 5) 
          panel.SetSizer(vbox) 
          self.Centre() 

          panel.Fit() 
          self.Show()

params={'xmod':np.zeros((3,3)),'specmod':np.zeros(2)}
app = wx.App() 
Mywin(None,  'staticboxsizer demo',params) 
app.MainLoop()
