import ROOT
from array import array
import sys
ROOT.gROOT.SetBatch(True)

def makeTGraphErrors(x,y,y_err,x_err=None,x_title='x',y_title='y',title='TGraph'):
        array_x = array('f')
        array_y = array('f')
        array_xerr = array('f')
        array_yerr = array('f')
        array_x.fromlist(x)
        array_y.fromlist(y)
        lenx = len(x)
        if x_err is None:
            array_xerr = array('f',lenx*[0.0])
        else:
            array_xerr.fromlist(x_err)
        array_yerr.fromlist(y_err)

        graph = ROOT.TGraphErrors(lenx,array_x,array_y,array_xerr,array_yerr)
        graph.SetTitle(title)
        graph.GetXaxis().SetTitle(x_title)
        graph.GetYaxis().SetTitle(y_title)
        return graph
