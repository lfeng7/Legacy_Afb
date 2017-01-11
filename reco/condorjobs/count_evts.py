"""
count total number of evts in given root files
"""
import ROOT

import sys

argv = sys.argv[1:]

counts = open('counts.txt','w')

towrite = ''
total = 0
for ifile in argv:
    tf = ROOT.TFile(ifile)
    ttree = tf.Get('selected')
    nevts = ttree.GetEntriesFast()
    towrite += '%s has %i entries \n'%(tf.GetName(),nevts)
    total += nevts
    tf.Close()
counts.write('Total num entries = %i\n'%total)
counts.write(towrite)
counts.close()

