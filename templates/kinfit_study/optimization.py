# take an txt files with results of efficiencies etc with different value of the cut

import os
import glob
import math
import ROOT
from ROOT import *
import sys
from optparse import OptionParser
from array import array	

ROOT.gROOT.Macro( os.path.expanduser( '~/rootlogon.C' ) )
# ROOT.gROOT.SetBatch(True)

parser = OptionParser()

parser.add_option('--file', metavar='F', type='string', action='store',
                  default='',
                  dest='file',
                  help='')
parser.add_option('--verbose', metavar='F', type='string', action='store',
                  default = 'no',
                  dest='verbose',
                  help='')

(options, args) = parser.parse_args()


def makeTGraph(x,y,x_title='x',y_title='y',title='TGraph',ROC=0):
	if ROC :
		pass
		# x.extend([0.,1.])
		# y.extend([0.,1.])
	array_x = array('f')
	array_y = array('f')
	array_x.fromlist(x)
	array_y.fromlist(y)
	graph = ROOT.TGraph(len(x),array_x,array_y)
	graph.SetTitle(title)
	graph.GetXaxis().SetTitle(x_title)
	graph.GetYaxis().SetTitle(y_title)
	return graph

def makeRefLine(num_points = 50,x_title='x',y_title='y',title='TGraph'):
	step = 1.0/num_points
	i = 0
	x,y = [],[]
	while i < 1 :
		x.append(i)
		y.append(i)
		i+=step
	return(makeTGraph(x,y,x_title='x',y_title='y',title='TGraph'))

def getColor(itype):
	color_types = ['signal','tt_bkg','singletop','wjets','bkg']
	color_colors = [2,4,6,3,ROOT.kRed+2]
	colors = [{'type':color_types[i],'color':color_colors[i]} for i in range(len(color_colors))]
	ientry = [ientry for ientry in colors if itype==ientry['type']][0]
	return ientry['color']

# x_y = [{'x_var':	,'x_title':	,'y_var':	,'y_title':	} 
def makePlottingPairs(x_var,y_var,ROC=0):
	return {'x_var':x_var,'x_title':x_var,'y_var':y_var,'y_title':y_var,'ROC':ROC}

############ Main program ############
f1 = open(options.file)

results = {}
results['cut_name'] = []
results['cut_value'] = []
results['zjets_efficiency'] = []
results['signal_efficiency'] = []
results['tt_bkg_efficiency'] = []
results['singletop_efficiency'] = [] 
results['wjets_efficiency'] = []
results['s_b'] = []
results['bkg_efficiency'] = []
results['total_fraction_left'] = []
results['s_b nominal'] = []

for line in f1:
	if '<' in line or '>' in line :
		line = line.split()
		if results['cut_name'] == []:
			results['cut_name'].append(line[1])
		results['cut_value'].append(float(line[-1]))

	if 'zjets' in line :
		line = line.split()
		results['zjets_efficiency'].append(float(line[-1]))

	if 'tt_bkg' in line :
		line = line.split()
		results['tt_bkg_efficiency'].append(float(line[-1]))

	if 'singletop' in line :
		line = line.split()
		results['singletop_efficiency'].append(float(line[-1]))

	if 'WJets' in line or 'wjets' in line :
		line = line.split()
		results['wjets_efficiency'].append(float(line[-1]))


	# only keep one number of s_b before cut
	if 'sig/bkg before' in line and results['s_b nominal'] == [] :
		line = line.split()
		results['s_b nominal'].append(float(line[-1]))

	if 'sig/bkg after' in line :
		line = line.split()
		results['s_b'].append(float(line[-1]))

	if 'sig efficiency' in line :
		line = line.split()
		results['signal_efficiency'].append(float(line[-1]))

	if 'bkg efficiency' in line :
		line = line.split()
		results['bkg_efficiency'].append(float(line[-1]))


	if 'fraction of total events left' in line :
		line = line.split()
		results['total_fraction_left'].append(float(line[-1]))

print 'Finish reading lines from txt file!'

# Making plots
vars_to_plot = ['signal_efficiency','bkg_efficiency','total_fraction_left','s_b','wjets_efficiency','singletop_efficiency']
cut_name = results['cut_name'][0]

x_y = [{'x_var':'cut_value','x_title':cut_name,'y_var':iy,'y_title':iy,'ROC':0} for iy in vars_to_plot]

x_y.append(makePlottingPairs('signal_efficiency','bkg_efficiency',1))
x_y.append(makePlottingPairs('signal_efficiency','wjets_efficiency',1))
x_y.append(makePlottingPairs('signal_efficiency','singletop_efficiency',1))
x_y.append(makePlottingPairs('signal_efficiency','tt_bkg_efficiency',1))

# Make plot dir
plotdir = 'plots/'
if not os.path.exists(plotdir):
    os.mkdir(plotdir)
    print 'Creating new dir '+plotdir

roc_plots = []

for item in x_y:
	x_var_name = item['x_var']
	x_var = results[x_var_name]
	x_title = item['x_title']
	y_var_name = item['y_var']
	y_var = results[y_var_name]
	y_title = item['y_title']
	title = y_title+' vs '+x_title
	is_ROC = item['ROC']

	graph = makeTGraph(x_var,y_var,x_title,y_title,title)
	graph.SetMarkerStyle(7)

	# Draw a straight line of slope one if it is ROC
	if is_ROC:
		roc_plots.append({'plot':graph,'title':y_title.split('_efficiency')[0]})
	else :
		c = ROOT.TCanvas()
		graph.Draw("Apl")
		c.SaveAs(plotdir+x_title+'_'+y_title+'.png')

# plot ROCs
canvas = ROOT.TCanvas()
# set legend
leg = ROOT.TLegend(0.2408027,0.6,0.5117057,0.8)
leg.SetFillColor(0)
leg.SetLineColor(0) 
leg.SetTextSize(0.03448276)
# create a multigraph
mg = ROOT.TMultiGraph()
# loop over all ROC
for i,item in enumerate(roc_plots):
	iroc = item['plot']
	ititle = item['title']
	icolor = getColor(ititle)
	iroc.SetMarkerColor(icolor)
	iroc.SetLineColor(icolor)
	iroc.SetTitle(ititle)
	iroc.SetFillStyle(0)
	mg.Add(iroc)
	l1 = leg.AddEntry(iroc,ititle,'p')
	l1.SetTextColor(icolor)
# config multigraph
mg.Draw('apl')
mg.SetTitle('ROC of %s'%cut_name)
mg.GetXaxis().SetTitle('signal_efficiency')
mg.GetYaxis().SetTitle('bkg_efficiency')
# Add a legend
# canvas.BuildLegend();
leg.Draw('same')
# Draw a reference line
ref_plot = ROOT.TF1('ref','x',0,1)
ref_plot.Draw('same')
# Plot and save ROC
canvas.SaveAs(plotdir+'ROC.png')


