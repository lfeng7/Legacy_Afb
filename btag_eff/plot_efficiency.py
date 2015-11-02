##############################################################################################################
##                                      Imports and Other Header Stuff                                      ##
##############################################################################################################
from utility import *
from ROOT import *
from optparse import OptionParser

parser = OptionParser()

parser.add_option('--inputfiles', metavar='F', type='string', action='store',
                  default = "",
                  dest='inputFiles',
                  help='Input files')

parser.add_option('--upload', metavar='F', type='string', action='store',
                  default = "no",
                  dest='upload',
                  help='If upload to webpage')

(options, args) = parser.parse_args()

argv = []

# Some prepends
plotsdir = 'plots/ejets/'

eta = [0.8,1.1,1.6,2.3]
pt = 75
n =40
pt_max = 1000
step = pt_max/n


def main() :
    # Get the file list with all input files.
    if options.inputFiles != '':
        allfiles = glob.glob( options.inputFiles )
        for ifile in allfiles: 
            tfile = ROOT.TFile(ifile)
            plot(tfile)
            tfile.Close()
            # Upload to webpage
        if options.upload == 'yes':
            Upload(plotsdir,'btagging_efficiency')
    else:
        print 'No correct input files given. Will do nothing.'


def get_btag_eff (tfile,pt,eta,jet_flavor):
    efficiency_b = tfile.Get('efficiency_b')
    efficiency_c = tfile.Get('efficiency_c')
    efficiency_udsg = tfile.Get('efficiency_udsg') 
    # x,y of TH2F of efficiency are pt and eta
    if jet_flavor == 5 :
        binx = efficiency_b.GetXaxis().FindBin(pt)
        biny = efficiency_b.GetYaxis().FindBin(eta)
        bin = efficiency_b.GetBin(binx,biny)
        return efficiency_b.GetBinContent(bin)
    elif jet_flavor == 4 :
        binx = efficiency_c.GetXaxis().FindBin(pt)
        biny = efficiency_c.GetYaxis().FindBin(eta)
        bin = efficiency_c.GetBin(binx,biny)
        return efficiency_c.GetBinContent(bin)  
    else :
        binx = efficiency_udsg.GetXaxis().FindBin(pt)
        biny = efficiency_udsg.GetYaxis().FindBin(eta)
        bin = efficiency_udsg.GetBin(binx,biny)
        return efficiency_udsg.GetBinContent(bin)



def plot(tfile) :
    # Setup input and output files
    print '\nProcessing',tfile.GetName()
    inputname = tfile.GetName().split('/')
    inputname = inputname[len(inputname)-1]
    inputname = inputname.split('.root')[0].split('_CSVM_bTaggingEfficiencyMap')[0]
    outputdir = plotsdir+inputname+'/'
    MakeDirectory(outputdir)    
    channel = inputname
    # Actualy program starts here
    title =[]
    file_names = []
    for i in range(len(eta)) :
        title.append(channel+' btagging efficiency eta = '+str(eta[i])) # +' flavor = '+str(jet_flavor)
        file_names += [outputdir+'btag_eff_eta_'+str(i)+'.png']
    x_title = 'pt (GeV)'
    y_title = 'btag_eff'
        
    x = array('d')
    flavors = [5,4,3]
    flavor_names = ['b','c','udsg']
    color_list = [2,4,7]

    # Make 2D plots of efficiency
    efficiency_b = tfile.Get('efficiency_b')
    efficiency_c = tfile.Get('efficiency_c')
    efficiency_udsg = tfile.Get('efficiency_udsg') 
    list_effs = [efficiency_b,efficiency_c,efficiency_udsg]
    for i in range(len(list_effs)):
        tmp_h = list_effs[i]
        tmp_c = ROOT.TCanvas('tmp_c_'+str(i))
        tmp_name = channel+' b-tagging '+tmp_h.GetName()     
        tmp_file = outputdir+'btagg_eff_'+tmp_h.GetName()+'.png'   
        tmp_h.SetTitle(tmp_name)
        tmp_h.SetStats(0)
        tmp_h.GetXaxis().SetTitle('pT (GeV)')
        tmp_h.GetYaxis().SetTitle('abs(eta)')
        tmp_h.Draw('colz')
        tmp_c.SaveAs(tmp_file)

    for i in range(n):
        pt_tmp = i*step
        x.append(pt_tmp)

    for eta_bin in range(len(eta)) :
        y = []
        for i in range(len(flavors)) :
            y_tmp = array('d')
            for j in range(n) :
                pt_tmp = j*step
                btag_eff = get_btag_eff(tfile,pt_tmp,eta[eta_bin],flavors[i])
                y_tmp.append(btag_eff)
            y.append(y_tmp)
        # print x
        # print y
        c1 = TCanvas("c_eta_"+str(eta_bin),"A Simple Graph ",200,10,700,500)
        c1.SetGrid();
        mg = TMultiGraph();
        leg = TLegend(0.7,0.7,0.9,0.9)
        for i in range(len(flavors)):
            # i=2
            gr = TGraph(n,x,y[i]);
            gr.SetLineColor(color_list[i]);
            gr.SetLineWidth(4);
            gr.SetMarkerColor(1);
            gr.SetMarkerStyle(21);
            #debug
            # print title[eta_bin]
            gr.SetTitle(title[eta_bin]);
            gr.GetXaxis().SetTitle(x_title);
            gr.GetYaxis().SetTitle(y_title);
            mg.Add(gr)
            leg.AddEntry(gr,flavor_names[i],'LP')

        mg.SetTitle(title[eta_bin]);    
        mg.Draw("ACP");
        leg.Draw()

        c1.SaveAs(file_names[eta_bin])

main()
