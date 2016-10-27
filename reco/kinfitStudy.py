# Take in ttree of reco files, keep some branches for kinfit result study

# from Legacy_Afb.Tools.ttbar_utility import *
# from Legacy_Afb.Tools.angles_tools import *
import glob
from optparse import OptionParser
import ROOT
# from angles_tools import *

# Job steering

# Input inputFiles to use. This is in "glob" format, so you can use wildcards.
# If you get a "cannot find file" type of error, be sure to use "\*" instead
# of "*" to make sure you don't confuse the shell. 

parser = OptionParser()

parser.add_option('--inputfiles', metavar='F', type='string', action='store',
                  default = "",
                  dest='inputFiles',
                  help='Input files')

parser.add_option('--evtsperjob', metavar='F', type='int', action='store',
                  default = 1000,
                  dest='evtsperjob',
                  help='number of events to run for each job')

parser.add_option('--evtstart', metavar='F', type='int', action='store',
                  default = 0,
                  dest='evtstart',
                  help='the evt to start')

parser.add_option('--verbose', metavar='F', type='string', action='store',
                  default = 'no',
                  dest='verbose',
                  help='If you want more information than you usually need.')

parser.add_option('--slim', metavar='F', type='string', action='store',
                  default = 'yes',
                  dest='slim',
                  help='If you want slimmed ttree with angles and stuff')

parser.add_option('--fakelep', metavar='F', type='string', action='store',
                  default = 'no',
                  dest='fakelep',
                  help='If run on selected events with fake lepton.')

parser.add_option('--lepisocut', metavar='F', type='float', action='store',
                  default = 0.15,
                  dest='lepisocut',
                  help='Lower bound for fake electron isolation.')

parser.add_option('--nbcut', metavar='F', type='int', action='store',
                  default = 2,
                  dest='nbcut',
                  help='Number of b-tagged jets cut')

parser.add_option('--nlepcut', metavar='F', type='int', action='store',
                  default = 1,
                  dest='nlepcut',
                  help='Number of selected leptons cut')

parser.add_option('--applytrigger', metavar='F', type='string', action='store',
                  default = 'yes',
                  dest='applytrigger',
                  help='If apply trigger on MC')

(options, args) = parser.parse_args()

argv = []

# Some preset constants
data_lumi = 19748 
csvm = 0.679 

def main():
    # Get the file list with all input files.
    if options.inputFiles != '':
        allfiles = glob.glob( options.inputFiles )
    else:
        allfiles = []

    timer = ROOT.TStopwatch()
    timer.Start()

    # Job splitting is done here
    for ifile in allfiles:
        sample_name = ifile.split('/')
        sample_name = sample_name[len(sample_name)-1].split('.root')[0]
        evt_start = options.evtstart
        evt_to_run = options.evtsperjob
        isFakeLep = options.fakelep
        tfile = ROOT.TFile(ifile)
        # Do make angles
        print 'Do stuff now!'
        print 'run options:',tfile.GetName(),sample_name,evt_start,evt_to_run,isFakeLep
        reco_message = makeAngles(tfile,sample_name,evt_start,evt_to_run,isFakeLep)   
        print reco_message

    print 'All done!'

    # Stop our timer
    timer.Stop()
    # Print out our timing information
    print '\n'
    rtime = timer.RealTime(); # Real time (or "wall time")
    ctime = timer.CpuTime(); # CPU time
    print("RealTime={0:6.2f} seconds, CpuTime={1:6.2f} seconds").format(rtime,ctime)



# Do reconstruction for the given file, process only a specific part of it, and output to a new file
# ex. tf1, 'Tbar_s','singletopbar' (for b-tagging correction purpose),2000,10000
def makeAngles(tfile,sample_name,evt_start=0,evt_to_run=1000,isFakeLep='no'):
    # Get ttree
    tmptree = tfile.Get('selected')
    print 'Start to process file',tfile.GetName()
    if evt_to_run>0 :
        evt_end = evt_start+evt_to_run
    else : # (evt_to_run = -1 indicates run all events)
        evt_end = -1 
    print 'Prcess from event',evt_start,'to event',evt_end
    # Check if the starting event is more than total number of evts in current files
    if evt_start > tmptree.GetEntries():
        return 'evt_start > total_num_events in input ttree. Will stop here.'
        
    # Make output file
    fout_name = sample_name+'_kinfit_'+str(evt_start)+'.root'
    fout = ROOT.TFile(fout_name,'recreate')
    # Make output ttree
    br_tokeep = ['gen_type','gen_is_ejets','N_btag','N_jets','final_chi2','kinfit_results','kinfit_results','N_combos','N_combo_errs']
    br_tokeep += ['final_errflags','reco_*']
    if options.slim == 'yes':
        tmptree.SetBranchStatus('*',0)
        for ibr in br_tokeep:
            tmptree.SetBranchStatus(ibr,1)
    newtree = tmptree.CloneTree(0)
    tmptree.SetBranchStatus('*',1)

    # Add new branches to the output tree
    vecs = []
    br_names = []
    # Angles and other template vars
    isMatched = ROOT.vector('unsigned int')()
    delR = ROOT.vector('float')()

    vecs += [isMatched,delR]
    br_names += ['isMatched','delR']
    # Add branches to the tree
    branches = zip(br_names,vecs)
    for ibr in branches:
        newtree.Branch(ibr[0],ibr[1])

    # Add cutflow diagrams
    h_cutflow = ROOT.TH1D('cutflow_extra',' cutflow extra;cuts;events',5,0.,5.)
    # h_cutflow.SetBit(ROOT.TH1.kCanRebin)
    h_cutflow.SetDirectory(fout)

    # Find out if the sample is TTbar MC
    ttbar_names = ['TT','signal']
    is_ttbar = 0
    for item in ttbar_names:
        if item in sample_name:
            is_ttbar = 1
    if is_ttbar == 1:
        print 'Is ttbar MC sample!'

    # Loop over entries
    n_evt = 0
    for iev in range(evt_start,tmptree.GetEntries()):
        # Reset all vector containers
        for ivec in vecs: ivec.clear()    

        # Progress report
        if iev%5000 == 0 : print 'processing event',iev 
        # Break at the given event number
        if iev == evt_end : print 'Finish processing. Quit at event',evt_end; break ;
        
        tmptree.GetEntry(iev)

        h_cutflow.Fill('no cut',1)
        # Apply trigger....
        if not tmptree.trigger[0] : continue
        h_cutflow.Fill('trigger',1)
        # Skip events that kinfit did not converge ( error flag == 4 )
        if not tmptree.final_errflags[0]==0 : continue
        h_cutflow.Fill('kinfit error',1)

        # Book the branches from input ttree
        # reco_p4 is a list of tlep_p4,thad_p4,wlep_p4,whad_p4
        reco_pt = tmptree.reco_pt
        reco_eta = tmptree.reco_eta
        reco_phi = tmptree.reco_phi
        reco_mass = tmptree.reco_mass      
        lep_charge = tmptree.lep_charge[0]
        # Make 4vec for reco t and W's
        reco_p4 = []
        for i in range(4):
            ip4 = ROOT.TLorentzVector()
            ip4.SetPtEtaPhiM(reco_pt[i],reco_eta[i],reco_phi[i],reco_mass[i])
            reco_p4.append(ip4)
        # Make 4vecs of gen tlep_p4,thad_p4,wlep_p4,whad_p4
        gen_pt = tmptree.gen_pt
        gen_eta = tmptree.gen_eta
        gen_phi = tmptree.gen_phi
        gen_mass = tmptree.gen_phi
        gen_side = tmptree.gen_side
        gen_pdgid = tmptree.gen_pdgid
        gen_index = [0,0,0,0]
        if tmptree.gen_type[0] == 'e_jets':
            for i in range(len(gen_side)):
                if gen_side[i] == 'lep' and abs(gen_pdgid[i])==6 : gen_index[0]=i 
                if gen_side[i] == 'had' and abs(gen_pdgid[i])==6 : gen_index[1]=i 
                if gen_side[i] == 'lep' and abs(gen_pdgid[i])==24 : gen_index[2]=i 
                if gen_side[i] == 'had' and abs(gen_pdgid[i])==24 : gen_index[3]=i 
            gen_p4 = []
            for j in range(4):
                ip4 = ROOT.TLorentzVector()
                i = gen_index[j]
                ip4.SetPtEtaPhiM(gen_pt[i],gen_eta[i],gen_phi[i],gen_mass[i])
                gen_p4.append(ip4)  
            # Do matching
            for i in range(4):
                delR_ = reco_p4[i].DeltaR(gen_p4[i])
                if delR_<0.2: isMatched.push_back(1)
                else : isMatched.push_back(0)
                delR.push_back(delR_)
        # Fill this entry
        newtree.Fill()

    # Write and close files
    fout.Write()
    fout.Close()
    tfile.Close()

    return 'Get angles finished!'

main()
