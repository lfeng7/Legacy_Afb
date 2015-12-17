# Take in the selection output ttree, do reco reconstruction via kinfit

from Legacy_Afb.Tools.ttbar_utility import *
from Legacy_Afb.Tools.angles_tools import *
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
                  default = 'no',
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
                  default = 1,
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
        print 'Do makeAngles now!'
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
    fout_name = sample_name+'_angles_'+str(evt_start)+'.root'
    fout = ROOT.TFile(fout_name,'recreate')
    # Make output ttree
    if options.slim == 'yes':
        tmptree.SetBranchStatus('*',0)
        tmptree.SetBranchStatus('final_chi2',1)
    newtree = tmptree.CloneTree(0)
    tmptree.SetBranchStatus('*',1)

    # Add new branches to the output tree
    vecs = []
    br_names = []
    # Angles and other template vars
    cos_theta = ROOT.vector('float')()
    xf = ROOT.vector('float')() 
    mtt = ROOT.vector('float')()

    cos_theta_mc = ROOT.vector('float')()
    xf_mc = ROOT.vector('float')() 
    mtt_mc = ROOT.vector('float')()
    init_type = ROOT.vector('string')() 

    vecs += [cos_theta,xf,mtt,cos_theta_mc,xf_mc,mtt_mc,init_type]
    br_names += ['cos_theta','xf','mtt','cos_theta_mc','xf_mc','mtt_mc','init_type']

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
        print 'Is ttbar MC sample! Will get true angles ,x_f and Mtt!'

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
        if not tmptree.trigger[0] and options.applytrigger == 'yes' : continue
        h_cutflow.Fill('trigger',1)
        # Skip events that kinfit did not converge ( error flag == 4 )
        if not tmptree.final_errflags[0]==0 : continue
        h_cutflow.Fill('kinfit error',1)

        # Book the branches from input ttree
        # reco_p4 is a list of tlep_p4,thad_p4,wlep_p4,whad_p4
        tlep_pt = tmptree.reco_pt[0]
        tlep_eta = tmptree.reco_eta[0]
        tlep_phi = tmptree.reco_phi[0]
        tlep_mass = tmptree.reco_mass[0]
        thad_pt = tmptree.reco_pt[1]
        thad_eta = tmptree.reco_eta[1]
        thad_phi = tmptree.reco_phi[1]
        thad_mass = tmptree.reco_mass[1]        
        lep_charge = tmptree.lep_charge[0]
        # Make 4vec of tlep and thad
        tlep_p4 = ROOT.TLorentzVector()
        thad_p4 = ROOT.TLorentzVector()
        tlep_p4.SetPtEtaPhiM(tlep_pt,tlep_eta,tlep_phi,tlep_mass)
        thad_p4.SetPtEtaPhiM(thad_pt,thad_eta,thad_phi,thad_mass)
        # Assign t and tbar according to the lep-charge
        if lep_charge >0 :
            t_p4 = tlep_p4.Clone()
            tbar_p4 = thad_p4.Clone()
        elif lep_charge <0 :
            t_p4 = thad_p4.Clone()
            tbar_p4 = tlep_p4.Clone()
        else :
            print 'lep_charge =',lep_charge,'will skip this event!'
            continue
        # Do angle calculation
        results = get_angles(t_p4,tbar_p4)
        xf.push_back(results[0])
        mtt.push_back(results[1])
        cos_theta.push_back(results[2])
    
        # Get true angles for ttbar MC
        # Set default value of true quantities. If the event is not qqbar->ttbar(j) events will pushback the default values.
        true_xf = -2.
        true_mtt = -2.
        true_cos_theta = -2.
        tmp_type = 'bkg'
        # Get new values for true qqbar->ttbar events
        if is_ttbar == 1:
            gen_pt = tmptree.gen_pt
            gen_eta = tmptree.gen_eta
            gen_phi = tmptree.gen_phi
            gen_mass = tmptree.gen_mass
            gen_pdgid = tmptree.gen_pdgid
            # Check if it is qq->ttbar
            if gen_pdgid[0]+gen_pdgid[1] ==0 : 
                tmp_type = 'qqbar'
            elif gen_pdgid[0] == 21 and gen_pdgid[1] == 21:
                tmp_type = 'gg'
            elif gen_pdgid[0] != 21 and gen_pdgid[1] != 21:
                tmp_type = 'q1q2'
            elif gen_pdgid[0] == 21 or gen_pdgid[1] == 21:
                tmp_type = 'qg'
            else : 
                tmp_type = 'unknown'

            # Make 4vecs for t,tbar,q,qbar
            true_t_p4 = ROOT.TLorentzVector()
            true_tbar_p4 = ROOT.TLorentzVector()
            true_q_p4 = ROOT.TLorentzVector()
            # Loop over gen pars to set p4
            q_pdgids = [1,2,3,4]
            q_ids,q_pzs = [],[]
            for i in range(gen_pt.size()):
                ipt = gen_pt[i]; ieta = gen_eta[i]; iphi = gen_phi[i]; imass = gen_mass[i];
                if gen_pdgid[i] in q_pdgids : 
                    true_q_p4.SetPtEtaPhiM(ipt,ieta,iphi,imass)
                    q_ids.append(gen_pdgid[i])
                    q_pzs.append(true_q_p4.Pz())
                if gen_pdgid[i] == 6 : true_t_p4.SetPtEtaPhiM(ipt,ieta,iphi,imass)
                if gen_pdgid[i] == -6 : true_tbar_p4.SetPtEtaPhiM(ipt,ieta,iphi,imass)
            # Get true cos,xf,mtt
            q_pz = 0.0
            if len(q_ids) == 1:
                q_pz = q_pzs[0]
            true_results = get_true_angles(true_t_p4,true_tbar_p4,q_pz)
            true_xf = true_results[0]
            true_mtt = true_results[1]
            true_cos_theta = true_results[2]           
            # push back true quantities
            xf_mc.push_back(true_xf)
            mtt_mc.push_back(true_mtt)
            cos_theta_mc.push_back(true_cos_theta)

        init_type.push_back(tmp_type)

        # Fill this entry
        if options.slim == 'yes' and tmp_type!='qqbar':continue
        newtree.Fill()

    # Write and close files
    fout.Write()
    fout.Close()
    tfile.Close()

    return 'Get angles finished!'

main()
