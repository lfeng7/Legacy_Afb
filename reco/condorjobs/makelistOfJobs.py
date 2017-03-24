import ROOT
import glob
import sys

# global constants
evtsperjob =  2000
postfix = '.root'

argv = sys.argv[1:]
if argv ==[]:
    print 'Usage: python makelistofjobs.py el/mu ct10/cteq dir/to/selected/files/'
    sys.exit(1)
lep_type = argv.pop(0)
pdf_name = argv.pop(0)
prepend = argv.pop(0)
print 'Dir of selected files: %s'%prepend
print 'Lep type:',lep_type
print 'Pdf: %s'%pdf_name

def make_joblist():
    fout = open('ana.listOfJobs','w')

    all_files = glob.glob(prepend)
    print all_files

    towrite = ''
    # Loop over files
    for ifile in all_files:
#        if ifile not in all_files: continue
        print 'making jobs for',ifile
        towrite += '\n# '+ifile+'\n'
        tmpfile = ROOT.TFile(ifile)
        nev = tmpfile.Get('selected').GetEntries()
        evtstart = 0
	print '%i events'%nev
        while evtstart<nev:
            towrite += 'python ./tardir/top_reco.py --lep_type %s --PDF_name %s --inputfiles %s --evtstart %i --evtsperjob %i \n'%(lep_type, pdf_name, ifile, evtstart, evtsperjob)
            evtstart += evtsperjob
        tmpfile.Close()

    fout.write(towrite)
    fout.close()

def LoadInputs():
    global flist 
    flist= []
       # 0,                 1         2          3         4                   5
    # (MC_sample_name, sample_type, Nevts_gen, x-sec, nevts_total_ntuple, btag_type)
    # Single Top
    flist.append(['T_s','singletop',259961,3.79,259176,'singletop'] )
    flist.append(['T_t','singletop',3758227,56.4,3748155,'singletop'] )
    flist.append(['T_tW','singletop',497658,11.1,495559,'singletop'])
    flist.append(['Tbar_s','singletop',139974, 1.76,139604,'singletopbar'])
    flist.append(['Tbar_t','singletop',1935072, 30.7,1930185,'singletopbar'])
    flist.append(['Tbar_tW','singletop',493460,11.1,491463,'singletopbar'])
    # Wjets
    flist.append(['W1JetsToLNu_TuneZ2Star_8TeV','wjets',23141598,6662.8,23038253,'wjets'])
    flist.append(['W2JetsToLNu_TuneZ2Star_8TeV','wjets',34044921,2159.2,33993463,'wjets'])
    flist.append(['W3JetsToLNu_TuneZ2Star_8TeV','wjets',15539503,640.4,15507852,'wjets'])
    flist.append(['W4JetsToLNu_TuneZ2Star_8TeV','wjets',13382803,246.0,13326400,'wjets'])
    # DYjets
    flist.append(['DY1JetsToLL_M','zjets',24045248,660.6,23802736,'zjets'])
    flist.append(['DY2JetsToLL_M','zjets',2352304,215.1,2345857,'zjets'])
    flist.append(['DY3JetsToLL_M','zjets',11015445,65.79,10655325,'zjets'])
    flist.append(['DY4JetsToLL_M','zjets',6402827,28.59,5843425,'zjets'])
    # QCD
    #flist.append(['QCD_Pt-15to3000','qcd',9991674,6662.6,9940092,'qcd'])
    # signal
    flist.append(['TT_CT10_TuneZ2star_8TeV','ttbar',21675970,245.9,21560109,'ttbar'])    
    flist.append(['TT_8TeV-mcatnlo','ttbar',21675970,245.9,21560109,'ttbar'])
    # data
    flist.append(['SingleEl_Run2012ABCD','data',19748])

make_joblist()
