import sys

argv = sys.argv[1:]

jobs = open('ana.listOfJobs','r')
new_jobs  = open('ana.listOfJobs_new','w')

# first get list of jobs that succeeds
white_list = []
for ifile in argv:
    #tmp/TT_8TeV-mcatnlo_selected_reco_80000_82000.root
#    print ifile
    white_list.append(int(ifile.split('/')[-1].split('_reco_')[-1].split('_')[0]))
white_list = set(white_list)
print 'number of files in white list = %i'%len(white_list)

# then modify the original job list by removing the jobs in whitelist
towrite = ''
for line in jobs: 
# python ./tardir/top_reco.py --inputfiles /uscms_data/d3/lfeng7/B2G_FW/CMSSW_7_2_0/src/Legacy_Afb/JHUntuple_selector/selected_files/v4_JEC/all/TT_8TeV-mcatnlo_selected.root --evtstart 0 --evtsperjob 2000
    if '--evtstart' not in line: continue
    evtstart = line.split('--evtstart')[-1].split('--')[0]
    evtstart = int(evtstart.strip())
    if evtstart not in white_list:
         towrite += line
new_jobs.write(towrite)
new_jobs.close()

    

