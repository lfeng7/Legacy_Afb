fout = open('ana.listOfJobs','w')
num_files = 119 
filesperjob =  10
startfile = 0
txtfile = 'W1JetsToLNu_TuneZ2Star_8TeV.txt'
sampletype = 'W1JetsToLNu_TuneZ2Star_8TeV'
mcordata = 'mc'
while startfile < num_files :
	toprint = 'python ./tardir/selection.py --txtfiles tardir/inputfiles/'+txtfile+' --makeplots no --mcordata '+mcordata
	toprint += ' --mctype wjets --maxevts -1 --type '+sampletype+' --grid yes --maxfiles '+str(filesperjob)+' --startfile '+str(startfile)+'\n'
	fout.write(toprint)
	startfile += filesperjob

fout.close()
