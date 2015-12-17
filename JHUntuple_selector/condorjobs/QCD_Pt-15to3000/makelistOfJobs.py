fout = open('ana.listOfJobs','w')
num_files = 41 
filesperjob = 5 
startfile = 0
txtfile = 'QCD_Pt-15to3000.txt'
sampletype = 'QCD_Pt-15to3000'
mcordata = 'mc'
fakelep = 'no'
while startfile < num_files :
	toprint = 'python ./tardir/selection.py --txtfiles tardir/inputfiles/'+txtfile+' --makeplots no --mcordata '+mcordata
	toprint += ' --fakelep '+fakelep+' --mctype qcd --maxevts -1 --type '+sampletype+' --grid yes --maxfiles '+str(filesperjob)+' --startfile '+str(startfile)+'\n'
	fout.write(toprint)
	startfile += filesperjob

fout.close()
