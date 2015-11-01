fout = open('ana.listOfJobs','w')
num_files = 188
filesperjob =  10
startfile = 0
txtfile = 'TT_CT10_TuneZ2star_8TeV.txt'
sampletype = 'TT_CT10_TuneZ2star_8TeV'
mcordata = 'mc'
while startfile < num_files :
	toprint = 'python ./tardir/selection.py --txtfiles tardir/inputfiles/'+txtfile+' --makeplots no --mcordata '+mcordata
	toprint += ' --mctype ttbar --maxevts -1 --type '+sampletype+' --grid yes --maxfiles '+str(filesperjob)+' --startfile '+str(startfile)+'\n'
	fout.write(toprint)
	startfile += filesperjob

fout.close()
