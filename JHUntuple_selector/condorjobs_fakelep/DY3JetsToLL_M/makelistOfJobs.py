fout = open('ana.listOfJobs','w')
num_files = 70
filesperjob =  10
startfile = 0
txtfile = 'DY3JetsToLL_M.txt'
sampletype = 'DY3JetsToLL_M'
mcordata = 'mc'
while startfile < num_files :
	toprint = 'python ./tardir/selection.py --txtfiles tardir/inputfiles/'+txtfile+' --makeplots no --mcordata '+mcordata
	toprint += ' --mctype zjets  --maxevts -1 --type '+sampletype+' --grid yes --maxfiles '+str(filesperjob)+' --startfile '+str(startfile)+'\n'
	fout.write(toprint)
	startfile += filesperjob

fout.close()
