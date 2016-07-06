fout = open('ana.listOfJobs','w')
num_files = 1
filesperjob =  10
startfile = 0
txtfile = 'Tbar_s.txt'
sampletype = 'Tbar_s'
mcordata = 'mc'
while startfile < num_files :
	toprint = 'python ./tardir/selection.py --txtfiles tardir/inputfiles/'+txtfile+' --makeplots no --mcordata '+mcordata
	toprint += '  --selection_type sideband --mctype singletopbar  --maxevts -1 --type '+sampletype+' --grid yes --maxfiles '+str(filesperjob)+' --startfile '+str(startfile)+'\n'
	fout.write(toprint)
	startfile += filesperjob

fout.close()
