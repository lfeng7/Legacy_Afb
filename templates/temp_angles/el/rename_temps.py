import glob
allfiles= glob.glob('*.root')
fout = open('do_rename.sh','w')
towrite = '' 
for ifile in allfiles:
    fname = ifile.split('_')
    postfix='.root'
    if 'gg' in ifile: postfix = '_gg.root'
    if 'qq' in ifile: postfix = '_qq.root'
    if 'bkg' in ifile: postfix = '_bkg.root'
    if 'Run' in ifile:
        prefix = 'angles_data_'
    else :
        prefix = 'angles_'
    fname = prefix+fname[0]+'_'+fname[1]+postfix
    towrite += 'mv '+ifile+'	'+fname+'\n'
fout.write(towrite)
fout.close()
