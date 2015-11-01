flist = []
# Single Top
flist.append(['T_s','singletop',259961,3.79,259176] )
flist.append(['T_t','singletop',3758227,56.4,3748155] )
flist.append(['T_tW','singletop',497658,11.1,495559])
flist.append(['Tbar_s','singletopbar',139974, 1.76,139604])
flist.append(['Tbar_t','singletopbar',1935072, 30.7,1930185])
flist.append(['Tbar_tW','singletopbar',493460,11.1,491463])
# Wjets
flist.append(['W1JetsToLNu_TuneZ2Star_8TeV','wjets',23141598,6662.8,23038253])
flist.append(['W2JetsToLNu_TuneZ2Star_8TeV','wjets',34044921,2159.2,33993463])
flist.append(['W3JetsToLNu_TuneZ2Star_8TeV','wjets',15539503,640.4,15507852])
flist.append(['W4JetsToLNu_TuneZ2Star_8TeV','wjets',13382803,246.0,13326400])
# DYjets
flist.append(['DY1JetsToLL_M','zjets',24045248,660.6,23802736])
flist.append(['DY2JetsToLL_M','zjets',2352304,215.1,2345857])
flist.append(['DY3JetsToLL_M','zjets',11015445,65.79,10655325])
flist.append(['DY4JetsToLL_M','zjets',6402827,28.59,5843425])
# signal
flist.append(['TT_CT10_TuneZ2star_8TeV','ttbar',21675970,245.9,21560109])


def MakeBtaggingEfficiency():
    # Get list of selected files for each type of events
    all_types = ['singletop','singletopbar','wjets','zjets','ttbar']    
    selected_files = []
    for itype in all_types:
        selected_files.append([ifile[0] for ifile in flist if ifile[1] == itype ])
    for item in selected_files: 
        for i in item: print i

MakeBtaggingEfficiency()
