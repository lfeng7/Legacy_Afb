import ROOT

def M3(jets):
    if not len(jets)>= 3: return 0
    jets_pt = list(ipar.pt() for ipar in jets)
    jets_pt.sort()
    return jets_pt[0]+jets_pt[1]+jets_pt[2]

def GetSampleColor(sample_type):
    if sample_type == 'singletop': return ROOT.kMagenta
    elif sample_type == 'wjets'    : return ROOT.kGreen-3
    elif sample_type == 'ttbar'    : return ROOT.kRed+1
    elif sample_type == 'zjets'    : return ROOT.kAzure-2
    else :
        return 0

def GetSampleType(color):
    if color == ROOT.kMagenta : return 'singletop'
    elif color == ROOT.kGreen-3 : return 'wjets'  
    elif color == ROOT.kRed+1 : return 'ttbar' 
    elif color == ROOT.kAzure-2 : return 'zjets'
    elif color == 0 : return 'data'
    else : return 'unknown'

# This will return a list contain the number of entries in each bin
def GetBinEntry(hist) :
    nbins = hist.GetSize()
    entry = []
    for i in range(1,nbins) : 
        entry.append(int(hist.GetArray()[i]))
    return entry
