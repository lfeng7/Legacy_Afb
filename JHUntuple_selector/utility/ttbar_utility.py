import ROOT

def M3(jets):
    if not len(jets)>= 3: return 0
    jets_pt = list(ipar.pt() for ipar in jets)
    jets_pt.sort()
    return jets_pt[0]+jets_pt[1]+jets_pt[2]

def sample_color(sample_type):
    if sample_type == 'singletop': return ROOT.kMagenta
    elif sample_type == 'wjets'    : return ROOT.kGreen-3
    elif sample_type == 'ttbar'    : return ROOT.kRed+1
    elif sample_type == 'zjets'    : return ROOT.kAzure-2
    else :
        return 0

