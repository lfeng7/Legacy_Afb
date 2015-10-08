import ROOT

def M3(jets):
    if not len(jets)>= 3: return 0
    jets_pt = list(ipar.pt() for ipar in jets)
    jets_pt.sort()
    return jets_pt[0]+jets_pt[1]+jets_pt[2]

def sample_color(sample_type):
    type_color = [('singletop',ROOT.kMagenta),('ttbar',ROOT.kRed+1),('zjets',ROOT.kAzure-2)]
    type_color.extend[('wjets',ROOT.kGreen-3)]    
    color = [ icolor[1] for icolor in type_color if icolor[0] == sample_type]
    if len(color)==1 :
        return color[0]
    else :
        return 0
        print 'the sample type is not in the color list!'

