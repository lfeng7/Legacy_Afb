def M3(jets):
    if not len(jets)>= 3: return 0
    jets_pt = list(ipar.pt() for ipar in jets)
    jets_pt.sort()
    return jets_pt[0]+jets_pt[1]+jets_pt[2]
    

