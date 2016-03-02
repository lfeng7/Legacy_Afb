# This file contains all necessary plotting helper modules
import ROOT
def overflow(hist):
# actual program
    h = hist.Clone()

    name  = h.GetName();
    title = h.GetTitle();
    nx    = h.GetNbinsX()+1;
    x1 = h.GetBinLowEdge(1);
    bw = h.GetBinWidth(nx);
    x2 = h.GetBinLowEdge(nx)+bw;

    # Book a temporary histogram having ab extra bin for overflows
    htmp = ROOT.TH1F(name+'_overflow', title, nx, x1, x2);

    # Fill the new hitogram including the extra bin for overflows
    # The nbinx+1 bin is the overflow bin
    for i in range(1,nx+1):
       # htmp.Fill(htmp.GetBinCenter(i), h.GetBinContent(i));
       htmp.SetBinContent(i,h.GetBinContent(i))

    # Fill the underflows
    # The 0 bin is underflow bin
    htmp.Fill(x1-1, h.GetBinContent(0));

    # Restore the number of entries
    htmp.SetEntries(h.GetEntries());

    return htmp