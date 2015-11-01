from fwlite_boilerplate import *

def makeplots(data,title="A Simple Graph Example",x_title="X title",y_title="Y title",out_fn="example.png"):
#Draw a simple graph
    c1 = ROOT.TCanvas("c1",title,200,10,700,500);
    c1.SetFillColor(42);
    c1.SetGrid();
    # unpack the 2d lists into x and y array
    import numpy as np
    x = np.asarray(list(x for x,y in data))
    y = np.asarray(list(y for x,y in data)) 
    
    gr = ROOT.TGraph(len(data),x,y);
    gr.SetLineColor(2);
    gr.SetLineWidth(4);
    gr.SetMarkerColor(4);
    gr.SetMarkerStyle(21);
    gr.SetTitle(title);
    gr.GetXaxis().SetTitle(x_title);
    gr.GetYaxis().SetTitle(y_title);
    gr.Draw("ACP");

    # TCanvas::Update() draws the frame, after which one can change it
    c1.Update();
    c1.GetFrame().SetFillColor(21);
    c1.GetFrame().SetBorderSize(12);
    c1.Modified();
    c1.SaveAs(out_fn)

def main():
    import numpy as np
    theta = np.linspace(0,np.pi,100)
    example_fc = list(
        (angle,np.sin(angle))
        for angle in theta)
    makeplots(example_fc)

if __name__ == "__main__":
    main()

