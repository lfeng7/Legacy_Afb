//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Dec 18 10:06:18 2013 by ROOT version 5.32/00
// from TTree output/output
// found on file: tagged_semilep_0.root
//////////////////////////////////////////////////////////

#ifndef ttbar_h
#define ttbar_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class ttbar {
	public :
	TTree          *fChain;   //!pointer to the analyzed TTree or TChain
	Int_t           fCurrent; //!current Tree number in a TChain
	
	// Declaration of leaf types
	Float_t         pt[8];
	Float_t         eta[8];
	Float_t         phi[8];
	Float_t         mass[8];
	Float_t         px[8];
	Float_t         py[8];
	Float_t         pz[8];
	Float_t         E[8];
	Float_t         PID[8];
	Float_t         is_leptonic_side[8];
	Float_t			bestFitParValues[6];
	Float_t			pileup_events;
	Float_t         finalChi[8];
	UInt_t          nbTags;
	Float_t         mc_pt[10];
	Float_t         mc_eta[10];
	Float_t         mc_phi[10];
	Float_t         mc_mass[10];
	Float_t         mc_px[10];
	Float_t         mc_py[10];
	Float_t         mc_pz[10];
	Float_t         mc_E[10];
	UInt_t          mc_was_matched[8];
	UInt_t			nValidJets;
	Double_t		Pdf_weights[100];
	Float_t			qScale;
	Float_t			parton_id[2];
	Float_t			parton_x[2];
	Int_t			motherParticles[2];
	Float_t 		weight_GJR_scale;
	Float_t 		weight_CT10_scale;
	Float_t 		weight_cteq_scale;
	Float_t 		weight_top_pT;
	Float_t 		weight_btag_eff;
	Float_t 		weight_btag_eff_err;
	Float_t 		weight_pileup;
	Float_t 		weight_tracking;
	Float_t 		weight_tracking_low;
	Float_t 		weight_tracking_hi;
	Float_t 		weight_lep_ID;
	Float_t 		weight_lep_ID_low;
	Float_t 		weight_lep_ID_hi;
	Float_t 		weight_lep_iso;
	Float_t 		weight_lep_iso_low;
	Float_t 		weight_lep_iso_hi;
	Float_t 		weight_trig_eff;
	Float_t 		weight_trig_eff_low;
	Float_t 		weight_trig_eff_hi;
	Float_t 		mc_pileup_events;
	
	// List of branches
	TBranch        *b_pt;   //!
	TBranch        *b_eta;   //!
	TBranch        *b_phi;   //!
	TBranch        *b_mass;   //!
	TBranch        *b_px;   //!
	TBranch        *b_py;   //!
	TBranch        *b_pz;   //!
	TBranch        *b_E;   //!
	TBranch        *b_PID;   //!
	TBranch        *b_is_leptonic_side;   //!
	TBranch		   *b_bestFitParValues;
	TBranch		   *b_pileup_events;
	TBranch        *b_finalChi;   //!
	TBranch        *b_nbTags;   //!
	TBranch        *b_mc_pt;   //!
	TBranch        *b_mc_eta;   //!
	TBranch        *b_mc_phi;   //!
	TBranch        *b_mc_mass;   //!
	TBranch        *b_mc_px;   //!
	TBranch        *b_mc_py;   //!
	TBranch        *b_mc_pz;   //!
	TBranch        *b_mc_E;   //!
	TBranch        *b_mc_was_matched;   //!
	TBranch		   *b_nValidJets;	//!
	TBranch		   *b_evt_weight;
	TBranch		   *b_Pdf_weights;
	TBranch		   *b_qScale;
	TBranch		   *b_parton_id;
	TBranch		   *b_parton_x;
	TBranch		   *b_motherParticles;
	TBranch 	   *b_weight_GJR_scale;
	TBranch 	   *b_weight_CT10_scale;
	TBranch 	   *b_weight_cteq_scale;
	TBranch 	   *b_weight_top_pT;
	TBranch 	   *b_weight_btag_eff;
	TBranch 	   *b_weight_btag_eff_err;
	TBranch 	   *b_weight_pileup;
	TBranch 	   *b_weight_tracking;
	TBranch 	   *b_weight_tracking_low;
	TBranch 	   *b_weight_tracking_hi;
	TBranch 	   *b_weight_lep_ID;
	TBranch 	   *b_weight_lep_ID_low;
	TBranch 	   *b_weight_lep_ID_hi;
	TBranch 	   *b_weight_lep_iso;
	TBranch 	   *b_weight_lep_iso_low;
	TBranch 	   *b_weight_lep_iso_hi;
	TBranch 	   *b_weight_trig_eff;
	TBranch 	   *b_weight_trig_eff_low;
	TBranch 	   *b_weight_trig_eff_hi;
	TBranch 	   *b_mc_pileup_events;
	
	ttbar(TTree *tree=0);
	ttbar(const char* input_filename);
	void SetOutputFilename(const char* outname);
	void SetOutputTreeName(const char* treename);
	virtual ~ttbar();
	virtual Int_t    Cut(Long64_t entry);
	virtual Int_t    GetEntry(Long64_t entry);
	virtual Long64_t LoadTree(Long64_t entry);
	virtual void     Init(TTree *tree);
	virtual void     Loop(int data_or_mc);
	virtual Bool_t   Notify();
	virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef ttbar_cxx
ttbar::ttbar(TTree *tree) : fChain(0) 
{
	// if parameter tree is not specified (or zero), connect the file
	// used to generate this class and read the Tree.
	if (tree == 0) {
		TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("TEST.root");
		if (!f || !f->IsOpen()) {
			f = new TFile("tagged_semilep_0.root");
		}
		f->GetObject("output",tree);
		
	}
	Init(tree);
}

ttbar::ttbar(const char* input_filename) : fChain(0) 
{
	TTree *tree;
	// connect the file to the tree in the filename
	string inputfile;
	TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("test.root");
	if (!f || !f->IsOpen()) {
		f = new TFile(input_filename);
	}
	f->GetObject("output",tree);
	
	Init(tree);
}

ttbar::~ttbar()
{
	if (!fChain) return;
	delete fChain->GetCurrentFile();
}

Int_t ttbar::GetEntry(Long64_t entry)
{
	// Read contents of entry.
	if (!fChain) return 0;
	return fChain->GetEntry(entry);
}
Long64_t ttbar::LoadTree(Long64_t entry)
{
	// Set the environment to read one entry
	if (!fChain) return -5;
	Long64_t centry = fChain->LoadTree(entry);
	if (centry < 0) return centry;
	if (fChain->GetTreeNumber() != fCurrent) {
		fCurrent = fChain->GetTreeNumber();
		Notify();
	}
	return centry;
}

void ttbar::Init(TTree *tree)
{
	// The Init() function is called when the selector needs to initialize
	// a new tree or chain. Typically here the branch addresses and branch
	// pointers of the tree will be set.
	// It is normally not necessary to make changes to the generated
	// code, but the routine can be extended by the user if needed.
	// Init() will be called many times when running on PROOF
	// (once per file to be processed).
	
	// Set branch addresses and branch pointers
	if (!tree) return;
	fChain = tree;
	fCurrent = -1;
	fChain->SetMakeClass(1);
	
	fChain->SetBranchAddress("pt", pt, &b_pt);
	fChain->SetBranchAddress("eta", eta, &b_eta);
	fChain->SetBranchAddress("phi", phi, &b_phi);
	fChain->SetBranchAddress("mass", mass, &b_mass);
	fChain->SetBranchAddress("px", px, &b_px);
	fChain->SetBranchAddress("py", py, &b_py);
	fChain->SetBranchAddress("pz", pz, &b_pz);
	fChain->SetBranchAddress("E", E, &b_E);
	fChain->SetBranchAddress("PID", PID, &b_PID);
	fChain->SetBranchAddress("is_leptonic_side", is_leptonic_side, &b_is_leptonic_side);
	fChain->SetBranchAddress("bestFitParValues",bestFitParValues,&b_bestFitParValues);
	fChain->SetBranchAddress("pileup_events",&pileup_events,&b_pileup_events);
	fChain->SetBranchAddress("finalChi", finalChi, &b_finalChi);
	fChain->SetBranchAddress("nbTags", &nbTags, &b_nbTags);
	fChain->SetBranchAddress("mc_pt", mc_pt, &b_mc_pt);
	fChain->SetBranchAddress("mc_eta", mc_eta, &b_mc_eta);
	fChain->SetBranchAddress("mc_phi", mc_phi, &b_mc_phi);
	fChain->SetBranchAddress("mc_mass", mc_mass, &b_mc_mass);
	fChain->SetBranchAddress("mc_px", mc_px, &b_mc_px);
	fChain->SetBranchAddress("mc_py", mc_py, &b_mc_py);
	fChain->SetBranchAddress("mc_pz", mc_pz, &b_mc_pz);
	fChain->SetBranchAddress("mc_E", mc_E, &b_mc_E);
	fChain->SetBranchAddress("mc_was_matched", mc_was_matched, &b_mc_was_matched);
	fChain->SetBranchAddress("nValidJets", &nValidJets, &b_nValidJets);
	fChain->SetBranchAddress("Pdf_weights", Pdf_weights, &b_Pdf_weights);
	fChain->SetBranchAddress("qScale", &qScale, &b_qScale);
	fChain->SetBranchAddress("parton_id", parton_id, &b_parton_id);
	fChain->SetBranchAddress("parton_x", parton_x, &b_parton_x);
	fChain->SetBranchAddress("motherParticles", motherParticles, &b_motherParticles);
	fChain->SetBranchAddress("weight_GJR_scale",&weight_GJR_scale,&b_weight_GJR_scale);
	fChain->SetBranchAddress("weight_CT10_scale",&weight_CT10_scale,&b_weight_CT10_scale);
	fChain->SetBranchAddress("weight_cteq_scale",&weight_cteq_scale,&b_weight_cteq_scale);
	fChain->SetBranchAddress("weight_top_pT",&weight_top_pT,&b_weight_top_pT);
	fChain->SetBranchAddress("weight_btag_eff",&weight_btag_eff,&b_weight_btag_eff);
	fChain->SetBranchAddress("weight_btag_eff_err",&weight_btag_eff_err,&b_weight_btag_eff_err);
	fChain->SetBranchAddress("weight_pileup",&weight_pileup,&b_weight_pileup);
	fChain->SetBranchAddress("weight_tracking",&weight_tracking,&b_weight_tracking);
	fChain->SetBranchAddress("weight_tracking_low",&weight_tracking_low,&b_weight_tracking_low);
	fChain->SetBranchAddress("weight_tracking_hi",&weight_tracking_hi,&b_weight_tracking_hi);
	fChain->SetBranchAddress("weight_lep_ID",&weight_lep_ID,&b_weight_lep_ID);
	fChain->SetBranchAddress("weight_lep_ID_low",&weight_lep_ID_low,&b_weight_lep_ID_low);
	fChain->SetBranchAddress("weight_lep_ID_hi",&weight_lep_ID_hi,&b_weight_lep_ID_hi);
	fChain->SetBranchAddress("weight_lep_iso",&weight_lep_iso,&b_weight_lep_iso);
	fChain->SetBranchAddress("weight_lep_iso_low",&weight_lep_iso_low,&b_weight_lep_iso_low);
	fChain->SetBranchAddress("weight_lep_iso_hi",&weight_lep_iso_hi,&b_weight_lep_iso_hi);
	fChain->SetBranchAddress("weight_trig_eff",&weight_trig_eff,&b_weight_trig_eff);
	fChain->SetBranchAddress("weight_trig_eff_low",&weight_trig_eff_low,&b_weight_trig_eff_low);
	fChain->SetBranchAddress("weight_trig_eff_hi",&weight_trig_eff_hi,&b_weight_trig_eff_hi);
	fChain->SetBranchAddress("mc_pileup_events",&mc_pileup_events,&b_mc_pileup_events);	
	Notify();
}

Bool_t ttbar::Notify()
{
	// The Notify() function is called when a new file is opened. This
	// can be either for a new TTree in a TChain or when when a new TTree
	// is started when using PROOF. It is normally not necessary to make changes
	// to the generated code, but the routine can be extended by the
	// user if needed. The return value is currently not used.
	
	return kTRUE;
}

void ttbar::Show(Long64_t entry)
{
	entry=entry;
	// Print contents of entry.
	// If entry is not specified, print current entry
	if (!fChain) return;
	fChain->Show(entry);
}
Int_t ttbar::Cut(Long64_t entry)
{
	entry=entry;
	// This function may be called from Loop.
	// returns  1 if entry is accepted.
	// returns -1 otherwise.
	return 1;
}
#endif // #ifdef ttbar_cxx
