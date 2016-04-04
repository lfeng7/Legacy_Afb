//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Sep 26 00:33:19 2013 by ROOT version 5.30/06
// from TTree angles_data/Recreate
// found on file: all_pp_z_tt~_data_angles.root
//////////////////////////////////////////////////////////

#ifndef angles_data_h
#define angles_data_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

class angles_data {
	public :
	TTree          *fChain;   //!pointer to the analyzed TTree or TChain
	Int_t           fCurrent; //!current Tree number in a TChain
	
	// Declaration of leaf types
	Float_t        ttbar_mass;
	Float_t        cos_theta_cs;
	Float_t        Feynman_x;
	Int_t			Q_l;
	Int_t			n_valid_jets;
	Float_t			ln_L;
	Int_t			n_bTags;

	// List of branches
	TBranch        *b_ttbar_mass;   //!
	TBranch        *b_cos_theta_cs;   //!
	TBranch        *b_Feynman_x;   //!
	TBranch        *b_Q_l;   //!
	TBranch        *b_n_valid_jets;   //!
	TBranch        *b_ln_L;   //!
	TBranch        *b_n_bTags;   //!
	
	angles_data(TTree *tree=0);
	angles_data(const char* input_filename);
	virtual ~angles_data();
	virtual Int_t    Cut(Long64_t entry);
	virtual Int_t    GetEntry(Long64_t entry);
	virtual Long64_t LoadTree(Long64_t entry);
	virtual void     Init(TTree *tree);
	virtual double   Loop(double Rqqbar, double Rbck, double RWJets, double xi, double delta, double Afb, int plot);
	virtual Bool_t   Notify();
	virtual void     Show(Long64_t entry = -1);
	void LoadHistogramsCombined();
};

#endif

#ifdef angles_data_cxx
angles_data::angles_data(TTree *tree)
{
	// if parameter tree is not specified (or zero), connect the file
	// used to generate this class and read the Tree.
	if (tree == 0) {
		TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("all_pp_z_tt~_data_angles.root");
		if (!f || !f->IsOpen()) {
			f = new TFile("all_pp_z_tt~_data_angles.root");
		}
		f->GetObject("angles_data",tree);
		
	}
	Init(tree);
}

angles_data::angles_data(const char* input_filename) : fChain(0) 
{
	TTree *tree;
	// connect the file to the tree in the filename
	string inputfile;
	TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("test.root");
	if (!f || !f->IsOpen()) {
		f = new TFile(input_filename);
	}
	f->GetObject("angles_data",tree);
	
	Init(tree);
}

angles_data::~angles_data()
{
	if (!fChain) return;
	delete fChain->GetCurrentFile();
}

Int_t angles_data::GetEntry(Long64_t entry)
{
	// Read contents of entry.
	if (!fChain) return 0;
	return fChain->GetEntry(entry);
}
Long64_t angles_data::LoadTree(Long64_t entry)
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

void angles_data::Init(TTree *tree)
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
	
	fChain->SetBranchAddress("ttbar_mass", &ttbar_mass, &b_ttbar_mass);
	fChain->SetBranchAddress("cos_theta_cs", &cos_theta_cs, &b_cos_theta_cs);
	fChain->SetBranchAddress("Feynman_x", &Feynman_x, &b_Feynman_x);
	fChain->SetBranchAddress("Q_l", &Q_l, &b_Q_l);
	fChain->SetBranchAddress("n_valid_jets", &n_valid_jets, &b_n_valid_jets);
	fChain->SetBranchAddress("lnL", &ln_L, &b_ln_L);
	fChain->SetBranchAddress("n_bTags",&n_bTags,&b_n_bTags); 


	Notify();
}

Bool_t angles_data::Notify()
{
	// The Notify() function is called when a new file is opened. This
	// can be either for a new TTree in a TChain or when when a new TTree
	// is started when using PROOF. It is normally not necessary to make changes
	// to the generated code, but the routine can be extended by the
	// user if needed. The return value is currently not used.
	
	return kTRUE;
}

void angles_data::Show(Long64_t entry)
{
	entry=entry;
	// Print contents of entry.
	// If entry is not specified, print current entry
	if (!fChain) return;
	fChain->Show(entry);
}
Int_t angles_data::Cut(Long64_t entry)
{
	entry=entry;
	// This function may be called from Loop.
	// returns  1 if entry is accepted.
	// returns -1 otherwise.
	return 1;
}
#endif // #ifdef angles_data_cxx
