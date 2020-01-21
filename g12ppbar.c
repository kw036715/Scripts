// Kevin Ward (original code by Nick Compton.)
// Last edited: May 22, 2019.

// Creates histograms of skimmed g12 data.
// Reads in a skim for two protons and one anti-proton, according to g12 procedures.

#include "TRint.h"
#include <string.h>
#include <iostream>
#include <stdio.h>
#include <fstream>
#include <TROOT.h>
#include <TChain.h>  
#include <TFile.h>
#include "TObject.h"
#include "TString.h"
#include "TCanvas.h"
#include "TGraph.h"
#include <TVector3.h>
#include "TLorentzVector.h"
#include "TGraphErrors.h"
#include "TSystem.h"
#include "TSpectrum.h"
#include "TH1I.h"
#include "TH3F.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TF1.h"
#include "TMath.h"
#include "TPaveStats.h"
#include "TVector3.h"
#include <stdlib.h>

using namespace std;

const int TRIGGER = 3;
 
int main(int __argc,char *__argv[]){

	char *outFileName = "PPPbar.root";
	extern int optind;   
	TROOT troot();  
	TFile outFile(outFileName,"recreate");
 

	// Declares parameters.
	Int_t nentries, j;
	Float_t egam, pxp1, pxp2, pxpb, pyp1, pyp2, pypb, pzp1, pzp2, pzpb;
	Float_t pbeta, pBeta, p2beta, p2Beta, pbbeta, pbBeta, p1delbeta, p2delbeta, pbdelbeta;
	Float_t vxp1, vyp1, vzp1, vxp2, vyp2, vzp2, vxpb, vypb, vzpb;

	Float_t MCegam, MCpxp1, MCpxp2, MCpxpb, MCpyp1, MCpyp2, MCpypb, MCpzp1, MCpzp2, MCpzpb;
	Float_t MCpbeta, MCpBeta, MCp2beta, MCp2Beta, MCpbbeta, MCpbBeta, MCp1delbeta, MCp2delbeta, MCpbdelbeta;
	Float_t MCvxp1, MCvyp1, MCvzp1, MCvxp2, MCvyp2, MCvzp2, MCvxpb, MCvypb, MCvzpb;

	Float_t     p1gpx, p1gpy, p1gpz, p1gpe;
	Float_t     p2gpx, p2gpy, p2gpz, p2gpE;
	Float_t     pbgpx, pbgpy, pbgpz, pbgpE;
	Float_t     pIgpx, pIgpy, pIgpz, pIgpE;

	Float_t xmm, xmm2;
	Float_t MCxmm, MCxmm2;

	Float_t xAve, yAve, zAve, pp1, pp2, ppb;
	Float_t MCxAve, MCyAve, MCzAve, MCpp1, MCpp2, MCppb;

	TLorentzVector P4p1, P4p2, P4pb, P4miss, P4target, P4pho, P4inc, P4all, P4prots, P4prMM;
	TLorentzVector P4MCp1, P4MCp2, P4MCpb, P4MCmiss, P4MCtarget, P4MCpho, P4MCinc, P4MCall, P4MCprots, P4MCprMM;
	TLorentzVector P4p1g, P4p2g, P4pbg;

	TVector3 V3p1, V3p2, V3pb;
	TVector3 V3MCp1, V3MCp2, V3MCpb;

    Float_t     boost, boostg, MCboost;

	Float_t     p1bcos, p2bcos, pbbcos;
	//Float_t     p1bphi, p2bphi, pbbphi;

	Float_t     MCp1bcos, MCp2bcos, MCpbbcos;
	//Float_t     MCp1bphi, MCp2bphi, MCpbbphi;

	//Float_t     MCsecp1, MCsecp2, MCsecpb;
	//Float_t     secp1, secp2, secpb;

	Float_t     p1bgcos, p2bgcos, pbbgcos;
	//Float_t     p1bgphi, p2bgphi, pbbgphi;

	Float_t     pb, pbe, ER1, ER2, ER3, Y, Ye, run, flux, runf, ER , cts, fl;
	Float_t     e1,e2,e3,e4,e5,e6,e7,e8,e9,e0;

	const Float_t xmpr = 0.938272; // Proton mass.

	const Float_t lum = 0.1 * 0.071 * 40 * 6e+23;

	/* ******************************************************************* */

	TH1F *Hegam = new TH1F("Hegam","Photon energy",200,0.0,6.0);
	TH1F *HMMsq = new TH1F("MMsq","MM^{2}(p p)",500,0.7,1.1);
	TH1F *HMMsqall = new TH1F("MMsqa","MM^{2}(p p pbar)",500,-0.05,0.05);
	TH1F *Hvz= new TH1F("Hvz","Multi-particle Vertex", 500,-120.0,0.0);

	TH2F *Hvxvy = new TH2F("VxVy","Y(multi) vs. X(multi) {vert}",100,-10,10,100,-10,10);

	TH1F *HMM = new TH1F("HMM","MM(p p pbar)",500,-0.25,0.25);
	TH1F *HMMb = new TH1F("MM base","MM(p p m-pbar)", 500, -0.5, 0.5);

	TH2F *Hp1beta = new TH2F("P1_betavsp","delBeta vs. Momentum (P1)",200,0.5,4.0,300,-0.15,0.15);
	TH2F *Hp2beta = new TH2F("P2_betavsp","delBeta vs. Momentum (P2)",200,0.5,4.0,300,-0.15,0.15);
	TH2F *Hpbbeta= new TH2F("Pb_betavsp","delBeta vs. Momentum (Pbar)",200,0.5,4.0,300,-0.15,0.15);	

	//TH1F *p1cos = new TH1F("P1_CosTheta","CosTheta(p)",400,-1.2,1.2);
	//TH1F *p2cos = new TH1F("P2_CosTheta","CosTheta(p)",400,-1.2,1.2);
	//TH1F *pbcos = new TH1F("Pb_CosTheta","CosTheta(pbar)",400,-1.2,1.2);

	//TH1F *p1phi = new TH1F("P1_Phi","Phi(p)",400,-3.2,3.2);
	//TH1F *p2phi = new TH1F("P2_Phi","Phi(p)",400,-3.2,3.2);
	//TH1F *pbphi = new TH1F("Pb_Phi","Phi(pbar)",400,-3.2,3.2);

	TH2F *p1vsegam = new TH2F("P1_EgamvsCosTheta","Egam(GeV) vs. CosTheta",10,4.5,5.5,20,-1,1);
	TH2F *p2vsegam = new TH2F("P2_EgamvsCosTheta","Egam(GeV) vs. CosTheta",10,4.5,5.5,20,-1,1);
	TH2F *pbvsegam = new TH2F("Pb_EgamvsCosTheta","Egam(GeV) vs. CosTheta",10,4.5,5.5,20,-1,1);

	/* ******************************************************************* */

	TH1F *HMCegam = new TH1F("HMCegam","MC Photon energy",200,0.0,6.0);
	TH1F *HMCMMsq = new TH1F("MCMMsq","MC MM^{2}(p p)",500,0.7,1.1);
	TH1F *HMCMMsqall = new TH1F("MCMMsqa","MC MM^{2}(p p pbar)",500,-1,1);
	TH1F *HMCvz= new TH1F("HMCvz","MC Multi-particle Vertex", 500,-120.0,0.0);

	TH2F *HMCvxvy = new TH2F("MCVxVy","MC Y(multi) vs. X(multi) {vert}",100,-10,10,100,-10,10);

	TH1F *HMCMM = new TH1F("MC2MM","MC MM(p p pbar)",500,0.84,1.05);
	TH1F *HMCMMb = new TH1F("MCMM base","MC MM(p p m-pbar)", 500, -0.5, 0.5);

	TH2F *HMCp1beta = new TH2F("MCP1_betavsp","MC delBeta vs. Momentum (P1)",200,0.5,4.0,300,-0.15,0.15);
	TH2F *HMCp2beta = new TH2F("MCP2_betavsp","MC delBeta vs. Momentum (P2)",200,0.5,4.0,300,-0.15,0.15);
	TH2F *HMCpbbeta= new TH2F("MCPb_betavsp","MC delBeta vs. Momentum (Pbar)",200,0.5,4.0,300,-0.15,0.15);

	//TH1F *MCp1cos = new TH1F("MCP1_CosTheta","MC CosTheta(p)",400,-1.2,1.2);
	//TH1F *MCp2cos = new TH1F("MCP2_CosTheta","MC CosTheta(p)",400,-1.2,1.2);
	//TH1F *MCpbcos = new TH1F("MCPb_CosTheta","MC CosTheta(pbar)",400,-1.2,1.2);

	//TH1F *MCp1phi = new TH1F("MCP1_Phi","MC Phi(p)",400,-3.2,3.2);
	//TH1F *MCp2phi = new TH1F("MCP2_Phi","MC Phi(p)",400,-3.2,3.2);
	//TH1F *MCpbphi = new TH1F("MCPb_Phi","MC Phi(pbar)",400,-3.2,3.2);

	TH2F *MCp1vsegam = new TH2F("MCP1_EgamvsCosTheta","MC Egam(GeV) vs. CosTheta",10,4.5,5.5,20,-1,1);
	TH2F *MCp2vsegam = new TH2F("MCP2_EgamvsCosTheta","MC Egam(GeV) vs. CosTheta",10,4.5,5.5,20,-1,1);
	TH2F *MCpbvsegam = new TH2F("MCPb_EgamvsCosTheta","MC Egam(GeV) vs. CosTheta",10,4.5,5.5,20,-1,1);

	/* ******************************************************************* */

	//TH1F *p1gcos = new TH1F("P1g_CosTheta","CosTheta(p)(Gen)",400,-1.2,1.2);
	//TH1F *p2gcos = new TH1F("P2g_CosTheta","CosTheta(p)(Gen)",400,-1.2,1.2);
	//TH1F *pbgcos = new TH1F("Pbg_CosTheta","CosTheta(pbar)(Gen)",400,-1.2,1.2);

	//TH1F *p1gphi = new TH1F("P1g_Phi","Phi(p)(Gen)",400,-3.2,3.2);
	//TH1F *p2gphi = new TH1F("P2g_Phi","Phi(p)(Gen)",400,-3.2,3.2);
	//TH1F *pbgphi = new TH1F("Pbg_Phi","Phi(pbar)(Gen)",400,-3.2,3.2);

	TH2F *p1gvsegam = new TH2F("P1g_EgamvsCosTheta","Egam(GeV) vs. CosTheta(Gen)",10,4.5,5.5,20,-1,1);
	TH2F *p2gvsegam = new TH2F("P2g_EgamvsCosTheta","Egam(GeV) vs. CosTheta(Gen)",10,4.5,5.5,20,-1,1);
	TH2F *pbgvsegam = new TH2F("Pbg_EgamvsCosTheta","Egam(GeV) vs. CosTheta(Gen)",10,4.5,5.5,20,-1,1);

	/* ******************************************************************* */

	TH2F *P1Accept = new TH2F("P1Acceptance","P1 Acceptance",10,4.5,5.5,20,-1,1);
	TH2F *P2Accept = new TH2F("P2Acceptance","P2 Acceptance",10,4.5,5.5,20,-1,1);
	TH2F *PbAccept = new TH2F("PbAcceptance","Pb Acceptance",10,4.5,5.5,20,-1,1);

	TH1F *PbAcc1 = new TH1F("Pb_Acceptance1","Egam(4.50-4.60 GeV), Acceptance(%) vs. CosTheta (Pbar)",20,-1,1);
	TH1F *PbAcc2 = new TH1F("Pb_Acceptance2","Egam(4.60-4.70 GeV), Acceptance(%) vs. CosTheta (Pbar)",20,-1,1);
	TH1F *PbAcc3 = new TH1F("Pb_Acceptance3","Egam(4.70-4.80 GeV), Acceptance(%) vs. CosTheta (Pbar)",20,-1,1);
	TH1F *PbAcc4 = new TH1F("Pb_Acceptance4","Egam(4.80-4.90 GeV), Acceptance(%) vs. CosTheta (Pbar)",20,-1,1);
	TH1F *PbAcc5 = new TH1F("Pb_Acceptance5","Egam(4.90-5.00 GeV), Acceptance(%) vs. CosTheta (Pbar)",20,-1,1);
	TH1F *PbAcc6 = new TH1F("Pb_Acceptance6","Egam(5.00-5.10 GeV), Acceptance(%) vs. CosTheta (Pbar)",20,-1,1);
	TH1F *PbAcc7 = new TH1F("Pb_Acceptance7","Egam(5.10-5.20 GeV), Acceptance(%) vs. CosTheta (Pbar)",20,-1,1);
	TH1F *PbAcc8 = new TH1F("Pb_Acceptance8","Egam(5.20-5.30 GeV), Acceptance(%) vs. CosTheta (Pbar)",20,-1,1);
	TH1F *PbAcc9 = new TH1F("Pb_Acceptance9","Egam(5.30-5.40 GeV), Acceptance(%) vs. CosTheta (Pbar)",20,-1,1);
	TH1F *PbAcc0 = new TH1F("Pb_Acceptance0","Egam(5.40-5.50 GeV), Acceptance(%) vs. CosTheta (Pbar)",20,-1,1);

	TH2F *P1Cross = new TH2F("P1Cross","P1 Cross",10,4.5,5.5,20,-1,1);
	TH2F *P2Cross = new TH2F("P2Cross","P2 Cross",10,4.5,5.5,20,-1,1);
	TH2F *PbCross = new TH2F("PbCross","Pb Cross",10,4.5,5.5,20,-1,1);

	TH1F *PbCross1 = new TH1F("Pb_Cross1","Egam(4.50-4.60 GeV), Cross Section(b) vs. CosTheta (Pbar)",20,-1,1);
	TH1F *PbCross2 = new TH1F("Pb_Cross2","Egam(4.60-4.70 GeV), Cross Section(b) vs. CosTheta (Pbar)",20,-1,1);
	TH1F *PbCross3 = new TH1F("Pb_Cross3","Egam(4.70-4.80 GeV), Cross Section(b) vs. CosTheta (Pbar)",20,-1,1);
	TH1F *PbCross4 = new TH1F("Pb_Cross4","Egam(4.80-4.90 GeV), Cross Section(b) vs. CosTheta (Pbar)",20,-1,1);
	TH1F *PbCross5 = new TH1F("Pb_Cross5","Egam(4.90-5.00 GeV), Cross Section(b) vs. CosTheta (Pbar)",20,-1,1);
	TH1F *PbCross6 = new TH1F("Pb_Cross6","Egam(5.00-5.10 GeV), Cross Section(b) vs. CosTheta (Pbar)",20,-1,1);
	TH1F *PbCross7 = new TH1F("Pb_Cross7","Egam(5.10-5.20 GeV), Cross Section(b) vs. CosTheta (Pbar)",20,-1,1);
	TH1F *PbCross8 = new TH1F("Pb_Cross8","Egam(5.20-5.30 GeV), Cross Section(b) vs. CosTheta (Pbar)",20,-1,1);
	TH1F *PbCross9 = new TH1F("Pb_Cross9","Egam(5.30-5.40 GeV), Cross Section(b) vs. CosTheta (Pbar)",20,-1,1);
	TH1F *PbCross0 = new TH1F("Pb_Cross0","Egam(5.40-5.50 GeV), Cross Section(b) vs. CosTheta (Pbar)",20,-1,1);

	TH1F *PbYield1 = new TH1F("Pb_Yield1","Egam(4.50-4.60 GeV), Counts vs. CosTheta (Pbar)",20,-1,1);
	TH1F *PbYield2 = new TH1F("Pb_Yield2","Egam(4.60-4.70 GeV), Counts vs. CosTheta (Pbar)",20,-1,1);
	TH1F *PbYield3 = new TH1F("Pb_Yield3","Egam(4.70-4.80 GeV), Counts vs. CosTheta (Pbar)",20,-1,1);
	TH1F *PbYield4 = new TH1F("Pb_Yield4","Egam(4.80-4.90 GeV), Counts vs. CosTheta (Pbar)",20,-1,1);
	TH1F *PbYield5 = new TH1F("Pb_Yield5","Egam(4.90-5.00 GeV), Counts vs. CosTheta (Pbar)",20,-1,1);
	TH1F *PbYield6 = new TH1F("Pb_Yield6","Egam(5.00-5.10 GeV), Counts vs. CosTheta (Pbar)",20,-1,1);
	TH1F *PbYield7 = new TH1F("Pb_Yield7","Egam(5.10-5.20 GeV), Counts vs. CosTheta (Pbar)",20,-1,1);
	TH1F *PbYield8 = new TH1F("Pb_Yield8","Egam(5.20-5.30 GeV), Counts vs. CosTheta (Pbar)",20,-1,1);
	TH1F *PbYield9 = new TH1F("Pb_Yield9","Egam(5.30-5.40 GeV), Counts vs. CosTheta (Pbar)",20,-1,1);
	TH1F *PbYield0 = new TH1F("Pb_Yield0","Egam(5.40-5.50 GeV), Counts vs. CosTheta (Pbar)",20,-1,1);

	TH2F *PbYield = new TH2F("Yield-Bac","YIELD",10,4.5,5.5,20,-1,1);
	TH2F *PbCro = new TH2F("CrossSec","CrossSec",10,4.5,5.5,20,-1,1);

	TH1F *PbCrossTot = new TH1F("Pb_CrossTot"," Total Cross Section Counts vs. Egam(GeV)(Pbar)",10,4.5,5.5);

	TH1F *Consist = new TH1F("Consistency", "Run Consistency",916,56401,57317);

	TH1F *Flux = new TH1F("Run_Flux", "Run Flux",916,56401,57317);
	TH1F *Counts = new TH1F("Run_Counts", "Run Counts",916,56401,57317);

	TH1F *cut1 = new TH1F("cut1","Min Energy Cut only",200,0.0,6.0);
	TH1F *cut1s = new TH1F("cut1s","Min Energy Cut seq",200,0.0,6.0);
	TH1F *cut2 = new TH1F("cut2","Vertex only",200,0.0,6.0);
	TH1F *cut2s = new TH1F("cut2s","Vertex seq",200,0.0,6.0);
	TH1F *cut3 = new TH1F("cut3","delBeta only",200,0.0,6.0);
	TH1F *cut3s = new TH1F("cut3s","delBeta seq",200,0.0,6.0);
	TH1F *cut4 = new TH1F("cut4","Missing Mass only",200,0.0,6.0);
	TH1F *cut4s = new TH1F("cut4s","Missing Mass seq",200,0.0,6.0);



	char hold[10];
	int x = 0;
	int yield;
	TH1F *Place[10][20];
	TH1 *Bac;

	for(int a = 1; a <= 10; a++){

		for(int b = 1; b <= 20; b++){

			sprintf(hold, "place%d", x);
			Place[a-1][b-1] = new TH1F(hold,"BackSub",200,0.84,1.00);
			x++;				

		}

	}

/* ***************************MAIN PROGRAM*************************** */
/* *******************************FLUX FILE************************** */
TFile *flFile = new TFile("Flux.root"); 
	TTree *pf = (TTree*)flFile->Get("T");
    pf->SetBranchAddress("run",&runf);
	pf->SetBranchAddress("flux",&flux);

	nentries = (Int_t)pf->GetEntries();

	for (j=0;j<nentries;j++){ 

		pf->GetEntry(j);

		Flux->Fill(runf,flux);
	}
	flFile->Close();


/* ***********************GENERATED EVENTS FILE********************** */
	TFile *genFile = new TFile("Events.root"); 
	TTree *pg = (TTree*)genFile->Get("T");
			
	pg->SetBranchAddress("p1px",&p1gpx);
	pg->SetBranchAddress("p1py",&p1gpy);
	pg->SetBranchAddress("p1pz",&p1gpz);
	pg->SetBranchAddress("p1pe",&p1gpe);
	pg->SetBranchAddress("p2px",&p2gpx);
	pg->SetBranchAddress("p2py",&p2gpy);
	pg->SetBranchAddress("p2pz",&p2gpz);
	pg->SetBranchAddress("p2pE",&p2gpE);
	pg->SetBranchAddress("pbpx",&pbgpx);
	pg->SetBranchAddress("pbpy",&pbgpy);
	pg->SetBranchAddress("pbpz",&pbgpz);
	pg->SetBranchAddress("pbpE",&pbgpE);
	pg->SetBranchAddress("pIpx",&pIgpx);
	pg->SetBranchAddress("pIpy",&pIgpy);
	pg->SetBranchAddress("pIpz",&pIgpz);
	pg->SetBranchAddress("pIpE",&pIgpE);

	nentries = (Int_t)pg->GetEntries();
			
	/* Loops over events. */
	for (j=0;j<nentries;j++){ 

		pg->GetEntry(j);

		boostg = -pIgpE/(pIgpE + xmpr);

		// Sets the 4-Vectors of the particles. 
		P4p1g.SetPxPyPzE(p1gpx,p1gpy,p1gpz,sqrt(p1gpx*p1gpx+p1gpy*p1gpy+p1gpz*p1gpz+xmpr*xmpr));
		P4p2g.SetPxPyPzE(p2gpx,p2gpy,p2gpz,sqrt(p2gpx*p2gpx+p2gpy*p2gpy+p2gpz*p2gpz+xmpr*xmpr));
		P4pbg.SetPxPyPzE(pbgpx,pbgpy,pbgpz,sqrt(pbgpx*pbgpx+pbgpy*pbgpy+pbgpz*pbgpz+xmpr*xmpr));

		P4p1g.Boost(0,0,boostg);
		P4p2g.Boost(0,0,boostg);
		P4pbg.Boost(0,0,boostg);

		p1bgcos = P4p1g.Vect().CosTheta();
		p2bgcos = P4p2g.Vect().CosTheta();
		pbbgcos = P4pbg.Vect().CosTheta();
	
		//p1bgphi = P4p1g.Vect().Phi();
		//p2bgphi = P4p2g.Vect().Phi();
		//pbbgphi = P4pbg.Vect().Phi();

		if(pIgpE > 4.5 && pIgpE < 5.5){

			//p1gcos->Fill(p1bgcos);
			//p2gcos->Fill(p2bgcos);
			//pbgcos->Fill(pbbgcos);

			//p1gphi->Fill(p1bgphi);
			//p2gphi->Fill(p2bgphi);
			//pbgphi->Fill(pbbgphi);

			p1gvsegam->Fill(pIgpE,p1bgcos);
			p2gvsegam->Fill(pIgpE,p2bgcos);
			pbgvsegam->Fill(pIgpE,pbbgcos);

		}

	}
	genFile->Close();

									   
/* ************************MONTE-CARLO FILE************************** */

	TFile *mcFile = new TFile("mc_ppmpbar_mass.root"); 
	TTree *pMC = (TTree*)mcFile->Get("ppbar");
			
	pMC->SetBranchAddress("ebeamcor",&MCegam);                                           
	pMC->SetBranchAddress("pxp1pcor", &MCpxp1);
	pMC->SetBranchAddress("pyp1pcor", &MCpyp1);
	pMC->SetBranchAddress("pzp1pcor", &MCpzp1);
	pMC->SetBranchAddress("pxp2pcor", &MCpxp2);
	pMC->SetBranchAddress("pyp2pcor", &MCpyp2);
	pMC->SetBranchAddress("pzp2pcor", &MCpzp2);
	pMC->SetBranchAddress("pxpbarpcor", &MCpxpb);
	pMC->SetBranchAddress("pypbarpcor", &MCpypb);
	pMC->SetBranchAddress("pzpbarpcor", &MCpzpb);	

	//pMC->SetBranchAddress("secp1", &MCsecp1);	
	//pMC->SetBranchAddress("secp2", &MCsecp2);	
	//pMC->SetBranchAddress("secpbar", &MCsecpb);

	pMC->SetBranchAddress("pbeta", &MCpbeta);
	pMC->SetBranchAddress("pBeta", &MCpBeta);
	pMC->SetBranchAddress("p2beta", &MCp2beta);
	pMC->SetBranchAddress("p2Beta", &MCp2Beta);
	pMC->SetBranchAddress("pbarbeta", &MCpbbeta);
	pMC->SetBranchAddress("pbarBeta", &MCpbBeta);		

	pMC->SetBranchAddress("vxp1", &MCvxp1);
	pMC->SetBranchAddress("vyp1", &MCvyp1);
	pMC->SetBranchAddress("vzp1", &MCvzp1);
	pMC->SetBranchAddress("vxp2", &MCvxp2);
	pMC->SetBranchAddress("vyp2", &MCvyp2);
	pMC->SetBranchAddress("vzp2", &MCvzp2);
	pMC->SetBranchAddress("vxpbar", &MCvxpb); 
	pMC->SetBranchAddress("vypbar", &MCvypb);
	pMC->SetBranchAddress("vzpbar", &MCvzpb);

	nentries = (Int_t)pMC->GetEntries();
			
	/* Loops over events. */
	for (j=0;j<nentries;j++){ 

		pMC->GetEntry(j);

		MCboost = -MCegam/(MCegam + xmpr);
		
		// Set the target and photon 4-Vectors.
		P4MCtarget.SetPxPyPzE(0.0,0.0,0.0,xmpr);		
		P4MCpho.SetPxPyPzE(0.0,0.0,MCegam,MCegam);
		P4MCinc =  P4MCtarget + P4MCpho;
       
		// Sets the 4-Vectors of the particles. 
		P4MCpb.SetPxPyPzE(MCpxpb,MCpypb,MCpzpb,sqrt(MCpxpb*MCpxpb+MCpypb*MCpypb+MCpzpb*MCpzpb+xmpr*xmpr));
		P4MCp1.SetPxPyPzE(MCpxp1,MCpyp1,MCpzp1,sqrt(MCpxp1*MCpxp1+MCpyp1*MCpyp1+MCpzp1*MCpzp1+xmpr*xmpr));
		P4MCp2.SetPxPyPzE(MCpxp2,MCpyp2,MCpzp2,sqrt(MCpxp2*MCpxp2+MCpyp2*MCpyp2+MCpzp2*MCpzp2+xmpr*xmpr));
	
		// Vertices of the event
		V3MCp1.SetXYZ(MCvxp1,MCvyp1,MCvzp1); 
		V3MCp2.SetXYZ(MCvxp2,MCvyp2,MCvzp2);
		V3MCpb.SetXYZ(MCvxpb,MCvypb,MCvzpb);
		MCxAve = (V3MCp1.X() + V3MCp2.X() + V3MCpb.X())/3;
		MCyAve = (V3MCp1.Y() + V3MCp2.Y() + V3MCpb.Y())/3;
		MCzAve = (V3MCp1.Z() + V3MCp2.Z() + V3MCpb.Z())/3;

		// Calculate various 4-vectors.
		P4MCall = P4MCp1 + P4MCp2 + P4MCpb;
		P4MCprots = P4MCp1 + P4MCp2;
		P4MCmiss = P4MCinc - P4MCall;
		P4MCprMM = P4MCinc - P4MCprots;

		// Calculate Momenta
		MCpp1 = sqrt(MCpxp1*MCpxp1+MCpyp1*MCpyp1+MCpzp1*MCpzp1);
		MCpp2 = sqrt(MCpxp2*MCpxp2+MCpyp2*MCpyp2+MCpzp2*MCpzp2);
		MCppb = sqrt(MCpxpb*MCpxpb+MCpypb*MCpypb+MCpzpb*MCpzpb);

		// Get masses.
		MCxmm  = P4MCmiss.M();
		MCxmm2 = P4MCprMM.M();	

		MCp1delbeta = (MCpBeta-MCpbeta);
		MCp2delbeta = (MCp2Beta-MCp2beta);
		MCpbdelbeta = (MCpbBeta-MCpbbeta);

		//Boost the vectors
		P4MCp1.Boost(0,0,MCboost);
		P4MCp2.Boost(0,0,MCboost);
		P4MCpb.Boost(0,0,MCboost);

		MCp1bcos = P4MCp1.Vect().CosTheta();
		MCp2bcos = P4MCp2.Vect().CosTheta();
		MCpbbcos = P4MCpb.Vect().CosTheta();
	
		//MCp1bphi = P4MCp1.Vect().Phi();
		//MCp2bphi = P4MCp2.Vect().Phi();
		//MCpbbphi = P4MCpb.Vect().Phi();


		Float_t Vradius = sqrt(MCxAve*MCxAve + MCyAve*MCyAve);

		HMCegam->Fill(MCegam);

		//Cut on z-vertex inside the target
		if(MCzAve > -108. && MCzAve < -72. && Vradius < 2.0) {
			HMCvxvy->Fill(MCxAve, MCyAve);

         		if(MCegam > 4.5 && MCegam < 5.5 && /*P4MCprMM.M2() > 0.7 && P4MCprMM.M2() < 1.1 &&*/
         		   MCp1delbeta < 0.04 && MCp1delbeta > -0.04 && MCp2delbeta < 0.04 && MCp2delbeta > -0.04 && MCpbdelbeta < 0.04 && MCpbdelbeta > -0.04 && P4MCmiss.M2() > -0.02 && P4MCmiss.M2()<0.02) {
			
				HMCMMb->Fill(MCxmm);
			
				HMCvz->Fill(MCzAve);

				HMCMM->Fill(MCxmm2);



				HMCMMsq->Fill(P4MCprMM.M2());
				HMCMMsqall->Fill(P4MCmiss.M2());

				HMCp1beta->Fill(MCpp1,(MCpBeta-MCpbeta));
				HMCp2beta->Fill(MCpp2,(MCp2Beta-MCp2beta));
				HMCpbbeta->Fill(MCppb,(MCpbBeta-MCpbbeta));

				//MCp1cos->Fill(MCp1bcos);
				//MCp2cos->Fill(MCp2bcos);
				//MCpbcos->Fill(MCpbbcos);

				//MCp1phi->Fill(MCp1bphi);
				//MCp2phi->Fill(MCp2bphi);
				//MCpbphi->Fill(MCpbbphi);

				MCp1vsegam->Fill(MCegam,MCp1bcos);
				MCp2vsegam->Fill(MCegam,MCp2bcos);
				MCpbvsegam->Fill(MCegam,MCpbbcos);

				P1Accept->Fill(MCegam,MCp1bcos);
				P2Accept->Fill(MCegam,MCp2bcos);
				PbAccept->Fill(MCegam,MCpbbcos);															

          		}//egam

		}//z-vertex

	}
	mcFile->Close();

/* ****************************DATA FILES**************************** */
	/* Loops over data files. */
	for(int n_arg = optind; n_arg < __argc; n_arg++){
   
		TFile inFile(__argv[n_arg]);
		
		/* Opens each input file. */        
		if(TTree *p1 = (TTree*)inFile.Get("ppbar")){
			
			p1->SetBranchAddress("ebeamcor",&egam); // Photon energy.
			p1->SetBranchAddress("pxp1pcor", &pxp1);
			p1->SetBranchAddress("pyp1pcor", &pyp1);
			p1->SetBranchAddress("pzp1pcor", &pzp1);
			p1->SetBranchAddress("pxp2pcor", &pxp2);
			p1->SetBranchAddress("pyp2pcor", &pyp2);
			p1->SetBranchAddress("pzp2pcor", &pzp2);
			p1->SetBranchAddress("pxpbarpcor", &pxpb);
			p1->SetBranchAddress("pypbarpcor", &pypb);
			p1->SetBranchAddress("pzpbarpcor", &pzpb);	

			//p1->SetBranchAddress("secp1", &secp1);
			//p1->SetBranchAddress("secp2", &secp2);
			//p1->SetBranchAddress("secpbar", secpb);

			p1->SetBranchAddress("run", &run);

			p1->SetBranchAddress("pbeta", &pbeta);
			p1->SetBranchAddress("pBeta", &pBeta);
			p1->SetBranchAddress("p2beta", &p2beta);
			p1->SetBranchAddress("p2Beta", &p2Beta);
			p1->SetBranchAddress("pbarbeta", &pbbeta);
			p1->SetBranchAddress("pbarBeta", &pbBeta);

			p1->SetBranchAddress("vxp1", &vxp1);
			p1->SetBranchAddress("vyp1", &vyp1);
			p1->SetBranchAddress("vzp1", &vzp1);
			p1->SetBranchAddress("vxp2", &vxp2);
			p1->SetBranchAddress("vyp2", &vyp2);
			p1->SetBranchAddress("vzp2", &vzp2);
			p1->SetBranchAddress("vxpbar", &vxpb); 
			p1->SetBranchAddress("vypbar", &vypb);
			p1->SetBranchAddress("vzpbar", &vzpb);

			nentries = (Int_t)p1->GetEntries();
			
			/* Loops over events. */
			for (j=0;j<nentries;j++){ 

				p1->GetEntry(j);

				boost = -egam/(egam + xmpr);
		
				// Set the target and photon 4-Vectors.
				P4target.SetPxPyPzE(0.0,0.0,0.0,xmpr);		
				P4pho.SetPxPyPzE(0.0,0.0,egam,egam);
				P4inc =  P4target + P4pho;
       
				// Sets the 4-Vectors of the particles. 
				P4pb.SetPxPyPzE(pxpb,pypb,pzpb,sqrt(pxpb*pxpb+pypb*pypb+pzpb*pzpb+xmpr*xmpr));
				P4p1.SetPxPyPzE(pxp1,pyp1,pzp1,sqrt(pxp1*pxp1+pyp1*pyp1+pzp1*pzp1+xmpr*xmpr));
				P4p2.SetPxPyPzE(pxp2,pyp2,pzp2,sqrt(pxp2*pxp2+pyp2*pyp2+pzp2*pzp2+xmpr*xmpr));
	
				// Vertices of the event
				V3p1.SetXYZ(vxp1,vyp1,vzp1); 
				V3p2.SetXYZ(vxp2,vyp2,vzp2);
				V3pb.SetXYZ(vxpb,vypb,vzpb);
				xAve = (V3p1.X() + V3p2.X() + V3pb.X())/3;
				yAve = (V3p1.Y() + V3p2.Y() + V3pb.Y())/3;
				zAve = (V3p1.Z() + V3p2.Z() + V3pb.Z())/3;

				// Calculate various 4-vectors.
				P4all = P4p1 + P4p2 + P4pb;
				P4prots = P4p1 + P4p2;
				P4miss = P4inc - P4all;
				P4prMM = P4inc - P4prots;

				// Calculate Momenta
				pp1 = sqrt(pxp1*pxp1+pyp1*pyp1+pzp1*pzp1);
				pp2 = sqrt(pxp2*pxp2+pyp2*pyp2+pzp2*pzp2);
				ppb = sqrt(pxpb*pxpb+pypb*pypb+pzpb*pzpb);

				// Get masses.
				xmm  = P4miss.M();
				xmm2 = P4prMM.M();	

				p1delbeta = (pBeta-pbeta);
				p2delbeta = (p2Beta-p2beta);
				pbdelbeta = (pbBeta-pbbeta);

				//Boost the vectors
				P4p1.Boost(0,0,boost);
				P4p2.Boost(0,0,boost);
				P4pb.Boost(0,0,boost);

				p1bcos = P4p1.Vect().CosTheta();
				p2bcos = P4p2.Vect().CosTheta();
				pbbcos = P4pb.Vect().CosTheta();
	
				//p1bphi = P4p1.Vect().Phi();
				//p2bphi = P4p2.Vect().Phi();
				//pbbphi = P4pb.Vect().Phi();

				//Fill Histograms
				Hegam->Fill(egam);
			
				HMMb->Fill(xmm);
			
				Hvz->Fill(zAve);
				Float_t Vradius = sqrt( xAve*xAve + yAve*yAve );

				if(egam > 4.5 && egam < 5.5) {cut1->Fill(egam);}
				if(zAve > -108. && zAve < -72. && Vradius < 2.0) {cut2->Fill(egam);}
				if(p1delbeta < 0.04 && p1delbeta > -0.04 && p2delbeta < 0.04 && p2delbeta > -0.04 && pbdelbeta < 0.04 && pbdelbeta > -0.045) {cut3->Fill(egam);}
				if(P4miss.M2() > -0.02 && P4miss.M2()<0.02) {cut4->Fill(egam);}
				if(egam > 4.5 && egam < 5.5){
					cut1s->Fill(egam);
					if(zAve > -108. && zAve < -72. && Vradius < 2.0){
						cut2s->Fill(egam);
						if(p1delbeta < 0.04 && p1delbeta > -0.04 && p2delbeta < 0.04 && p2delbeta > -0.04 && pbdelbeta < 0.04 && pbdelbeta > -0.045){
							cut3s->Fill(egam);
							if(P4miss.M2() > -0.02 && P4miss.M2()<0.02){
								cut4s->Fill(egam);
							}
						}
					}
				}


				

				//Cut on z-vertex inside the target
				if(zAve > -108. && zAve < -72. && Vradius < 2.0) {
					Hvxvy->Fill(xAve, yAve);

         				if(egam > 4.5 && egam < 5.5 && /*P4prMM.M2() > 0.7 && P4prMM.M2() < 1.1 &&*/
         				   p1delbeta < 0.04 && p1delbeta > -0.04 && p2delbeta < 0.04 && p2delbeta > -0.04 && pbdelbeta < 0.04 && pbdelbeta > -0.04 && P4miss.M2() > -0.02 && P4miss.M2()<0.02) {

						HMM->Fill(xmm);
						HMMsq->Fill(P4prMM.M2());
						HMMsqall->Fill(P4miss.M2());

						Hp1beta->Fill(pp1,(pBeta-pbeta));
						Hp2beta->Fill(pp2,(p2Beta-p2beta));
						Hpbbeta->Fill(ppb,(pbBeta-pbbeta));

						Consist->Fill(run);
						Counts->Fill(run);

						//p1cos->Fill(p1bcos);
						//p2cos->Fill(p2bcos);
						//pbcos->Fill(pbbcos);

						//p1phi->Fill(p1bphi);
						//p2phi->Fill(p2bphi);
						//pbphi->Fill(pbbphi);

						for(int a = 0; a < 10; a++){

							for(int b = 0; b < 20; b++){

								if(egam > (4.5 + ((double)a/10)) && egam <= (4.5 + (((double)a+1)/10)) && pbbcos > ((0.1 * (double)b)-1) && pbbcos <= ((0.1 * ((double)b+1)))-1){
									Place[a][b] -> Fill(xmm2);
								}

							}

						}

						p1vsegam->Fill(egam,p1bcos);
						p2vsegam->Fill(egam,p2bcos);
						pbvsegam->Fill(egam,pbbcos);

						P1Cross->Fill(egam,p1bcos);
						P2Cross->Fill(egam,p2bcos);
						PbCross->Fill(egam,pbbcos);	

          				}//egam

				}//z-vertex
	
			}// entries
		}//if TTree *nbar = (TTree*)
	
	}// files

	for(int a = 0; a < 10; a++){

		for(int b = 0; b < 20; b++){

		Bac = Place[a][b] -> ShowBackground(50);
		Place[a][b] -> Add(Bac,-1);
		Place[a][b] -> SetMinimum(0);

		yield = Place[a][b] -> Integral(75,160);
		PbYield -> SetBinContent(a+1,b+1,yield);
		PbCro -> SetBinContent(a+1,b+1,yield);
		}

	}


	P1Accept->Divide(p1gvsegam);
	P2Accept->Divide(p2gvsegam);
	PbAccept->Divide(pbgvsegam);

	P1Cross->Divide(P1Accept);
	P2Cross->Divide(P2Accept);
	PbCross->Divide(PbAccept);

	PbCro->Divide(PbAccept);

	
	Consist->Divide(Flux);

	for(int i = 1; i < 917; i++){

		cts = Counts->GetBinContent(i);
		fl = Flux->GetBinContent(i);

		ER = (cts / (fl+1)) * sqrt((1/(cts+1))+(1/(fl+1)));

		Consist->SetBinError(i,ER);
		Counts->SetBinError(i,sqrt(cts));
		Flux->SetBinError(i,sqrt(fl));

	}


	/* Trigger 1 */
	if(TRIGGER == 1){
		e1 = lum * 5.05176143094e+11;
		e2 = lum * 4.73176495078e+11;
		e3 = lum * 4.17117938873e+11;
		e4 = lum * 4.62978840232e+11;
		e5 = lum * 4.70693544999e+11;
		e6 = lum * 4.57069088076e+11;
		e7 = lum * 4.58548090432e+11;
		e8 = lum * 4.34902850828e+11;
		e9 = lum * 4.46131165432e+11;
		e0 = lum * 2.29742689709e+11;
	}

	/* Trigger 2 */
	if(TRIGGER == 2){
		e1 = lum * 1.71787301675e+12;
		e2 = lum * 1.60822833343e+12;
		e3 = lum * 1.4085899412e+12;
		e4 = lum * 1.56631796846e+12;
		e5 = lum * 1.59688388124e+12;
		e6 = lum * 1.54712464412e+12;
		e7 = lum * 1.56201234707e+12;
		e8 = lum * 1.28316654195e+12;
		e9 = lum * 1.52079767527e+12;
		e0 = lum * 7.79553704428e+11;
	}
	// Total 
	if(TRIGGER == 3){
		e1 = lum * 2.22304915985e+12;
		e2 = lum * 2.08140482851e+12;
		e3 = lum * 1.82570788007e+12;
		e4 = lum * 2.0292968087e+12;
		e5 = lum * 2.06757742624e+12;
		e6 = lum * 2.0041937322e+12;
		e7 = lum * 2.0205604375e+12;
		e8 = lum * 1.71806939278e+12;
		e9 = lum * 1.96692884071e+12;
		e0 = lum * 1.00929639414e+12;
	}

	for(int a = 0; a < 10; a++){
		for(int b = 0; b < 20; b++){

			pb = PbAccept->GetBinContent(a,b);
			Y = PbCross->GetBinContent(a,b);

			ER1 = MCpbvsegam->GetBinContent(a,b);
			ER2 = pbgvsegam->GetBinContent(a,b);
			ER3 = pbvsegam->GetBinContent(a,b);

			if(ER1 != 0 && ER2 != 0 && ER3 !=  0){

				pbe = (ER1/ER2)*sqrt(1/ER1+1/ER2);
				Ye = ER3/(ER1/ER2)*sqrt(1/ER3+(pbe/pb)*(pbe/pb));

			}
			else{

				pbe = 0;
				Ye = 0;

			}

			if(a==1){
				PbYield1->SetBinContent(b,ER3);
				PbYield1->SetBinError(b,sqrt(ER3));
				PbAcc1->SetBinContent(b, pb);
				PbAcc1->SetBinError(b, pbe);
				PbCross1->SetBinContent(b,Y/e1);
				PbCross1->SetBinError(b,Ye/e1);

			}
			if(a==2){
				PbYield2->SetBinContent(b,ER3);
				PbYield2->SetBinError(b,sqrt(ER3));
				PbAcc2->SetBinContent(b, pb);
				PbAcc2->SetBinError(b, pbe);
				PbCross2->SetBinContent(b,Y/e2);
				PbCross2->SetBinError(b,Ye/e2);
			}
			if(a==3){
				PbYield3->SetBinContent(b,ER3);
				PbYield3->SetBinError(b,sqrt(ER3));
				PbAcc3->SetBinContent(b, pb);
				PbAcc3->SetBinError(b, pbe);
				PbCross3->SetBinContent(b,Y/e3);
				PbCross3->SetBinError(b,Ye/e3);
			}
			if(a==4){
				PbYield4->SetBinContent(b,ER3);
				PbYield4->SetBinError(b,sqrt(ER3));
				PbAcc4->SetBinContent(b, pb);
				PbAcc4->SetBinError(b, pbe);
				PbCross4->SetBinContent(b,Y/e4);
				PbCross4->SetBinError(b,Ye/e4);
			}
			if(a==5){
				PbYield5->SetBinContent(b,ER3);
				PbYield5->SetBinError(b,sqrt(ER3));
				PbAcc5->SetBinContent(b, pb);
				PbAcc5->SetBinError(b, pbe);
				PbCross5->SetBinContent(b,Y/e5);
				PbCross5->SetBinError(b,Ye/e5);
			}
			if(a==6){
				PbYield6->SetBinContent(b,ER3);
				PbYield6->SetBinError(b,sqrt(ER3));
				PbAcc6->SetBinContent(b, pb);
				PbAcc6->SetBinError(b, pbe);
				PbCross6->SetBinContent(b,Y/e6);
				PbCross6->SetBinError(b,Ye/e6);
			}
			if(a==7){
				PbYield7->SetBinContent(b,ER3);
				PbYield7->SetBinError(b,sqrt(ER3));
				PbAcc7->SetBinContent(b, pb);
				PbAcc7->SetBinError(b, pbe);
				PbCross7->SetBinContent(b,Y/e7);
				PbCross7->SetBinError(b,Ye/e7);
			}
			if(a==8){
				PbYield8->SetBinContent(b,ER3);
				PbYield8->SetBinError(b,sqrt(ER3));
				PbAcc8->SetBinContent(b, pb);
				PbAcc8->SetBinError(b, pbe);
				PbCross8->SetBinContent(b,Y/e8);
				PbCross8->SetBinError(b,Ye/e8);
			}
			if(a==9){
				PbYield9->SetBinContent(b,ER3);
				PbYield9->SetBinError(b,sqrt(ER3));
				PbAcc9->SetBinContent(b, pb);
				PbAcc9->SetBinError(b, pbe);
				PbCross9->SetBinContent(b,Y/e9);
				PbCross9->SetBinError(b,Ye/e9);
			}
			if(a==10){
				PbYield0->SetBinContent(b,ER3);
				PbYield0->SetBinError(b,sqrt(ER3));
				PbAcc0->SetBinContent(b, pb);
				PbAcc0->SetBinError(b, pbe);
				PbCross0->SetBinContent(b,Y/e0);
				PbCross0->SetBinError(b,Ye/e0);
			}
																					
		}
	}

	
	PbCrossTot->SetBinContent(0, PbCross1->Integral(0,20));
	PbCrossTot->SetBinContent(1, PbCross2->Integral(0,20));
	PbCrossTot->SetBinContent(2, PbCross3->Integral(0,20));
	PbCrossTot->SetBinContent(3, PbCross4->Integral(0,20));
	PbCrossTot->SetBinContent(4, PbCross5->Integral(0,20));
	PbCrossTot->SetBinContent(5, PbCross6->Integral(0,20));
	PbCrossTot->SetBinContent(6, PbCross7->Integral(0,20));
	PbCrossTot->SetBinContent(7, PbCross8->Integral(0,20));
	PbCrossTot->SetBinContent(8, PbCross9->Integral(0,20));
	PbCrossTot->SetBinContent(9, PbCross0->Integral(0,20));


	PbCross1->SetMinimum(0);
	PbCross2->SetMinimum(0);
	PbCross3->SetMinimum(0);
	PbCross4->SetMinimum(0);
	PbCross5->SetMinimum(0);
	PbCross6->SetMinimum(0);
	PbCross7->SetMinimum(0);
	PbCross8->SetMinimum(0);
	PbCross9->SetMinimum(0);
	PbCross0->SetMinimum(0);

	Consist->SetMinimum(0);
	Consist->SetMaximum(6e-9);

	PbCross1->SetBinContent(1,0);
	PbCross2->SetBinContent(1,0);
	PbCross3->SetBinContent(1,0);
	PbCross4->SetBinContent(1,0);
	PbCross5->SetBinContent(1,0);
	PbCross6->SetBinContent(1,0);
	PbCross7->SetBinContent(1,0);
	PbCross8->SetBinContent(1,0);
	PbCross9->SetBinContent(1,0);
	PbCross0->SetBinContent(1,0);
	PbCross1->SetBinError(1,0);
	PbCross2->SetBinError(1,0);
	PbCross3->SetBinError(1,0);
	PbCross4->SetBinError(1,0);
	PbCross5->SetBinError(1,0);
	PbCross6->SetBinError(1,0);
	PbCross7->SetBinError(1,0);
	PbCross8->SetBinError(1,0);
	PbCross9->SetBinError(1,0);
	PbCross0->SetBinError(1,0);

outFile.Write(); // Write to the output file.
outFile.Close(); // Close the output file.

}





