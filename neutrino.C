// example code to run Bulk Graviton->ZZ->ZlepZhad selections on electron-channel

#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <TString.h>
#include <map>
#include <TH2F.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TLegend.h>
#include "untuplizer.h"
#include <TClonesArray.h>
#include <TLorentzVector.h>
#include "setNCUStyle.C"


using namespace std;
void neutrino(std::string inputFile) {

  //get TTree from file ...
  TreeReader data(inputFile.data());

  TH2F *h_leadJet_neuP4 = new TH2F("","",100,0,3600,100,0,860);
  TH2F *h_leadJet_nNeu  = new TH2F("","",36,0,3600,10,0,10);
  TH2F *h_leadPR_neuP4  = new TH2F("","",100,0,3600,100,0,860);
  TH2F *h_leadPR_nNeu   = new TH2F("","",100,0,3600,100,0,860);
  TH2F *h_sublJet_neuP4 = new TH2F("","",100,0,3600,100,0,860);
  TH2F *h_sublJet_nNeu  = new TH2F("","",100,0,3600,100,0,860);
  TH2F *h_sublPR_neuP4  = new TH2F("","",100,0,3600,100,0,860);
  TH2F *h_sublPR_nNeu   = new TH2F("","",100,0,3600,100,0,860);

  //Event loop
  for(Long64_t jEntry=0; jEntry<data.GetEntriesFast() ;jEntry++){

    if (jEntry % 50000 == 0)
      fprintf(stderr, "Processing event %lli of %lli\n", jEntry + 1, data.GetEntriesFast());

    data.GetEntry(jEntry);
    
    //0. has a good vertex
    Int_t nVtx        = data.GetInt("nVtx");

    int nFATJet         = data.GetInt("FATnJet");
    const int nJets=nFATJet;
    TClonesArray* genParP4    = (TClonesArray*) data.GetPtrTObject("genParP4");
    TClonesArray* fatjetP4    = (TClonesArray*) data.GetPtrTObject("FATjetP4");
    Float_t*  fatjetCISVV2    = data.GetPtrFloat("FATjetCISVV2");
    Float_t*  fatjetSDmass    = data.GetPtrFloat("FATjetSDmass");
    Float_t*  fatjetPRmass    = data.GetPtrFloat("FATjetPRmass");
    Int_t*    nSubSoftDropJet = data.GetPtrInt("FATnSubSDJet");
    Int_t     nGenPar         = data.GetInt("nGenPar");
    Int_t*    genParId        = data.GetPtrInt("genParId");
    Int_t*    genParSt        = data.GetPtrInt("genParSt");
    Int_t*    genMomParId     = data.GetPtrInt("genMomParId");
    Int_t*    genDa1          = data.GetPtrInt("genDa1");
    Int_t*    genDa2          = data.GetPtrInt("genDa2");
    vector<float>   *subjetSDCSV =  data.GetPtrVectorFloat("FATsubjetSDCSV", nFATJet);
    vector<float>   *subjetSDPx  =  data.GetPtrVectorFloat("FATsubjetSDPx", nFATJet);
    vector<float>   *subjetSDPy  =  data.GetPtrVectorFloat("FATsubjetSDPy", nFATJet);
    vector<float>   *subjetSDPz  =  data.GetPtrVectorFloat("FATsubjetSDPz", nFATJet);
    vector<bool>    &passFatJetLooseID = *((vector<bool>*) data.GetPtr("FATjetPassIDLoose"));
    
    vector<int> fatty;
    for(int ij=0; ij<nJets; ij++)
      {
    	
     	TLorentzVector* thisJet = (TLorentzVector*)fatjetP4->At(ij);
    	if(thisJet->Pt()<30)continue;
    	if(fabs(thisJet->Eta())>2.5)continue;
	if(fatjetSDmass[ij]<95 || fatjetSDmass[ij]>145)continue;
	if(!passFatJetLooseID[ij])continue;

    	if(fatjetCISVV2[ij] < 0.605)continue;
    	if(fatjetCISVV2[ij] > 1)continue;

	fatty.push_back(ij);
      }
    

    vector<int> pruned;
    for(int ij=0; ij<nJets; ij++)
      {

        TLorentzVector* thisJet = (TLorentzVector*)fatjetP4->At(ij);
        if(thisJet->Pt()<300)continue;
        if(fabs(thisJet->Eta())>2.5)continue;
	if(!passFatJetLooseID[ij])continue;

        if(fatjetCISVV2[ij] < 0.605)continue;
        if(fatjetCISVV2[ij] > 1)continue;

        pruned.push_back(ij);
      }
    
    if(fatty.size()<2)continue;
    if(pruned.size()<2)continue;

    Int_t lead = fatty[0]; //first leading jet
    Int_t subl = fatty[1]; //second leading jet
    Int_t leadPR = pruned[0];
    Int_t sublPR = pruned[1];

    TLorentzVector *leadjet = (TLorentzVector*)fatjetP4->At(lead);
    TLorentzVector *subljet = (TLorentzVector*)fatjetP4->At(subl);
    TLorentzVector *PRleadjet = (TLorentzVector*)fatjetP4->At(leadPR);
    TLorentzVector *PRsubljet = (TLorentzVector*)fatjetP4->At(sublPR);
    
    int n_neutrinos=0;
    TLorentzVector neutrinos_p4;
    
    for(int ig=0; ig < nGenPar; ig++){
      
      int pid =abs(genParId[ig]);
      int mom_pid = abs(genMomParId[ig]);
      if(genParSt[ig]!=1)continue;
      
      if(pid!=12 && pid!=14 && pid!=16)continue;
    
      // look for daughters of charm and b hadrons
      // bhadrons can decay to l nu + charm hadrons 
      // charm hadrons could further decay to l nu + other hadrons

      int da1_pid = abs(genDa1[ig]);
      int da2_pid = abs(genDa2[ig]);
       
      if(da1_pid == pid && (da2_pid<500 && da2_pid>400)) mom_pid < 500; 
      if(da1_pid == pid && (da2_pid<600 && da2_pid>99 )) mom_pid > 500;
      
      if(mom_pid < 400 || mom_pid > 600)continue;
      
      TLorentzVector* thisGen = (TLorentzVector*)genParP4->At(ig);
      if(thisGen->DeltaR(*leadjet)>0.8)continue;
      //if(thisGen->DeltaR(*subljet)>0.8)continue;
      //if(thisGen->DeltaR(*PRleadjet)>0.8)continue;
      //if(thisGen->DeltaR(*PRsubljet)>0.8)continue;
      n_neutrinos++;
      neutrinos_p4 += *thisGen; 
      
    }
    
    h_leadJet_neuP4->Fill(leadjet->Pt(),neutrinos_p4.Pt()); 
    h_leadJet_nNeu->Fill(leadjet->Pt(),n_neutrinos);    
    //h_leadPR_neuP4->Fill(fatjetPRmass[leadPR],neutrinosp4.Pt());
    //h_leadPR_nNeu->Fill(fatjetPRmass[leadPR],n_neutrinos);
    //h_sublJet_neuP4->Fill(subljet->Pt(),neutrinos_p4.Pt());
    //h_sublJet_nNeu->Fill(subljet->Pt(),n_neutrinos); 
    ///h_sublPR_neuP4->Fill(fatjetPRmass[sublPR],neutrinosp4.Pt());
    //h_sublPR_nNeu->Fill(fatjetPRmass[sublPR],n_neutrinos);
  

  } // end of loop over entries
  
  setNCUStyle();
  TCanvas* c1 = new TCanvas("c1","c1",0,0,800,600);
  c1->cd();
  h_leadJet_neuP4->Draw("colz");
  TCanvas* c2 = new TCanvas("c2","c2",0,0,800,600);
  c2->cd();
  h_leadJet_nNeu->Draw("colz");
}
