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
#include <TLatex.h>
#include <TAxis.h>
#include <TPad.h>
#include <TStyle.h>
#include "untuplizer.h"
#include <TClonesArray.h>
#include <TLorentzVector.h>
#include "setNCUStyle.C"
//#include "additional_style.C"

using namespace std;
void neutrino(std::string inputFile) {

  //get TTree from file ...
  TreeReader data(inputFile.data());

  TH2F *h_leadJet_neuP4 = new TH2F("","",75,0,3600,45,0,900);
  TH2F *h_leadJet_nNeu  = new TH2F("","",36,0,3600,17,0,17);
  TH2F *h_sublJet_neuP4 = new TH2F("","",70,0,3000,45,0,900);
  TH2F *h_sublJet_nNeu  = new TH2F("","",30,0,3000,17,0,17);
  TH2F *h_leadPR_neuP4  = new TH2F("","",50,40,200,45,0,900);//
  TH2F *h_leadPR_nNeu   = new TH2F("","",25,0,500,17,0,17);//
  TH2F *h_sublPR_neuP4  = new TH2F("","",50,40,200,45,0,900);//
  TH2F *h_sublPR_nNeu   = new TH2F("","",25,0,500,17,0,17);//

  //Event loop
  for(Long64_t jEntry=0; jEntry<data.GetEntriesFast() ;jEntry++){

    if (jEntry % 50000 == 0)
      fprintf(stderr, "Processing event %lli of %lli\n", jEntry + 1, data.GetEntriesFast());

    data.GetEntry(jEntry);
    
    int nFATJet         = data.GetInt("FATnJet");
    const int nJets=nFATJet;
    TClonesArray* genParP4    = (TClonesArray*) data.GetPtrTObject("genParP4");
    TClonesArray* fatjetP4    = (TClonesArray*) data.GetPtrTObject("FATjetP4");
    Float_t*  fatjetCISVV2    = data.GetPtrFloat("FATjetCISVV2");
    Float_t*  fatjetPRmass    = data.GetPtrFloat("FATjetPRmass");
    Int_t     nGenPar         = data.GetInt("nGenPar");
    Int_t*    genParId        = data.GetPtrInt("genParId");
    Int_t*    genParSt        = data.GetPtrInt("genParSt");
    Int_t*    genMomParId     = data.GetPtrInt("genMomParId");
    vector<bool>    &passFatJetLooseID = *((vector<bool>*) data.GetPtr("FATjetPassIDLoose"));

    vector<int> fatjet;
    for(int ij=0; ij<nJets; ij++)
      {

        TLorentzVector* thisJet = (TLorentzVector*)fatjetP4->At(ij);
        if(thisJet->Pt()<300)continue;
        if(fabs(thisJet->Eta())>2.5)continue;
	if(!passFatJetLooseID[ij])continue;
	if(fatjetPRmass[ij]<40)continue;
        if(fatjetCISVV2[ij] < 0.605)continue;
        if(fatjetCISVV2[ij] > 1)continue;

        fatjet.push_back(ij);
      }
    

    if(fatjet.size()<2)continue;

    Int_t lead = fatjet[0];
    Int_t subl = fatjet[1];

    TLorentzVector *leadjet = (TLorentzVector*)fatjetP4->At(lead);
    TLorentzVector *subljet = (TLorentzVector*)fatjetP4->At(subl);
    
    int n_neutrinos_lead=0;
    TLorentzVector neutrinos_lead_p4;

    int n_neutrinos_subl=0;
    TLorentzVector neutrinos_subl_p4;    

    for(int ig=0; ig < nGenPar; ig++){
      
      int pid =abs(genParId[ig]);
      int mom_pid = abs(genMomParId[ig]);
      if(genParSt[ig]!=1)continue;
      
      if(pid!=12 && pid!=14 && pid!=16)continue;
    
      // look for daughters of charm and b hadrons
      // bhadrons can decay to l nu + charm hadrons 
      // charm hadrons could further decay to l nu + other hadrons
      
      if(mom_pid < 400 || mom_pid > 600)continue;
      
      TLorentzVector* thisGen = (TLorentzVector*)genParP4->At(ig);
      
      if(thisGen->DeltaR(*leadjet)<0.8) { 
	n_neutrinos_lead++;
	neutrinos_lead_p4 += *thisGen; 
      }
      else if(thisGen->DeltaR(*subljet)<0.8) { 
	n_neutrinos_subl++;
	neutrinos_subl_p4 += *thisGen;
      }
      else continue;
      
    }

    if(n_neutrinos_lead > 0) {
      h_leadJet_neuP4->Fill(leadjet->Pt(),neutrinos_lead_p4.Pt());
      h_leadPR_neuP4->Fill(fatjetPRmass[lead],neutrinos_lead_p4.Pt());
    } 
    h_leadJet_nNeu->Fill(leadjet->Pt(),n_neutrinos_lead);
    h_leadPR_nNeu->Fill(fatjetPRmass[lead],n_neutrinos_lead);

    if(n_neutrinos_subl > 0) {
      h_sublJet_neuP4->Fill(subljet->Pt(),neutrinos_subl_p4.Pt());
      h_sublPR_neuP4->Fill(fatjetPRmass[subl],neutrinos_subl_p4.Pt());
    }
    h_sublJet_nNeu->Fill(subljet->Pt(),n_neutrinos_subl);
    h_sublPR_nNeu->Fill(fatjetPRmass[subl],n_neutrinos_subl);
    

  } // end of loop over entries
  
  setNCUStyle(true);

  TAxis *axis1 = h_leadJet_neuP4->GetYaxis();
  TAxis *axis2 = h_sublJet_neuP4->GetYaxis();
  TAxis *axis3 = h_leadPR_neuP4->GetYaxis();
  TAxis *axis4 = h_sublPR_neuP4->GetYaxis();

  TCanvas* c1 = new TCanvas("c1","c1",0,0,600,600);
  c1->cd();
  h_leadJet_neuP4->Draw("colz");
  h_leadJet_neuP4->SetStats(0);
  h_leadJet_neuP4->SetXTitle(" Leading Jet p_{T} [GeV]");
  h_leadJet_neuP4->SetYTitle("Neutrino p_{T} [GeV]");
  axis1->SetTitleOffset(1.4);

  TCanvas* c2 = new TCanvas("c2","c2",0,0,600,600);
  c2->cd();
  h_leadJet_nNeu->Draw("colz");
  h_leadJet_nNeu->SetStats(0);
  h_leadJet_nNeu->SetXTitle("Leading Jet p_{T} [GeV]");
  h_leadJet_nNeu->SetYTitle("No. of neutrino");

  TCanvas* c3 =  new TCanvas("c3","c3",0,0,600,600);
  c3->cd();
  h_sublJet_neuP4->Draw("colz");
  h_sublJet_neuP4->SetStats(0);
  h_sublJet_neuP4->SetXTitle("Subleading Jet p_{T} [GeV]");
  h_sublJet_neuP4->SetYTitle("Neutrino p_{T} [GeV]");
  axis2->SetTitleOffset(1.4);

  TCanvas* c4 = new TCanvas("c4","c4",0,0,600,600);
  c4->cd();
  h_sublJet_nNeu->Draw("colz");
  h_sublJet_nNeu->SetStats(0);
  h_sublJet_nNeu->SetXTitle("Subleading Jet p_{T} [GeV]");
  h_sublJet_nNeu->SetYTitle("No. of neutrino");

  TCanvas* c5 = new TCanvas("c5","c5",0,0,600,600);
  c5->cd();
  h_leadPR_neuP4->Draw("colz");
  h_leadPR_neuP4->SetStats(0);
  h_leadPR_neuP4->SetXTitle("Leading Jet M_{Pruned} [GeV]");
  h_leadPR_neuP4->SetYTitle("Neutrino p_{T} [GeV]");
  axis3->SetTitleOffset(1.4);

  TCanvas* c6 = new TCanvas("c6","c6",0,0,600,600);
  c6->cd();
  h_leadPR_nNeu->Draw("colz");
  h_leadPR_nNeu->SetStats(0);
  h_leadPR_nNeu->SetXTitle("Leading Jet M_{Pruned} [GeV]");
  h_leadPR_nNeu->SetYTitle("No. of neutrino");

  TCanvas* c7 = new TCanvas("c7","c7",0,0,600,600);                                                                                                         
  c7->cd();                                                                                                                                                 
  h_sublPR_neuP4->Draw("colz");                                                                                                                             
  h_sublPR_neuP4->SetStats(0);                                                                                                                              
  h_sublPR_neuP4->SetXTitle("Subleading Jet M_{Pruned} [GeV]");                                                                                
  h_sublPR_neuP4->SetYTitle("Neutrino p_{T} [GeV]");
  axis4->SetTitleOffset(1.4);                                                        
  
  TCanvas* c8 = new TCanvas("c8","c8",0,0,600,600);                                                                                                         
  c8->cd();                                                                                                                                                 
  h_sublPR_nNeu->Draw("colz");                                                                                                                              
  h_sublPR_nNeu->SetStats(0);                                                                                                                               
  h_sublPR_nNeu->SetXTitle("Subleading Jet M_{Pruned} [GeV]");                                                  
  h_sublPR_nNeu->SetYTitle("No. of neutrino");
    
}
