#include <Riostream.h>
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TH1F.h"
#include "TCanvas.h"

void ReadTree(){
//definizione istogrammi da stampare
  TH1F *hist1 = new TH1F("Z_m", "Z_m", 160, -16.5, 14.5);
  TH1F *hist2 = new TH1F("Meson_m", "Meson_m", 160, -0.15, 0.17);
  
  //Apertura file di input
  TFile hfile("ZMesonGamma_output.root"); 

  //Lettura TDirectory
  TDirectoryFile *dir = (TDirectoryFile*)hfile.Get("ZMesonGamma");
  
  //definizione variabili 
  int genZ_ID;
  float genZ_pT;
  float genZ_eta;
  float genZ_phi;
  float genZ_E;
  float genZ_mass;

  int genGamma_ID;
  float genGamma_pT;
  float genGamma_eta;
  float genGamma_phi;
  float genGamma_E;

  int genMeson_ID;
  float genMeson_pT;
  float genMeson_eta;
  float genMeson_phi;
  float genMeson_E;
  float genMeson_mass;

  int genTrackminus_ID;
  float genTrackminus_pT;
  float genTrackminus_eta;
  float genTrackminus_phi;
  float genTrackminus_E;

  int genTrackplus_ID;
  float genTrackplus_pT;
  float genTrackplus_eta;
  float genTrackplus_phi;
  float genTrackplus_E;

  float m_Z_true   = 91.1876;
  float m_Rho_true = 0.77;
  float m_Phi_true = 1.019;
  
  //Lettura TTree  e branch
  TTree *tree = (TTree*)dir->Get("mytree");
  TBranch *b1=tree->GetBranch("genZ_ID");
  TBranch *b2=tree->GetBranch("genZ_pT");
  TBranch *b3=tree->GetBranch("genZ_eta");
  TBranch *b4=tree->GetBranch("genZ_phi");
  TBranch *b5=tree->GetBranch("genZ_E");
  TBranch *b6=tree->GetBranch("genZ_mass");
  TBranch *b7=tree->GetBranch("genGamma_ID");
  TBranch *b8=tree->GetBranch("genGamma_pT");
  TBranch *b9=tree->GetBranch("genGamma_eta");
  TBranch *b10=tree->GetBranch("genGamma_phi");
  TBranch *b11=tree->GetBranch("genGamma_E");
  TBranch *b12=tree->GetBranch("genMeson_ID");
  TBranch *b13=tree->GetBranch("genMeson_pT");
  TBranch *b14=tree->GetBranch("genMeson_eta");
  TBranch *b15=tree->GetBranch("genMeson_phi");
  TBranch *b16=tree->GetBranch("genMeson_E");
  TBranch *b17=tree->GetBranch("genMeson_mass");
  TBranch *b18=tree->GetBranch("genTrackminus_ID");
  TBranch *b19=tree->GetBranch("genTrackminus_pT");
  TBranch *b20=tree->GetBranch("genTrackminus_eta");
  TBranch *b21=tree->GetBranch("genTrackminus_phi");
  TBranch *b22=tree->GetBranch("genTrackminus_E");
  TBranch *b23=tree->GetBranch("genTrackplus_ID");
  TBranch *b24=tree->GetBranch("genTrackplus_pT");
  TBranch *b25=tree->GetBranch("genTrackplus_eta");
  TBranch *b26=tree->GetBranch("genTrackplus_phi");
  TBranch *b27=tree->GetBranch("genTrackplus_E");  

  // Definizione degli indirizzi per la lettura dei dati sul ttree
  b1->SetAddress(&genZ_ID);
  b2->SetAddress(&genZ_pT);
  b3->SetAddress(&genZ_eta);
  b4->SetAddress(&genZ_phi);
  b5->SetAddress(&genZ_E);
  b6->SetAddress(&genZ_mass);
  b7->SetAddress(&genGamma_ID);
  b8->SetAddress(&genGamma_pT);
  b9->SetAddress(&genGamma_eta);
  b10->SetAddress(&genGamma_phi);
  b11->SetAddress(&genGamma_E);
  b12->SetAddress(&genMeson_ID);
  b13->SetAddress(&genMeson_pT);
  b14->SetAddress(&genMeson_eta);
  b15->SetAddress(&genMeson_phi);
  b16->SetAddress(&genMeson_E);
  b17->SetAddress(&genMeson_mass);
  b18->SetAddress(&genTrackminus_ID);
  b19->SetAddress(&genTrackminus_pT);
  b20->SetAddress(&genTrackminus_eta);
  b21->SetAddress(&genTrackminus_phi);
  b22->SetAddress(&genTrackminus_E);
  b23->SetAddress(&genTrackplus_ID);
  b24->SetAddress(&genTrackplus_pT);
  b25->SetAddress(&genTrackplus_eta);
  b26->SetAddress(&genTrackplus_phi);
  b27->SetAddress(&genTrackplus_E);

   // loop sugli ingressi nel TTree
  for(int ev=0; ev<tree->GetEntries(); ev++){
    tree->GetEvent(ev);
    hist1->Fill(genZ_mass - m_Z_true);
    hist2->Fill(genMeson_mass - m_Phi_true);
    }
  
  hfile.Close();
  
  hist1->Draw();
  new TCanvas();
  hist2->Draw();  
  
}