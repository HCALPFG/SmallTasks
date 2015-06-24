// ------------------------------------------------------------------------------------
//  ROOT macro that produces average RecHit energy from PFG ntuples
//
//  Author : Jae Hyeok Yoo (jae.hyeok.yoo@cern.ch)
//  Written on 12 June 2015 
// ------------------------------------------------------------------------------------
//  
// Pre-requisite :
//
//   You should have the PFG ntuple for the Run from which you want to do a measurement. 
//   Instruction on how to make PFG ntuples can be found here : FIXME link here 
//
//   You should have "Fig" directory for plots 
//
// Usage : 
//
//   $ root -b  
//   root> .L HCALResponseCheck.C++ 
//   root> HCALResponseCheck("PFGtuples_local/*root")
//    
// -----------------------------------------------------------------------------------
// 

// 
// Indices of channels in the subdetectors 
//
//  HBHE -----------------------
//      IEta    = -29 - 29
//      IPhi    =  1 - 72 
//      Depth   =  1,2,3 
//  HO -------------------------
//      IEta    = -15 - 15 
//      IPhi    =  1 - 72 
//      Depth   =  4 
//  HF ------------------------
//      IEta    = -41--29, 29-41
//      IPhi    =  3,5,7,9,...,25 
//      Depth   =  1,2 
//

#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <iomanip> // for setw()
#include <algorithm> 

#include "TROOT.h"
#include "TF1.h"
#include "TMath.h"
#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TDirectory.h"
#include "TBranch.h"
#include "TString.h"
#include "TStyle.h"
#include "TInterpreter.h"
#include "TStyle.h"

// In order to use vector of vectors : vector<vector<data type> >
// ACLiC makes dictionary for this
// [ref] http://root.cern.ch/phpBB3/viewtopic.php?f=3&t=10236&p=44117#p44117
#ifdef __MAKECINT__
#pragma link C++ class std::vector < std::vector<int> >+;
#pragma link C++ class std::vector < std::vector<float> >+;
#endif

using namespace std;

bool DRAWPLOTS  = false;  // draw plots or not (make "Fig" directory first before turning this on)
bool VERBOSE    = false;  // print out mean +/- sigma for each channel or not

//
// h2 cosmetics
//
void h2cosmetic(TH2F* &h2, char* title, TString Xvar="", TString Yvar="", TString Zvar="Events/bin")
{
    h2->SetTitle(title);
    h2->SetXTitle(Xvar);
    h2->SetYTitle(Yvar);
    h2->SetZTitle(Zvar);
    h2->SetStats(0);
}


const char* GetDetName(int Subdet) 
{ 
    const char* DetName;
    if(Subdet==1) DetName = "HB"; 
    if(Subdet==2) DetName = "HE"; 
    if(Subdet==3) DetName = "HO"; 
    if(Subdet==4) DetName = "HF"; 
    return DetName;
}

//
void HCALResponseCheckSubdet(TString rootfile="../HcalNtuples_239995_Viktor.root", TString SubDet="HB", int option=2) 
{ 

    //gInterpreter->ExecuteMacro("~/macros/JaeStyle.C");
    //gStyle->SetOptStat(0); 

    cout << "[HCAL response check] Running option " << option << " for " << SubDet << endl; 

    // fit pannel display option
    gStyle->SetOptFit(1011);

    //
    // Get the tree from the PFG ntuple 
    //
    TChain *ch = new TChain("hcalTupleTree/tree");
    ch->Add(rootfile);

    //
    // Set up branch address
    //

    // event,ls and run  
    UInt_t   event_ = 0;
    ch->SetBranchAddress("event", &event_);
    UInt_t   ls_ = 0;
    ch->SetBranchAddress("ls", &ls_);
    UInt_t   run_ = 0;
    ch->SetBranchAddress("run", &run_);
    
    // HBHE
    vector<int>   *HBHEDigiRawID_ = 0;
    ch->SetBranchAddress("HBHEDigiRawID", &HBHEDigiRawID_);
    vector<int>   *HBHEDigiSubdet_ = 0;
    ch->SetBranchAddress("HBHEDigiSubdet", &HBHEDigiSubdet_);
    vector<int>   *HBHEDigiIEta_ = 0;
    ch->SetBranchAddress("HBHEDigiIEta", &HBHEDigiIEta_);
    vector<int>   *HBHEDigiIPhi_ = 0;
    ch->SetBranchAddress("HBHEDigiIPhi", &HBHEDigiIPhi_);
    vector<int>   *HBHEDigiDepth_ = 0;
    ch->SetBranchAddress("HBHEDigiDepth", &HBHEDigiDepth_);
    vector<float> *HBHEDigiRecEnergy_ = 0;
    ch->SetBranchAddress("HBHEDigiRecEnergy", &HBHEDigiRecEnergy_);
//    vector<vector<int> >   *HBHEDigiCapID_ = 0;
//    ch->SetBranchAddress("HBHEDigiCapID", &HBHEDigiCapID_);
//    vector<vector<float> >   *HBHEDigiNomFC_ = 0; // linearlized ADC count
//    ch->SetBranchAddress("HBHEDigiNomFC", &HBHEDigiNomFC_);
//    vector<vector<float> >   *HBHEDigiADC_ = 0; // unlinearlized ADC count
//    ch->SetBranchAddress("HBHEDigiADC", &HBHEDigiADC_);
    

    // 
    // Define histograms for each channel 
    //  - One channel has 4 capacitors, so there are four plots per channel
    //  - Unlearized ADC count goes from 0 to 127, so there are 128 bins 
    //    and the range is from -0.5 to 127.5 
    //    
    
    // number of indices in eta, phi, depth
    int nieta = 83;
    int niphi = 72;
    int ndepth = 4;
    TH1F *h1_RecEnergy[nieta][niphi][ndepth][4]; // the last dimention is capid 
    for(int ieta=0; ieta<nieta; ieta++) 
    { 
        for(int iphi=0; iphi<niphi; iphi++) 
        {
            for(int idepth=0; idepth<ndepth; idepth++) 
            { 
                for(int icap=0; icap<4; icap++)  
                {
                    h1_RecEnergy[ieta][iphi][idepth][icap] = new 
                    TH1F( Form("h1_RecEnergy_ieta%s_iphi%i_depth%i_cap%i", (ieta>=41?Form("%i",ieta-41):Form("m%i",41-ieta)), (iphi+1), (idepth+1), icap),
                          Form("h1_RecEnergy_ieta%s_iphi%i_depth%i_cap%i", (ieta>=41?Form("%i",ieta-41):Form("m%i",41-ieta)), (iphi+1), (idepth+1), icap),
                          60, -10, 50);
                    h1_RecEnergy[ieta][iphi][idepth][icap]->Sumw2();
                }
            }
        }
    }
    TH1F *h1_RecEnergy_avg[nieta][ndepth]; // RecEnergy avg for a given ieta 
    for(int ieta=0; ieta<nieta; ieta++) 
    {
        for(int idepth=0; idepth<ndepth; idepth++) 
        { 
            h1_RecEnergy_avg[ieta][idepth] = new 
                TH1F( Form("h1_RecEnergy_avg_ieta%s_depth%i",(ieta>=41?Form("%i",ieta-41):Form("m%i",41-ieta)),(idepth+1)),
                      Form("h1_RecEnergy_avg_ieta%s_depth%i",(ieta>=41?Form("%i",ieta-41):Form("m%i",41-ieta)),(idepth+1)),
                      60, -10, 50);
            h1_RecEnergy_avg[ieta][idepth]->Sumw2();
        }
    }

    TH2F *h2[ndepth]; 
    TH2F *h2phiavg[ndepth]; 
    TH2F *h2diff[ndepth]; 
    TH2F *h2fracdiff[ndepth]; 
    for(int idepth=0; idepth<ndepth; idepth++) 
    { 
        h2[idepth]          = new TH2F(Form("h2_depth%i",(idepth+1)),   Form("h2_depth%i",(idepth+1)),      83,-41.5,41.5,72,0.5,72.5);
        h2phiavg[idepth]    = new TH2F(Form("h2phiavg_depth%i",(idepth+1)),   Form("h2phiavg_depth%i",(idepth+1)),      83,-41.5,41.5,72,0.5,72.5);
        h2diff[idepth]      = new TH2F(Form("h2diff_depth%i",(idepth+1)),Form("h2diff_depth%i",(idepth+1)),   83,-41.5,41.5,72,0.5,72.5);
        h2fracdiff[idepth]      = new TH2F(Form("h2fracdiff_depth%i",(idepth+1)),Form("h2fracdiff_depth%i",(idepth+1)),   83,-41.5,41.5,72,0.5,72.5);
        h2[idepth]->Sumw2();
        h2phiavg[idepth]->Sumw2();
        h2diff[idepth]->Sumw2();
        h2fracdiff[idepth]->Sumw2();
    }

    //
    // Define and initialize arrays to be used to make text file 
    //
    float RecEnergy_mean[nieta][niphi][ndepth][4];    
    float RecEnergy_avg_mean[nieta][niphi][ndepth][4];   
    int DetId[nieta][niphi][ndepth];            // Id for channel : it is decimal in the ntuple but to be converted into Heximal  
    int Subdet[nieta][niphi][ndepth];           // Id for subdetectors : HB=1, HE=2, HO=3, HF=4
    for(int ieta=0; ieta<nieta; ieta++) 
    { 
        for(int iphi=0; iphi<niphi; iphi++) 
        {
            for(int idepth=0; idepth<ndepth; idepth++) 
            { 
                DetId[ieta][iphi][idepth] = -999.; 
                Subdet[ieta][iphi][idepth] = -999.; 
                for(int icap=0; icap<4; icap++)  
                {
                    RecEnergy_mean[ieta][iphi][idepth][icap] = -999.; 
                    RecEnergy_avg_mean[ieta][iphi][idepth][icap] = -999.; 
                }
            }
        }
    }

    //
    // Loop over entries
    //
    unsigned int nentries = (Int_t)ch->GetEntries();
    //nentries = 20000; // FIXME 
    cout << "[HCAL Pedestal table maker] The number of entries is: " << nentries << endl;

    // main event loop
    for(unsigned int ievent = 0; ievent<nentries; ievent++) 
    {
        ch->GetEntry(ievent); 
        
        // Progress indicator 
        if(ievent%100==0) cout << "[HCAL Pedestal table maker] Processed " << ievent << " out of " << nentries << " events" << endl; 

        // Fill HBHE
        if(SubDet=="HB" || SubDet=="HE" || SubDet=="HBHE") 
        { 
            for(unsigned int i=0; i<HBHEDigiSubdet_->size(); i++) 
            {
                if(SubDet=="HB" && HBHEDigiSubdet_->at(i)!=1) continue;
                if(SubDet=="HE" && HBHEDigiSubdet_->at(i)!=2) continue;
                
                int ieta =  HBHEDigiIEta_->at(i);
                int iphi =  HBHEDigiIPhi_->at(i);
                int idepth =  HBHEDigiDepth_->at(i);

                DetId[ieta+41][iphi-1][idepth-1] = HBHEDigiRawID_->at(i);  
                Subdet[ieta+41][iphi-1][idepth-1] = HBHEDigiSubdet_->at(i);  


                if(HBHEDigiRecEnergy_->at(i) < 5) continue;

                //h1_RecEnergy[ieta+41][iphi-1][idepth-1][0]->Fill(TMath::Min(HBHEDigiRecEnergy_->at(i),49.999)); 
                //h1_RecEnergy_avg[ieta+41][idepth-1]->Fill(TMath::Min(HBHEDigiRecEnergy_->at(i),49.999)); 
                h1_RecEnergy[ieta+41][iphi-1][idepth-1][0]->Fill(HBHEDigiRecEnergy_->at(i)); 
                h1_RecEnergy_avg[ieta+41][idepth-1]->Fill(HBHEDigiRecEnergy_->at(i)); 
            }
        } 

    } //for(unsigned int ievent = 0; ievent<nentries; ievent++) 
   
    // 
    // Extract mean and sigma 
    // 
    cout << endl; 
    cout << " ........................................................................................  " << endl; 
    cout << " ........................... Extraction of mean and sigma ...............................  " << endl; 
    cout << " ........................................................................................  " << endl; 
    cout << endl; 

    for(int ieta=0; ieta<nieta; ieta++) 
    { 
        for(int iphi=0; iphi<niphi; iphi++) 
        {
            for(int idepth=0; idepth<ndepth; idepth++) 
            { 
                if( h1_RecEnergy[ieta][iphi][idepth][0]->Integral()==0 ) continue;
                if( Subdet[ieta][iphi][idepth]==-999. ) continue; 

                if(VERBOSE) 
                { 
                    cout << "[HCAL Pedestal table maker] For ieta, iphi, depth, icap = ";
                    cout << (ieta-41) <<  ", " << (iphi+1) << ", " << (idepth+1) << ", " << 0 << endl;
                    cout << "[HCAL Pedestal table maker]   pedestal = " << h1_RecEnergy[ieta][iphi][idepth][0]->GetMean() << " +/- " 
                        << h1_RecEnergy[ieta][iphi][idepth][0]->GetRMS() << endl;  
                } 

                RecEnergy_mean[ieta][iphi][idepth][0]       = h1_RecEnergy[ieta][iphi][idepth][0]->GetMean();   
                RecEnergy_avg_mean[ieta][iphi][idepth][0]   = h1_RecEnergy_avg[ieta][idepth]->GetMean();

                h2[idepth]->SetBinContent(       ieta+1, iphi+1, RecEnergy_mean[ieta][iphi][idepth][0]);
                h2phiavg[idepth]->SetBinContent( ieta+1, iphi+1, RecEnergy_avg_mean[ieta][iphi][idepth][0]);
                h2diff[idepth]->SetBinContent(   ieta+1, iphi+1, (RecEnergy_mean[ieta][iphi][idepth][0]-RecEnergy_avg_mean[ieta][iphi][idepth][0]));
                h2fracdiff[idepth]->SetBinContent(   ieta+1, iphi+1, (RecEnergy_mean[ieta][iphi][idepth][0]-RecEnergy_avg_mean[ieta][iphi][idepth][0])/RecEnergy_avg_mean[ieta][iphi][idepth][0]);
            }
        }
    }

    // 
    // Drawing : pedestal distribution per channel 
    // 
    if(DRAWPLOTS)
    {
        cout << endl; 
        cout << " ........................................................................................  " << endl; 
        cout << " ..................................... Drawing ..........................................  " << endl; 
        cout << " ........................................................................................  " << endl; 
        cout << endl; 

        for(int ieta=0; ieta<nieta; ieta++) 
        { 
            for(int iphi=0; iphi<niphi; iphi++) 
            {
                for(int idepth=0; idepth<ndepth; idepth++) 
                { 

                    if(TMath::Abs(RecEnergy_mean[ieta][iphi][idepth][0]-RecEnergy_avg_mean[ieta][iphi][idepth][0])<1.) continue;

                    cout << ieta-41 << " " 
                         << iphi+1 << " " 
                         << idepth+1 << " " 
                         << RecEnergy_mean[ieta][iphi][idepth][0]-RecEnergy_avg_mean[ieta][iphi][idepth][0] 
                         << endl; 

                    if( h1_RecEnergy[ieta][iphi][idepth][0]->Integral()==0 ) continue;
                    if( Subdet[ieta][iphi][idepth]==-999. ) continue; 

                    // Canvas for each channel
                    TCanvas *c = new TCanvas("c", "c", 800, 400); 
                    c->Divide(2,1);  

                    c->cd(1); h1_RecEnergy[ieta][iphi][idepth][0]->Draw(); 
                    c->cd(2); h1_RecEnergy_avg[ieta][idepth]->Draw(); 

                    c->Print(Form("Fig/RecEnergy_ieta%s_iphi%i_depth%i_%s_option%i.C",(ieta>=41?Form("%i",ieta-41):Form("m%i",41-ieta)),(iphi+1),(idepth+1),GetDetName(Subdet[ieta][iphi][idepth]),option)); 
                    c->Print(Form("Fig/RecEnergy_ieta%s_iphi%i_depth%i_%s_option%i.pdf",(ieta>=41?Form("%i",ieta-41):Form("m%i",41-ieta)),(iphi+1),(idepth+1),GetDetName(Subdet[ieta][iphi][idepth]),option)); 
                    delete c;
                }
            }
        }
    } 

    // Always draw 2D plots
    for(int idepth=0; idepth<ndepth; idepth++) 
    { 
        h2cosmetic(h2[idepth], "Average rechit energy", "ieta", "iphi", "Events/bin");
        h2cosmetic(h2phiavg[idepth], "Average rechit energy averaged in iphi", "ieta", "iphi", "Events/bin");
        h2cosmetic(h2diff[idepth], "rechit energy difference ", "ieta", "iphi", "Events/bin");
        h2cosmetic(h2fracdiff[idepth], "Fractional rechit energy w.r.t. iphi-averaged rechit energy  ", "ieta", "iphi", "Events/bin");

        TCanvas *c2d = new TCanvas("c2d", "c2d", 1200, 800);  
        c2d->Divide(2,2);
        c2d->cd(1); 
        h2[idepth]->Draw("colz"); 
        c2d->cd(2); 
        h2phiavg[idepth]->Draw("colz"); 
        c2d->cd(3); 
        h2diff[idepth]->Draw("colz"); 
        c2d->cd(4); 
        h2fracdiff[idepth]->Draw("colz"); 
        c2d->Print(Form("Fig/RecEnergy_depth%i.pdf",(idepth+1))); 

    } 

}

//
// Main function
//
void HCALResponseCheck(TString rootfile="PFGtuples_local/*root")
{
    HCALResponseCheckSubdet(rootfile, "HBHE", 0);
}
