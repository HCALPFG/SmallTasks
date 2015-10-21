#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <iomanip> // for setw()
#include <algorithm>    // std::min

#include "TROOT.h"
#include "TChain.h"
#include "TF1.h"
#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TDirectory.h"
#include "TBranch.h"
#include "TString.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TProfile.h"

// In order to use vector of vectors : vector<vector<data type> >
// ACLiC makes dictionary for this
// [ref] http://root.cern.ch/phpBB3/viewtopic.php?f=3&t=10236&p=44117#p44117
#ifdef __MAKECINT__
#pragma link C++ class std::vector < std::vector<float> >+;
#endif

int Nrun=45;
int run[45]={254231,254232,254790,254852,254879,254906,254907,254914,  // 2015C
             256630,256673,256674,256675,256676,256677,256801,256842,
             256843,256866,256867,256868,256869,256926,256941,257461,
             257531,257599,257613,257614,257645,257682,257722,257723,
             257735,257751,257804,257805,257816,257819,257968,257969,
             258129,258136,258157,258158,258159};

//int Ethres[4]={5,1000,1001,1002};
//int HistColor[4]={kBlack,kBlack,kBlack,kBlack};
//int Ethres[4]={5,10,15,20};
int Ethres[4]={50,100000,1000001,1000002};
int HistColor[4]={kBlack,kRed,kBlue,kGreen};
bool DoNorm=false;

float Min(float a, float b) { return a <= b ? a : b; }
float Max(float a, float b) { return a >= b ? a : b; }

//
// h1 cosmetics
//
void h1cosmetic(TH1F* &h1, /*char* title, */int linecolor=kBlack, int linewidth=1, int fillcolor=0/*, TString var=""*/)
{
    h1->SetLineColor(linecolor);
    h1->SetLineWidth(linewidth);
    h1->SetMarkerColor(linecolor);
    h1->SetFillColor(fillcolor);
    //h1->SetTitle(title);
    //h1->SetXTitle(var);
    h1->SetStats(0);
    h1->SetMinimum(0.1);
}

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

//
//
//
void HFTimingOne(int TStoCheck = 2, int TSadjacent = 1, int IETA=999, int IPHI=999, int DEPTH=999)
{ 
    TChain ch("hcalTupleTree/tree");
    ch.Add("~/public/JetHT_Run2015C-v1_RAW/*.root"); 
    ch.Add("~/public/JetHT_Run2015D-v1_RAW/*.root"); 

    // 
    // Branches 
    // 
    int   run_ = 0;
    ch.SetBranchAddress("run", &run_);
    int   ls_ = 0;
    ch.SetBranchAddress("ls", &ls_);
    int   bx_ = 0;
    ch.SetBranchAddress("bx", &bx_);
    vector<int>   *HFDigiSubdet_ = 0;
    ch.SetBranchAddress("HFDigiSubdet", &HFDigiSubdet_);
    vector<int>   *HFDigiIEta_ = 0;
    ch.SetBranchAddress("HFDigiIEta", &HFDigiIEta_);
    vector<int>   *HFDigiIPhi_ = 0;
    ch.SetBranchAddress("HFDigiIPhi", &HFDigiIPhi_);
    vector<int>   *HFDigiDepth_ = 0;
    ch.SetBranchAddress("HFDigiDepth", &HFDigiDepth_);
    vector<float>   *HFDigiRecEnergy_ = 0;
    ch.SetBranchAddress("HFDigiRecEnergy", &HFDigiRecEnergy_);
    vector<vector<float> >   *HFDigiFC_ = 0;
    ch.SetBranchAddress("HFDigiFC", &HFDigiFC_);  // Note that FC is after ped subtraction while NomFC is not 
    vector<vector<float> >   *HFDigiPedFC_ = 0;
    ch.SetBranchAddress("HFDigiPedFC", &HFDigiPedFC_); 
    vector<vector<float> >   *HFDigiNomFC_ = 0;
    ch.SetBranchAddress("HFDigiNomFC", &HFDigiNomFC_); 

    //
    // histograms
    //
    TH1F *h1[Nrun][4], *h2[Nrun][4], *h12[Nrun][4], *h2over12[Nrun][4], *ht[Nrun][4], *h1over2[Nrun][4], *havgtime[Nrun][4]; 
    TProfile *h2profile[Nrun][4];

    for(int irun=0; irun<Nrun; irun++) 
    {
        for(int i=0; i<4; i++)
        {
            int Ethreslow=Ethres[i];
            int Ethreshigh;
            if(i==3) Ethreshigh=999;
            else Ethreshigh=Ethres[i+1];

            h1[irun][i]        = new TH1F(Form("h1_run%i_E%iTo%i", run[irun],         Ethreslow,Ethreshigh),  Form("TS%i",TSadjacent),                                 50, 0, 300);
            h2[irun][i]        = new TH1F(Form("h2_run%i_E%iTo%i", run[irun],         Ethreslow,Ethreshigh),  Form("TS%i",TStoCheck),                                  50, 0, 300);
            h12[irun][i]       = new TH1F(Form("h12_run%i_E%iTo%i", run[irun],        Ethreslow,Ethreshigh),  Form("TS%i+TS%i",TSadjacent,TStoCheck),                  50, 0, 600);
            h2over12[irun][i]  = new TH1F(Form("h2over12_run%i_E%iTo%i", run[irun],    Ethreslow,Ethreshigh),  Form("TS%i/(TS%i+TS%i)",TStoCheck,TSadjacent,TStoCheck), 20, 0, 1);
            h1over2[irun][i]   = new TH1F(Form("h1over2_run%i_E%iTo%i", run[irun],    Ethreslow,Ethreshigh),  Form("TS%i/TS%i",TSadjacent,TStoCheck),                  20, 0, 5);
            havgtime[irun][i]  = new TH1F(Form("havgtime_run%i_E%iTo%i", run[irun],   Ethreslow,Ethreshigh),  "Energy-avg timinig (in the unit of TS)",                30, 0, 3);
            h2profile[irun][i]  = new TProfile(Form("h2profile_run%i_E%iTo%i", run[irun],   Ethreslow,Ethreshigh),  "h2profile",        1000,0,1000,0,1);
        
            h1[irun][i]->Sumw2();
            h2[irun][i]->Sumw2();
            h12[irun][i]->Sumw2();
            h2over12[irun][i]->Sumw2();
            h1over2[irun][i]->Sumw2();
            havgtime[irun][i]->Sumw2();

            h1cosmetic(h1[irun][i],       HistColor[i], 2, 0); 
            h1cosmetic(h2[irun][i],       HistColor[i], 2, 0); 
            h1cosmetic(h12[irun][i],      HistColor[i], 2, 0); 
            h1cosmetic(h2over12[irun][i], HistColor[i], 2, 0); 
            h1cosmetic(h1over2[irun][i],  HistColor[i], 2, 0); 
            h1cosmetic(havgtime[irun][i], HistColor[i], 2, 0); 
        } 
    } 
   


    //
    // main event loop
    //
    unsigned int nentries = (Int_t)ch.GetEntries();
    cout << "[HF Timing] The number of entries is: " << nentries << endl;
    for(unsigned int ievent = 0; ievent<nentries; ievent++)
    {
        ch.GetEntry(ievent); 
    
        // Status
        if((ievent%100000)==0) cout << "[HF Timing] Event: " << ievent << " / " << nentries << "(" << (int)((float)ievent/(float)nentries*100) << "%)" << endl;

        // Choose isolated bunch crossing
        if( run_==254231                    && bx_!=895 && bx_!=1780 && bx_!=2674   ) continue;
        if( run_==254232                    && bx_!=895 && bx_!=1780 && bx_!=2674   ) continue;
        if( run_==254790                    && bx_!=1   && bx_!=61   && bx_!=141    ) continue;
        if( run_==254852                    && bx_!=39  && bx_!=91   && bx_!=141    ) continue;
        if((run_>=254879 && run_<=254914)   && bx_!=39  && bx_!=91                  ) continue;
        if( run_>=256630                    && bx_!=39                              ) continue;

        // loop over channels
        for(unsigned int ich=0; ich<HFDigiSubdet_->size(); ich++)
        { 
            // Selected only interesting channel
            if( IETA!=999 && !(HFDigiIEta_->at(ich)==IETA && HFDigiIPhi_->at(ich)==IPHI && HFDigiDepth_->at(ich)==DEPTH)) continue; 

            // Fill histograms 
            for(int i=0; i<4; i++) // for different E or Q cuts
            { 
                if(i<3 && ((HFDigiFC_->at(ich).at(TSadjacent)+HFDigiFC_->at(ich).at(TStoCheck))<Ethres[i] || 
                          (HFDigiFC_->at(ich).at(TSadjacent)+HFDigiFC_->at(ich).at(TStoCheck))>Ethres[i+1]) ) continue;
                if(i==3 && ((HFDigiFC_->at(ich).at(TSadjacent)+HFDigiFC_->at(ich).at(TStoCheck))<Ethres[i]) ) continue;

                if(0) // DEBUG
                {
                    cout << Form("[HF Timing] Subdet=%i IEta=%i IPhi=%i Depth=%i RecEnergy=%.3f\n",
                            HFDigiSubdet_->at(ich), 
                            HFDigiIEta_->at(ich),    
                            HFDigiIPhi_->at(ich),    
                            HFDigiDepth_->at(ich),   
                            HFDigiRecEnergy_->at(ich) ); 
                } 
                
                //
                // Filling histogram
                //
                
                // Get irun   
                int ithisrun=-1;
                for(int irun=0; irun<Nrun; irun++) if(run_==run[irun]) ithisrun = irun; 
                if(ithisrun==-1) { cout << "Run number does not match!!!" << endl; continue; }

                h1[ithisrun][i]->Fill(Min(HFDigiFC_->at(ich).at(TSadjacent),299.999));
                h2[ithisrun][i]->Fill(Min(HFDigiFC_->at(ich).at(TStoCheck),299.999)); 
                h12[ithisrun][i]->Fill(Min(HFDigiFC_->at(ich).at(TSadjacent)+HFDigiFC_->at(ich).at(TStoCheck),599.999));
                h1over2[ithisrun][i]->Fill(Min(HFDigiFC_->at(ich).at(TSadjacent)/HFDigiFC_->at(ich).at(TStoCheck),4.999));
                h2over12[ithisrun][i]->Fill(Min(HFDigiFC_->at(ich).at(TStoCheck)/(HFDigiFC_->at(ich).at(TSadjacent)+HFDigiFC_->at(ich).at(TStoCheck)),9.999));
                h2profile[ithisrun][i]->Fill(ls_, HFDigiFC_->at(ich).at(TStoCheck)/(HFDigiFC_->at(ich).at(TSadjacent)+HFDigiFC_->at(ich).at(TStoCheck))); 
                

                //  
                // energy average timing (considering only TS1 and TS2) 
                //  
                float avgtime=1; 
                float num =0;
                float deno =0;
                for(int icap=1; icap<3; icap++)
                { 
                    num = num + icap*(HFDigiFC_->at(ich).at(icap));
                    deno = deno + HFDigiFC_->at(ich).at(icap);
                }
                avgtime = num/deno;
                havgtime[ithisrun][i]->Fill(avgtime);  
            }

        } //for(unsigned int ich=0; ich<HFDigiSubdet_->size(); ich++)
    } //for(unsigned int ievent = 0; ievent<nentries; ievent++)


    //
    // Draw 
    //
    if(DoNorm)
    {
        for(int irun=0; irun<Nrun; irun++) 
        { 
            for(int i=0; i<4; i++) 
            { 
                h1[irun][i]->Scale(1./h1[irun][i]->Integral());
                h2[irun][i]->Scale(1./h2[irun][i]->Integral());
                h12[irun][i]->Scale(1./h12[irun][i]->Integral());
                h1over2[irun][i]->Scale(1./h1over2[irun][i]->Integral());
                h2over12[irun][i]->Scale(1./h2over12[irun][i]->Integral());
                havgtime[irun][i]->Scale(1./havgtime[irun][i]->Integral());
            }
        }
    }

    TLatex *tex_RUN = new TLatex(0.4,0.85,Form("Run number = %i", run_));
    tex_RUN->SetNDC();
    tex_RUN->SetTextSize(0.04);
    tex_RUN->SetLineWidth(2);
    tex_RUN->SetTextAlign(12);
    
    TLatex *tex_ch = new TLatex(0.4,0.8,Form("ieta=%i iphi=%i depth=%i",IETA,IPHI,DEPTH));
    tex_ch->SetNDC();
    tex_ch->SetTextSize(0.04);
    tex_ch->SetLineWidth(2);
    tex_ch->SetTextAlign(12);

    for(int irun=0; irun<Nrun; irun++) 
    {  

       // if(h1[irun][0]->Integral()!=0) continue;
        
        //
        // Legend 
        //
        TLegend *l1 = new TLegend(0.4, 0.5, 0.9, 0.75);
        l1->SetBorderSize(0);
        l1->SetFillColor(0);
        l1->SetFillStyle(0);
        l1->SetFillColor(kWhite);
        l1->SetLineColor(kWhite);
        l1->SetShadowColor(kWhite);
        l1->AddEntry(h1[irun][0],        Form("%i < E_{rechit} < %i GeV ",Ethres[0],Ethres[1]),    "lp");
        l1->AddEntry(h1[irun][1],        Form("%i < E_{rechit} < %i GeV ",Ethres[1],Ethres[2]),    "lp");
        l1->AddEntry(h1[irun][2],        Form("%i < E_{rechit} < %i GeV ",Ethres[2],Ethres[3]),    "lp");
        l1->AddEntry(h1[irun][3],        Form("E_{rechit} > %i GeV",Ethres[3]),    "lp");
        
        //
        // Canvas 
        //
        TCanvas *c = new TCanvas("c", "c", 1200, 800);
        c->Divide(3,2);
        c->cd(1);
        h1[irun][0]->SetMaximum(h1[irun][0]->GetMaximum()*2);
        for(int i=0; i<4; i++) h1[irun][i]->Draw(Form("E %s",i==0?"":"SAME"));  
        tex_RUN->Draw("SAME");
        tex_ch->Draw("SAME");
        l1->Draw("same");
        c->cd(2);
        h2[irun][0]->SetMaximum(h2[irun][0]->GetMaximum()*2);
        for(int i=0; i<4; i++) h2[irun][i]->Draw(Form("E %s",i==0?"":"SAME"));
        tex_RUN->Draw("SAME");
        tex_ch->Draw("SAME");
        c->cd(3);
        h12[irun][0]->SetMaximum(h12[irun][0]->GetMaximum()*2);
        for(int i=0; i<4; i++) h12[irun][i]->Draw(Form("E %s",i==0?"":"SAME"));
        tex_RUN->Draw("SAME");
        tex_ch->Draw("SAME");
        c->cd(4);
        c->cd(4)->SetLogy(1);
        h1over2[irun][0]->SetMaximum(h1over2[irun][0]->GetMaximum()*50);
        for(int i=0; i<4; i++) h1over2[irun][i]->Draw(Form("E %s",i==0?"":"SAME"));
        tex_RUN->Draw("SAME");
        tex_ch->Draw("SAME");
        c->cd(5);
        c->cd(5)->SetLogy(0);
        h2over12[irun][0]->SetMaximum(h2over12[irun][0]->GetMaximum()*2);
        for(int i=0; i<4; i++) h2over12[irun][i]->Draw(Form("E %s",i==0?"":"SAME"));
        tex_RUN->Draw("SAME");
        tex_ch->Draw("SAME");
        c->cd(6);
        havgtime[irun][0]->SetMaximum(havgtime[irun][0]->GetMaximum()*2);
        for(int i=0; i<4; i++) havgtime[irun][i]->Draw(Form("E %s",i==0?"":"SAME"));
        tex_RUN->Draw("SAME");
        tex_ch->Draw("SAME");  
        

        cout << Form("RUN=%i IETA=%i IPHI=%i DEPTH=%i",run[irun],IETA,IPHI,DEPTH) << endl;
        cout << "Peak : " << havgtime[irun][0]->GetMean() << " +/- " << havgtime[irun][0]->GetRMS()<< endl;

        c->Print(Form("/afs/cern.ch/user/j/jaehyeok/www/HCAL/HFTiming/ForTimeX/RUN%i_IETA%i_IPHI%i_DEPTH%i.pdf",run[irun],IETA,IPHI,DEPTH));
        //c->Print(Form("/afs/cern.ch/user/j/jaehyeok/www/HCAL/HFTiming/ForTimeX/RUN%i_IETA%i_IPHI%i_DEPTH%i.C",run[irun],IETA,IPHI,DEPTH));

        
        TCanvas *cprofile = new TCanvas("cprofile", "cprofile", 600, 400);
        cprofile->cd(1);
        h2profile[irun][0]->Draw();
        cprofile->Print(Form("/afs/cern.ch/user/j/jaehyeok/www/HCAL/HFTiming/ForTimeX/Profile_RUN%i_IETA%i_IPHI%i_DEPTH%i.pdf",run[irun],IETA,IPHI,DEPTH));
    }

    //
    // Make a summary plot 
    //
    TH1F *hsummary = new TH1F("hsummary",  "hsummary", Nrun, 0, Nrun); 
    for(int irun=0; irun<Nrun; irun++) 
    {   
        hsummary->SetBinContent(irun+1,h2over12[irun][0]->GetMean());
        hsummary->SetBinError(irun+1,h2over12[irun][0]->GetMeanError());  
        //hsummary->SetBinError(irun+1,0);  
        hsummary->GetXaxis()->SetBinLabel(irun+1,Form("%i",run[irun]));  
        hsummary->GetXaxis()->SetLabelSize(0.07);  
        hsummary->SetStats(0);  
        hsummary->SetTitle("");  
    } 
    TCanvas *csum = new TCanvas("csum", "csum", 800, 400);
    csum->cd(1);
    hsummary->SetMinimum(0);
    hsummary->SetMaximum(1);
    hsummary->SetMarkerStyle(20);
    hsummary->SetMarkerSize(1);
    hsummary->Draw("ep");
    csum->Print(Form("/afs/cern.ch/user/j/jaehyeok/www/HCAL/HFTiming/ForTimeX/Summary_IETA%i_IPHI%i_DEPTH%i.pdf",IETA,IPHI,DEPTH));
    csum->Print(Form("/afs/cern.ch/user/j/jaehyeok/www/HCAL/HFTiming/ForTimeX/Summary_IETA%i_IPHI%i_DEPTH%i.C",IETA,IPHI,DEPTH));

    //
    // clean 
    //
    for(int irun=0; irun<Nrun; irun++) 
    {
        for(int i=0; i<4; i++)
        {
            delete h1[irun][i]; 
            delete h2[irun][i]; 
            delete h12[irun][i]; 
            delete h1over2[irun][i]; 
            delete h2over12[irun][i]; 
            delete havgtime[irun][i]; 
            delete h2profile[irun][i]; 
        } 
    }
    delete hsummary; 
}

//
//
//
void HFTiming() 
{ 

    HFTimingOne(2, 1, 41, 3, 2);
    HFTimingOne(2, 1, -41, 3, 2);
  
}
