// Requires using output from github rep: RohanBhandari/Misc/findFillScheme.py

#include <algorithm> 
#include <ctime>
#include <fstream>
#include <iomanip> // for setw()
#include <iostream>
#include <map>
#include <sstream>
#include <vector>

#include "TBranch.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TDirectory.h"
#include "TError.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TMath.h"
#include "TProfile.h"
#include "TROOT.h"
#include "TString.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TTree.h"

#include "FillScheme.h"
#include "FillMaps.h"
#include "HFTiming.h"

TString dataset = "ZeroBias";
TString which_bx = "all"; //Choose "all","iso", "first", "middle", "last", or "noniso"
TString fc_thresh = "100";
TString runset = "Run2016E";
bool doPUVeto = true;
TString pu_fcthresh; //Is set below if doPUVeto true

TString filldir="fillschemes/";
TString outdir ="/afs/cern.ch/user/r/rbhandar/www/hcal/hftiming/2016/"+dataset+"/"+which_bx+"_fc"+fc_thresh+"/"+runset+"/";

TChain* ch = new TChain("hcalTupleTree/tree");

int Ethres[4]={fc_thresh.Atoi(),100000,100001,100002};

// Macro
void HFTiming(){ 
  //Set Up
  gROOT->SetBatch(true);  //Don't show canvases
  gErrorIgnoreLevel=1001; //Don't print message about canvas being saved

  if(dataset=="ZeroBias"){
    if(runset=="Run2016B"){
      ch->Add("~/Work/public/hcaltuples/2016/ZeroBias_Run2016B-v1_RAW_DCS_272023_273146/*.root");
      ch->Add("~/Work/public/hcaltuples/2016/ZeroBias_Run2016B-v2_RAW_DCS_273150_273730/*.root");    
      ch->Add("~/Work/public/hcaltuples/2016/ZeroBias_Run2016B-v2_RAW_DCS_274094_274250/*.root");    
      ch->Add("~/Work/public/hcaltuples/2016/ZeroBias_Run2016B-v2_RAW_DCS_274282_274443/*.root");    
      ch->Add("~/Work/public/hcaltuples/2016/ZeroBias_Run2016B-v2_RAW_DCS_274954_275127/*.root");    
      ch->Add("~/Work/public/hcaltuples/2016/ZeroBias_Run2016B-v2_RAW_DCS_275282_275376/*.root");    
    }
    else if(runset=="Run2016C"){
      ch->Add("~/Work/public/hcaltuples/2016/ZeroBias_Run2016C-v2_RAW_DCS_275656_276283/*.root");    
    }
    else if(runset=="Run2016D"){
      ch->Add("~/Work/public/hcaltuples/2016/ZeroBias_Run2016D-v2_RAW_DCS_276315_276582/*.root");    
      ch->Add("~/Work/public/hcaltuples/2016/ZeroBias_Run2016D-v2_RAW_DCS_276583_276811/*.root");    
    }
    else if(runset=="Run2016E"){
      ch->Add("~/Work/public/hcaltuples/2016/ZeroBias_Run2016E-v2_RAW_DCS_276831_276870/*.root");    
    }
  }
  else if(dataset=="JetHT"){
    if(runset=="Run2016B"){
      ch->Add("~/Work/public/hcaltuples/2016/JetHT_Run2016B-v1_RAW_DCS_272023_273146/*.root");
      ch->Add("~/Work/public/hcaltuples/2016/JetHT_Run2016B-v2_RAW_DCS_273150_273730/*.root");    
    }
  }

  if(doPUVeto){
    pu_fcthresh = "10";
    outdir.ReplaceAll("2016/","2016/Pileup/");
    outdir.ReplaceAll("_fc"+fc_thresh+"/","_fc"+fc_thresh+"_pu"+pu_fcthresh+"/");
  }

  int dir = gSystem->mkdir(outdir,true);
  if(dir==0){
    cout<<"\n[HF Timing] Configuration: Created directory: "<<outdir<<"\n"<<endl;
    gSystem->CopyFile("/afs/cern.ch/user/r/rbhandar/www/index.php",outdir+"index.php");
  }
  else if(dir==-1) 
    cout<<"\n[HF Timing] Saving plots to "<<outdir<<"\n"<<endl;

  //Analysis
  map<unsigned int, int> fillmap = getFillMap(runset);

  vector<unsigned int> run;
  for(map<unsigned int,int>::iterator imap=fillmap.begin(); imap!=fillmap.end(); ++imap) //Fill run vector with analysis runs
    run.push_back((*imap).first);

  vector<vector<bool> > fillSchemes = getFillSchemes(fillmap, run, filldir);
  map<unsigned int, vector<int> > selectedBXMap = selectBXs(run, fillSchemes, which_bx);
  printSelectedBXs(true, fillmap, selectedBXMap, outdir);
  
  HFTimingOne(selectedBXMap, 2, 1, 41, 3, 2);
  HFTimingOne(selectedBXMap, 2, 1, -41, 3, 2);
  
  cout<<"\n[HF Timing PU] Plots available at "<<outdir.Copy().ReplaceAll("/afs/cern.ch/user/r/rbhandar/www/","rbhandar.web.cern.ch/rbhandar/")<<"\n"<<endl;
}

// Main Looper
void HFTimingOne(map<unsigned int, vector<int> > selectedBXs, int TStoCheck = 2, int TSadjacent = 1, int IETA=999, int IPHI=999, int DEPTH=999){ 
  time_t begtime, endtime;
  time(&begtime);

  cout<<"\n[HF Timing] Analyzing channel "<<IETA<<":"<<IPHI<<":"<<DEPTH<<"..."<<endl;
  
  init_chain(ch);

  
  //Full run vector with analysis runs (same as done above)
  vector<unsigned int> run;
  for(map<unsigned int,vector<int> >::iterator imap=selectedBXs.begin(); imap!=selectedBXs.end(); ++imap)
    run.push_back((*imap).first);
  const unsigned int Nrun=selectedBXs.size();

  // histograms
  int HistColor[4]={kBlack,kRed,kBlue,kGreen};
  TH1F *h1[Nrun][4], *h2[Nrun][4], *h12[Nrun][4], *h2over12[Nrun][4], *ht[Nrun][4], *h1over2[Nrun][4], *havgtime[Nrun][4]; 
  TProfile *h2profile[Nrun][4], *hshape[Nrun][4];

  for(unsigned int irun=0; irun<Nrun; irun++){
    for(int i=0; i<4; i++){
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
      hshape[irun][i]  = new TProfile(Form("hshape_run%i_E%iTo%i", run[irun],   Ethreslow,Ethreshigh),  "hshape",        10,0,10,0,500);

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
   
  // main event loop
  unsigned int skiprun = 9999;
  unsigned int nentries = (Int_t)ch->GetEntries();
  for(unsigned int ievent = 0; ievent<nentries; ievent++){
    ch->GetEntry(ievent); 
    
    // Status
    if((ievent<100&&ievent%10==0) || (ievent<1000&&ievent%100==0) || (ievent<10000&&ievent%1000==0) || 
       (ievent<100000&&ievent%10000==0) || (ievent%100000==0) || (ievent+1==nentries)){
      if (isatty(1)) {
	printf("\r[HF Timing] Processing Events: %i / %i (%i%%)",ievent,nentries,static_cast<int>((ievent+1)*100./nentries));
	fflush(stdout);
	if(ievent+1==nentries) printf("\n");
      }
    }        
    
    //Check if bx of interest
    if(!isSelectedBX(selectedBXs[run_], bx_)) continue;
    
    // loop over channels
    for(unsigned int ich=0; ich<HFDigiSubdet_->size(); ich++){ 
      // Selected only interesting channel
      if( IETA!=999 && !(HFDigiIEta_->at(ich)==IETA && HFDigiIPhi_->at(ich)==IPHI && HFDigiDepth_->at(ich)==DEPTH)) continue; 

      // Fill histograms 
      for(int i=0; i<4; i++){  // for different E or Q cuts
	if(i<3 && ((HFDigiFC_->at(ich).at(TSadjacent)+HFDigiFC_->at(ich).at(TStoCheck))<Ethres[i] || 
		   (HFDigiFC_->at(ich).at(TSadjacent)+HFDigiFC_->at(ich).at(TStoCheck))>Ethres[i+1])) continue;
	if(i==3 && ((HFDigiFC_->at(ich).at(TSadjacent)+HFDigiFC_->at(ich).at(TStoCheck))<Ethres[i])) continue;
	if(doPUVeto)
	  if((HFDigiFC_->at(ich).at(0)+HFDigiFC_->at(ich).at(3))>=pu_fcthresh.Atoi()) continue;
	
	//Make sure this run is part of the run vector above
	int ithisrun=-1;
	for(unsigned int irun=0; irun<Nrun; irun++) 
	  if(run_==run[irun]) ithisrun = irun;

	if(ithisrun==-1) {
	  if(run_!=skiprun){
	    cout<<"\n\e[31m[HF Timing]\e[0m WARNING: Run = "<<run_<<" is being skipped!"<<endl;
	    skiprun = run_;
	    continue;
	  }
	  else		    
	    continue;
	}

	//Filling Histograms
	h1[ithisrun][i]->Fill(Min(HFDigiFC_->at(ich).at(TSadjacent),299.999));
	h2[ithisrun][i]->Fill(Min(HFDigiFC_->at(ich).at(TStoCheck),299.999)); 
	h12[ithisrun][i]->Fill(Min(HFDigiFC_->at(ich).at(TSadjacent)+HFDigiFC_->at(ich).at(TStoCheck),599.999));
	h1over2[ithisrun][i]->Fill(Min(HFDigiFC_->at(ich).at(TSadjacent)/HFDigiFC_->at(ich).at(TStoCheck),4.999));
	h2over12[ithisrun][i]->Fill(Min(HFDigiFC_->at(ich).at(TStoCheck)/(HFDigiFC_->at(ich).at(TSadjacent)+HFDigiFC_->at(ich).at(TStoCheck)),9.999));
	h2profile[ithisrun][i]->Fill(ls_, HFDigiFC_->at(ich).at(TStoCheck)/(HFDigiFC_->at(ich).at(TSadjacent)+HFDigiFC_->at(ich).at(TStoCheck))); 
        for(unsigned int iTS=0; iTS<HFDigiFC_->at(ich).size(); iTS++)
          hshape[ithisrun][i]->Fill(iTS,HFDigiFC_->at(ich).at(iTS),HFDigiFC_->at(ich).at(iTS));                

	// energy average timing (considering only TS1 and TS2) 
	float avgtime=1; 
	float num =0;
	float deno =0;
	for(int icap=1; icap<3; icap++){ 
	  num = num + icap*(HFDigiFC_->at(ich).at(icap));
	  deno = deno + HFDigiFC_->at(ich).at(icap);
	}
	avgtime = num/deno;
	havgtime[ithisrun][i]->Fill(avgtime);  
      }

    } //for(unsigned int ich=0; ich<HFDigiSubdet_->size(); ich++)
  } //for(unsigned int ievent = 0; ievent<nentries; ievent++)

  
  // Draw 
  for(unsigned int irun=0; irun<Nrun; irun++){ 
    for(int i=0; i<4; i++){ 
      h1[irun][i]->Scale(1./h1[irun][i]->Integral());
      h2[irun][i]->Scale(1./h2[irun][i]->Integral());
      h12[irun][i]->Scale(1./h12[irun][i]->Integral());
      h1over2[irun][i]->Scale(1./h1over2[irun][i]->Integral());
      h2over12[irun][i]->Scale(1./h2over12[irun][i]->Integral());
      havgtime[irun][i]->Scale(1./havgtime[irun][i]->Integral());
      hshape[irun][i]->Scale(1./hshape[irun][i]->Integral());
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

  for(unsigned int irun=0; irun<Nrun; irun++){  

    // Legend 
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
        
    // Canvas 
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

    c->Print(Form(outdir+"/RUN%i_IETA%i_IPHI%i_DEPTH%i.png",run[irun],IETA,IPHI,DEPTH));
    delete c;
        
    TCanvas *cprofile = new TCanvas("cprofile", "cprofile", 600, 400);
    cprofile->cd(1);
    h2profile[irun][0]->Draw();
    if(h2profile[irun][0]->GetEntries()<10) cout<<"\e[31m[HF Timing]\e[0m WARNING: Low stats run = "<<run[irun]<<endl;
    cprofile->Print(Form(outdir+"/Profile_RUN%i_IETA%i_IPHI%i_DEPTH%i.png",run[irun],IETA,IPHI,DEPTH));
	
    delete cprofile;

    TCanvas *cshape = new TCanvas("cshape", "cshape", 600, 400);
    cshape->cd(1);
    hshape[irun][0]->Draw();
    cshape->Print(Form(outdir+"/SHAPE_RUN%i_IETA%i_IPHI%i_DEPTH%i.png",run[irun],IETA,IPHI,DEPTH));

    delete cshape;
  }

  // Make a summary plot 
  TH1F *hsummary = new TH1F("hsummary",  "hsummary", Nrun, 0, Nrun); 
  for(unsigned int irun=0; irun<Nrun; irun++){   
    if(h2over12[irun][0]->GetEntries()>=10){
      hsummary->SetBinContent(irun+1,h2over12[irun][0]->GetMean());
      hsummary->SetBinError(irun+1,h2over12[irun][0]->GetMeanError());     
    }
    hsummary->GetXaxis()->SetBinLabel(irun+1,Form("%i",(run[irun]-270000)));  
    hsummary->GetXaxis()->SetLabelSize(0.06);  
    hsummary->GetXaxis()->LabelsOption("v");
    hsummary->GetYaxis()->SetTitle("<Q2/(Q1+Q2)>");
    hsummary->GetYaxis()->SetTitleSize(0.06);
    hsummary->GetYaxis()->SetTitleOffset(0.8);
    hsummary->GetYaxis()->SetLabelSize(0.06);  
    hsummary->SetStats(0);
    hsummary->SetTitle(Form("iEta=%i, iPhi=%i, Depth=%i",IETA,IPHI,DEPTH));    
  } 
  gStyle->SetPadTickY(2);
  TCanvas *csum = new TCanvas("csum", "csum", 1600, 400);
  csum->cd(1);
  csum->SetGridy(1);
  hsummary->SetMinimum(0.2);
  hsummary->SetMaximum(0.8);
  hsummary->SetMarkerStyle(20);
  hsummary->SetMarkerSize(1);
  hsummary->Draw("e0");
  csum->Print(Form(outdir+"/Summary_IETA%i_IPHI%i_DEPTH%i.png",IETA,IPHI,DEPTH));
  csum->Print(Form(outdir+"/Summary_IETA%i_IPHI%i_DEPTH%i.root",IETA,IPHI,DEPTH));
  
  delete csum;

  makeTextFile(h2over12, IETA, run, outdir);
  
  // clean 
  for(unsigned int irun=0; irun<Nrun; irun++) {
    for(int i=0; i<4; i++){
      delete h1[irun][i]; 
      delete h2[irun][i]; 
      delete h12[irun][i]; 
      delete h1over2[irun][i]; 
      delete h2over12[irun][i]; 
      delete havgtime[irun][i]; 
      delete h2profile[irun][i]; 
      delete hshape[irun][i];
    } 
  }
  delete hsummary; 
  
  time(&endtime);
  int seconds = difftime(endtime, begtime);
  float hertz = nentries; hertz /= seconds;
  cout<<"[HF Timing] Processed channel "<<IETA<<":"<<IPHI<<":"<<DEPTH<<" in "<<seconds<<" seconds ("<<hoursMinSec(seconds)<<") for "<<nentries
      <<" events -> "<<roundNumber(hertz,1,1000)<<" kHz, "<<roundNumber(1000,2,hertz)<<" ms per event"<<endl<<endl;
}
