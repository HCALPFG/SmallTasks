#!/usr/bin/env python
import argparse
import sys

# Batch mode hack
sys.argv.append("-b-")
import ROOT
ROOT.gROOT.SetBatch(True)
sys.argv.remove("-b-")

# Import the Fill class from LHCFillScheme (https://github.com/RohanBhandari/LHCFillScheme/)
sys.path.insert(0,"/home/users/rbhandar/hcal/LHCFillScheme/")
from python.Fill import Fill

# Declare the fills and their corresponding runs
def set_fills():

    folder = "/home/users/rbhandar/hcal/LHCFillScheme/bunchpatterns/"

    fills = []
    fills.append(Fill(folder+"bunchpattern5822.txt", range(296642,296643)))
    fills.append(Fill(folder+"bunchpattern5883.txt", range(297670,297671)))
    fills.append(Fill(folder+"bunchpattern5950.txt", range(299061,299062)))    
    fills.append(Fill(folder+"bunchpattern6019.txt", range(300087,300088)))

    return fills

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Plot the TDC timing of the HF')
    parser.add_argument('--firstRun', type=int, default=1,      help='The first run to plot')
    parser.add_argument('--lastRun',  type=int, default=999999, help='The last run to plot')
    parser.add_argument('--firstFill', type=int, help='The first fill to plot')
    parser.add_argument('--lastFill',  type=int, help='The last fill to plot')
    args = parser.parse_args()

    # Create chain of pfg-tuples
    chain = ROOT.TChain("hcalTupleTree/tree")
    chain.Add("Run2017A/*5822*.root")
    chain.Add("Run2017B/*5883*.root")
    chain.Add("Run2017B/*5950*.root")
    chain.Add("Run2017C/*6019*.root")

    # Declare the fills to look at
    fills = set_fills()

    # Get the first/last run for the first/last fill to be used
    if args.firstFill:
        for fill in fills:
            if fill.fillN==args.firstFill:
                args.firstRun = fill.runs[0] 
    if args.lastFill:
        for fill in fills:
            if fill.fillN==args.lastFill:
                args.lastRun = fill.runs[-1]

    # Declare the reference run. Currently set to the last phase adjustment: 
    # https://indico.cern.ch/event/645440/contributions/2626379/attachments/1478691/2291882/170619_hcal_ops_cms_week_hf_hep17_phase_settings.pdf
    ref_run = 296642

    # Create dictionary between run number and a list of TProfile2D hists (one for each depth) that are filled with the TDC value
    hf_timing2D = {run : [ROOT.TProfile2D("","Run {}".format(run),84,-42,42,72,0,72)]*4 for fill in fills for run in fill.runs}

    nevents = -1
    print("Looping over {} events".format(chain.GetEntries()))
    # Loop over events
    for i,event in enumerate(chain):

        if i % 10000 == 0: print("Event #{}".format(i))
        if i == nevents: break
    
        if not (args.firstRun <= event.run <= args.lastRun): continue

        for fill in fills:
            # Select the appropriate fill and use only isolatedbxs
            if event.run not in fill.runs: continue
            # Use only one isolated bx per fill to speed up the script (assuming high stat runs)
            if event.bx != fill.isoBXs[0]: continue

            # Loop over channels
            for chan in range(len(event.QIE10DigiDepth)):
                # Require at least 50 fC in the channel to suppress noise
                if event.QIE10DigiFC[chan][1] <= 50: continue
                # Don't fill with TDC special codes
                if event.QIE10DigiLETDC[chan][1] >= 60: continue
                # Fill the appropriate channel with the TDC information (TS1 is SOI)
                hf_timing2D[event.run][event.QIE10DigiDepth[chan]-1].Fill(event.QIE10DigiIEta[chan], 
                                                                          event.QIE10DigiIPhi[chan], 
                                                                          event.QIE10DigiLETDC[chan][1]/2 + 25)

    print("Making plots")
    # Plot histograms for all for depths on the same pad and save
    ROOT.gStyle.SetOptStat(0)
    for fill in fills:
        for run in fill.runs:

            # Check against empty runs
            if sum([hf_timing2D[run][idepth].Integral() for idepth in range(4)]) == 0: continue

            print("Run #{}".format(run))

            # Save a 2D plot of all channels/depths for a run
            can = ROOT.TCanvas("c_run{}".format(run),"Run {}".format(run),1200,1200)
            can.Divide(2,2)
            can.cd(1)
            hf_timing2D[run][0].Draw("colz")
            can.cd(2)
            hf_timing2D[run][1].Draw("colz")
            can.cd(3)
            hf_timing2D[run][2].Draw("colz")
            can.cd(4)
            hf_timing2D[run][3].Draw("colz")

            for idepth in range(4):
                hf_timing2D[run][idepth].SetMinimum(30)
                hf_timing2D[run][idepth].SetMaximum(40)
            
            can.SaveAs("plots/run{}.pdf".format(run))

            # Compare against the reference run and save the plots
            for ix in range(1, hf_timing2D[ref_run][0].GetNbinsX()+1):
                for iy in range(1, hf_timing2D[ref_run][0].GetNbinsY()+1):
                    for idepth in range(4):
                        hf_timing2D[run][idepth].SetBinContent(ix, iy, hf_timing2D[run][idepth].GetBinContent(ix,iy)
                                                               - hf_timing2D[ref_run][idepth].GetBinContent(ix,iy))

            can2 = ROOT.TCanvas("c2_run{}".format(run),"Run {}".format(run),1200,1200)
            can2.Divide(2,2)
            can2.cd(1)
            hf_timing2D[run][0].Draw("colz")
            can2.cd(2)
            hf_timing2D[run][1].Draw("colz")
            can2.cd(3)
            hf_timing2D[run][2].Draw("colz")
            can2.cd(4)
            hf_timing2D[run][3].Draw("colz")
            
            for idepth in range(4):
                hf_timing2D[run][idepth].SetMinimum(-5)
                hf_timing2D[run][idepth].SetMaximum(5)

            can2.SaveAs("plots/run{}_comparison.pdf".format(run))
