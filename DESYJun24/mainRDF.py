import ROOT as R
import numpy as np
import pandas as pd
import time
from tqdm import tqdm
import os
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument("--col", nargs="*", default=0)
parser.add_argument("--row", nargs="*", default=0)
parser.add_argument("--board", default=3)

args = parser.parse_args()
COL    = args.col if type(args.col)==list else [args.col]
ROW    = args.row if type(args.row)==list else [args.row]
BOARD  = int(args.board)

COL = int(COL[0])
ROW = int(ROW[0])

toDeclare = [
    "tools/declarations.C"
]
toRedefine = {}
toDefine = {
    "isGoodEvent": "isGoodEvent(board)",
    #"acal":   "getacal(hacal,col,row)",
    #"toa_ps": "1000*(12.5-toa*3.125/acal)",
    #"tot_ps": "1000*((2.*tot-std::floorf(tot/32.))*3.125/acal)",
}

# ------------------------------------------
#        Define additional variables
# ------------------------------------------
def declareVariables(RDF):
    for path in toDeclare:
        R.gInterpreter.Declare(open(path,'r').read())
    for var_ in toDefine.keys():
        RDF = RDF.Define(var_, toDefine[var_])
    for var_ in toRedefine.keys():
        RDF = RDF.Redefine(var_, toRedefine[var_])
    return RDF

if __name__=='__main__':

    files = []
    for i in range(4):
        fileIn = f"rootFiles/loop_{i}.root"
        files.append(fileIn)
    print(files)
    RDF = R.RDataFrame("data_tree", files)
    RDF = declareVariables(RDF)
    RDF = RDF.Filter("isGoodEvent==1")

    R.gROOT.ProcessLine('.L ./tdrstyle.C')
    R.gROOT.SetBatch(1)
    R.setTDRStyle()

    # Fill histograms of avergae cal per pixel and cal vs time
    histmodel = R.RDF.TH1DModel(f"h_cal",f";CAL code;N events",1024,0.,1024.)
    RDF = RDF.Define(f"cal_B{BOARD}C{COL}R{ROW}",f"cal[(board=={BOARD})&&(col=={COL})&&(row=={ROW})]")
    start = time.time()
    h_temp = RDF.Histo1D(histmodel,f"cal_B{BOARD}C{COL}R{ROW}").GetPtr()
    end = time.time()
    print(f"Time elapsed: {end-start} s")
    f = R.TF1(f"fitfn","gaus")
    if h_temp.GetEntries() != 0: 
        maxval = h_temp.GetMaximumBin()
        f.SetParameter(1,maxval)
        h_temp.Fit(f"fitfn","q","",maxval-3,maxval+3);
        acal_temp = f.GetParameter(1)
    else:
        maxval = 160
        acal_temp = 0
    print(f"AV_CAL['B{BOARD}C{COL}R{ROW}'] = {acal_temp}")

    # Plot cal hist
    c = R.TCanvas(f"c_cal_B{BOARD}C{COL}R{ROW}","",800,800)
    c.cd()
    c.SetLogy()
    h_temp.SetMinimum(0.1)
    h_temp.SetFillColor(R.kGreen-7)
    h_temp.SetFillColor(R.kGreen-10)
    h_temp.GetXaxis().SetRangeUser(maxval-5,maxval+5)
    f.SetLineWidth(2)
    f.SetLineColor(R.kRed)
    h_temp.SetMinimum(0.01)
    h_temp.Draw()
    l = R.TLine(f.GetParameter(1),0.,f.GetParameter(1),h_temp.GetMaximum())
    l.Draw("L,SAME")
    c.SaveAs(f"h_cal/h_cal_B{BOARD}C{COL}R{ROW}.png")
