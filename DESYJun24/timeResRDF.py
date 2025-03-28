import ROOT as R
import numpy as np
import pandas as pd
import time
from tqdm import tqdm
import os
from argparse import ArgumentParser
www = '/eos/r/rlopezru/www/MTD-ETL/DESYJun24/'

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
    "tools/declarations.C",
    "tools/aCalDict.C",
]
toRedefine = {}
toDefine = {
    "isGoodEvent": "isGoodEvent(board)",
    "acal":        "getAcal(board,col,row)",
    "toa_ps":      "toaToPs(toa,acal)",
    "tot_ps":      "totToPs(tot,acal)",
    "abscaldiff":  "abs(acal-cal)",
    "caldiff":     "acal-cal",
    "toa12":   "toa_ps[1]-toa_ps[2]",
    "toa23":   "toa_ps[2]-toa_ps[3]",
    "toa31":   "toa_ps[3]-toa_ps[1]",
    "deltatoa1": "(toa_ps[2]+toa_ps[3])/2-toa_ps[1]",
    "deltatoa2": "(toa_ps[3]+toa_ps[1])/2-toa_ps[2]",
    "deltatoa3": "(toa_ps[1]+toa_ps[2])/2-toa_ps[3]",
    "tot_ps1":   "tot_ps[1]",
    "tot_ps2":   "tot_ps[2]",
    "tot_ps3":   "tot_ps[3]",
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
    for i in range(1):
        fileIn = f"rootFiles/loop_{i}.root"
        files.append(fileIn)
    print(files)
    RDF = R.RDataFrame("data_tree", files)
    RDF = declareVariables(RDF)
    RDF = RDF.Filter("isGoodEvent==1")

    R.gROOT.ProcessLine('.L ./tdrstyle.C')
    R.gROOT.SetBatch(1)
    R.setTDRStyle()
    R.gStyle.SetPaintTextFormat("1.f")

    # Plot average cal per board
    '''for i in range(4):
        # Compute average cal
        h_acal_temp = R.TH2F(f"h_cal_B{i}", f"Board {i};col;row", 16, -0., 16., 16, -0., 16.)
        for col in range(16):
            for row in range(16):
                if R.getAcal(i,col,row) < 10.: continue
                h_acal_temp.Fill(col,row,R.getAcal(i,col,row))
        c = R.TCanvas(f"c_acal_B{i}","",800,800)
        c.cd()
        h_acal_temp.SetMaximum(h_acal_temp.GetMaximum()+5)
        h_acal_temp.SetMinimum(h_acal_temp.GetBinContent(15,15)-5)
        h_acal_temp.Draw("COLZ,TEXT")
        c.SaveAs(f"h_acal_B{i}.png")'''

    # Filter events with bad cal (|cal-acal|>4)
    RDF = RDF.Filter("passCalDiff(caldiff,3.)")
    '''c = R.TCanvas("c_caldiff","",800,800)
    c.cd()
    c.SetLogy()
    histmodel = R.RDF.TH1DModel("h_caldiff",";|av.CAL - CAL|;N events",25,-12.,13.)
    h = RDF.Histo1D(histmodel,"caldiff")
    h.SetFillColor(R.kRed-7)
    h.SetMinimum(0.1)
    h.Draw("HIST")
    c.SaveAs("caldiff.png")'''

    # Get only interesting pixel
    RDF = RDF.Filter(f"isInPixel(board,col,row,{BOARD},{COL},{ROW})")

    def getTimeRes(RDF, name):
        if RDF.Count() == 0: return [0.,0.,0.]
        c1 = R.TCanvas("c1", "#Delta TOA", 1800, 600)
        c1.Divide(3,1)
        funct = []
        hists = []
        for i in range(3):
            j = ((i+1) % 3) + 1
            i = i + 1
            fg = R.TF1(f"fg{i}","gaus")
            print(f"fg{i}")
            print(f"toa{i}{j}")
            xmean = RDF.Mean(f"toa{i}{j}").GetValue()
            xstd  = RDF.StdDev(f"toa{i}{j}").GetValue()
            print(xmean, xstd)
            xmin = xmean - 4*xstd
            xmax = xmean + 4*xstd
            histmodel = R.RDF.TH1DModel("h",";#Delta TOA_{"+str(i)+str(j)+"} [ps]; N events",50,xmin,xmax)
            htemp = RDF.Histo1D(histmodel, f"toa{i}{j}")
            print(htemp.GetEntries())
            hists.append(htemp)
            funct.append(fg)
            c1.cd(i)
            htemp.Fit(f"fg{i}","Wq","")
            htemp.SetMaximum(htemp.GetMaximum()*1.6)
            htemp.SetFillColor(R.kAzure+3)
            htemp.Draw()
        c1.SaveAs(name)
        s = []
        for fg in funct: s.append(fg.GetParameter(2))
        sigma1 = np.sqrt((s[2]**2+s[0]**2-s[1]**2)/2)
        sigma2 = np.sqrt((s[0]**2+s[1]**2-s[2]**2)/2)
        sigma3 = np.sqrt((s[1]**2+s[2]**2-s[0]**2)/2)
        sigmas = [sigma1,sigma2,sigma3]
        return sigmas

    def applyTWC(RDF):
        print('- Applying TWC...')
        f1 = R.TF1("f1","pol2",0.,8000.)
        f2 = R.TF1("f2","pol2",0.,8000.)
        f3 = R.TF1("f3","pol2",0.,8000.)
        h1 = R.RDF.TProfile1DModel("h1_prof",";TOT [ps];#Delta TOA [ps]",20,1500.,5500.)
        h2 = R.RDF.TProfile1DModel("h2_prof",";TOT [ps];#Delta TOA [ps]",20,3000.,5500.)
        h3 = R.RDF.TProfile1DModel("h3_prof",";TOT [ps];#Delta TOA [ps]",20,2000.,6500.)
        hmodels = [h1,h2,h3]
        hprof = []
        for j in range(3):
            hprof.append(RDF.Profile1D(hmodels[j],f"tot_ps{j+1}",f"deltatoa{j+1}"))
    
        ctwc = R.TCanvas("ctwc", "", 1800, 600)
        ctwc.Divide(3,1)
        for i in range(3):
            ctwc.cd(i+1)
            hprof[i].Draw("PE")
            hprof[i].Fit(f"f{i+1}","Wq")
        ctwc.SaveAs("TWC_profiles.png")
    
        f1pars = [f1.GetParameter(0),f1.GetParameter(1),f1.GetParameter(2)]
        f2pars = [f2.GetParameter(0),f2.GetParameter(1),f2.GetParameter(2)]
        f3pars = [f3.GetParameter(0),f3.GetParameter(1),f3.GetParameter(2)]
        RDF = RDF.Define("toa_ps1",f"toa_ps[1] + {f1pars[0]} + {f1pars[1]}*tot_ps1 + {f1pars[2]}*tot_ps1*tot_ps1")
        RDF = RDF.Define("toa_ps2",f"toa_ps[2] + {f2pars[0]} + {f2pars[1]}*tot_ps2 + {f2pars[2]}*tot_ps2*tot_ps2")
        RDF = RDF.Define("toa_ps3",f"toa_ps[3] + {f3pars[0]} + {f3pars[1]}*tot_ps3 + {f3pars[2]}*tot_ps3*tot_ps3")
        RDF = RDF.Redefine("toa12","toa_ps1-toa_ps2")
        RDF = RDF.Redefine("toa23","toa_ps2-toa_ps3")
        RDF = RDF.Redefine("toa31","toa_ps3-toa_ps1")
        return RDF

    # Get time resolution before TWC
    sigmasBefore = getTimeRes(RDF, "timeResBefore.png")
    print(sigmasBefore)

    # Compute TWC
    RDF = applyTWC(RDF)

    # Get time resolution after TWC 
    sigmasAfter = getTimeRes(RDF, "timeResAfter.png")
    print(sigmasAfter)
