import ROOT as R
import numpy as np
import pandas as pd
import time
from tqdm import tqdm
import os

''' 
Require events to have exactly one hit in each board
And to have CALcode values that are close to the average CALcode
'''
def filterEvents(event, pixel):
    col = pixel[0]
    row = pixel[1]
    if not len(event) == 4: return False
    if not len(event['board']) == len(np.unique(event['board'])): return False
    #if not (list(event['col'])[3]==col and list(event['row'])[3]==row): return False
    return True

def draw1Dhist(hist, opt="HIST"):
    c = R.TCanvas(f"c_{hist.GetName()}","")
    c.cd()
    c.SetLogy(0)
    hist.Draw(opt)
    return c

def drawProfile2D(hcalprof, zmin, zmax):
    c = R.TCanvas(f"c_{hcalprof.GetName()}","")
    c.cd()
    R.gStyle.SetPaintTextFormat("0.f")
    hcalprof.SetMinimum(zmin)
    hcalprof.SetMaximum(zmax)
    hcalprof.Draw("COLZ,TEXT")
    return c

def getAvCAL(cal, board, col, row):
    if len(cal) == 0: return 0
    c = R.TCanvas(f"cal_B{board}C{col}R{row}")
    c.cd()
    c.SetLogy()
    h = R.TH1F(f"hB{board}C{col}R{row}", "", 1024, 0., 1024.) # CAL ranges between 0 and 1023
    h.FillN(len(cal), cal, np.ones(len(cal)))
    maxval = h.GetMaximumBin()
    f = R.TF1(f"fB{board}C{col}R{row}","gaus")
    f.SetParameter(1,maxval)
    h.Fit(f"fB{board}C{col}R{row}","q","",maxval-3,maxval+3);
    h.SetFillColor(R.kGreen-10)
    h.GetXaxis().SetRangeUser(maxval-50,maxval+50)
    f.SetLineWidth(2)
    f.SetLineColor(R.kRed)
    h.SetMinimum(0.01)
    h.Draw()
    l = R.TLine(f.GetParameter(1),0.,f.GetParameter(1),h.GetMaximum())
    l.Draw("L,SAME")
    c.SaveAs(f"hists_cal/cal_B{board}C{col}R{row}.png")
    return f.GetParameter(1)

def plotMulti1D(c, hists, names, colors):
    c = R.TCanvas(f"{c}","",1600,800)
    c.cd()
    l = R.TLegend(0.6,0.2,0.8,0.4)
    for i,h in enumerate(hists):
        if i==0: h.Draw("HIST")
        else: h.Draw("HIST,SAME")
        h.SetLineColor(colors[i])
        h.SetLineWidth(2)
        l.AddEntry(h, names[i],"L")
    return c, l

def getTimeRes(tot, toa, name):
    c1 = R.TCanvas("c1", "'#Delta TOA", 1800, 600)
    c1.Divide(3,1)
    fg1 = R.TF1("fg1","gaus")
    fg2 = R.TF1("fg2","gaus")
    fg3 = R.TF1("fg3","gaus")
    hists = []
    for i in range(3):
        j = (i+1) % 3
        xmin = np.mean(toa[i]-toa[j]) - 4*np.std(toa[i]-toa[j])
        xmax = np.mean(toa[i]-toa[j]) + 4*np.std(toa[i]-toa[j])
        htemp = R.TH1F(f"h{i+1}",";#Delta TOA_{"+str(i+1)+str(j+1)+"} [ps]; N events",100,xmin,xmax)
        htemp.FillN(len(toa[i]), toa[i]-toa[j], np.ones(len(toa[i])))
        hists.append(htemp)
        c1.cd(i+1)
        htemp.Fit(f"fg{i+1}","q","")
        htemp.SetMaximum(htemp.GetMaximum()*1.6)
        htemp.SetFillColor(R.kAzure+3)
        htemp.Draw()
    c1.SaveAs(name)
    sigma12 = fg1.GetParameter(2)
    sigma23 = fg2.GetParameter(2)
    sigma13 = fg3.GetParameter(2)
    sigma1 = np.sqrt((sigma12**2+sigma13**2-sigma23**2)/2)
    sigma2 = np.sqrt((sigma12**2+sigma23**2-sigma13**2)/2)
    sigma3 = np.sqrt((sigma23**2+sigma13**2-sigma12**2)/2)
    return [sigma1,sigma2,sigma3]

def compute_time_res(df, board, col, row):
    print(f'- Computing time resolution for board {board} and pixel ({col},{row})')
    tracks_df = df.groupby('evt').filter(lambda x: len(x)==4 and list(x['col'])[board]==col and list(x['row'])[board]==row)
    if len(tracks_df)<2000: return [0,0,0]
    tot = []
    toa = []
    for i in range(3):
        bdata = tracks_df
        tot.append(np.array(bdata[bdata['board']==i+1]['tot_ps']))
        toa.append(np.array(bdata[bdata['board']==i+1]['toa_ps']))

    print('- Computing deltaTOA...')
    deltaTOA = []
    for i in range(3):
        j = (i+1) % 3
        k = (i+2) % 3
        deltaTOA.append((toa[j]+toa[k])/2-toa[i])

    h1 = R.TH2F(f"h1_{board}{col}{row}",";TOT_{1} [ps];#Delta TOA_{1} [ps]",75,1500,5500.,75,-3000,3000)
    h2 = R.TH2F(f"h2_{board}{col}{row}",";TOT_{2} [ps];#Delta TOA_{2} [ps]",75,3000,5500.,75,-3000,3000)
    h3 = R.TH2F(f"h3_{board}{col}{row}",";TOT_{3} [ps];#Delta TOA_{3} [ps]",75,2000,6500.,75,-3000,3000)
    h = [h1,h2,h3]
    for j in range(3):
        for i in range(len(tot[j])):
            h[j].Fill(tot[j][i],deltaTOA[j][i])

    print('- Time resolution before TWC')
    sigma = getTimeRes(tot,toa,f"timeres_{board}{col}{row}_beforeTWC.png")
    #print(f"sigma1 = {sigma[0]}")
    #print(f"sigma2 = {sigma[1]}")
    print(f"sigma3 = {sigma[2]}")

    print('- Applying TWC...')
    f1 = R.TF1("f1","pol2",0.,8000.)
    f2 = R.TF1("f2","pol2",0.,8000.)
    f3 = R.TF1("f3","pol2",0.,8000.)
    h1 = R.TProfile("h1_prof_{board}{col}{row}",";TOT [ps];#Delta TOA [ps]",20,1500.,5500.)
    h2 = R.TProfile("h2_prof_{board}{col}{row}",";TOT [ps];#Delta TOA [ps]",20,3000.,5500.)
    h3 = R.TProfile("h3_prof_{board}{col}{row}",";TOT [ps];#Delta TOA [ps]",20,2000.,6500.)
    h = [h1,h2,h3]
    for j in range(3):
        for i in range(len(tot[j])):
            h[j].Fill(tot[j][i],deltaTOA[j][i])

    ctwc = R.TCanvas("ctwc", "", 1800, 600)
    ctwc.Divide(3,1)
    for i in range(3):
        ctwc.cd(i+1)
        h[i].Draw("PE")
        h[i].Fit(f"f{i+1}","Wq")
    ctwc.SaveAs("TWC_profiles.png")

    f1pars = [f1.GetParameter(0),f1.GetParameter(1),f1.GetParameter(2)]
    f2pars = [f2.GetParameter(0),f2.GetParameter(1),f2.GetParameter(2)]
    f3pars = [f3.GetParameter(0),f3.GetParameter(1),f3.GetParameter(2)]
    ftot1 = f1pars[0] + f1pars[1]*tot[0] + f1pars[2]*tot[0]*tot[0]
    ftot2 = f2pars[0] + f2pars[1]*tot[1] + f2pars[2]*tot[1]*tot[1]
    ftot3 = f3pars[0] + f3pars[1]*tot[2] + f3pars[2]*tot[2]*tot[2]
    corrtoa1 = toa[0] + ftot1
    corrtoa2 = toa[1] + ftot2
    corrtoa3 = toa[2] + ftot3
    corrtoa = [corrtoa1,corrtoa2,corrtoa3]

    print('- Time resolution after TWC')
    sigma = getTimeRes(tot,corrtoa,f"timeres_{board}{col}{row}_afterTWC.png")
    #print(f"sigma1 = {sigma[0]}")
    #print(f"sigma2 = {sigma[1]}")
    print(f"sigma3 = {sigma[2]}")
    return sigma

if __name__=='__main__':
    # Read dataset
    dfs = []
    end_evt = 0
    for i in range(1):
        dftemp = pd.read_feather(f'datos/loop_{i}.feather')
        dftemp['evt'] += end_evt
        end_evt = max(dftemp['evt'])+1
        print(dftemp)
        dfs.append(dftemp)
    df = pd.concat(dfs)
    print(df)

    # Filter events
    print("Filtering events...")
    start = time.time()
    COL = 11
    ROW = 9
    df = df.groupby('evt').filter(filterEvents, pixel=[COL, ROW])
    end = time.time()
    print(f"Elapsed time: {end-start:.3f} s")

    R.gROOT.ProcessLine('.L ./tdrstyle.C')
    R.gROOT.SetBatch(1)
    R.setTDRStyle()

    # Fill histograms of avergae cal per pixel and cal vs time
    h_acal = []
    h_cal_time = []
    h_cal = []
    for i in range(4):
        # Compute average cal
        h_acal_temp = R.TH2F(f"h_cal_B{i}", f"Board {i};col;row", 16, -0., 16., 16, -0., 16.)
        for col in tqdm(range(16)):
            for row in range(16):
                dft = df[df['board']==i]
                dft = dft[dft['col']==col]
                dft = dft[dft['row']==row]
                cal = np.array(dft['cal'], dtype=float)
                acal = getAvCAL(cal, i, col, row)
                h_acal_temp.Fill(col, row, acal)
    
                event = np.array(dft['evt'], dtype=float)
                if len(event) == 0:
                    h_cal_time_temp = R.TH2F(f"h_cal_time_B{i}C{col}R{row}", f"B{i}C{col}R{row};Event;calCODE",
                                         1, 0., 1., 7, int(acal)-3, int(acal)+4)
                    h_cal_time.append(h_cal_time_temp)
                    h_cal_temp = R.TH1F(f"h_cal_B{i}C{col}R{row}", f"B{i}C{col}R{row};calCODE",
                                         7, int(acal)-3, int(acal)+4)
                    h_cal.append(h_cal_temp)
                else:
                    h_cal_time_temp = R.TH2F(f"h_cal_time_B{i}C{col}R{row}", f"B{i}C{col}R{row};Event;calCODE",
                                         len(event), min(event), max(event)+1, 7, int(acal)-3, int(acal)+4)
                    h_cal_time_temp.FillN(len(event), event, cal, np.ones(len(event)))
                    h_cal_time.append(h_cal_time_temp)
                    h_cal_temp = R.TH1F(f"h_cal_B{i}C{col}R{row}", f"B{i}C{col}R{row};calCODE",
                                         7, int(acal)-3, int(acal)+4)
                    h_cal_temp.FillN(len(event), cal, np.ones(len(event)))
                    h_cal.append(h_cal_temp)
        h_acal.append(h_acal_temp)

    # Produce plots
    '''if not os.path.exists('hists_acal'): os.makedirs('hists_acal')
    if not os.path.exists('hists_cal_evt'): os.makedirs('hists_cal_evt')
    if not os.path.exists('hists_cal'): os.makedirs('hists_cal')

    R.gStyle.SetPaintTextFormat("0.f")
    for i in range(4):
        cprof = R.TCanvas(f"cprof_B{i}","",800,800)
        cprof.cd()
        h_acal[i].Draw("COLZ,TEXT")
        h_acal[i].SetMinimum(h_acal[i].GetBinContent(16,14)-5)
        h_acal[i].SetMaximum(h_acal[i].GetBinContent(3,2)+5)
        cprof.SaveAs(f'hists_acal/acal_board_{i}.png')

    for board in range(4):
        for col in range(16):
            for row in range(16):
                cprof = R.TCanvas(f"c_B{board}C{col}R{row}","",800,800)
                cprof.cd()
                h_cal_time[32*board+16*col+row].Draw("COLZ")
                cprof.SaveAs(f'hists_cal_evt/cal_B{board}C{col}R{row}.png')
                #cprof = R.TCanvas(f"c_B{board}C{col}R{row}","",800,800)
                #cprof.cd()
                #h_cal[32*board+16*col+row].SetFillColor(R.kGreen-10)
                #h_cal[32*board+16*col+row].Draw("COLZ")
                #l = R.TLine(h_acal[board].GetBinContent(col,row),0.,h_acal[board].GetBinContent(col,row),h_cal[32*board+16*col+row].GetMaximum())
                #l.SetLineWidth(2)
                #l.Draw("SAME")
                #cprof.SaveAs(f'hists_cal/cal_B{board}C{col}R{row}.png')'''

    # Compute acal column
    print('- Adding average cal column')
    df['acal'] = df.apply(lambda x: h_acal[x['board']].GetBinContent(int(x['col']),int(x['row'])), axis=1)

    cdiff = R.TCanvas("cdiff", "", 800, 800)
    cdiff.cd()
    cdiff.SetLogy()
    h = R.TH1F("hdiff",";aCAL - CAL;N events",40,-20,20)
    diff = np.array(df[df['board']==3][df['col']==3][df['row']==11]['acal']-df[df['board']==3][df['col']==3][df['row']==11]['cal'], dtype=float)
    h.FillN(len(diff), diff, np.ones(len(diff)))
    h.Fit("gauss")
    h.SetFillColor(R.kRed-10)
    h.Draw("HIST")
    cdiff.SaveAs('caldiff.png')

    # Remove events with bad cal
    print('- Remove events with bad cal')
    df = df[abs(df['acal']-df['cal'])<4]
    #df = df.groupby('evt').filter(lambda x: len(x)==4 and list(x['col'])[3]==COL and list(x['row'])[3]==ROW)

    # Filter events
    print("Filtering events...")
    start = time.time()
    df = df.groupby('evt').filter(filterEvents, pixel=[COL,ROW])
    end = time.time()
    print(f"Elapsed time: {end-start:.3f} s")

    df['toa_ps'] = 1000*(12.5-df['toa']*3.125/df['acal'])
    df['tot_ps'] = 1000*((2.*df['tot']-np.floor(df['tot']/32.))*3.125/df['acal'])

    # Plot ToA and ToT (all pixels)
    NBINS = int(12500./50.)
    colors = [R.kBlack,
              R.TColor.GetColor("#44AF69"),
              R.TColor.GetColor("#F8333C"),
              R.TColor.GetColor("#5056FF")]

    leg = []
    toa, htoa = [], []
    for i in range(4):
        toa.append(np.array(df[df['board']==3]['toa_ps'], dtype=float))
        htoa.append(R.TH1F(f"htoa{i}",";TOA [ps];N events", NBINS, 0., 12500.))
        htoa[i].FillN(len(toa[i]), toa[i], np.ones(len(toa[i])))
        leg.append(f"Board {i}")
    ctoa, ltoa = plotMulti1D("ctoa", htoa, leg, colors)
    ctoa.cd()
    ltoa.Draw()
    leg = []
    tot, htot = [], []
    for i in range(4):
        tot.append(np.array(df[df['board']==3]['tot_ps'], dtype=float))
        htot.append(R.TH1F(f"htot{i}",";TOT [ps];N events", NBINS, 0., 12500.))
        htot[i].FillN(len(tot[i]), tot[i], np.ones(len(tot[i])))
        leg.append(f"Board {i}")
    ctot, ltot = plotMulti1D("ctot", htot, leg, colors)
    ctot.cd()
    ltot.Draw()
    ctoa.SaveAs("h_toa_Average_pm4.png")
    ctot.SaveAs("h_tot_Average_pm4.png")

    '''timeres = []
    h_timeres = R.TH2F("h_timeres", f"Board 3;col;row", 16, -0., 16., 16, -0., 16.)
    for col in range(16):
        for row in range(16):
            sigmas = compute_time_res(df, 3, col, row)
            h_timeres.Fill(col,row,sigmas[2])
            timeres.append(sigmas[2])
    c = R.TCanvas("c","",800,800)
    c.cd()
    R.gStyle.SetPaintTextFormat("1.f")
    h_timeres.Draw("COLZ,TEXT")
    c.SaveAs("timeRes_board2.png")
    print(timeres)'''
