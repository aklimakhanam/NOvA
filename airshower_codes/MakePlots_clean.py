#!/usr/bin/env python
# coding: utf-8

# ## Imports

# In[1]:
import os
import ROOT
import glob
import tqdm
from tqdm import tqdm
from datetime import datetime
import pytz
import iy10_ultrashower_utilities as shower
import healpy as hp
import numpy as np
import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord
from calendar import month_abbr
# ## Define ROOT Histograms

# In[2]:

def Preliminary2():
    prelim = ROOT.TLatex(0.935, 0.95, "NOvA Preliminary")
    prelim.SetTextColor(ROOT.kBlue)
    prelim.SetNDC()           # use normalized device coordinates
    prelim.SetTextSize(2.0/30.0)
    prelim.SetTextAlign(32)
    prelim.Draw()

Multiplicity = ROOT.TH1D("Multiplicity", "Multiplicity distribution (opposite through-goers); Tracks per event; Counts", 200, 0, 400)
TrackDensity = ROOT.TH1D("TrackDensity", "Track density distribution (all tracks); Track density (m^{-2})", 50, 0, 1)
XsecAltAz = ROOT.TH2D("XsecAltAz", "Cross section vs Altitude, Azimuth; Azimuth; Altitude", 360, 0, 360, 90, 0, 90)
VA = ROOT.TH1D("VA","Vertical asymmetry distribution (YZ view); Ratio; Counts", 25, -0.5, 0.5)
DM = ROOT.TH2D("DM","Density vs Multiplicity distribution; Multiplicity; Density", 200, 0, 400, 50, 0, 1)

Az_BU = ROOT.TH1D("Az_BU","Azimuth | Non Uniform showers(Compass axis);Azimuth Angle;Events",72,0,360)
Az_alt_l20 = ROOT.TH1D("Az_GU_alt_l20","Azimuth | Uniform showers (Alt < 20, Compass axis);Azimuth Angle;Events",72,0,360)
Az_LD = ROOT.TH1D("Az_LD","Azimuth | Uniform showers (LD, Compass axis);Azimuth Angle;Events",72,0,360)
Az_U = ROOT.TH1D("Az_U","Azimuth | Uniform showers(Compass axis);Azimuth Angle;Events",72,0,360)
Az_MD = ROOT.TH1D("Az_HD","Azimuth | Uniform showers(MD, Compass axis);Azimuth Angle;Events",72,0,360)
Az_HD = ROOT.TH1D("Az_HD", "Azimuth | Uniform Showers(HD, Compass axis);Azimuth Angle;Events",72,0,360)

Az_BU_DA = ROOT.TH1D("Az_BU","Azimuth | Non Uniform showers(Detector Axis);Azimuth Angle;Events",72,0,360)
Az_alt_l20_DA = ROOT.TH1D("Az_GU_alt_l20","Azimuth | Uniform showers (Alt < 20,Detector Axis);Azimuth Angle;Events",72,0,360)
Az_LD_DA = ROOT.TH1D("Az_LD","Azimuth | Uniform showers (LD, Detector Axis);Azimuth Angle;Events",72,0,360)
Az_U_DA = ROOT.TH1D("Az_U_DA","Azimuth | Uniform showers(Detector Axis);Azimuth Angle;Events",72,0,360)
Az_MD_DA = ROOT.TH1D("Az_MD","Azimuth | Uniform showers (MD, Detector Axis);Azimuth Angle;Events",72,0,360)
Az_HD_DA = ROOT.TH1D("Az_HD", "Azimuth | Uniform Showers (HD, Detector Axis);Azimuth Angle;Events",72,0,360)

Alt_BU = ROOT.TH1D("Alt_BU","Altitude | Non Uniform showers;Altitide Angle;Events",18,0,90)
Alt_alt_l20 = ROOT.TH1D("Az_GU_alt_l20","Altitude | Uniform showers (Alt < 20);Altitide Angle;Events",18,0,90)
Alt_LD = ROOT.TH1D("Alt_LD","Altitude | Uniform showers (LD);Altitide Angle;Events",18,0,90)
Alt_MD = ROOT.TH1D("Alt_HD","Altitude | Uniform showers (MD);Altitide Angle;Events",18,0,90)
Alt_HD = ROOT.TH1D("Alt_HD", "Altitude | Uniform Showers (HD);Altitide Angle;Events",18,0,90)

RA_BU = ROOT.TH1D("RA_BU","Right Ascension | Non Uniform showers;Right Ascension (deg);Events", 24, -180, 180)
RA_alt_l20 = ROOT.TH1D("RA_GU_alt_l20","Right Ascension | Non Uniform showers;Right Ascension (deg);Events", 24, -180, 180)
RA_LD = ROOT.TH1D("RA_LD"," ;Right Ascension (deg);Events", 24, -180, 180)
RA_MD = ROOT.TH1D("RA_MD"," ;Right Ascension (deg);Events", 24, -180, 180)
RA_HD = ROOT.TH1D("RA_HD"," ;Right Ascension (deg);Events", 24, -180, 180)

hEvts_Month = ROOT.TH1D("hEvts_Month","Events by Month | All;Month;Events",108,0,108)
hEvtsHD_Month = ROOT.TH1D("hEvtsHD_Month","Events by Month | High Density;Month;Events",108,0,108)
hEvtsMD_Month = ROOT.TH1D("hEvtsMD_Month","Events by Month | Medium Density;Month;Events",108,0,108)
hEvtsLD_Month = ROOT.TH1D("hEvtsLD_Month","Events by Month | Low Density;Month;Events",108,0,108)
hEvtsNU_Month = ROOT.TH1D("hEvtsNU_Month","Events by Month | Non-Uniform;Month;Events",108,0,108)
hLiveTime_Month = ROOT.TH1D("hLiveTime_Month","Live Time by Month;Month;Live Time (seconds)",108,0,108)

hEvts_Wrap = ROOT.TH1D("hEvts_Wrap","Events by Month | All;Month;Events",12,0,12)
hEvtsHD_Wrap = ROOT.TH1D("hEvtsHD_Wrap","Events by Month | High Density;Month;Events",12,0,12)
hEvtsMD_Wrap = ROOT.TH1D("hEvtsMD_Wrap","Events by Month | Medium Density;Month;Events",12,0,12)
hEvtsLD_Wrap = ROOT.TH1D("hEvtsLD_Wrap","Events by Month | Low Density;Month;Events",12,0,12)
hEvtsNU_Wrap = ROOT.TH1D("hEvtsNU_Wrap","Events by Month | Non-Uniform;Month;Events",12,0,12)
hLiveTime_Wrap = ROOT.TH1D("hLiveTime_Wrap","Live Time by Month;Month;Live Time (seconds)",12,0,12)

hEvts_Month.Sumw2()
hEvtsHD_Month.Sumw2()
hEvtsMD_Month.Sumw2()
hEvtsLD_Month.Sumw2()
hEvtsNU_Month.Sumw2()

hSRDuration = ROOT.TH1D("hSRDuration","Subrun Duration;Duration (s);Subruns",300,0,300)

# ## Loop through processed (text) files

# In[3]:


# Store Subrun Event Multiplicty Info
SubRunInfo = {}
AltAz_BU = []
AltAz_alt_l20 = []
AltAz_U = []
AltAz_LM = []
AltAz_HM = [] 
AltAz_L = []
AltAz_M = []
AltAz_H = []
RABU = []
DECBU = []
RAl20 = []
DECl20 = []
RALD = []
DECLD = []
RAMD = []
DECMD = []
RAHD = []
DECHD = []
DetL = 60
DetW = 15.6
DetH = 15.6
for i in range(0,90):
    for j in range(0,360):
        alti = np.radians(i)
        azi = np.radians(j)
        axsec = DetL * DetH * abs( np.cos(alti) * np.sin(azi) ) + DetW * DetH * abs( np.cos(alti) * np.cos(azi) ) + DetL * DetW * abs( np.sin(alti) )
        XsecAltAz.Fill(j,i,axsec)

# Loop through files
file_path = glob.glob("/home/aklima/MSAnalysis/Codes/ntuples_text/*_SubrunInfo.txt")
for n in range(len(file_path)):
#for n in range(0,500):
    base_name = os.path.basename(file_path[n])
    name_split = os.path.splitext(base_name)
    file_base_name = name_split[0]
    file_base_name = file_base_name[0:file_base_name.find("_SubrunInfo")]
    print("Parsing file "+file_base_name)

    si_file_path = f"/home/aklima/MSAnalysis/Codes/ntuples_text/{file_base_name}_SubrunInfo.txt"
    if os.path.getsize(si_file_path) == 0:
        #print(f"SubrunInfo file {si_file_path} is empty. Skipping...")
        continue
    run_si, subrun_si, start_si, stop_si = np.loadtxt(f"/home/aklima/MSAnalysis/Codes/ntuples_text/{file_base_name}_SubrunInfo.txt", unpack=True, ndmin=2)

    tz = pytz.timezone('America/Chicago')
    
    for sr in tqdm( range(len(run_si)), desc="Parsing Subruns...", ascii=False, ncols=75 ):
        start = start_si[sr]
        stop  = stop_si[sr]
        time = ( stop + start ) / 2

        duration = stop - start
        hSRDuration.Fill( duration )
        if duration < 150: continue

        date = datetime.utcfromtimestamp(time)
        datelocal = tz.localize(date)
        year = datelocal.year
        month = datelocal.month
        hLiveTime_Month.Fill( (year-2015)*12 + month-1, duration )
        hLiveTime_Wrap.Fill( month-1, duration )
        
        SubRunInfo[ ( int(run_si[sr]), int(subrun_si[sr]) ) ] = {"time":time, "duration":duration, "nNU":0, "nLD":0, "nMD":0, "nHD":0}

    ei_file_path = f"/home/aklima/MSAnalysis/Codes/ntuples_text/{file_base_name}_EventInfo.txt"
    if not os.path.exists(ei_file_path) or os.path.getsize(ei_file_path) == 0:
        #print(f"EventInfo file {ei_file_path} is empty. Skipping...")
        continue
    run_ei, subrun_ei, event, start_ei, stop_ei, ocxz, ocyz, asymxz, asymyz, anglexz, angleyz, sigmaxz, sigmayz, altitude, azimuth, ntracks, ntracks_n, ntracks_a, xsec, trkden, vasym, timestamp = np.loadtxt(f"/home/aklima/MSAnalysis/Codes/ntuples_text/{file_base_name}_EventInfo.txt", unpack=True, ndmin=2)
    
    # Loop through events
    for evt in tqdm( range(len(run_ei)), desc="Parsing Events....", ascii=False, ncols=75 ):

        start = start_ei[evt]
        stop  = stop_ei[evt]
        time = ( stop + start ) / 2
        duration = stop - start
        if (duration < 150): continue

        SubRunInfo[ ( int(run_ei[evt]), int(subrun_ei[evt]) ) ]["time"] = time
        SubRunInfo[ ( int(run_ei[evt]), int(subrun_ei[evt]) ) ]["duration"] = duration

        evtTime = shower.ProcessTimestamp(int(timestamp[evt]))
        date = datetime.utcfromtimestamp(evtTime.utc.value)
        datelocal = tz.localize(date)
        year = datelocal.year
        month = datelocal.month
        hEvts_Month.Fill( (year-2015)*12 + month-1 )
        hEvts_Wrap.Fill( month-1 )
        
        
        alt, az = shower.GetAltAz(anglexz[evt],angleyz[evt])
        #if (alt < 20): continue
        if ( anglexz[evt] == 0 and angleyz[evt] == 0 ): continue
        
        Multiplicity.Fill(ntracks_n[evt])         # Opposite face through goers
        TrackDensity.Fill(trkden[evt])            # All tracks
        DM.Fill(ntracks_a[evt],trkden[evt])       # All tracks
        VA.Fill(vasym[evt])

        # Uniform Shower. Determine Sky Coordinates and Such
        skyCoord = shower.MakeSkyCoord(alt, az, evtTime)
        ra = skyCoord.icrs.ra.degree
        RA = skyCoord.icrs.ra.degree
        #print('ra:', ra)
        if ra > 180 : ra = ra - 360
        dec = skyCoord.icrs.dec.degree
        DEC = 90 - skyCoord.icrs.dec.degree
            
        if(np.sqrt(asymxz[evt]**2 + asymyz[evt]**2) > 0.08 or vasym[evt] < -0.05 or vasym[evt] > 0.2):
            if(alt < 20):
                Az_alt_l20.Fill(az)
                Az_alt_l20_DA.Fill(anglexz[evt])
                Alt_alt_l20.Fill(alt)
                RA_alt_l20.Fill(ra)
                AltAz_alt_l20.append((az,alt))
                #AltAz_GU_alt_l20.append((anglexz[evt],alt))
                RAl20.append(RA)
                DECl20.append(DEC)
            else:
                Az_BU.Fill(az)
                Az_BU_DA.Fill(anglexz[evt])
                Alt_BU.Fill(alt)
                RA_BU.Fill(ra)
                AltAz_BU.append((az,alt))
                #AltAz_BU.append((anglexz[evt],alt))
                RABU.append(RA)
                DECBU.append(DEC)
                hEvtsNU_Month.Fill( (year-2015)*12 + month-1 )
                hEvtsNU_Wrap.Fill( month-1 )
                SubRunInfo[ ( int(run_ei[evt]), int(subrun_ei[evt]) ) ]["nNU"] += 1
        #elif((asymxz[evt]**2 + asymyz[evt]**2) <= 0.08 and vasym[evt] >= -0.05 and vasym[evt] <= 0.2):
        else:
            Az_U.Fill(az)
            Az_U_DA.Fill(anglexz[evt])
            if(alt < 20):
                Az_alt_l20.Fill(az)
                Az_alt_l20_DA.Fill(anglexz[evt])
                Alt_alt_l20.Fill(alt)
                RA_alt_l20.Fill(ra)
                AltAz_alt_l20.append((az,alt))
                #AltAz_GU_alt_l20.append((anglexz[evt],alt))
                RAl20.append(RA)
                DECl20.append(DEC)
            if (alt > 20 and trkden[evt] <= 0.2):
                AltAz_L.append((az,alt))
                #AltAz_L.append((anglexz[evt],alt))
                Az_LD.Fill(az)
                Az_LD_DA.Fill(anglexz[evt])
                Alt_LD.Fill(alt)
                RA_LD.Fill(ra)
                RALD.append(RA)
                DECLD.append(DEC)
                hEvtsLD_Month.Fill( (year-2015)*12 + month-1 )
                hEvtsLD_Wrap.Fill( month-1 )
                SubRunInfo[ ( int(run_ei[evt]), int(subrun_ei[evt]) ) ]["nLD"] += 1
            if (alt > 20 and trkden[evt] > 0.2 and trkden[evt] <= 0.3):
                AltAz_M.append((az,alt))
                #AltAz_M.append((anglexz[evt],alt))
                Az_MD.Fill(az)
                Az_MD_DA.Fill(anglexz[evt])
                Alt_MD.Fill(alt) 
                RA_MD.Fill(ra)
                RAMD.append(RA)
                DECMD.append(DEC)
                hEvtsMD_Month.Fill( (year-2015)*12 + month-1 )
                hEvtsMD_Wrap.Fill( month-1 )
                SubRunInfo[ ( int(run_ei[evt]), int(subrun_ei[evt]) ) ]["nMD"] += 1
            if (alt > 20 and trkden[evt] > 0.3):
                AltAz_H.append((az,alt))
                #AltAz_H.append((anglexz[evt],alt))
                Az_HD.Fill(az)
                Az_HD_DA.Fill(anglexz[evt])
                Alt_HD.Fill(alt)
                RA_HD.Fill(ra)
                RAHD.append(RA)
                DECHD.append(DEC)
                hEvtsHD_Month.Fill( (year-2015)*12 + month-1 )
                hEvtsHD_Wrap.Fill( month-1 )
                SubRunInfo[ ( int(run_ei[evt]), int(subrun_ei[evt]) ) ]["nHD"] += 1
            if(alt > 20 and ntracks_n[evt] <=120):
                AltAz_LM.append((az,alt))
            if(alt > 20 and ntracks_n[evt] > 120):
                AltAz_HM.append((az,alt))

# ## Histograms

# In[5]:

# Azimuth distributions in subplots
# Create a canvas
#canvas0 = ROOT.TCanvas("canvas0", "canvas0", 1200, 800)
canvas0 = ROOT.TCanvas("canvas0", "canvas0", 600, 800)

# Divide the canvas into a 2x2 grid
#canvas0.Divide(2, 3)
canvas0.Divide(1, 3)

# Set line colors
#Az_BU.SetLineColor(ROOT.kOrange+2)
#Az_alt_l20.SetLineColor(ROOT.kRed)
Az_LD.SetLineColor(ROOT.kGreen+2)
Az_MD.SetLineColor(ROOT.kBlue)
Az_HD.SetLineColor(ROOT.kMagenta)

# Set titles
#Az_BU.SetTitle("Non-Uniform Showers")
#Az_alt_l20.SetTitle("All Showers (Alt < 20)")
Az_LD.SetTitle("Low Density Uniform Showers")
Az_MD.SetTitle("Medium Density Uniform Showers")
Az_HD.SetTitle("High Density Uniform Showers")

# Draw each histogram in a separate pad
#canvas0.cd(1)
#Az_alt_l20.Draw("hist")
#canvas0.cd(2)
#Az_BU.Draw("hist")
#canvas0.cd(3)
canvas0.cd(1)
Az_LD.Draw("hist")
canvas0.Update()  # Ensure the statistics box is created
statBox = Az_LD.FindObject("stats")
if statBox:
    statBox.SetX1NDC(0.6)  # Left boundary (near 0.7 to center it)
    statBox.SetX2NDC(0.75)  # Right boundary
    statBox.SetY1NDC(0.7)  # Bottom boundary
    statBox.SetY2NDC(0.9)  # Top boundary
    statBox.Draw()
#canvas0.cd(4)
canvas0.cd(2)
Az_MD.Draw("hist")
#canvas0.cd(5)
canvas0.cd(3)
Az_HD.Draw("hist")
canvas0.Draw()
# Save the canvas as an image
canvas0.SaveAs("./Plots/Azimuth_dist_uniform_selections.png")

# Uniform azimuth distribution
canvas3 = ROOT.TCanvas("canvas3", "canvas3", 800, 600)
Az_U.SetTitle("Uniform showers")
#Az_U_DA.SetTitle("Uniform showers (Detector Axis)")
Az_U.SetLineColor(ROOT.kRed)
Az_U.Draw("hist")
#Az_U_DA.Draw('same')
#Az_U_DA.SetLineColor(ROOT.kBlack)
canvas3.Draw()
canvas3.Update()
statBox = Az_U.FindObject("stats")
if statBox:
    statBox.SetX1NDC(0.55)  # Left boundary (near 0.7 to center it)
    statBox.SetX2NDC(0.75)  # Right boundary
    statBox.Draw()
canvas3.SaveAs("./Plots/UniformAzDist.png")

# Altitude distributions in subplots
# Create a canvas
canvas2 = ROOT.TCanvas("canvas2", "canvas2", 1200, 800)

# Divide the canvas into a 2x3 grid
canvas2.Divide(2, 3)

# Set line colors
Alt_BU.SetLineColor(ROOT.kOrange+7)
Alt_alt_l20.SetLineColor(ROOT.kRed)
Alt_LD.SetLineColor(ROOT.kGreen+2)
Alt_MD.SetLineColor(ROOT.kBlue)
Alt_HD.SetLineColor(ROOT.kMagenta)

# Set titles
Alt_BU.SetTitle("Non-Uniform Showers")
Alt_alt_l20.SetTitle("Uniform Showers (Alt < 20)")
Alt_LD.SetTitle("Low Density Uniform Showers")
Alt_MD.SetTitle("Medium Density Uniform Showers")
Alt_HD.SetTitle("High Density Uniform Showers")

# Draw each histogram in a separate pad
canvas2.cd(1)
Alt_BU.Draw("hist")
canvas2.cd(2)
Alt_alt_l20.Draw("hist")
canvas2.cd(3)
Alt_LD.Draw("hist")
canvas2.Update()  # Ensure the statistics box is created
statBox = Alt_LD.FindObject("stats")
if statBox:
    statBox.SetX1NDC(0.7)  # Left boundary (near 0.7 to center it)
    statBox.SetX2NDC(0.9)  # Right boundary
    #statBox.SetY1NDC(0.8)  # Bottom boundary
    #statBox.SetY2NDC(0.9)  # Top boundary
    statBox.Draw()
canvas2.cd(4)
Alt_MD.Draw("hist")
canvas2.cd(5)
Alt_HD.Draw("hist")
canvas2.Draw()
# Save the canvas as an image
canvas2.SaveAs("./Plots/Altitude_dist.png")


# RA distributions in one plot

#RA_BU.Scale(1.0 / RA_BU.Integral())
#RA_GU_alt_l20.Scale(1.0 / RA_GU_alt_l20.Integral())
#RA_LD.Scale(1.0 / RA_LD.Integral())
#RA_MD.Scale(1.0 / RA_MD.Integral())
#RA_HD.Scale(1.0 / RA_HD.Integral())

canvas1 = ROOT.TCanvas("canvas1","canvas1",800,600)
canvas1.cd()
# --- Style for axes and overall appearance ---
ROOT.gStyle.SetLineWidth(2)  # thicker axis/frame
canvas1.SetLeftMargin(0.12)
canvas1.SetBottomMargin(0.12)
canvas1.SetRightMargin(0.08)
canvas1.SetTopMargin(0.12)

#RA_BU.SetLineColor(ROOT.kOrange+7)
#RA_alt_l20.SetLineColor(ROOT.kRed)
#RA_LD.SetLineColor(ROOT.kGreen+2)
RA_LD.SetMarkerStyle(20)
#RA_LD.SetMarkerSize(1.2)
RA_LD.SetMarkerColor(ROOT.kGreen+2)
RA_LD.SetLineColor(ROOT.kGreen+2)
RA_LD.Draw("ep")
RA_MD.SetMarkerStyle(20)
#RA_MD.SetMarkerSize(1.2)
RA_MD.SetMarkerColor(ROOT.kBlue)
RA_MD.SetLineColor(ROOT.kBlue)
RA_MD.Draw("same ep")
RA_HD.SetMarkerStyle(20)
#RA_HD.SetMarkerSize(1.2)
RA_HD.SetMarkerColor(ROOT.kMagenta)
RA_HD.SetLineColor(ROOT.kMagenta)
RA_HD.Draw("same ep")
RA_LD.SetStats(0)
RA_MD.SetStats(0)
RA_HD.SetStats(0)
# Y-axis range and title
RA_LD.GetYaxis().SetRangeUser(0,2100)
#RA_LD.SetTitle("RA Distribution")
RA_LD.GetXaxis().SetTitle("Right Ascension (degrees)")
RA_LD.GetYaxis().SetTitle("Counts")
#RA_BU.Draw("hist")
#RA_alt_l20.Draw("same")
#RA_LD.Draw("hist")

for axis in [RA_LD.GetXaxis(), RA_LD.GetYaxis()]:
    axis.SetLabelFont(62)       # bold labels
    axis.SetLabelSize(0.05)     # bigger labels
    axis.SetTitleFont(62)       # bold title
    axis.SetTitleSize(0.05)     # bigger title
    axis.SetTitleOffset(1.1)    # adjust spacing

# Combine sum (clone)
RA_sum = RA_LD.Clone("RA_sum")
RA_sum.Add(RA_MD)
RA_sum.Add(RA_HD)
RA_sum.SetMarkerStyle(20)
#RA_sum.SetMarkerSize(1.2)
RA_sum.SetMarkerColor(ROOT.kBlack)
RA_sum.SetLineColor(ROOT.kBlack)
RA_sum.SetStats(0)
RA_sum.Draw("same ep")

# Legend
legend = ROOT.TLegend(0.65, 0.75, 0.88, 0.88)
legend.SetBorderSize(0)
legend.SetFillStyle(0)
legend.SetTextFont(62)   # bold legend font
legend.SetTextSize(0.035)
legend.AddEntry(RA_LD, "Low density", "ep")
legend.AddEntry(RA_MD, "Medium density", "ep")
legend.AddEntry(RA_HD, "High density", "ep")
legend.AddEntry(RA_sum, "Sum", "ep")
legend.Draw()

canvas1.cd()
prelim = ROOT.TLatex()
prelim.SetNDC()                # normalized coordinates
prelim.SetTextColor(ROOT.kBlue)
prelim.SetTextSize(0.07)       # adjust size
prelim.SetTextAlign(32)        # right-top alignment
prelim.DrawLatex(0.93, 0.95, "NOvA Preliminary")

canvas1.Update()
#canvas1.Draw()
canvas1.SaveAs("./Plots/RA_dist.png")
canvas1.SaveAs("./Plots/RA_dist.eps")
canvas1.SaveAs("./Plots/RA_dist.pdf")
canvas1.SaveAs("./Plots/RA_dist.root")

canvasRate = ROOT.TCanvas("canvasRate","Muon Rate vs RA",800,600)
canvasRate.cd()

# Making percentage histogram

def make_deviation_hist(hist):
    # compute mean bin content (average number of events per bin)
    mean_content = hist.GetEntries() / hist.GetNbinsX()
    dev_hist = hist.Clone(hist.GetName() + "_dev")
    dev_hist.SetTitle(";Right Ascension (degrees);Deviation from mean (%)")
    for i in range(1, dev_hist.GetNbinsX()+1):
        bin_content = hist.GetBinContent(i)
        if mean_content != 0:
            new_content = 100.0 * (bin_content - mean_content) / mean_content
            #new_content = bin_content / mean_content
        else:
            new_content = 0
        dev_hist.SetBinContent(i, new_content)
        # scale errors into percentage as well
        bin_err = hist.GetBinError(i)
        dev_hist.SetBinError(i, 100.0 * bin_err / mean_content if mean_content != 0 else 0)
        #dev_hist.SetBinError(i, bin_err / mean_content if mean_content != 0 else 0)
    return dev_hist

# make transformed histograms
RA_LD_dev  = make_deviation_hist(RA_LD)
RA_MD_dev  = make_deviation_hist(RA_MD)
RA_HD_dev  = make_deviation_hist(RA_HD)
RA_sum_dev = make_deviation_hist(RA_sum)

# Style histograms
for hist, color in zip([RA_LD_dev, RA_MD_dev, RA_HD_dev, RA_sum_dev],
        [ROOT.kGreen+2, ROOT.kBlue, ROOT.kMagenta, ROOT.kBlack]):
    hist.SetMarkerStyle(20)
    hist.SetMarkerColor(color)
    hist.SetLineColor(color)
    hist.SetLineWidth(2)
    hist.SetStats(0)

# Axis titles and range
RA_LD_dev.GetXaxis().SetTitle("Right Ascension (degrees)")
RA_LD_dev.GetYaxis().SetTitle("(#Delta N / <N>) %")
RA_LD_dev.GetYaxis().SetRangeUser(-50, 50)

# Draw histograms
RA_LD_dev.Draw("EP")
RA_MD_dev.Draw("same EP")
RA_HD_dev.Draw("same EP")
RA_sum_dev.Draw("same EP")

# Legend
legend = ROOT.TLegend(0.65, 0.75, 0.88, 0.88)
legend.SetBorderSize(0)
legend.SetFillStyle(0)
legend.SetTextFont(62)
legend.SetTextSize(0.035)
legend.AddEntry(RA_LD_dev, "Low density", "ep")
legend.AddEntry(RA_MD_dev, "Medium density", "ep")
legend.AddEntry(RA_HD_dev, "High density", "ep")
legend.AddEntry(RA_sum_dev, "Sum", "ep")
legend.Draw()

# NOvA Preliminary text
prelim = ROOT.TLatex()
prelim.SetNDC()
prelim.SetTextColor(ROOT.kBlue)
prelim.SetTextSize(0.07)
prelim.SetTextAlign(32)   # right-top
prelim.DrawLatex(0.93, 0.95, "NOvA Preliminary")

canvasRate.Update()
canvasRate.SaveAs("./Plots/RA_EventPrcnt.png")
canvasRate.SaveAs("./Plots/RA_EventPrcnt.eps")
canvasRate.SaveAs("./Plots/RA_EventPrcnt.pdf")
canvasRate.SaveAs("./Plots/RA_EventPrcnt.root")

# New: Individual canvases
# -------------------------
def draw_hist_on_canvas(hist, title, fname):
    canvas = ROOT.TCanvas("c_" + hist.GetName(), title, 800, 600)
    hist.Draw("EP")

    # Axis titles and range
    hist.GetXaxis().SetTitle("Right Ascension (degrees)")
    hist.GetYaxis().SetTitle("(#Delta N / <N>) %")
    hist.GetYaxis().SetRangeUser(-50, 50)

    # NOvA Preliminary text
    prelim = ROOT.TLatex()
    prelim.SetNDC()
    prelim.SetTextColor(ROOT.kBlue)
    prelim.SetTextSize(0.07)
    prelim.SetTextAlign(32)
    prelim.DrawLatex(0.93, 0.95, "NOvA Preliminary")

    canvas.Update()
    canvas.SaveAs(f"./Plots/{fname}.png")
    canvas.SaveAs(f"./Plots/{fname}.eps")
    canvas.SaveAs(f"./Plots/{fname}.pdf")
    canvas.SaveAs(f"./Plots/{fname}.root")

# Call for each histogram
draw_hist_on_canvas(RA_LD_dev,  "Muon Rate vs RA (Low density)",   "RA_EventPrcnt_LD")
draw_hist_on_canvas(RA_MD_dev,  "Muon Rate vs RA (Medium density)","RA_EventPrcnt_MD")
draw_hist_on_canvas(RA_HD_dev,  "Muon Rate vs RA (High density)",  "RA_EventPrcnt_HD")
draw_hist_on_canvas(RA_sum_dev, "Muon Rate vs RA (Sum)",           "RA_EventPrcnt_Sum")

# Multiplicity distribution
canvasA = ROOT.TCanvas("canvasA","canvasA",800,600)
Multiplicity.Draw("hist")  
canvasA.Draw()
canvasA.SaveAs("./Plots/multiplicity.png")

# Vertical asymmetry
canvasB = ROOT.TCanvas("canvasB","canvasB",800,600)
VA.Draw("hist")
canvasB.Draw()
canvasB.SaveAs("./Plots/va.png")

# Track density
canvasC = ROOT.TCanvas("canvasC","canvasC",800,600)
TrackDensity.Draw("hist")  
canvasC.Draw()
canvasC.SaveAs("./Plots/trkden.png")

# Track density vs Multiplicity
canvasD = ROOT.TCanvas("canvasD", "canvasD", 800, 600)
DM.Draw('colz')
canvasD.Draw()
canvasD.SaveAs("./Plots/DensityVsMultiplicity.png")

# Cross sectional area vs Alt, Az
canvasE = ROOT.TCanvas("canvasE", "canvasE", 900, 600)
XsecAltAz.Draw('colz')
canvasE.Draw()
canvasE.SaveAs("./Plots/XsecAltAz.png")

# Polar Plots for multiplicity
fig, axs = plt.subplots(2, 2, subplot_kw=dict(polar=True), figsize=(12, 12))  # Create a 2x2 grid of polar subplots

# Plotting AltAzBU
for azbu, altbu in AltAz_BU:
    theta0 = np.radians(azbu)
    r0 = altbu
    axs[0, 0].plot(theta0, r0, marker='o', color='orange', markersize=2)
axs[0, 0].set_title('Non-Uniform Showers')

# Plotting AltAzl20
for azl20, altl20 in AltAz_alt_l20:
    theta1 = np.radians(azl20)
    r1 = altl20
    axs[0, 1].plot(theta1, r1, marker='o', color='green', markersize=2)
axs[0, 1].set_title('Uniform Showers (Alt < 20)')

# Plotting AltAzLM
for azlm, altlm in AltAz_LM:
    theta2 = np.radians(azlm)
    r2 = altlm
    axs[1, 0].plot(theta2, r2, marker='o', color='blue', markersize=2)
axs[1, 0].set_title('LM Showers')

# Plotting AltAzHM
for azhm, althm in AltAz_HM:
    theta3 = np.radians(azhm)
    r3 = althm
    axs[1, 1].plot(theta3, r3, marker='o', color='magenta', markersize=2)
axs[1, 1].set_title('HM Showers')

for subplot in axs.flat:
    subplot.set_theta_zero_location('N')  # Set the azimuth zero to the north
    subplot.set_theta_direction(-1)  # Set the azimuth direction to clockwise
    subplot.set_rlabel_position(60)  # Set radial label position
    subplot.set_rlim(91 , 1)
    subplot.grid(True)
    
plt.savefig("./Plots/AltitudeAzimuth.png")
plt.close()


# Polar plots for density
#fig, axs = plt.subplots(2, 3, subplot_kw=dict(polar=True), figsize=(12, 12))


# Plotting AltAzl20
#for azl20, altl20 in AltAz_alt_l20:
#    theta1 = np.radians(azl20)
#    r1 = altl20
#    axs[0,0].plot(theta1, r1, marker='o', color='red', markersize=2)
#axs[0,0].set_title('All Showers (Alt < 20)')

# Plotting AltAzBU
#for azbu, altbu in AltAz_BU:
#    theta0 = np.radians(azbu)
#    r0 = altbu
#    axs[0,1].plot(theta0, r0, marker='o', color='orange', markersize=2)
#axs[0,1].set_title('Non-Uniform Showers')

# Plotting AltAz_L
#for azl, altl in AltAz_L:
#    theta0 = np.radians(azl)
#    r0 = altl
#    axs[1,0].plot(theta0, r0, marker='o', color='green', markersize=2)
#axs[1,0].set_title('Low density showers')

# Plotting AltAz_M
#for azm, altm in AltAz_M:
#    theta1 = np.radians(azm)
#    r1 = altm
#    axs[1,1].plot(theta1, r1, marker='o', color='blue', markersize=2)
#axs[1,1].set_title('Medium density showers')

# Plotting AltAz_H
#for azh, alth in AltAz_H:
#    theta2 = np.radians(azh)
#    r2 = alth
#    axs[1,2].plot(theta2, r2, marker='o', color='magenta', markersize=2)
#axs[1,2].set_title('High density showers')

#for row in axs:
#    for subplot in row:
#        # Skip the empty subplot at axs[0, 2]
#       if subplot == axs[0, 2]:  
#            subplot.axis('off')  # Turn off the axis for the empty subplot
#            continue
        #if subplot is not None:  # Skip the empty subplot
#        subplot.set_theta_zero_location('N')  # Set the azimuth zero to the north
#        subplot.set_theta_direction(-1)  # Set the azimuth direction to clockwise
#        subplot.set_rlabel_position(60)  # Set radial label position
#        subplot.set_rlim(91, 1)  # Invert the radial axis
#        subplot.grid(True)
#plt.savefig("./Plots/AltitudeAzimuth_density.png", bbox_inches='tight')
#plt.close()

# Polar plots for density (only two rejected categories)
fig, axs = plt.subplots(1, 2, subplot_kw=dict(polar=True), figsize=(12, 6))

# Plotting AltAzl20 (All Showers with Alt < 20)
for azl20, altl20 in AltAz_alt_l20:
    theta1 = np.radians(azl20)
    r1 = altl20
    axs[0].plot(theta1, r1, marker='o', color='red', markersize=2)
axs[0].set_title('All Showers (Alt < 20)')

# Plotting AltAzBU (Non-Uniform Showers)
for azbu, altbu in AltAz_BU:
    theta0 = np.radians(azbu)
    r0 = altbu
    axs[1].plot(theta0, r0, marker='o', color='orange', markersize=2)
axs[1].set_title('Non-Uniform Showers')

# Adjust polar subplot settings
for subplot in axs:
    subplot.set_theta_zero_location('N')   # Azimuth zero at North
    subplot.set_theta_direction(-1)        # Clockwise
    subplot.set_rlabel_position(60)        # Radial label position
    subplot.set_rlim(91, 1)                # Invert radial axis
    subplot.grid(True)

plt.tight_layout()
plt.savefig("./Plots/AltitudeAzimuth_rejected.png", bbox_inches='tight')
plt.close()

# Polar plots for density
fig, axs = plt.subplots(1, 3, subplot_kw=dict(polar=True), figsize=(18, 6))

# Plotting AltAz_L
for azl, altl in AltAz_L:
    theta0 = np.radians(azl)
    r0 = altl
    axs[0].plot(theta0, r0, marker='o', color='green', markersize=2)
axs[0].set_title('Low density uniform Showers')

# Plotting AltAz_M (Medium density showers) on second column
for azm, altm in AltAz_M:
    theta1 = np.radians(azm)
    r1 = altm
    axs[1].plot(theta1, r1, marker='o', color='blue', markersize=2)
axs[1].set_title('Medium density uniform showers')

# Plotting AltAz_H (High density showers) on third column
for azh, alth in AltAz_H:
    theta2 = np.radians(azh)
    r2 = alth
    axs[2].plot(theta2, r2, marker='o', color='magenta', markersize=2)
axs[2].set_title('High density uniform showers')

# Adjust polar subplot settings
for subplot in axs:
    subplot.set_theta_zero_location('N')   # Azimuth zero at North
    subplot.set_theta_direction(-1)        # Clockwise
    subplot.set_rlabel_position(65)        # Radial label position
    subplot.set_rlim(91, 1)                # Invert radial axis
    subplot.grid(True)
    subplot.title.set_fontsize(18)
    subplot.title.set_fontweight('bold')
    subplot.tick_params(labelsize=12)  # larger tick labels
    subplot.xaxis.set_tick_params(labelsize=12, width=1.5)  # azimuth ticks
    subplot.yaxis.set_tick_params(labelsize=12, width=1.5)  # radial ticks
    
    # Thicker gridlines
    for gridline in subplot.yaxis.get_gridlines() + subplot.xaxis.get_gridlines():
        gridline.set_linewidth(1.2)  # increase grid line thickness
        gridline.set_alpha(0.7)   

plt.tight_layout()
plt.savefig("./Plots/AltitudeAzimuth_density_uniform.png", bbox_inches='tight')
plt.savefig("./Plots/AltitudeAzimuth_density_uniform.eps", bbox_inches='tight')
plt.savefig("./Plots/AltitudeAzimuth_density_uniform.pdf", bbox_inches='tight')
plt.close()


# Polar plot for AltAz_H
fig = plt.figure(figsize=(8, 8))
ax = plt.subplot(polar=True)

# Plotting AltAzH
for azh, alth in AltAz_H:
    theta = np.radians(azh)  # Convert azimuth to radians
    r = alth  # Altitude
    ax.plot(theta, r, marker='o', color='magenta', markersize=2)

# Customize the plot
ax.set_title('High density Showers', fontsize=16)
ax.set_theta_zero_location('N')  # Set azimuth zero to North
ax.set_theta_direction(-1)  # Set azimuth direction clockwise
ax.set_rlabel_position(60)  # Position radial labels
ax.set_rlim(91, 1)  # Invert radial axis
ax.grid(True)

# Save and display the plot
plt.savefig("./Plots/HighDensityAltitudeAzimuth.png")
plt.close()

# Skymaps
fig = plt.figure(figsize=(12, 12))

nside = 32

# BU
pix_bu = hp.ang2pix(nside, np.radians(DECBU), np.radians(RABU))
map_bu = np.bincount(pix_bu, minlength=hp.nside2npix(nside))
hp.mollview(map_bu, title='Non-uniform Shower', cmap='plasma', sub=(3, 2, 1))
hp.graticule()

# GU alt < 20
pix_gu = hp.ang2pix(nside, np.radians(DECl20), np.radians(RAl20))
map_gu = np.bincount(pix_gu, minlength=hp.nside2npix(nside))
hp.mollview(map_gu, title='All Shower (Alt < 20)', cmap='plasma', sub=(3, 2, 2))
hp.graticule()

# LD 
pix_ld = hp.ang2pix(nside, np.radians(DECLD), np.radians(RALD))
map_ld = np.bincount(pix_ld, minlength=hp.nside2npix(nside))
hp.mollview(map_ld, title='Uniform Shower (LD)', cmap='plasma', sub=(3, 2, 3))
hp.graticule()

# MD
pix_md = hp.ang2pix(nside, np.radians(DECMD), np.radians(RAMD))
map_md = np.bincount(pix_md, minlength=hp.nside2npix(nside))
hp.mollview(map_md, title='Uniform Shower (MD)', cmap='plasma', sub=(3, 2, 4))
hp.graticule()

# HD 
pix_hd = hp.ang2pix(nside, np.radians(DECHD), np.radians(RAHD))
map_hd = np.bincount(pix_hd, minlength=hp.nside2npix(nside))
hp.mollview(map_hd, title='Uniform Shower (HD)', cmap='plasma', sub=(3, 2, 5))
hp.graticule()

plt.savefig("./Plots/skymaps.png")
plt.close()


# Time series plots

labelsM = []
for i in range(108):
    month = (10 + i - 1) % 12 + 1
    year = 2015 + (10 + i - 1) // 12
    label = f"{month_abbr[month]}-{year}"
    labelsM.append(label)


month_abbr = ["Jan", "Feb", "Mar", "Apr", "May", "Jun", 
              "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]


hRateDist_NU = ROOT.TH1D("hRateDist_NU","Shower Rates;Rate (h^{-1})",100,0,100)
hRateDist_LD = ROOT.TH1D("hRateDist_LD","Shower Rates;Rate (h^{-1})",100,0,100)
hRateDist_MD = ROOT.TH1D("hRateDist_MD","Shower Rates;Rate (h^{-1})",100,0,100)
hRateDist_HD = ROOT.TH1D("hRateDist_HD","Shower Rates;Rate (h^{-1})",100,0,100)

hRateDist_NU.Reset()
hRateDist_LD.Reset()
hRateDist_MD.Reset()
hRateDist_HD.Reset()

for sr in SubRunInfo.values():
    rate_NU = sr["nNU"]/(sr["duration"]/3600.)
    rate_LD = sr["nLD"]/(sr["duration"]/3600.)
    rate_MD = sr["nMD"]/(sr["duration"]/3600.)
    rate_HD = sr["nHD"]/(sr["duration"]/3600.)
    time = sr["time"]
    hRateDist_NU.Fill(rate_NU)
    hRateDist_LD.Fill(rate_LD)
    hRateDist_MD.Fill(rate_MD)
    hRateDist_HD.Fill(rate_HD)

canvasT = ROOT.TCanvas("canvasT","canvasT",800,600)
hRateDist_NU.SetLineColor(ROOT.kRed)
hRateDist_LD.SetLineColor(ROOT.kGreen+2)
hRateDist_MD.SetLineColor(ROOT.kBlue)
hRateDist_HD.SetLineColor(ROOT.kMagenta)

hRateDist_NU.Draw("hist")
hRateDist_LD.Draw("same hist")
hRateDist_MD.Draw("same hist")
hRateDist_HD.Draw("same hist")
canvasT.SetLogy()

canvasT.Draw()
canvasT.SaveAs("./Plots/hRateDist.png")
######################################################
canvasT1 = ROOT.TCanvas("canvasT1","canvasT1",800,600)
hEvts_Month.SetLineColor(ROOT.kBlack)
hEvtsNU_Month.SetLineColor(ROOT.kRed)
hEvtsLD_Month.SetLineColor(ROOT.kGreen+2)
hEvtsMD_Month.SetLineColor(ROOT.kBlue)
hEvtsHD_Month.SetLineColor(ROOT.kMagenta)

hEvts_Month.Draw("hist E0")
hEvtsNU_Month.Draw("same hist E0")
hEvtsLD_Month.Draw("same hist E0")
hEvtsMD_Month.Draw("same hist E0")
hEvtsHD_Month.Draw("same hist E0")

canvasT1.Draw()
canvasT1.SaveAs("./Plots/hEvts_Month.png")
######################################################
canvasT2 = ROOT.TCanvas("canvasT2","canvasT2",800,600)
hEvts_Wrap.SetLineColor(ROOT.kBlack)
hEvtsNU_Wrap.SetLineColor(ROOT.kRed)
hEvtsLD_Wrap.SetLineColor(ROOT.kGreen+2)
hEvtsMD_Wrap.SetLineColor(ROOT.kBlue)
hEvtsHD_Wrap.SetLineColor(ROOT.kMagenta)

hEvts_Wrap.GetYaxis().SetRangeUser(0,6700)

hEvts_Wrap.Draw("hist E0")
hEvtsNU_Wrap.Draw("same hist E0")
hEvtsLD_Wrap.Draw("same hist E0")
hEvtsMD_Wrap.Draw("same hist E0")
hEvtsHD_Wrap.Draw("same hist E0")

canvasT2.Draw()
canvasT2.SaveAs("./Plots/hEvts_Wrap.png")
######################################################
canvasT3 = ROOT.TCanvas("canvasT3","canvasT3",800,600)
hLiveTime_Month.SetLineColor(ROOT.kBlack)

for i in range(1,97):
    hLiveTime_Month.SetBinError(i,0)

hLiveTime_Month.Draw("hist")
for i in range(108):
    label = labelsM[i] if i % 6 == 0 else ""
    hLiveTime_Month.GetXaxis().SetBinLabel(i + 1, label)

canvasT3.Draw()
canvasT3.SaveAs("./Plots/hLiveTime_Month.png")
######################################################
canvasT4 = ROOT.TCanvas("canvasT4","canvasT4",800,600)

#hRatio_Month = hEvts_Month.Clone()
#hRatio_Month.Divide(hLiveTime_Month*(1/3600.))
#hRatio_Month.SetTitle("Event Rate (per hour);Month;Rate (h^{-1})")

#hRatioNU_Month = hEvtsNU_Month.Clone()
#hRatioNU_Month.Divide(hLiveTime_Month*(1/3600.))
#hRatioNU_Month.SetTitle("Event Rate (per hour);Month;Rate (h^{-1})")
#hRatioNU_Month.SetLineColor(ROOT.kRed)

hRatioLD_Month = hEvtsLD_Month.Clone()
hRatioLD_Month.Divide(hLiveTime_Month*(1/3600.))
hRatioLD_Month.SetTitle("Event Rate (per hour);Month;Rate (h^{-1})")
hRatioLD_Month.SetLineColor(ROOT.kGreen+2)

hRatioMD_Month = hEvtsMD_Month.Clone()
hRatioMD_Month.Divide(hLiveTime_Month*(1/3600.))
hRatioMD_Month.SetTitle("Event Rate (per hour);Month;Rate (h^{-1})")
hRatioMD_Month.SetLineColor(ROOT.kBlue)

hRatioHD_Month = hEvtsHD_Month.Clone()
hRatioHD_Month.Divide(hLiveTime_Month*(1/3600.))
hRatioHD_Month.SetTitle("Event Rate (per hour);Month;Rate (h^{-1})")
hRatioHD_Month.SetLineColor(ROOT.kMagenta)

hRationU_Month = hRatioLD_Month.Clone()
hRationU_Month.Add(hRatioMD_Month)
hRationU_Month.Add(hRatioHD_Month)
hRationU_Month.SetLineColor(ROOT.kBlack)

hRatioLD_Month.GetYaxis().SetRangeUser(0,3.5)
#hRatio_Month.Draw("E0")
#hRatioNU_Month.Draw("same E0")
hRatioLD_Month.Draw("E0")
hRatioMD_Month.Draw("same E0")
hRatioHD_Month.Draw("same E0")
hRationU_Month.Draw("same E0")

for i in range(108):
    label = labelsM[i] if i % 6 == 0 else ""
    hRatioLD_Month.GetXaxis().SetBinLabel(i + 1, label)

canvasT4.Draw()
canvasT4.SaveAs("./Plots/hRatio_Month.png")
######################################################
canvasT5 = ROOT.TCanvas("canvasT5", "canvasT5", 800, 600)

# Clone the original month histograms and rebin them
hRatio_Year = hEvts_Month.Clone()
hRatio_Year.Divide(hLiveTime_Month*(1/3600.))
hRatio_Year.Rebin(12)
hRatio_Year.SetTitle("Event Rate (per year);Year;Rate (h^{-1})")

hRatioNU_Year = hEvtsNU_Month.Clone()
hRatioNU_Year.Divide(hLiveTime_Month*(1/3600.))
hRatioNU_Year.Rebin(12)
hRatioNU_Year.SetLineColor(ROOT.kRed)

hRatioLD_Year = hEvtsLD_Month.Clone()
hRatioLD_Year.Divide(hLiveTime_Month*(1/3600.))
hRatioLD_Year.Rebin(12)
hRatioLD_Year.SetLineColor(ROOT.kGreen+2)

hRatioMD_Year = hEvtsMD_Month.Clone()
hRatioMD_Year.Divide(hLiveTime_Month*(1/3600.))
hRatioMD_Year.Rebin(12)
hRatioMD_Year.SetLineColor(ROOT.kBlue)

hRatioHD_Year = hEvtsHD_Month.Clone()
hRatioHD_Year.Divide(hLiveTime_Month*(1/3600.))
hRatioHD_Year.Rebin(12)
hRatioHD_Year.SetLineColor(ROOT.kMagenta)

# Set Y-axis range
hRatio_Year.GetYaxis().SetRangeUser(0, 140)

# Draw the histograms
hRatio_Year.Draw("hist E0")
hRatioNU_Year.Draw("same hist E0")
hRatioLD_Year.Draw("same hist E0")
hRatioMD_Year.Draw("same hist E0")
hRatioHD_Year.Draw("same hist E0")

# Add year labels to the x-axis (every 12 bins)
labels = []
for i in range(9):  # For 9 years (since you have 108 months)
    year = 2015 + i
    label = str(year)
    labels.append(label)

for i in range(9):  # Set the labels for each year
    hRatio_Year.GetXaxis().SetBinLabel(i + 1, labels[i])

# Draw the canvas
canvasT5.Draw()
canvasT5.SaveAs("./Plots/hRatio_Year.png")
######################################################
canvasT6 = ROOT.TCanvas("canvasT6","canvasT6",800,600)
canvasT6.SetLeftMargin(0.12)
canvasT6.SetBottomMargin(0.12)
canvasT6.SetRightMargin(0.10)
canvasT6.SetTopMargin(0.12)

#hRatio_Wrap = hEvts_Wrap.Clone()
#hRatio_Wrap.Divide(hLiveTime_Wrap*(1/3600.))
#hRatio_Wrap.SetTitle("Event Rate (per hour);Wrap;Rate (h^{-1})")

#hRatioNU_Wrap = hEvtsNU_Wrap.Clone()
#hRatioNU_Wrap.Divide(hLiveTime_Wrap*(1/3600.))
#hRatioNU_Wrap.SetTitle("Event Rate (per hour);Wrap;Rate (h^{-1})")
#hRatioNU_Wrap.SetLineColor(ROOT.kRed)

hRatioLD_Wrap = hEvtsLD_Wrap.Clone()
hRatioLD_Wrap.Divide(hLiveTime_Wrap*(1/3600.))
hRatioLD_Wrap.SetTitle(" ;Month;Event rate (h^{-1})")
hRatioLD_Wrap.SetMarkerStyle(20)
#hRatioLD_Wrap.SetMarkerSize(1.2)
hRatioLD_Wrap.SetMarkerColor(ROOT.kGreen+2)
hRatioLD_Wrap.SetLineColor(ROOT.kGreen+2)
hRatioLD_Wrap.SetStats(0)

hRatioMD_Wrap = hEvtsMD_Wrap.Clone()
hRatioMD_Wrap.Divide(hLiveTime_Wrap*(1/3600.))
hRatioMD_Wrap.SetTitle(" ;Month;Event rate (h^{-1})")
hRatioMD_Wrap.SetMarkerStyle(20)
#hRatioMD_Wrap.SetMarkerSize(1.2)
hRatioMD_Wrap.SetMarkerColor(ROOT.kBlue)
hRatioMD_Wrap.SetLineColor(ROOT.kBlue)
hRatioMD_Wrap.SetStats(0)

hRatioHD_Wrap = hEvtsHD_Wrap.Clone()
hRatioHD_Wrap.Divide(hLiveTime_Wrap*(1/3600.))
hRatioHD_Wrap.SetTitle(" ;Month;Event rate (h^{-1})")
hRatioHD_Wrap.SetMarkerStyle(20)
#hRatioHD_Wrap.SetMarkerSize(1.2)
hRatioHD_Wrap.SetMarkerColor(ROOT.kMagenta)
hRatioHD_Wrap.SetLineColor(ROOT.kMagenta)
hRatioHD_Wrap.SetStats(0)

hRatioSum_Wrap = hRatioLD_Wrap.Clone()
hRatioSum_Wrap.Add(hRatioMD_Wrap)
hRatioSum_Wrap.Add(hRatioHD_Wrap)
hRatioSum_Wrap.SetMarkerStyle(20)
#hRatioSum_Wrap.SetMarkerSize(1.2)
hRatioSum_Wrap.SetMarkerColor(ROOT.kBlack)
hRatioSum_Wrap.SetLineColor(ROOT.kBlack)
hRatioSum_Wrap.SetStats(0)

# --- Axis style ---
for axis in [hRatioLD_Wrap.GetXaxis(), hRatioLD_Wrap.GetYaxis()]:
    axis.SetLabelFont(62)   # bold
    axis.SetLabelSize(0.05)
    axis.SetTitleFont(62)   # bold
    axis.SetTitleSize(0.05)
    axis.SetTitleOffset(1.2)


hRatioLD_Wrap.GetYaxis().SetRangeUser(0,2.5)
#hRatio_Wrap.Draw("hist E0")
#hRatioNU_Wrap.Draw("same hist E0")
hRatioLD_Wrap.Draw("ep")
hRatioMD_Wrap.Draw("same ep")
hRatioHD_Wrap.Draw("same ep")
hRatioSum_Wrap.Draw("same ep")

for i in range(12):
    hRatioLD_Wrap.GetXaxis().SetBinLabel(i + 1, month_abbr[i])

canvasT6.cd()

# Legend
legendT = ROOT.TLegend(0.65, 0.75, 0.88, 0.88)
legendT.SetBorderSize(0)
legendT.SetFillStyle(0)
legendT.SetTextFont(62)   # bold legend font
legendT.SetTextSize(0.035)
legendT.AddEntry(hRatioLD_Wrap, "Low density", "ep")
legendT.AddEntry(hRatioMD_Wrap, "Medium density", "ep")
legendT.AddEntry(hRatioHD_Wrap, "High density", "ep")
legendT.AddEntry(hRatioSum_Wrap, "Sum", "ep")
legendT.Draw()

prelim = ROOT.TLatex()
prelim.SetNDC()                # normalized coordinates
prelim.SetTextColor(ROOT.kBlue)
prelim.SetTextSize(0.07)       # adjust size
prelim.SetTextAlign(32)        # right-top alignment
prelim.DrawLatex(0.93, 0.95, "NOvA Preliminary")

canvasT6.Update()

#canvasT6.Draw()

canvasT6.SaveAs("./Plots/hRatio_wrap.png")
canvasT6.SaveAs("./Plots/hRatio_wrap.eps")
canvasT6.SaveAs("./Plots/hRatio_wrap.pdf")
canvasT6.SaveAs("./Plots/hRatio_wrap.root")

# In[16]:

canvasT7 = ROOT.TCanvas("canvas","canvas",800,600)
hSRDuration.Draw()
canvasT7.Draw()
canvasT7.SaveAs("./Plots/hSRDuration.png")

outfile = TFile("AirshowerHistograms.root", RECREATE)
canvasT6.Write("Timeseries")
canvasRate.Write("MuonRatevsRA")
canvas1.Write("RADist")
outfile.Close()
