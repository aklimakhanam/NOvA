#!/usr/bin/env python
# coding: utf-8

# In[1]:


# Primary imports and constants
import ROOT
import numpy as np
import matplotlib.pyplot as plot
import datetime
import pandas as pd
import os
import glob
import math
import io
from PIL import Image, ImageOps
import iy10_ultrashower_utilities as shower
import astropy.units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation, AltAz 
px = 1/plot.rcParams['figure.dpi']  # pixel in inches

# Functions
def ProcessTimestamp( evtTimestamp ):
    evtNanosec = evtTimestamp & 0xFFFFFFFF
    evtTimestamp = (evtTimestamp >> 32) & 0xFFFFFFFF
    evtTime = datetime.datetime.fromtimestamp(evtTimestamp + evtNanosec*1e-9)
    return evtTime

# Custom rectangular markers
from matplotlib.path import Path
verts = [
   (0., 0.),  # left, bottom
   (0., 3.9),  # left, top
   (6, 3.9),  # right, top
   (6, 0.),  # right, bottom
   (0., 0.),  # back to left, bottom
]
codes = [
    Path.MOVETO, #begin drawing
    Path.LINETO, #straight line
    Path.LINETO,
    Path.LINETO,
    Path.CLOSEPOLY, #close shape. This is not required for this shape but is "good form"
]
recMarker = Path(verts, codes)

ratio = ROOT.TH1D("TBR","Top-Bottom symmetry of events; Ratio; Counts", 50, -0.5, 0.5)
tracks = ROOT.TH1D("nTracks","Tracks distribution; nTracks; Counts", 150, 0, 300)
tracksr = ROOT.TH1D("nTracksR","Tracks distribution (rejected); nTracks (rejected); Counts", 150, 0, 300)

#Reading subrun_tree and event_tree
file_path = glob.glob("/home/aklima/MSAnalysis/Codes/RootFiles/*.root")
#file_path = glob.glob("/home/dwwhitti/nova/airshower/N23-07-05/*.root")
print(len(file_path))
for n in range(0,500): #Run this batchwise because the kernel dies after 58 files, create a self controlled loop that runs a batch of files each time
    base_name = os.path.basename(file_path[n])
    #print(base_name)
    name_split = os.path.splitext(base_name)
    #print(name_split)
    file_base_name = name_split[0]
    print(file_base_name, n)
    inFile = ROOT.TFile.Open(file_path[n],"READ")
    eventTree = inFile.Get("ultrashowerimage/Events")
    subrunTree = inFile.Get("ultrashowerimage/SubrunInfo")
    AsymXZ = []
    AsymYZ = []
    occupancyXZ = []
    occupancyYZ = []
    SubrunInfo = []
    Subrun_start = []
    Subrun_stop = []
    Run_number = []
    Subrun_number = []
    Run = []
    Subrun = []
    Event = []
    S_start = []
    S_stop = []
    AngleXZ = []
    StDevXZ = []
    AngleYZ = []
    StDevYZ = []
    Alt = []
    Az = []
    R = []
    nTracks = []
    nTracks_new = []
    nTracks_all = []
    TrkDen = []
    A_xsec = []
    Evt_timestamp = []
    Evt_timestamp = []
    CombinedList = []
    # Calculating number of tracks
    #angle, mean, stdev = np.loadtxt("mean_stdev.txt", unpack=True)
    angle, mean = np.loadtxt("mean_stdev1.txt", unpack=True)
    azs, alts, nhits = np.loadtxt("AzAltHits.txt", unpack=True)
    azall, altall, nhitsall = np.loadtxt("AzAltHits_all.txt", unpack=True)

    for index in range(subrunTree.GetEntries()):
        subrunTree.GetEntry(index)
        run = getattr(subrunTree, "run")
        Run_number.append(run)
        subrun = getattr(subrunTree, "subrun")
        Subrun_number.append(subrun)
        srStart = getattr(subrunTree,"srStart")
        Subrun_start.append(srStart)
        srStop = getattr(subrunTree, "srStop")
        Subrun_stop.append(srStop)
        SubrunInfo.append((run,subrun,srStart,srStop))

    np.array(SubrunInfo)
    header_sinfo = ['Run','Subrun','Subrun_start','Subrun_stop']
    file_name = f"/home/aklima/MSAnalysis/Codes/TextAll/{file_base_name}_SubrunInfo.txt"
    np.savetxt(file_name, SubrunInfo, delimiter='\t', fmt='%d', header='\t'.join(header_sinfo), comments='#')

    for index in range(eventTree.GetEntries()):
        eventTree.GetEntry(index)
        hitsXZ = getattr(eventTree,"hitsXZ")
        hitsYZ = getattr(eventTree,"hitsYZ")
        run = getattr(eventTree,"run")
        subrun = getattr(eventTree,"subrun")
        event = getattr(eventTree,"event")
        timestamp = getattr(eventTree,"timestamp")
        
        #unique_subruns = set()
        for i in range(len(Subrun_number)):
            if Subrun_number[i]==subrun and Run_number[i]==run:
                #unique_subruns.add(Subrun_number[i])
                S_start.append(Subrun_start[i])
                S_stop.append(Subrun_stop[i])
            else:
                continue
        
        Run.append(run)
        Subrun.append(subrun)
        Event.append(event)
        
        arrayXZ = np.stack( ( hitsXZ.GetX(), hitsXZ.GetY() ), axis=1 )
        arrayYZ = np.stack( ( hitsYZ.GetX(), hitsYZ.GetY() ), axis=1 )
        
        nHitsY = len([i for i in arrayYZ[:,0] if i>0])
        nHitsYUP = len([i for i in arrayYZ[:,1] if i>0])
        nHitsYDN = len([i for i in arrayYZ[:,1] if i<0])
        
        print(f"Run {run}, Subrun {subrun}, Event {event}")
        print(f"nHitsYUP: {nHitsYUP}, nHitsYDN: {nHitsYDN}")

        if (nHitsYUP + nHitsYDN) != 0:
            r = (nHitsYUP - nHitsYDN) / (nHitsYUP + nHitsYDN)
            print(f"ratio: {r}")
            ratio.Fill(r)
            R.append(r)
        else:
            print("Warning: nHitsYUP + nHitsYDN is zero, skipping ratio calculation")
            R.append(0)
        
        # Test Everything by Drawing
        plot.rcParams['figure.max_open_warning'] = 100
        plot.figure(figsize=(3000*px,800*px))
        plot.scatter(arrayXZ[:,0],arrayXZ[:,1],s=1.4,marker=recMarker,c='black')
        plot.ylim(-800,800)
        plot.xlim(0,6000)
        plot.axis('off')
        imgbuf_XZ = io.BytesIO()
        plot.savefig(imgbuf_XZ, format='png', bbox_inches='tight') #for phi (0 to 360), estimate from horizontal
        #plot.savefig(f"/home/aklima/MSAnalysis/Codes/TextAll/{file_base_name}_{event}_xz.png", bbox_inches='tight')
        im = Image.open(imgbuf_XZ)
        im.show(title="XZ")

        print("Run %i Subrun %i Event %i, %s"%(run,subrun,event,ProcessTimestamp(timestamp).strftime('%Y-%m-%d %H:%M:%S.%f')))

        plot.figure(figsize=(3000*px,800*px))
        plot.scatter(arrayYZ[:,0],arrayYZ[:,1],s=1.4,marker=recMarker,c='black')
        plot.ylim(-800,800)
        plot.xlim(0,6000)
        plot.axis('off')
        imgbuf_YZ = io.BytesIO() #opening image buffer
        plot.savefig(imgbuf_YZ, format='png', bbox_inches='tight') #for theta (0 to 180), estimate from vertical
        #plot.savefig(f"/home/aklima/MSAnalysis/Codes/TextAll/{file_base_name}_{event}_yz.png", bbox_inches='tight')
        im = Image.open(imgbuf_YZ)
        im.show(title="YZ")

        imageXZ = Image.open(imgbuf_XZ)
        imageYZ = Image.open(imgbuf_YZ)
        occupyXZ = len([i for i in arrayXZ[:,0] if i>0])/172000
        occupyYZ = len([i for i in arrayYZ[:,0] if i>0])/172000
        occupancyXZ.append(occupyXZ)
        occupancyYZ.append(occupyYZ)
        asymXZ = shower.ImageAsym(imageXZ)
        AsymXZ.append(asymXZ)
        asymYZ = shower.ImageAsym(imageYZ)
        AsymYZ.append(asymYZ)

        print ("Asymmetry XZ = ", asymXZ)
        print ("Asymmetry YZ = ", asymYZ)
        print ("Occupancy XZ = ", occupyXZ)
        print ("Occupancy YZ = ", occupyYZ)

        AnglesDet = (shower.GetAnglesDet (imageXZ, imageYZ))
        print(AnglesDet)

        AngleXZ.append(AnglesDet[0])
        StDevXZ.append(AnglesDet[1])
        AngleYZ.append(AnglesDet[2])
        StDevYZ.append(AnglesDet[3])
        Evt_timestamp.append(timestamp)
        
        if (nHitsYUP + nHitsYDN) != 0:
            rt = (nHitsYUP - nHitsYDN) / (nHitsYUP + nHitsYDN)
        else:
            rt = 0
        
        DetL = 60
        DetW = 15.6
        DetH = 15.6
        if (AnglesDet[0] != 0 or AnglesDet[2] != 0):
            alt, az = shower.GetAltAz(AnglesDet[0],AnglesDet[2])
            alti = np.radians(alt)
            azi = np.radians(AnglesDet[0])
            axsec = DetL * DetH * abs( np.cos(alti) * np.sin(azi) ) + DetW * DetH * abs( np.cos(alti) * np.cos(azi) ) + DetL * DetW * abs( np.sin(alti) )
            Alt.append(alt)
            Az.append(az)
            A_xsec.append(axsec)
            print("A_xsec: ", axsec)
            for i in range(0,len(azs), 18):
                if (az > azs[i]-5 and az <= azs[i]):
                    for j in range(i, i+18):
                        if (alt > alts[j]-5 and alt <= alts[j]):
                            if nhits[j] != 0:
                                nT1 = nHitsY / nhits[j]
                                print(f"nT1: {nT1}, Altitude: {alt}, Azimuth: {az}")
                                nTracks_new.append(nT1)
                            else:
                                nTracks_new.append(0)
                                
                            if nhitsall[j] != 0:
                                nT2 = nHitsY / nhitsall[j]
                                trkden = nT2 / axsec
                                print(f"nT2: {nT2}, Altitude: {alt}, Azimuth: {az}")
                                nTracks_all.append(nT2)
                                TrkDen.append(trkden)
                            else:
                                nTracks_all.append(0)
                                TrkDen.append(0)
        else:
            nTracks_new.append(0)
            nTracks_all.append(0)
            TrkDen.append(0)
            A_xsec.append(0)
            Alt.append(0)
            Az.append(0)
            
            
        if AnglesDet[2] < 0:
            ang = 180 + AnglesDet[2]
            for i in range(len(angle)):
                if ang >= angle[i] and ang <= angle[i+1]:
                    nT = nHitsY / mean[i]
                    print("nT",nT)
                    nTracks.append(nT)
                    if rt >= -0.05 and rt <= 0.2 and asymYZ <= 0.08:
                        tracks.Fill(nT)
                    else: 
                        tracksr.Fill(nT)    
                    print("angle1: ",ang)
        if AnglesDet[2] == 0:
            nTracks.append(0)
        if AnglesDet[2] > 0:
            ang = 180 - AnglesDet[2]
            print("angle2: ",ang)
            for i in range(len(angle)):
                if ang >= angle[i] and ang <= angle[i+1]:
                    nT = nHitsY / mean[i]
                    print("nT",nT)
                    nTracks.append(nT)
                    if rt >= -0.05 and rt <= 0.2 and asymYZ <= 0.08:
                        tracks.Fill(nT)
                    else: 
                        tracksr.Fill(nT)    
                    print("angle1: ",ang)

    print(len(Run), len(Subrun), len(Event), len(S_start), len(S_stop), len(occupancyXZ), len(occupancyYZ), len(AsymXZ), len(AsymYZ), len(AngleXZ), len(AngleYZ), len(StDevXZ), len(StDevYZ),len(Alt),len(Az),len(nTracks),len(nTracks_new),len(nTracks_all),len(A_xsec),len(TrkDen),len(R),len(Evt_timestamp))
    for i in range(len(Run)):
        CombinedList.append((Run[i], Subrun[i], Event[i], S_start[i], S_stop[i], occupancyXZ[i], occupancyYZ[i], AsymXZ[i], AsymYZ[i], AngleXZ[i], AngleYZ[i], StDevXZ[i], StDevYZ[i], Alt[i], Az[i], nTracks[i], nTracks_new[i], nTracks_all[i], A_xsec[i], TrkDen[i], R[i], Evt_timestamp[i]))

    #print(len(AsymXZ), len(AsymYZ))
    if CombinedList and len(CombinedList) > 0:
        headers=['Run','SRun','Evt','SR_start','SR_stop', 'Oc_XZ', 'Oc_YZ','AsymXZ','AsymYZ','AngXZ','AngYZ', 'StdXZ', 'StdYZ', 'Alt', 'Az', 'nTY_o', 'nTY_n', 'nT_a', 'XSec', 'TrkDen', 'TBSym', 'Timestamp'] #, 'OccupancyXZ', 'OccupancyYZ'
        formats=['%d','%d','%d','%d','%d','%.3f','%.3f','%.3f','%.3f','%.3f','%.3f','%.5f','%.5f','%.2f','%.2f','%.0f','%.0f','%.0f','%.2f','%.2f','%.2f','%d']
        output_file_name = f"/home/aklima/MSAnalysis/Codes/TextAll/{file_base_name}_EventInfo.txt"
        np.savetxt(output_file_name, CombinedList, delimiter='\t',fmt=formats, header='\t'.join(headers), comments='#')
    else:
        output_file_name = f"/home/aklima/MSAnalysis/Codes/TextAll/{file_base_name}_EventInfo.txt" 
        with open(output_file_name, "w"):
            pass
       
    print(n)
    plot.close('all') #closing all plots and canvas
    imgbuf_XZ.close() #closing image buffer
    imgbuf_YZ.close()
    

canvas0 = ROOT.TCanvas("canvas0","canvas0",800,600)
ratio.Draw("") 
canvas0.Draw()
canvas0.SaveAs("canvas0_plot.png")

canvas1 = ROOT.TCanvas("canvas1","canvas1",800,600)
tracksr.Draw("") 
canvas1.Draw()
canvas1.SaveAs("canvas1_plot.png")

canvas2 = ROOT.TCanvas("canvas2","canvas2",800,600)
tracks.Draw("") 
canvas2.Draw()
canvas2.SaveAs("canvas2_plot.png")


ROOT.gApplication.Run(True)

# In[ ]:




