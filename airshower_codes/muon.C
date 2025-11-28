/* HOWTO: change the strings infile and outfile point to based on your needs.
   This macro will TChain together infiles with a wildcard.

   There are more variables in the TTree (see the analyzer module for
   what they are) than are in this example.  Some are "per event"
   things, most are "per track" so repeated *ntracks times per event.
   Examples of histogramming both are below */

/* NOTE - this is a root6 macro. */

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TTimeStamp.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"
#include "TProfile2D.h"

using namespace std;

void Preliminary2()
  {
    TLatex* prelim = new TLatex(.935, .95, "NOvA Preliminary");
    prelim->SetTextColor(kBlue);
    prelim->SetNDC();
    prelim->SetTextSize(2/30.);
    prelim->SetTextAlign(32);
    prelim->Draw();
  }

void muon() {

  // use TChain just to get list of files in directory with a wildcard
  TChain chain;
  std::string infilename;
  infilename = "/exp/nova/data/users/aklima/textfiles/cosmicfilter_r*.root";
  chain.Add(infilename.c_str());
  
  TObjArray* filelist = chain.GetListOfFiles(); 
  TIter iter(filelist); 
  
  bool firstEvent = true;
  
  // make histograms
  TH1D *ntrackHisto = new TH1D("ntrackHisto","number of tracks per event",200,0.0,200.0);
  
  TH1D *altdist = new TH1D("altdist", "Altitude distribution", 18, 0, 90);
  
  TH1D *azdist = new TH1D("azdist", "Azimuth distribution", 72, 0, 360);
  
  TH1D *trkxzdist = new TH1D("trkxzdist", "TrkXZ distribution", 360, 0, 360);
  
  TH1D *trkyzdist = new TH1D("trkyzdist", "TrkYZ distribution", 180, -90, 90);
  
  TH1D *dirx = new TH1D("dirx", "DirX distribution", 100, -2, 2);
  
  TH1D *diry = new TH1D("diry", "DirY distribution", 100, -2, 2);
  
  TH1D *dirz = new TH1D("dirz", "DirZ distribution", 100, -2, 2);

  TH1F *subrunHisto = new TH1F("subrunHisto","subrun",64,0.0,63.0);
  
  TH2F *hitstrackYZHisto = new TH2F("hitstrackYZHisto", "Number of hits per track vs. track angles (YZ view)",40,100,180,72,0.0,720.0);
  
  TProfile2D *hitstrackYZAltAz = new TProfile2D("hitstrackYZAltAz"," ",72,0,360,18,0,90);
  
  //TH2F *hMeanHits = new TH2F("hMeanHits","Average number of hits per track vs. Altitude,Azimuth",360,0,360,90,0,90);
  
  //std::ofstream outf("dir.txt");
  //std::ofstream outf1("altaz.txt");
  

  // Loop over all files in chain element list until done selecting those 
  TChainElement* element = 0; 
  while ((element = dynamic_cast<TChainElement*>(iter.Next())) ) 
    { 
      cout<<" new input file "<<element->GetTitle()<<endl; 
       
      TFile* CondenseFile = new TFile(element->GetTitle()); 
      if ( !CondenseFile->IsOpen() || CondenseFile->IsZombie() )  
        { 
          cerr << "Warning! Failed to open file " << CondenseFile->GetName() << endl; 
          delete CondenseFile; 
          continue; 
        } 

      // make a tree reader instance: get the "cosmicTree" tree out of CondenseFile
      // This is the root6 magic.  It's gross in root5
      TTreeReader muonReader("cosmicntuple/cosmicTree", CondenseFile);

      // Let's look at number of tracks.  Things where there are one value
      // per event are accessed like this
      TTreeReaderValue<Int_t> ntrack(muonReader, "ntrack");

      // Let's look at run, which is a UInt_t
      TTreeReaderValue<UInt_t> run(muonReader, "run");
      UInt_t lastRun;
      
      // Let's look at subrun, which is a UChar_t
      TTreeReaderValue<UChar_t> subrun(muonReader, "subrun");
      UInt_t lastSubrun;

      // Let's look at event, which is a UInt_t
      TTreeReaderValue<UInt_t> event(muonReader, "event");
      UInt_t lastEvent;

      // Get a TTimeStamp ready to play with
      TTimeStamp *ts = new TTimeStamp();
      TTreeReaderValue<Int_t> tsHigh(muonReader, "tsHigh");
      TTreeReaderValue<Int_t> tsLow(muonReader, "tsLow");
      TTimeStamp lastts;

      // and the length of all the tracks
      // Things where there's an array per event (in our case, could
      // be a lot of tracks) are accessed like this
      TTreeReaderArray<Float_t> totalLength(muonReader, "totalLength");

      // or the zenith angle
      TTreeReaderArray<Float_t> startDirZ(muonReader, "startDirZ");
      
      // YZ view angles
      TTreeReaderArray<Float_t> startDirY(muonReader, "startDirY");
      
      TTreeReaderArray<Float_t> startDirX(muonReader, "startDirX");
      
      // YZ view hits
      TTreeReaderArray<unsigned short> nYhits(muonReader, "nYhits");
      
      // positions
      TTreeReaderArray<Float_t> startPosZ(muonReader, "startPosZ");
      TTreeReaderArray<Float_t> stopPosZ(muonReader, "stopPosZ");
      TTreeReaderArray<Float_t> startPosY(muonReader, "startPosY");
      TTreeReaderArray<Float_t> stopPosY(muonReader, "stopPosY");
      TTreeReaderArray<Float_t> startPosX(muonReader, "startPosX");
      TTreeReaderArray<Float_t> stopPosX(muonReader, "stopPosX");
      
      
      
      // Loop over all entries of the TTree, one event at a time
      while (muonReader.Next()) {
	// TTreeReaderValue's are accessed with the "*" dereference
	ntrackHisto->Fill((Float_t)*ntrack);
	subrunHisto->Fill((Float_t)*subrun);
	// decode timestamp
	ts->SetSec(*tsHigh);
	ts->SetNanoSec(*tsLow);
	if (firstEvent) {
	  firstEvent = false;
	  cout << "First Event: Run " << *run << " subrun " << (UInt_t)*subrun
	       << " event " << *event << endl
	       << ts->AsString("") << endl;
	}
	
		       
	// Loop over the tracks in the event
	for ( Int_t i = 0; i < *ntrack; i++ ) {
	
	  /* if((startPosY[i] > 750  && stopPosY[i] < -750) || (startPosZ[i] < 50 && stopPosZ[i] > 5950) || (startPosZ[i] > 5950 && stopPosZ[i] < 50) || 
	       (startPosY[i] > 750 && stopPosZ[i] < 50) || (startPosY[i] > 750 && stopPosZ[i] > 5950) || (startPosZ[i] < 50 && stopPosY[i] < -750) || (startPosZ[i] > 5950 && stopPosY[i] < -750)) {
	       
	   (startPosY[i] > 750  && stopPosY[i] < -750) || (startPosZ[i] < 50 && stopPosZ[i] > 5950) || (startPosZ[i] > 5950 && stopPosZ[i] < 50) || 
	        (startPosY[i] > 750 && stopPosZ[i] > 5950)
	       */
	 
        //if ((startPosY[i] > 750  && stopPosY[i] < -750) )      {            
	      
	      hitstrackYZHisto->Fill(acos(startDirY[i])*180.0/M_PI,nYhits[i]);
	      double trkyz;
	      double trkxz;
	      double az;
	      double alt;
	      double FDrot = -27.857222; // degrees from North (along +z axis)
              alt = (M_PI/2 - acos(-startDirY[i])) * 180.0/M_PI;
              az = atan2(startDirX[i],startDirZ[i]) * 180.0/M_PI;
              az = az + FDrot;
              if (az < 0){ az = az + 360;}
              dirx->Fill(startDirX[i]);
              diry->Fill(startDirY[i]);
              dirz->Fill(startDirZ[i]);
              //outf << startDirX[i] << "\t" << startDirY[i] << "\t" << startDirZ[i] << "\t" << trkxz << "\t" << trkyz << "\n"; 
              //outf1 << alt << "\t" << az << "\n";
              trkxzdist->Fill(trkxz);
              trkyzdist->Fill(trkyz);
	      altdist->Fill(alt);
	      azdist->Fill(az);
	      hitstrackYZAltAz->Fill(az,alt,nYhits[i]);
	      
	      //}
	  
	}
	  lastRun = *run;
	  lastSubrun = (UInt_t)*subrun;
	  lastEvent = *event;
	  ts->Copy(lastts);
      }
      cout << "Last Event: Run " << lastRun << " subrun " << lastSubrun
	   << " event " << lastEvent << endl
	   << lastts.AsString("") << endl;
      firstEvent = true;
      // free input tree stuff, close input file
      if (CondenseFile) { delete CondenseFile; CondenseFile = 0; } // also destructs input file tree
      
      cout << "# of tracks: " << *ntrack << "\n";
      Int_t totalLengthSize = totalLength.GetSize(); 
      Int_t startDirZSize = startDirZ.GetSize();
      Int_t startDirYSize = startDirY.GetSize();
      Int_t nYhitsSize = nYhits.GetSize();
      Int_t startPosZSize = startPosZ.GetSize();
      Int_t stopPosZSize = stopPosZ.GetSize();
      Int_t startPosYSize = startPosY.GetSize();
      Int_t stopPosYSize = stopPosY.GetSize();
      Int_t startPosXSize = startPosX.GetSize();
      Int_t stopPosXSize = stopPosX.GetSize();
  
      cout << "length of totalLength array is: " << totalLengthSize << "\n";
      cout << "length of startDirZ array is: " << startDirZSize << "\n";
      cout << "length of startDirY array is: " << startDirYSize << "\n";
      cout << "length of nYhits array is: " << nYhitsSize << "\n";
      cout << "length of startPosZ array is: " << startPosZSize << "\n";
      cout << "length of stopPosZ array is: " << stopPosZSize << "\n";
      cout << "length of startPosY array is: " << startPosYSize << "\n";
      cout << "length of stopPosY array is: " << stopPosYSize << "\n";
      cout << "length of startPosX array is: " << startPosXSize << "\n";
      cout << "length of stopPosX array is: " << stopPosXSize << "\n";
    }
 
 std::ofstream outfile2("AzAltHits_all.txt"); 
 for (int binX = 1; binX <= hitstrackYZAltAz->GetNbinsX(); binX++) {
    for (int binY = 1; binY <= hitstrackYZAltAz->GetNbinsY(); binY++) {
        double content = hitstrackYZAltAz->GetBinContent(binX, binY);
        double azimuth = hitstrackYZAltAz->GetXaxis()->GetBinCenter(binX);
        double altitude = hitstrackYZAltAz->GetYaxis()->GetBinCenter(binY);
        std::cout << "Azimuth: " << azimuth + 2.5 << ", Altitude: " << altitude + 2.5 << ", Bin Content: " << content << std::endl;
        outfile2 << azimuth + 2.5  << '\t' << altitude + 2.5  << '\t' << content << '\n';
    }
 }

  
  // Create a canvas
  TCanvas* canvas = new TCanvas("canvas", "2D Histogram with Vertical Slice Projection", 800, 600);

  // Draw the 2D histogram
  hitstrackYZHisto->Draw("COLZ");
  gStyle->SetStatX(0.9); // X position
  gStyle->SetStatY(0.9); // Y position
  canvas->SetLogz();
  
  TCanvas* canvas2 = new TCanvas("canvas2", "2D Histogram", 800, 600);
  // enlarge right margin for Z axis
  canvas2->SetRightMargin(0.16); 
  canvas2->SetBottomMargin(0.14);

  // Draw the 2D histogram
  hitstrackYZAltAz->SetStats(kFALSE);

  hitstrackYZAltAz->GetXaxis()->SetTitle("Azimuth (degrees)");
  hitstrackYZAltAz->GetXaxis()->SetTitleOffset(1.1);
  hitstrackYZAltAz->GetYaxis()->SetTitle("Altitude (degrees)");
  hitstrackYZAltAz->GetYaxis()->SetTitleOffset(0.9);
  hitstrackYZAltAz->GetZaxis()->SetTitle("Avg. number of hits/track");
  hitstrackYZAltAz->GetZaxis()->SetTitleOffset(1.1);

  // Axis label size and font
  hitstrackYZAltAz->GetXaxis()->SetLabelSize(0.05);
  hitstrackYZAltAz->GetYaxis()->SetLabelSize(0.05);
  hitstrackYZAltAz->GetZaxis()->SetLabelSize(0.05);

  hitstrackYZAltAz->GetXaxis()->SetTitleSize(0.05);
  hitstrackYZAltAz->GetYaxis()->SetTitleSize(0.05);
  hitstrackYZAltAz->GetZaxis()->SetTitleSize(0.05);

  hitstrackYZAltAz->GetXaxis()->SetLabelFont(62);
  hitstrackYZAltAz->GetYaxis()->SetLabelFont(62);
  hitstrackYZAltAz->GetZaxis()->SetLabelFont(62);

  hitstrackYZAltAz->GetXaxis()->SetTitleFont(62);
  hitstrackYZAltAz->GetYaxis()->SetTitleFont(62);
  hitstrackYZAltAz->GetZaxis()->SetTitleFont(62);

  // Ticks and frame
  hitstrackYZAltAz->GetXaxis()->SetTickLength(0.03);
  hitstrackYZAltAz->GetYaxis()->SetTickLength(0.03);
  gStyle->SetLineWidth(2);  // thicker axis lines

  hitstrackYZAltAz->Draw("COLZ");
  //gStyle->SetStatX(0.9); // X position
  //gStyle->SetStatY(0.9); // Y position
  //canvas2->SetLogz();
  
  //canvas2->SetGrid();
  canvas2->Update();
  canvas2->cd();           // make sure we are drawing on canvas2
  Preliminary2();  
  canvas2->Draw(); 
  canvas2->SaveAs("hitstrackYZAltAz.png");
  
  TCanvas* canvas3 = new TCanvas("canvas3", "Alt distribution", 800, 600);
  altdist->Draw();
  gStyle->SetStatX(0.9); // X position
  gStyle->SetStatY(0.9); // Y position
  canvas3->Update();
  canvas3->Draw(); 
  canvas3->SaveAs("AltDist.png");
  
  TCanvas* canvas4 = new TCanvas("canvas4", "Az distribution", 800, 600);
  azdist->Draw();
  gStyle->SetStatX(0.9); // X position
  gStyle->SetStatY(0.9); // Y position
  canvas4->Update();
  canvas4->Draw(); 
  canvas4->SaveAs("AzDist.png");
  
  TCanvas* canvas5 = new TCanvas("canvas5", "TrkXZ distribution", 800, 600);
  trkxzdist->Draw();
  gStyle->SetStatX(0.9); // X position
  gStyle->SetStatY(0.9); // Y position
  canvas5->Update();
  canvas5->Draw(); 
  canvas5->SaveAs("trkxzDist.png");
  
  TCanvas* canvas6 = new TCanvas("canvas6", "TrkYZ distribution", 800, 600);
  trkyzdist->Draw();
  gStyle->SetStatX(0.9); // X position
  gStyle->SetStatY(0.9); // Y position
  canvas6->Update();
  canvas6->Draw(); 
  canvas6->SaveAs("trkyzDist.png");
  
  TCanvas* canvas7 = new TCanvas("canvas7", "DirX distribution", 800, 600);
  dirx->Draw();
  gStyle->SetStatX(0.9); // X position
  gStyle->SetStatY(0.9); // Y position
  canvas7->Update();
  canvas7->Draw(); 
  canvas7->SaveAs("DirX.png");
  
  TCanvas* canvas8 = new TCanvas("canvas8", "DirY distribution", 800, 600);
  diry->Draw();
  gStyle->SetStatX(0.9); // X position
  gStyle->SetStatY(0.9); // Y position
  canvas8->Update();
  canvas8->Draw(); 
  canvas8->SaveAs("DirY.png");
  
  TCanvas* canvas9 = new TCanvas("canvas9", "DirZ distribution", 800, 600);
  dirz->Draw();
  gStyle->SetStatX(0.9); // X position
  gStyle->SetStatY(0.9); // Y position
  canvas9->Update();
  canvas9->Draw(); 
  canvas9->SaveAs("DirZ.png");
   
  

 // Get the number of bins in x direction of the 2D histogram
 int nbinsx = hitstrackYZHisto->GetXaxis()->GetNbins(); 
 
 
 // Create a TFile to save the histograms and fits
 TFile* outputFile = new TFile("output.root", "RECREATE");
 
 // Store mean and standard deviation
 std::ofstream outFile("mean_stdev.txt");
 std::ofstream outFile1("mean_stdev1.txt");
 
 
 for(int bin_x = 1; bin_x <= nbinsx; bin_x++) {
    
    int modeBin = 0;
    double mean = 0.0;
    double stdev = 0.0;

    TH1D* projY = hitstrackYZHisto->ProjectionY("", bin_x, bin_x);

    double avg = projY->GetMean();
    
    
    modeBin = projY->GetMaximumBin();

    // Fit a Gaussian around the mode
    TF1* gaussian = new TF1("gaussian", "gaus", (modeBin-5)*10 , (modeBin+5)*10);
    projY->Fit(gaussian, "QNR");

    
    // Get the mean of the Gaussian fit
    mean = gaussian->GetParameter(1);
    stdev = gaussian->GetParameter(2);
    
    
    // Write mean and standard deviation to the text file
    //if(bin_x >= 45){
    if(mean < 0)  {mean = 0;}
    outFile << bin_x*2 +100 << "\t" << mean << "\t" << stdev << "\n";
    outFile1 << bin_x*2 +100<< "\t" << avg <<"\n";
    //}

    // Plot the mean on the 2D histogram
    canvas->cd();
    TMarker *marker = new TMarker(hitstrackYZHisto->GetXaxis()->GetBinCenter(bin_x), mean, 20);
    marker->SetMarkerSize(1);
    marker->SetMarkerColor(kRed);
    marker->Draw();
    
    TMarker *marker1 = new TMarker(hitstrackYZHisto->GetXaxis()->GetBinCenter(bin_x), avg, 20);
    marker1->SetMarkerSize(1);
    marker1->SetMarkerColor(kGreen);
    marker1->Draw();
    
    TString gaussianName = Form("gaussian_%d", bin_x);
    gaussian->SetName(gaussianName);
    gaussian->Write();

    TString projYName = Form("projection_%d", bin_x);
    projY->SetName(projYName);
    projY->Write();
    
    // Clean up
    delete gaussian;
   
 }
 
  // Close the output file
  outputFile->Close();
  delete outputFile;
  
  // Update the canvas
  canvas->SetGrid();
  canvas->Update();
  canvas->Draw();
  canvas->SaveAs("hitstrackYZHisto.png");
  
  
  // save some histograms to a root file to play with later
  
  TFile* outfile = new TFile("Read_muons-wt.root","RECREATE");
  outfile->cd();
  hitstrackYZHisto->Write();
  canvas->Write("canvas_with_projections");
  outfile->Close();
  
  
}
