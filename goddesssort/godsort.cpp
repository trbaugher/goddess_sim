#include "TTree.h"
#include "TFile.h"
#include "TRandom3.h"
#include "TH2F.h"
#include "TH1F.h"

#include "GodHit.h"

#include <fstream>
#include <iostream>

using namespace std;

int main(int argc, char** argv)
{
  if (argc == 1) {
    cout << "Usage: goddesssort <filename.out>\n";
    cout << "Sorts simulation output.\n";
    cout << "travis.baugher@gmail.com\n";
    return 0;
  }

  GodEvent theEvent;

  TTree* GodTree = new TTree("GodTree", "GodTree");

  GodTree->Branch("nHit", &theEvent.nHits, "theEventnHit/I");

  GodTree->Branch("GS.nHits", &theEvent.theGS.nHits, "theEventtheGSnHits/I");
  GodTree->Branch("GS.ID",     theEvent.theGS.ID,    "theEventtheGSID[theEventtheGSnHits]/I");
  GodTree->Branch("GS.E",      theEvent.theGS.E,     "theEventtheGSE[theEventtheGSnHits]/D");
  GodTree->Branch("GS.T",      theEvent.theGS.T,     "theEventtheGST[theEventtheGSnHits]/D");
  GodTree->Branch("GS.X",      theEvent.theGS.X,     "theEventtheGSX[theEventtheGSnHits]/D");
  GodTree->Branch("GS.Y",      theEvent.theGS.Y,     "theEventtheGSY[theEventtheGSnHits]/D");
  GodTree->Branch("GS.Z",      theEvent.theGS.Z,     "theEventtheGSZ[theEventtheGSnHits]/D");

  GodTree->Branch("BGO.nHits", &theEvent.theBGO.nHits, "theEventtheBGOnHits/I");
  GodTree->Branch("BGO.ID",     theEvent.theBGO.ID,    "theEventtheBGOID[theEventtheBGOnHits]/I");
  GodTree->Branch("BGO.E",      theEvent.theBGO.E,     "theEventtheBGOE[theEventtheBGOnHits]/D");
  GodTree->Branch("BGO.T",      theEvent.theBGO.T,     "theEventtheBGOT[theEventtheBGOnHits]/D");
  GodTree->Branch("BGO.X",      theEvent.theBGO.X,     "theEventtheBGOX[theEventtheBGOnHits]/D");
  GodTree->Branch("BGO.Y",      theEvent.theBGO.Y,     "theEventtheBGOY[theEventtheBGOnHits]/D");
  GodTree->Branch("BGO.Z",      theEvent.theBGO.Z,     "theEventtheBGOZ[theEventtheBGOnHits]/D");

  ifstream file;
  file.open(argv[1]);

  // expecting data format: 
  // Event# #hits
  // detID edep x y x

  int NHits = 0;
  int EventID = 0;
  int DetID = 0;
  double Edep = 0.0; // exact energy from output file
	double E = 0.0; // Energy with detector resolution applied
  double T = 0.0;
  double X = 0.0;
  double Y = 0.0;
  double Z = 0.0;
  double sigma = 1;
  int EventCount = 0;

  int i = 0;

  // random number generator for detector resolution
  TRandom3 random = TRandom3();

  while (file >> EventID >> NHits) {
    //cout << "Event ID " << EventID << " # hits " << NHits << endl;
    if (EventCount++ % 1000 == 0) {
      cout << "Event " << EventCount-1 << endl;
    }
    theEvent.Reset();
    theEvent.nHits = NHits;
    for (i = 0; i < NHits; i++) {
      file >> DetID >> Edep >> T >> X >> Y >> Z;
      //cout << DetID << " " << Edep*1000 << " " << T << " " << X << " " <<  Y << " " << Z << endl;;
      if ((DetID >= 401) && (DetID < 511)) { // GammaSphere
        theEvent.theBGO.ID[theEvent.theBGO.nHits] = DetID;
        theEvent.theBGO.E[theEvent.theBGO.nHits] = Edep*1000; // convert to keV
        theEvent.theBGO.T[theEvent.theBGO.nHits] = T;
        theEvent.theBGO.X[theEvent.theBGO.nHits] = X;
        theEvent.theBGO.Y[theEvent.theBGO.nHits] = Y;
        theEvent.theBGO.Z[theEvent.theBGO.nHits] = Z;
        theEvent.theBGO.nHits++;
      }
      if ((DetID >= 1) && (DetID < 111)) { // GammaSphere
        theEvent.theGS.ID[theEvent.theGS.nHits] = DetID;
        // apply detecror resolution
        sigma = 2.5;
				E = random.Gaus(Edep*1000, sigma);
        theEvent.theGS.E[theEvent.theGS.nHits] = E;
        theEvent.theGS.T[theEvent.theGS.nHits] = T;
        theEvent.theGS.X[theEvent.theGS.nHits] = X;
        theEvent.theGS.Y[theEvent.theGS.nHits] = Y;
        theEvent.theGS.Z[theEvent.theGS.nHits] = Z;
        theEvent.theGS.nHits++;
      }
    }

    GodTree->Fill();
  }

  cout << "Processed " << EventCount << " events\n";

  GodTree->Print();

  string rFileName = argv[1];
  rFileName+=".root";
  TFile f(rFileName.c_str(), "recreate");
  GodTree->Write();
  f.Close();

  file.close();
  return 0;
}
