#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <tuple>
#include <cstdlib>
#include <math.h>       
#include <cmath>        
#include <stdlib.h>  
#include <stdio.h>
#include <string.h>
#include "TVector3.h"
#include "TROOT.h"
#include "TMath.h"
#include "TF1.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TProfile3D.h"
#include "TProfile2D.h"
#include "TProfile.h"
#include "TGraph.h"
#include "TBrowser.h"
#include "TRandom.h"
#include <sstream>
#include <algorithm>
#include <random>
using namespace std;

struct Proton{
    
  Float_t x12,y12,z12;
  Float_t p1x,p1y,p1z;
   
  Float_t x21,y21,z21; 
  Float_t p2x,p2y,p2z;
  
  Float_t wepl;
};

int main(int argc, char *argv[]){
  // Load Proton data
  Proton Point;
  //Prepare in file	
  char rootf[10] = ".root";
  char inf[50] = "../outputFiles/yourFileName_"; //*** make sure you are using either single or double-sided projections
  strcat(inf,argv[1]);
  strcat(inf,rootf);   
  TFile* phaseFile = new TFile(inf,"update"); //Associate any new objects to the proper working file
 
  TTree* tt = (TTree*)phaseFile->Get("filteredTree"); //Filtered hull data from running filterSS.cc!

  tt->SetBranchAddress("x12_h",&Point.x12);
  tt->SetBranchAddress("y12_h",&Point.y12);
  tt->SetBranchAddress("z12_h",&Point.z12);
  
  tt->SetBranchAddress("p1x_h",&Point.p1x);
  tt->SetBranchAddress("p1y_h",&Point.p1y);
  tt->SetBranchAddress("p1z_h",&Point.p1z);  
  
  tt->SetBranchAddress("x21_h",&Point.x21);
  tt->SetBranchAddress("y21_h",&Point.y21);
  tt->SetBranchAddress("z21_h",&Point.z21);
  
  tt->SetBranchAddress("p2x_h",&Point.p2x);
  tt->SetBranchAddress("p2y_h",&Point.p2y);
  tt->SetBranchAddress("p2z_h",&Point.p2z);
  
  tt->SetBranchAddress("weplReal_h",&Point.wepl);
  
  int NEntries = tt->GetEntries(); 
  
//--------------------------------
// Defining the header variables for the binary
//-------------------------------
  float projection_angle = atof(argv[1]); //projection angle as integer
  float beam_energy = 230.; //Proton initial energy
  char magic_number[]="PCTD";
  char PHANTOM_NAME[]="LinePair"; //name of the phantom
  char PREPARED_BY[]="Jarle R S"; //Author
  int version_id = 0;
  int current_time = time(NULL);
  int phantom_name_size = sizeof(PHANTOM_NAME);
  int data_source_size = sizeof(inf); 
  int prepared_by_size = sizeof(PREPARED_BY);

  Float_t t[4], v[4], u[4], wepl, angle;
//------------------------------
// Writing the header
//------------------------------  
  printf("Writing headers to file...\n");
  ofstream outfile(Form("../outputSS/projection_%03.0f.bin",projection_angle), ios::binary); //*** Set folder
  //Check File Type
  outfile.write(magic_number, 4);
  outfile.write(reinterpret_cast<char*>(&version_id), sizeof(int));
  outfile.write(reinterpret_cast<char*>(&NEntries), sizeof(int));  
  outfile.write(reinterpret_cast<char*>(&projection_angle), sizeof(float));
  outfile.write(reinterpret_cast<char*>(&beam_energy), sizeof(float));
  outfile.write(reinterpret_cast<char*>(&current_time), sizeof(int));//NOT USED
  outfile.write(reinterpret_cast<char*>(&current_time), sizeof(int));//NOT USED (Just for consistency check)
  outfile.write(reinterpret_cast<char*>(&phantom_name_size), sizeof(int));
  outfile.write(PHANTOM_NAME, phantom_name_size);
  outfile.write(reinterpret_cast<char*>(&data_source_size), sizeof(int));
  outfile.write(inf, data_source_size); //Name of the input file for consitency checks
  outfile.write(reinterpret_cast<char*>(&prepared_by_size), sizeof(int));
  outfile.write(PREPARED_BY, prepared_by_size);

  int data_size = (NEntries) * sizeof(float);
  
  // Allocate memory
  float* t_in1_h         = (float*) malloc(data_size);
  float* t_in2_h         = (float*) malloc(data_size);
  float* t_out1_h        = (float*) malloc(data_size);
  float* t_out2_h        = (float*) malloc(data_size);
  float* u_in1_h         = (float*) malloc(data_size);
  float* u_in2_h         = (float*) malloc(data_size);
  float* u_out1_h        = (float*) malloc(data_size);
  float* u_out2_h        = (float*) malloc(data_size);
  float* v_in1_h         = (float*) malloc(data_size);
  float* v_in2_h         = (float*) malloc(data_size);
  float* v_out1_h        = (float*) malloc(data_size);
  float* v_out2_h        = (float*) malloc(data_size);
  float* WEP_h           = (float*) malloc(data_size);

  // Initialize memory
  memset(t_in1_h, 0, data_size);
  memset(t_in2_h, 0, data_size);
  memset(t_out1_h, 0, data_size);
  memset(t_out2_h, 0, data_size);
  memset(v_in1_h, 0, data_size);
  memset(v_in2_h, 0, data_size);
  memset(v_out1_h, 0, data_size);
  memset(v_out2_h, 0, data_size);
  memset(u_in1_h, 0, data_size);
  memset(u_in2_h, 0, data_size);
  memset(u_out1_h, 0, data_size);
  memset(u_out2_h, 0, data_size);
  memset(WEP_h, 0, data_size);
  
  // All positions in mm! 
  float distance_tracker  = 50.0;
  
  int NPerformed = 0; 
  for(int i=0;i<NEntries;i++){//Loop over all protons again
	  if(i%100000 == 0) cout<<i<<" "<<NPerformed<<endl;
    tt->GetEntry(i);  
    
		  //t=transveral
      t_in1_h[i] = Point.x12-(Point.p1x*distance_tracker);  
      t_in2_h[i] = Point.x12;
      t_out1_h[i] = Point.x21;
      t_out2_h[i] = Point.x21+(Point.p2x*distance_tracker);;
      
      //v=vertical 
      v_in1_h[i]  = Point.z12-(Point.p1z*distance_tracker); 
      v_in2_h[i]  = Point.z12;
      v_out1_h[i] = Point.z21;
      v_out2_h[i] = Point.z21+(Point.p2z*distance_tracker);
        
      //u=distal
      u_in1_h[i]  = Point.y12-distance_tracker; 
      u_in2_h[i]  = Point.y12;
      u_out1_h[i] = Point.y21;
      u_out2_h[i] = Point.y21+distance_tracker;
      WEP_h[i]    = Point.wepl;
	    NPerformed++;
  }
  phaseFile->Close();
    
  outfile.write((char*)t_in1_h, data_size);
  outfile.write((char*)t_in2_h, data_size);
  outfile.write((char*)t_out1_h, data_size);
  outfile.write((char*)t_out2_h, data_size);
  outfile.write((char*)v_in1_h, data_size);
  outfile.write((char*)v_in2_h, data_size);
  outfile.write((char*)v_out1_h, data_size);
  outfile.write((char*)v_out2_h, data_size);
  outfile.write((char*)u_in1_h, data_size);
  outfile.write((char*)u_in2_h, data_size);
  outfile.write((char*)u_out1_h, data_size);
  outfile.write((char*)u_out2_h, data_size);
  outfile.write((char*)WEP_h, data_size);
  outfile.close();
  cout << "Binary file written" << endl;  
  return 0; 
}
