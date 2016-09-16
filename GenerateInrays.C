// oblancog. 2016.07.xx NO purpose tracking
// this code generates the initial coordinates of the particles

// For the particle generation
#include "TRandom2.h"//period=10**26, fast, 32 bits for the state
// To check if the courant-snyder invariant 
#include "TMath.h"
// To read/write into a file
#include <iostream>
#include "TFile.h"
#include <fstream>
#include <stdio.h>
#include <string>

using namespace std;

int GenerateInrays () {
  //  gROOT->Reset();
  int N =1000;//number of particles to generate
  int debug=1; // if debug 1, print goes to stdout
  int madxtrac=1;//track in MAD-X and PTC_MAD-X

  // Energy distribution 
  // 0=Uniform -EnergySpread/2 to EnergySpread/2
  // 1=Gaussian(0,sigma=EnergySpread)
  int Edistr=0;

  //Relativistic factors
  double_t gammar = 1;
  double_t betar  = 1;

  //Beam geometrical emittances
  double_t ex = 1;
  double_t ey = 1;
  double_t et = 1;
  // Phase space at input
  double_t gammax = 0.5;
  double_t betax = 2;
  double_t alfax = -0.5;
  double_t gammay = 1;
  double_t betay = 1;
  double_t alfay = 0;
  double_t Energyspread = 1;
  double_t etax=1;//off-momentum function
  double_t etapx=1;//off-momentum function
  double_t etay=0;//off-momentum function
  double_t etapy=0;//off-momentum function
  // Linear beam sizes from emittances and twiss
  double_t sigmax0 = 0;
  double_t sigmapx0 = 0;
  double_t sigmay0 = 0;
  double_t sigmapy0 = 0;
  double_t sigmas0 = 0;
  double_t sigmad0 = 0;

  TRandom2 *xrnd = new TRandom2(0);//any different seed will do
  TRandom2 *pxrnd = new TRandom2(1);
  TRandom2 *yrnd = new TRandom2(0);//any different seed will do
  TRandom2 *pyrnd = new TRandom2(1);
  TRandom2 *Espreadrnd = new TRandom2(2);
  int i = 0; //particle counter
  double_t xbeta, pxbeta, ybeta, pybeta, Esprd = 0;
  double_t xgausvalue;
  double_t ygausvalue;
  double_t tgausvalue;
  int gauslimit=1 ;

  //Output values
  double_t ux=0;
  double_t uy=0;
  double_t upx=0;
  double_t upy=0;
  double_t upt=0;


  ofstream mydebug;
  ofstream mymadxtrac;
  ifstream beta0in;
  string beta0string;
  ifstream beam0in;
  string beam0string;

  int nlines;
  char beta0line[80];
  char beam0line[80];

  char madx00[20];
  char madx01[20];
  char madx02[20];
  char madx03[20];
  char madx04[20];


  // Read twiss info
  beta0in.open("beta0.txt");
  if (beta0in == 0) {
    // if we cannot open the file, 
    // print an error message and return immediatly
    printf("Error: cannot open beta0.txt!\n");
    return 1;
  }
  cout << "  ... reading file beta0.txt (twiss params at input)"<<endl;
  beta0in >>  madx00 >> madx01 >> madx02 >> madx03;
  while(!beta0in.eof()){
    beta0in >> madx00 >> madx01 >> madx02 >> madx03 >> madx04;
    //    cout << madx01<<endl;
    if (!strcmp(madx01,"betx")) {cout<<"    "<<madx01<<" "<<madx04<<endl;betax=atof(madx04);}
    if (!strcmp(madx01,"alfx")) {cout<<"    "<<madx01<<" "<<madx04<<endl;alfax=atof(madx04);}
    if (!strcmp(madx01,"bety")) {cout<<"    "<<madx01<<" "<<madx04<<endl;betay=atof(madx04);}
    if (!strcmp(madx01,"alfy")) {cout<<"    "<<madx01<<" "<<madx04<<endl;alfay=atof(madx04);}
    if (!strcmp(madx01,"dx")) {cout<<"    "<<madx01<<" "<<madx04<<endl;etax=atof(madx04);}
    if (!strcmp(madx01,"dpx")) {cout<<"    "<<madx01<<" "<<madx04<<endl;etapx=atof(madx04);}
  }
  cout << "    ... all others ignored.";
  cout << "  beta0.txt read."<<endl;
  beam0in.close();
  cout << "  Calculating gamma[xy]..."<<endl;
  gammax = (1 + alfax*alfax )/betax;
  gammay = (1 + alfay*alfay )/betay;
  cout <<"    gammax "<<gammax<<endl;
  cout <<"    gammay "<<gammay<<endl;

  // Reading beam info
  beam0in.open("beam0.txt");
  if (beam0in == 0) {
    // if we cannot open the file, 
    // print an error message and return immediatly
    printf("Error: cannot open beam0.txt!\n");
    return 1;
  }
  cout << "  ... reading file beam0.txt (beam params)"<<endl;
  beam0in >>  madx00 >> madx01 >> madx02 >> madx03;
  while(!beam0in.eof()){
    beam0in >> madx01;// >> madx02 >> madx03 >> madx04;
    //    cout << madx01<<endl;
    //    cout << madx01<<endl;
    if (!strcmp(madx01,"ex")){
      beam0in >> madx02;
      beam0in >> madx03;
      beam0in >> madx04;
      cout<<"    "<<madx01<<" "<<madx04<<endl;
      ex=atof(madx04);}
    if (!strcmp(madx01,"ey")){
      beam0in >> madx02;
      beam0in >> madx03;
      beam0in >> madx04;
      cout<<"    "<<madx01<<" "<<madx04<<endl;
      ey=atof(madx04);}
    if (!strcmp(madx01,"sige")){
      beam0in >> madx02;
      beam0in >> madx03;
      beam0in >> madx04;
      cout<<"    "<<madx01<<" "<<madx04<<endl;
      Energyspread=atof(madx04);}
  }
  cout << "    ... all others ignored.";
  cout << "    beam0.txt read."<<endl;
  beam0in.close();



  if (debug) mydebug.open ("debug.txt");
  if (madxtrac) mymadxtrac.open ("madxInrays.madx");
  mymadxtrac << "! GenerateInrays. orblancog. 2016.07.28\n";  
  mymadxtrac << "! Dummy file generated in root\n";


  // Calculate sigma0
  sigmax0  = TMath::Sqrt(ex*betax);
  sigmapx0 = TMath::Sqrt(ex/betax);
  sigmay0  = TMath::Sqrt(ey*betay);
  sigmapy0 = TMath::Sqrt(ey/betay);
  sigmas0  = TMath::Sqrt(et*0);
  //  sigmad0  = TMath::Sqrt(et*0);

  while (i<N){
    xbeta = xrnd->Gaus(0,sigmax0);
    //x
    pxbeta = pxrnd->Gaus(0,sigmapx0);
    xgausvalue = (gammax*xbeta*xbeta+2*alfax*xbeta*pxbeta+betax*pxbeta*pxbeta)/ex;
    //y
    ybeta = yrnd->Gaus(0,sigmay0);
    pybeta = pyrnd->Gaus(0,sigmapy0);
    ygausvalue = (gammay*ybeta*ybeta+2*alfay*ybeta*pybeta+betay*pybeta*pybeta)/ey;
    //d
    if (Edistr){
      upt = Espreadrnd->Gaus(0,1.0/betar*Energyspread);//\Delta E/(Pc) = 1/\beta_r * \Delta E/E
      tgausvalue = upt/(1.0/betar*Energyspread);
    }else{
      upt = Espreadrnd->Uniform(-0.5/betar*Energyspread,0.5/betar*Energyspread);//\Delta E/(Pc) = 1/\beta_r * \Delta E/E
      tgausvalue = 0;//dummy value to pass the gaus limit check
    }
    if (xgausvalue<gauslimit && ygausvalue<gauslimit && tgausvalue<gauslimit){
      i++;
      ux = xbeta + etax*Esprd;
      upx = pxbeta + etapx*Esprd;
      uy = ybeta + etay*Esprd;
      upy = pybeta + etapy*Esprd;
      if (debug) mydebug <<ux<<'\t'<<upx<<'\t'<<uy<<'\t'<<upy<<"\t0\t"<<upt<<endl ;
      if (madxtrac) mymadxtrac << "start, x="<<ux<<",px="<<upx<<",y="<<uy<<",py="<<upy<<",t=0,pt="<<upt<<";\n";
    }
  }
  if (debug) mydebug.close();
  mymadxtrac.close();
  
  cout << "  " << i << " rays generated. All OK. Adios !";
  return 0;
};
