// oblancog. 2016.07.xx NO purpose tracking
// this code generates the initial coordinates of the particles
// reads the show,beam; and a show,betaXXX; output from madx

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

int GenerateInrays (const char * k, int N) {
  cout << "  Using flag : "<< k << endl;
  TString * myflname = new TString(k);
  //  gROOT->Reset();
  //  int N =10000;//number of particles to generate
  int debug=1; // if debug 1, print goes to stdout
  int tracktable=1; // if tracktable 1, print goes to madx track table format
  int madxtrac=1;//track in 1=MAD-X and 2=PTC_MAD-X

  // Energy distribution 
  // 0=Uniform -EnergySpread/2 to EnergySpread/2
  // 1=Gaussian(0,sigma=EnergySpread)
  int Edistr=1;

  // Limit of the gaussian
  int gauslimit=6;


  //Relativistic factors
  double_t gammar = 1;
  double_t betar  = 1;

  //Beam geometrical emittances
  double_t ebeam = 1;
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
  double_t offsetx0 = 0;
  double_t offsety0 = 0;

  TRandom2 *xrnd = new TRandom2(0);//any different seed will do
  TRandom2 *pxrnd = new TRandom2(1);
  TRandom2 *yrnd = new TRandom2(2);//any different seed will do
  TRandom2 *pyrnd = new TRandom2(3);
  TRandom2 *Espreadrnd = new TRandom2(4);
  TRandom2 *bunchlrnd = new TRandom2(5);
  int i = 0; //particle counter
  double_t xbeta, pxbeta, ybeta, pybeta, Esprd = 0;
  double_t xgausvalue;
  double_t ygausvalue;
  double_t tgausvalue;

  //Output values
  double_t ux=0;
  double_t uy=0;
  double_t upx=0;
  double_t upy=0;
  double_t ups=0;
  double_t upd=0;

  ofstream mydebug;
  ofstream mytracktable;
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
  TString * betafl = new TString("beta");
  betafl->Append(k);
  betafl->Append(".txt");
  beta0in.open(betafl->Data());
  if (beta0in == 0) {
    // if we cannot open the file, 
    // print an error message and return immediatly
    cout << "Error: cannot open "<<betafl->Data()<<" !"<<endl;
    //   printf("Error: cannot open betaXXXX.txt!\n");
    return 1;
  }
  cout << "  ... reading file "<<betafl->Data()<<" (twiss params at input)"<<endl;
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
    if (!strcmp(madx01,"dy")) {cout<<"    "<<madx01<<" "<<madx04<<endl;etay=atof(madx04);}
    if (!strcmp(madx01,"dpy")) {cout<<"    "<<madx01<<" "<<madx04<<endl;etapy=atof(madx04);}
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
    if (!strcmp(madx01,"energy")){
      beam0in >> madx02;
      beam0in >> madx03;
      beam0in >> madx04;
      cout<<"    "<<madx01<<" "<<madx04<<endl;
      ebeam=atof(madx04);}
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
    if (!strcmp(madx01,"sigt")){
      beam0in >> madx02;
      beam0in >> madx03;
      beam0in >> madx04;
      cout<<"    "<<madx01<<" "<<madx04<<endl;
      sigmas0=atof(madx04);}
    if (!strcmp(madx01,"et")){
      beam0in >> madx02;
      beam0in >> madx03;
      beam0in >> madx04;
      cout<<"    "<<madx01<<" "<<madx04<<endl;
      sigmas0=atof(madx04);}
  }
  cout << "    ... all others ignored.";
  cout << "    beam0.txt read."<<endl;
  beam0in.close();

  if (debug) mydebug.open ("debug.txt");
  if (tracktable) mytracktable.open ("trackTABLE");
  if (madxtrac) mymadxtrac.open ("madxInrays.madx");
  mymadxtrac << "! GenerateInrays. orblancog. 2016.07.28\n";  
  mymadxtrac << "! Dummy file generated in root\n";

  // Calculate sigma0
  sigmax0  = TMath::Sqrt(ex*betax);
  sigmapx0 = TMath::Sqrt(ex/betax);
  sigmay0  = TMath::Sqrt(ey*betay);
  sigmapy0 = TMath::Sqrt(ey/betay);

  //  cout << sigmax0 << sigmay0;
  //  sigmas0  = TMath::Sqrt(et*0);//already assigned when reading beam0
  sigmad0  = 1.0/betar*Energyspread;//TMath::Sqrt(et*0);//\Delta E/(Pc) = 1/\beta_r * \Delta E/E
  Double_t xbetarnd,pxbetarnd;
  Double_t ybetarnd,pybetarnd;
  Double_t tbetarnd,ptbetarnd;

  while (i<N){
    //x
    xbetarnd   =  xrnd->Gaus(0,1);
    pxbetarnd  =  pxrnd->Gaus(0,1);
    xbeta      =  TMath::Sqrt(ex*betax)*xbetarnd;
    pxbeta     =  TMath::Sqrt(ex/betax)*pxbetarnd - alfax*TMath::Sqrt(ex/betax)*xbetarnd;
    xgausvalue = (gammax*xbeta*xbeta+2*alfax*xbeta*pxbeta+betax*pxbeta*pxbeta)/(TMath::TwoPi()/2*ex);
    //y
    ybetarnd   = yrnd->Gaus(0,1);
    pybetarnd  = pyrnd->Gaus(0,1);
    ybeta      =  TMath::Sqrt(ey*betay)*ybetarnd;
    pybeta     =  TMath::Sqrt(ey/betay)*pybetarnd - alfay*TMath::Sqrt(ey/betay)*ybetarnd;    
    ygausvalue = (gammay*ybeta*ybeta+2*alfay*ybeta*pybeta+betay*pybeta*pybeta)/(TMath::TwoPi()/2*ey);
    //d
    ups = bunchlrnd->Gaus(0,sigmas0);
    if (Edistr){
      upd = Espreadrnd->Gaus(0,sigmad0);
      tgausvalue = upd*upd/(sigmad0*sigmad0)/(TMath::TwoPi()/2) + ups*ups/(sigmas0*sigmas0)/(TMath::TwoPi()/2);
    }else{
      upd = Espreadrnd->Uniform(-0.5*sigmad0,0.5*sigmad0);//\Delta E/(Pc) = 1/\beta_r * \Delta E/E
      tgausvalue = ups*ups/(sigmas0*sigmas0);
    }
      
    if (xgausvalue<gauslimit && ygausvalue<gauslimit && tgausvalue<gauslimit){
      i++;
      ux  = xbeta  + etax *upd   + 0e-3;
      upx = pxbeta + etapx*upd   - 0e-3;
      uy  = ybeta  + etay *upd   + 0e-3;
      upy = pybeta + etapy*upd   + 0e-3;
      upd = upd                  - 0e-3;
      if (debug) mydebug <<ux<<'\t'<<upx<<'\t'<<uy<<'\t'<<upy<<"\t"<<ups<<"\t"<<upd<<endl ;
      if (madxtrac==1) mymadxtrac << "start, x="<<ux<<",px="<<upx<<",y="<<uy<<",py="<<upy<<",t="<<ups<<",pt="<<upd<<";\n";
      if (madxtrac==2) mymadxtrac << "ptc_start, x="<<ux<<",px="<<upx<<",y="<<uy<<",py="<<upy<<",t="<<ups<<",pt="<<upd<<";\n";
      //      if (madxtrac==2) mymadxtrac << "ptc_start, x="<<ux<<",px="<<upx<<",y="<<uy<<",py="<<upy<<",t="<<ups<<",pt="<<"0"<<";\n";
      if (tracktable==1) mytracktable <<' '<<i<<" 0 "<<ux<<' '<<upx<<' '<<uy<<' '<<upy<<' '<<ups<<' '<<upd<<" 0 "<<ebeam<<endl;//NUMBER,TURN,X,PX,Y,PY,T,PT,S,E
    }
  }
  if (tracktable) mytracktable.close();
  if (debug) mydebug.close();
  mymadxtrac.close();
  
  cout << "  " << i << " rays generated. All OK. Adios !";
  return 0;
};
