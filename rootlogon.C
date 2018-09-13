#ifndef __CINT__
#include "TROOT.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TCanvas.h"

#include <iostream>
#include <string>
#include <vector>
#include <cstdlib>
#endif

//#include <sys/types.h>
//#include <unistd.h>

void AddIncludePath(std::string dir);
void AddLinkedLibs();
void LoadModule(std::string dir);

void rootlogon()
{

  //  system(Form("xterm -e gdb -p %d",getpid()));
  //  cout << getpid() << endl;
  //  system("export SAMURAI_CATANA_BETA=0.8"); 
  char* c = getenv("TARTSYS");
  if(!c){
    std::cout << "set environment variable \"TARTSYS\"" << std::endl;
    std::cout << "quit ROOT" << std::endl;      
    gROOT->ProcessLine(".q");
  }
  std::string install_dir(c);

  AddIncludePath(install_dir);
  AddLinkedLibs();
  LoadModule(install_dir);

  std::cout << "\n Manual of ANAROOT \n  http://be.nucl.ap.titech.ac.jp/~kondo/moin/moin.cgi/ANAROOT/Manual\n" << std::endl;

  //Base Style
  //  gROOT->SetStyle("Plain");
  gROOT->SetStyle("Modern");
  //  gROOT->SetStyle("Classic");

  //Force Style
  gStyle->SetHistFillColor(7);
  gStyle->SetHistFillStyle(3002);
  gStyle->SetHistLineColor(kBlue);
  gStyle->SetFuncColor(kRed);
  gStyle->SetFrameLineWidth(2);
  gStyle->SetPadGridX(1);
  gStyle->SetPadGridY(1);
  gStyle->SetCanvasColor(0);
  gStyle->SetTitleFillColor(0);
  gStyle->SetTitleStyle(0);
  gStyle->SetStatColor(0);
  gStyle->SetStatStyle(0);
  gStyle->SetStatX(0.9);  
  gStyle->SetStatY(0.9);  
  gStyle->SetPalette(1);
  gStyle->SetOptLogz(1);
  //  gStyle->SetOptTitle(0);
  gStyle->SetOptFit(1);
  gStyle->SetOptStat(1111111);
  gStyle->SetPadBorderMode(1);
  //  gStyle->SetOptDate(1);
  gStyle->SetOptDate(0);

  // gStyle->SetLabelFont(132,"XYZ");
  // gStyle->SetTitleFont(132,"XYZ");
  // gStyle->SetTitleFont(132,"");
  // gStyle->SetTextFont(132);
  // gStyle->SetStatFont(132);
  gStyle->SetLabelFont(132,"XYZ");
  gStyle->SetTitleFont(132,"XYZ");
  gStyle->SetTitleFont(132,"");
  gStyle->SetTextFont(132);
  gStyle->SetStatFont(132);
  gStyle->SetTitleSize(0.045);
  gStyle->SetLabelSize(0.04,"XYZ");
  gStyle->SetCanvasDefW(500);
  gStyle->SetCanvasDefH(500);
  
  //gStyle->SetCanvasSize(500,500);
  
  gStyle->SetPaperSize(30,16);
  const Int_t NRGBs = 5;
  const Int_t NCont = 99;
  
  Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
  Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
  Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
  Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
  TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);   
  gStyle->SetNumberContours(NCont);

}

void AddIncludePath(std::string install_dir)
{
  std::vector<std::string> include;
  include.push_back("-I"+install_dir+"/include");
  //  include.push_back("`xml2-config --cflags`");

  std::vector<std::string>::iterator it = include.begin();
  while(it != include.end()){
    gSystem->AddIncludePath((*it).c_str());
    std::cout << "add include path : " << *it << std::endl;
    ++it;
  }
}

void AddLinkedLibs()
{
//  std::vector<std::string> include;
//  include.push_back("`xml2-config --libs`");

//  std::vector<std::string>::iterator it = include.begin();
//  while(it != include.end()){
//    gSystem->AddLinkedLibs((*it).c_str());
//    std::cout << "add linked libs : " << *it << std::endl;
//    ++it;
//  }
}

void LoadModule(std::string install_dir)
{
  std::vector<std::string> modules;
  modules.push_back("libXMLParser.so");
  modules.push_back(install_dir+"/lib/"+"libanaroot.so"); // load at once
  modules.push_back("/usr/lib64/libpthread.so"); // for multithreading
  modules.push_back("$ROOTSYS/lib/libThread.so"); // for multithreading

//  modules.push_back(install_dir+"/lib/"+"libananadeko.so"); // load each modules one by one
//  modules.push_back(install_dir+"/lib/"+"libanacore.so");
//  modules.push_back(install_dir+"/lib/"+"libanabrips.so");
//  modules.push_back(install_dir+"/lib/"+"libanadali.so");
//  modules.push_back(install_dir+"/lib/"+"libanasamurai.so");
//  modules.push_back(install_dir+"/lib/"+"libanaanaloop.so");
//  modules.push_back(install_dir+"/lib/"+"libanaespri.so");

  std::vector<std::string>::iterator it = modules.begin();
  while(it != modules.end()){
    std::cout << "reading " << *it << std::endl;
    if(gSystem->Load((*it).c_str()) < 0){
      std::cout << "cannot read in " << *it << std::endl;      
    }
    ++it;
  }
}

