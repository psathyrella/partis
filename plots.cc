#include <iostream>
#include <fstream>
#include <sstream>
#include <cassert>
#include <iomanip>
#include <algorithm>
#include <math.h>
#include <vector>
#include <map>
#include "TH1F.h"
#include "TLine.h"
#include "TCanvas.h"
#include "TString.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TLegend.h"

#include "HistUtils.h"
using namespace std;

//----------------------------------------------------------------------------------------
void plots() {
  gStyle->SetOptStat(0);
  TCanvas c1("c1","",700,600);
  TString naivety("M");
  map<TString,int> string_map;
  TString human("A");
  float xmin_force(-330),xmax_force(60);
  int n_bins(60),n_dummy_columns(2);
  TH1F hsame = make_hist("output/test.txt", "same",      "double", "", string_map, n_bins, xmin_force, xmax_force, n_dummy_columns, 1);
  TH1F hdiff = make_hist("output/test.txt", "different", "double", "", string_map, n_bins, xmin_force, xmax_force, n_dummy_columns, 0);

  double xmin(min(hsame.GetBinLowEdge(1), hdiff.GetBinLowEdge(1)));
  double xmax(max(hsame.GetXaxis()->GetBinUpEdge(hsame.GetNbinsX()), hdiff.GetXaxis()->GetBinUpEdge(hdiff.GetNbinsX())));
  xmax = 1.1 * xmax;
  // xmin = -150;
  // xmax = 0.3;
  TH1F hframe("hframe", "", hsame.GetNbinsX(), xmin, xmax);
  hframe.SetMaximum(1.13*(max(hsame.GetMaximum(), hdiff.GetMaximum())));
  hframe.SetTitle(";Bayes factor;frequency");
  hframe.Draw("txt");
  TLegend leg(0.16, 0.75, 0.7, 0.9);
  leg.SetFillColor(0);
  leg.SetBorderSize(0);
  leg.AddEntry(&hsame, "same rearr. event", "l");
  leg.AddEntry(&hdiff, "different", "l");
  leg.Draw();

  hdiff.SetLineColor(kRed+1);
  hsame.SetLineColor(kBlue);
  hdiff.SetMarkerSize(0);
  hsame.SetMarkerSize(0);
  hdiff.SetLineWidth(4);
  hdiff.Draw("ehist same");
  hsame.Draw("ehist same");

  TLine line(xmin, 0, xmax, 0);
  line.SetLineWidth(4);
  line.Draw();
  
  // c1.SetLogx();
  // c1.SetLogy(log.Contains("y"));
  // TString plotdir("/var/www/sharing/dralph/work/plotting/human-beings/" + human + "/M/" + var + "/plots");
  // gSystem->mkdir(plotdir, true);
  // c1.SaveAs(plotdir + "/" + region + "-" + imatch + ".png");
  // c1.SaveAs("output/" + human + "-pairscores.png");
  c1.SaveAs("/tmp/dralph/test.svg");
}
