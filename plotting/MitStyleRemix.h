#ifndef MITSTYLEREMIX
#define MITSTYLEREMIX

#include <TCanvas.h>
#include <TH1.h>
#include <TPad.h>

void     MitStyleRemix();
TCanvas* MakeCanvas   (const char* name, const char *title, int dX = 500, int dY = 500);
void     InitSubPad   (TPad* pad, int i);
void     InitHist     (TH1 *hist, const char *xtit, const char *ytit  = "Number of Entries",
		       EColor color = kBlack);
void     SetStyle     ();

#endif
