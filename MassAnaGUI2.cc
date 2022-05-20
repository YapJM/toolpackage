// Mainframe macro generated from application: /home/xian/rootbuild/root_v6_22_06/bin/root.exe
// By ROOT version 6.22/06 on 2020-12-31 09:54:37

#ifndef ROOT_TGDockableFrame
#include "TGDockableFrame.h"
#endif
#ifndef ROOT_TGMdiDecorFrame
#include "TGMdiDecorFrame.h"
#endif
#ifndef ROOT_TG3DLine
#include "TG3DLine.h"
#endif
#ifndef ROOT_TGMdiFrame
#include "TGMdiFrame.h"
#endif
#ifndef ROOT_TGMdiMainFrame
#include "TGMdiMainFrame.h"
#endif
#ifndef ROOT_TGMdiMenu
#include "TGMdiMenu.h"
#endif
#ifndef ROOT_TGColorDialog
#include "TGColorDialog.h"
#endif
#ifndef ROOT_TGListBox
#include "TGListBox.h"
#endif
#ifndef ROOT_TGNumberEntry
#include "TGNumberEntry.h"
#endif
#ifndef ROOT_TGScrollBar
#include "TGScrollBar.h"
#endif
#ifndef ROOT_TGComboBox
#include "TGComboBox.h"
#endif
#ifndef ROOT_TGuiBldHintsEditor
#include "TGuiBldHintsEditor.h"
#endif
#ifndef ROOT_TRootBrowser
#include "TRootBrowser.h"
#endif
#ifndef ROOT_TGuiBldNameFrame
#include "TGuiBldNameFrame.h"
#endif
#ifndef ROOT_TGFrame
#include "TGFrame.h"
#endif
#ifndef ROOT_TGMenu
#include "TGMenu.h"
#endif
#ifndef ROOT_TGFileDialog
#include "TGFileDialog.h"
#endif
#ifndef ROOT_TGShutter
#include "TGShutter.h"
#endif
#ifndef ROOT_TGButtonGroup
#include "TGButtonGroup.h"
#endif
#ifndef ROOT_TGCommandPlugin
#include "TGCommandPlugin.h"
#endif
#ifndef ROOT_TGCanvas
#include "TGCanvas.h"
#endif
#ifndef ROOT_TGFSContainer
#include "TGFSContainer.h"
#endif
#ifndef ROOT_TGuiBldEditor
#include "TGuiBldEditor.h"
#endif
#ifndef ROOT_TGColorSelect
#include "TGColorSelect.h"
#endif
#ifndef ROOT_TGTextEdit
#include "TGTextEdit.h"
#endif
#ifndef ROOT_TGButton
#include "TGButton.h"
#endif
#ifndef ROOT_TRootContextMenu
#include "TRootContextMenu.h"
#endif
#ifndef ROOT_TGFSComboBox
#include "TGFSComboBox.h"
#endif
#ifndef ROOT_TGLabel
#include "TGLabel.h"
#endif
#ifndef ROOT_TGView
#include "TGView.h"
#endif
#ifndef ROOT_TGProgressBar
#include "TGProgressBar.h"
#endif
#ifndef ROOT_TGMsgBox
#include "TGMsgBox.h"
#endif
#ifndef ROOT_TRootGuiBuilder
#include "TRootGuiBuilder.h"
#endif
#ifndef ROOT_TGFileBrowser
#include "TGFileBrowser.h"
#endif
#ifndef ROOT_TGTab
#include "TGTab.h"
#endif
#ifndef ROOT_TGListView
#include "TGListView.h"
#endif
#ifndef ROOT_TGSplitter
#include "TGSplitter.h"
#endif
#ifndef ROOT_TGTextEditor
#include "TGTextEditor.h"
#endif
#ifndef ROOT_TRootCanvas
#include "TRootCanvas.h"
#endif
#ifndef ROOT_TGStatusBar
#include "TGStatusBar.h"
#endif
#ifndef ROOT_TGListTree
#include "TGListTree.h"
#endif
#ifndef ROOT_TGuiBldGeometryFrame
#include "TGuiBldGeometryFrame.h"
#endif
#ifndef ROOT_TGToolTip
#include "TGToolTip.h"
#endif
#ifndef ROOT_TGToolBar
#include "TGToolBar.h"
#endif
#ifndef ROOT_TRootEmbeddedCanvas
#include "TRootEmbeddedCanvas.h"
#endif
#ifndef ROOT_TCanvas
#include "TCanvas.h"
#endif
#ifndef ROOT_TGuiBldDragManager
#include "TGuiBldDragManager.h"
#endif
#ifndef ROOT_TGHtmlBrowser
#include "TGHtmlBrowser.h"
#endif

#include "Riostream.h"


#include "TROOT.h"
#include <vector>
#include <thread>
#include "DriftCorrect2.C"
#include "TCanvas.h"

using namespace std;


class MassAnaGUI : public TGMainFrame {

private:
   TGMainFrame *fMainFrame2609;
   TGCompositeFrame *fMainFrame1933, *fMainFrame1817;
   TGHorizontalFrame *MassAnalysis_v1;
   TGFont *ufont;
   TGGC   *uGC;
   TGTextEntry *fTextEntryPATH ,*fTextEntryFilename;
   TGLabel *fLabelPATH , *fLabelFilename;
   TGTextButton     *fTextButtonDecoder, *fTextButtonLoadTree, *fTextButtonExit;
   TGRadioButton *fTextButtonOriginal , *fTextButtonCorrected;
   //TGHProgressBar *fHProgressBar1546;

   // drift correct part
   TGGroupFrame *fGroupFrameDC;
   TGNumberEntry *fNumberEntryRefCento;
   TGNumberEntry *fNumberEntryNbins;
   TGNumberEntry *fNumberEntryEveAkill;
   TGNumberEntry *fNumberEntryRefSigma;
   TGNumberEntry *fNumberEntryHisto_halfwidth;
   
   TGLabel *fLabelRefCenter;  // share
   TGLabel *fLabelnbins;		
   TGLabel *fLabelEveAkill;
   TGLabel *fLabelRefSigma;
   TGLabel *fLabelHalf_histoWidth; // share
   TGTextButton *fTextButtonOK;

     //T0 entry;
   TGLabel *fLabelT0;
   TGNumberEntry *fNumberEntryT0;
//thread* t;
   TCanvas* cc;
   bool runcorrect;
   void runthread(string _PAHT,string _filename,double _centro, int _nbins, int _event_akill, double _ref_sigma,double _Half_hiswidth,double __intTo,TGProgressBar* _bar);
   void exc(string _PAHT="",string _filename="",double _centro=0, int _nbins=0, int _event_akill=0, double _ref_sigma=0,double _Half_hiswidth=0,double __intTo=0,TGProgressBar* _bar=NULL);


public:
   MassAnaGUI();
   virtual ~MassAnaGUI();
   // slots
	void ActionSetDriftON();
	void ActionSetDriftOFF();
	void ActivateDecoder();
	void ActivateDrifCorrect();
	void CloseWindow();

	TGHProgressBar *fHProgressBar1546;

   ClassDef(MassAnaGUI, 0)
};


MassAnaGUI :: MassAnaGUI():cc(NULL),runcorrect(false)
{

   // main frame
	
   fMainFrame2609 = new TGMainFrame(gClient->GetRoot(),10,10,kMainFrame | kVerticalFrame);
   fMainFrame2609->SetCleanup(kDeepCleanup);
   fMainFrame2609->SetName("fMainFrame2609");

  // composite frame
   fMainFrame1933 = new TGCompositeFrame(fMainFrame2609,781,452,kVerticalFrame);
   fMainFrame1933->SetName("fMainFrame1933");
   fMainFrame1933->SetLayoutBroken(kTRUE);


// composite frame
   fMainFrame1817 = new TGCompositeFrame(fMainFrame1933,876,539,kVerticalFrame);
   fMainFrame1817->SetName("fMainFrame1817");

   // horizontal frame
   MassAnalysis_v1 = new TGHorizontalFrame(fMainFrame1817,874,537,kVerticalFrame);
   MassAnalysis_v1->SetName("MassAnalysis_v1");
   MassAnalysis_v1->SetLayoutBroken(kTRUE);

   //TGFont *ufont;         // will reflect user font changes
   ufont = gClient->GetFont("-*-helvetica-medium-r-*-*-12-*-*-*-*-*-iso8859-1");

   //TGGC   *uGC;           // will reflect user GC changes
   // PATH input bar
   GCValues_t valEntryPATH;
   valEntryPATH.fMask = kGCForeground | kGCBackground | kGCFillStyle | kGCFont | kGCGraphicsExposures;
   gClient->GetColorByName("#000000",valEntryPATH.fForeground);
   gClient->GetColorByName("#e8e8e8",valEntryPATH.fBackground);
   valEntryPATH.fFillStyle = kFillSolid;
   valEntryPATH.fFont = ufont->GetFontHandle();
   valEntryPATH.fGraphicsExposures = kFALSE;
   uGC = gClient->GetGC(&valEntryPATH, kTRUE);
   fTextEntryPATH = new TGTextEntry(MassAnalysis_v1, new TGTextBuffer(14),-1,uGC->GetGC(),ufont->GetFontStruct(),kSunkenFrame | kOwnBackground);
   fTextEntryPATH->SetMaxLength(4096);
   fTextEntryPATH->SetAlignment(kTextLeft);
   fTextEntryPATH->SetText("../");
   fTextEntryPATH->Resize(296,fTextEntryPATH->GetDefaultHeight());
   MassAnalysis_v1->AddFrame(fTextEntryPATH, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
   fTextEntryPATH->MoveResize(120,8,296,20);

   // Text label of "PATH" & "filename"
   fLabelPATH = new TGLabel(MassAnalysis_v1,"PATH");
   fLabelPATH->SetTextJustify(36);
   fLabelPATH->SetMargins(0,0,0,0);
   fLabelPATH->SetWrapLength(-1);
   MassAnalysis_v1->AddFrame(fLabelPATH, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
   fLabelPATH->MoveResize(0,0,112,32);
   fLabelFilename = new TGLabel(MassAnalysis_v1,"File Name");
   fLabelFilename->SetTextJustify(36);
   fLabelFilename->SetMargins(0,0,0,0);
   fLabelFilename->SetWrapLength(-1);
   MassAnalysis_v1->AddFrame(fLabelFilename, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
   fLabelFilename->MoveResize(8,32,96,32);

   ufont = gClient->GetFont("-*-helvetica-medium-r-*-*-12-*-*-*-*-*-iso8859-1");

   // Filename input bar
   GCValues_t valEntryFilename;
   valEntryFilename.fMask = kGCForeground | kGCBackground | kGCFillStyle | kGCFont | kGCGraphicsExposures;
   gClient->GetColorByName("#000000",valEntryFilename.fForeground);
   gClient->GetColorByName("#e8e8e8",valEntryFilename.fBackground);
   valEntryFilename.fFillStyle = kFillSolid;
   valEntryFilename.fFont = ufont->GetFontHandle();
   valEntryFilename.fGraphicsExposures = kFALSE;
   uGC = gClient->GetGC(&valEntryFilename, kTRUE);
   fTextEntryFilename = new TGTextEntry(MassAnalysis_v1, new TGTextBuffer(14),-1,uGC->GetGC(),ufont->GetFontStruct(),kSunkenFrame | kOwnBackground);
   fTextEntryFilename->SetMaxLength(4096);
   fTextEntryFilename->SetAlignment(kTextLeft);
   fTextEntryFilename->SetText("Input filename without .lst or .root");
   fTextEntryFilename->Resize(296,fTextEntryFilename->GetDefaultHeight());
   MassAnalysis_v1->AddFrame(fTextEntryFilename, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
   fTextEntryFilename->MoveResize(120,32,296,20);

   // "decorder" and "load tree" button
   fTextButtonDecoder = new TGTextButton(MassAnalysis_v1,"Decoder",-1,TGTextButton::GetDefaultGC()(),TGTextButton::GetDefaultFontStruct(),kRaisedFrame);
   fTextButtonDecoder->SetTextJustify(36);
   fTextButtonDecoder->SetMargins(0,0,0,0);
   fTextButtonDecoder->SetWrapLength(-1);
   fTextButtonDecoder->Resize(120,24);
   fTextButtonDecoder->Connect("Clicked()", "MassAnaGUI", this, "ActivateDecoder()");
   MassAnalysis_v1->AddFrame(fTextButtonDecoder, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
   fTextButtonDecoder->MoveResize(440,8,120,24);
   fTextButtonLoadTree = new TGTextButton(MassAnalysis_v1,"Load Tree",-1,TGTextButton::GetDefaultGC()(),TGTextButton::GetDefaultFontStruct(),kRaisedFrame);
   fTextButtonLoadTree->SetTextJustify(36);
   fTextButtonLoadTree->SetMargins(0,0,0,0);
   fTextButtonLoadTree->SetWrapLength(-1);
   fTextButtonLoadTree->Resize(120,24);
   MassAnalysis_v1->AddFrame(fTextButtonLoadTree, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
   fTextButtonLoadTree->MoveResize(440,32,120,24);

   // "original" & "correct" button
   fTextButtonOriginal = new TGRadioButton(MassAnalysis_v1,"Run ON");
   fTextButtonOriginal->SetTextJustify(36);
   fTextButtonOriginal->SetMargins(0,0,0,0);
   fTextButtonOriginal->SetWrapLength(-1);
   fTextButtonOriginal->Connect("Clicked()", "MassAnaGUI", this, "ActionSetDriftON()");
   MassAnalysis_v1->AddFrame(fTextButtonOriginal, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
   fTextButtonOriginal->MoveResize(576,32,113,24);

   fTextButtonCorrected = new TGRadioButton(MassAnalysis_v1,"Run OFF");
   fTextButtonCorrected->SetState(kButtonDown);
   fTextButtonCorrected->SetTextJustify(36);
   fTextButtonCorrected->SetMargins(0,0,0,0);
   fTextButtonCorrected->SetWrapLength(-1);
   fTextButtonCorrected->Connect("Clicked()", "MassAnaGUI", this, "ActionSetDriftOFF()");
   MassAnalysis_v1->AddFrame(fTextButtonCorrected, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
   fTextButtonCorrected->MoveResize(696,32,113,24);

   ufont = gClient->GetFont("-*-helvetica-medium-r-*-*-12-*-*-*-*-*-iso8859-1");

   //T0 input
   fNumberEntryT0 = new TGNumberEntry(MassAnalysis_v1, (Double_t) 0,3,-1,(TGNumberFormat::EStyle) 5);
   fNumberEntryT0->SetName("fNumberEntryT0");
   MassAnalysis_v1->AddFrame(fNumberEntryT0, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
   fNumberEntryT0->MoveResize(735,66,40,20);
   fLabelT0= new TGLabel(MassAnalysis_v1,"T0 [ns]");
   fLabelT0->SetTextJustify(36);
   fLabelT0->SetMargins(0,0,0,0);
   fLabelT0->SetWrapLength(-1);
   MassAnalysis_v1->AddFrame(fLabelT0, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
   fLabelT0->MoveResize(688,64,40,24);

   // Exit button
   GCValues_t valButtonExit;
   valButtonExit.fMask = kGCForeground | kGCBackground | kGCFillStyle | kGCFont | kGCGraphicsExposures;
   gClient->GetColorByName("#ff0000",valButtonExit.fForeground);
   gClient->GetColorByName("#e8e8e8",valButtonExit.fBackground);
   valButtonExit.fFillStyle = kFillSolid;
   valButtonExit.fFont = ufont->GetFontHandle();
   valButtonExit.fGraphicsExposures = kFALSE;
   uGC = gClient->GetGC(&valButtonExit, kTRUE);
   fTextButtonExit = new TGTextButton(MassAnalysis_v1,"Exit",-1,uGC->GetGC(),ufont->GetFontStruct(),kRaisedFrame);
   fTextButtonExit->SetTextJustify(36);
   fTextButtonExit->SetMargins(0,0,0,0);
   fTextButtonExit->SetWrapLength(-1);
   fTextButtonExit->Resize(99,24);
   fTextButtonExit->Connect("Clicked()", "MassAnaGUI", this, "CloseWindow()");
   MassAnalysis_v1->AddFrame(fTextButtonExit, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
   fTextButtonExit->MoveResize(680,6,99,24);

   // ProgressBar
   TGProgressBar::EBarType ProgressBarType = TGProgressBar::EBarType::kStandard;
   fHProgressBar1546 = new TGHProgressBar(MassAnalysis_v1,ProgressBarType,296);  
   fHProgressBar1546->SetName("fHProgressBar1546");
   fHProgressBar1546->SetFillType(TGProgressBar::kBlockFill);

   ULong_t ucolor;        // will reflect user color changes
   gClient->GetColorByName("#ffffff",ucolor);
   fHProgressBar1546->SetBackgroundColor(ucolor);
   fHProgressBar1546->SetPosition(5.);  // Set Bar "Percentage"!!!!!!!!!!
   fHProgressBar1546->SetBarColor("#0000ff");
   MassAnalysis_v1->AddFrame(fHProgressBar1546, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
   fHProgressBar1546->MoveResize(120,56,296,17);


ufont = gClient->GetFont("-*-helvetica-medium-r-*-*-12-*-*-*-*-*-iso8859-1");

   // graphics context changes
   GCValues_t valpFrame1634;
   valpFrame1634.fMask = kGCForeground | kGCBackground | kGCFillStyle | kGCFont | kGCGraphicsExposures;
   gClient->GetColorByName("#000000",valpFrame1634.fForeground);
   gClient->GetColorByName("#e8e8e8",valpFrame1634.fBackground);
   valpFrame1634.fFillStyle = kFillSolid;
   valpFrame1634.fFont = ufont->GetFontHandle();
   valpFrame1634.fGraphicsExposures = kFALSE;
   uGC = gClient->GetGC(&valpFrame1634, kTRUE);

   // "Drift Correct" group frame and number for input
   fGroupFrameDC = new TGGroupFrame(MassAnalysis_v1,"Drift Correct",kVerticalFrame,uGC->GetGC());
   fGroupFrameDC->SetLayoutBroken(kTRUE);
   fNumberEntryRefCento = new TGNumberEntry(fGroupFrameDC, (Double_t) 1e+07,6,-1,(TGNumberFormat::EStyle) 5);
   fGroupFrameDC->AddFrame(fNumberEntryRefCento, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
   fNumberEntryRefCento->MoveResize(99,16,72,20);
   fNumberEntryNbins = new TGNumberEntry(fGroupFrameDC, (Double_t) 200,6,-1,(TGNumberFormat::EStyle) 5);
   fGroupFrameDC->AddFrame(fNumberEntryNbins, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
   fNumberEntryNbins->MoveResize(232,16,58,20);
   fNumberEntryEveAkill = new TGNumberEntry(fGroupFrameDC, (Double_t) 250,6,-1,(TGNumberFormat::EStyle) 5);
   fGroupFrameDC->AddFrame(fNumberEntryEveAkill, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY));
   fNumberEntryEveAkill->MoveResize(368,16,56,20);
   fNumberEntryRefSigma = new TGNumberEntry(fGroupFrameDC, (Double_t) 10,6,-1,(TGNumberFormat::EStyle) 5);
   fGroupFrameDC->AddFrame(fNumberEntryRefSigma, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY));
   fNumberEntryRefSigma->MoveResize(504,16,56,20);
   fNumberEntryHisto_halfwidth = new TGNumberEntry(fGroupFrameDC, (Double_t) 150,6,-1,(TGNumberFormat::EStyle) 5);
   fGroupFrameDC->AddFrame(fNumberEntryHisto_halfwidth, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY));
   fNumberEntryHisto_halfwidth->MoveResize(658,15,56,20);


   ufont = gClient->GetFont("-*-helvetica-medium-r-*-*-12-*-*-*-*-*-iso8859-1");

   // graphics context changes
   GCValues_t vall2224;
   vall2224.fMask = kGCForeground | kGCBackground | kGCFillStyle | kGCFont | kGCGraphicsExposures;
   gClient->GetColorByName("#000000",vall2224.fForeground);
   gClient->GetColorByName("#e8e8e8",vall2224.fBackground);
   vall2224.fFillStyle = kFillSolid;
   vall2224.fFont = ufont->GetFontHandle();
   vall2224.fGraphicsExposures = kFALSE;
   uGC = gClient->GetGC(&vall2224, kTRUE);

   //Drift correct group label and "OK" buttorn
   fLabelRefCenter = new TGLabel(fGroupFrameDC,"Ref. center[ns]",uGC->GetGC());
   fLabelRefCenter->SetTextJustify(36);
   fLabelRefCenter->SetMargins(0,0,0,0);
   fLabelRefCenter->SetWrapLength(-1);
   fGroupFrameDC->AddFrame(fLabelRefCenter, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
   fLabelRefCenter->MoveResize(8,16,80,24);

   fLabelnbins = new TGLabel(fGroupFrameDC,"nbins");
   fLabelnbins->SetTextJustify(36);
   fLabelnbins->SetMargins(0,0,0,0);
   fLabelnbins->SetWrapLength(-1);
   fGroupFrameDC->AddFrame(fLabelnbins, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
   fLabelnbins->MoveResize(174,20,48,16);

   fLabelEveAkill = new TGLabel(fGroupFrameDC,"EveAkill");
   fLabelEveAkill->SetTextJustify(36);
   fLabelEveAkill->SetMargins(0,0,0,0);
   fLabelEveAkill->SetWrapLength(-1);
   fGroupFrameDC->AddFrame(fLabelEveAkill, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
   fLabelEveAkill->MoveResize(312,16,48,24);

   fLabelRefSigma = new TGLabel(fGroupFrameDC,"Sigma");
   fLabelRefSigma->SetTextJustify(36);
   fLabelRefSigma->SetMargins(0,0,0,0);
   fLabelRefSigma->SetWrapLength(-1);
   fGroupFrameDC->AddFrame(fLabelRefSigma, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
   fLabelRefSigma->MoveResize(447,15,56,24);

   fLabelHalf_histoWidth = new TGLabel(fGroupFrameDC,"Half_histoWidth");
   fLabelHalf_histoWidth->SetTextJustify(36);
   fLabelHalf_histoWidth->SetMargins(0,0,0,0);
   fLabelHalf_histoWidth->SetWrapLength(-1);
   fGroupFrameDC->AddFrame(fLabelHalf_histoWidth, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
   fLabelHalf_histoWidth->MoveResize(561,16,96,22);

   fTextButtonOK = new TGTextButton(fGroupFrameDC,"OK",-1,TGTextButton::GetDefaultGC()(),TGTextButton::GetDefaultFontStruct(),kRaisedFrame);
   fTextButtonOK->SetTextJustify(36);
   fTextButtonOK->SetMargins(0,0,0,0);
   fTextButtonOK->SetWrapLength(-1);
   fTextButtonOK->Resize(38,22);
   fTextButtonOK->Connect("Clicked()", "MassAnaGUI", this, "ActivateDrifCorrect()");
   fGroupFrameDC->AddFrame(fTextButtonOK, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
   fTextButtonOK->MoveResize(724,16,38,22);




   // main frame setting....

   fGroupFrameDC->SetLayoutManager(new TGVerticalLayout(fGroupFrameDC));
   fGroupFrameDC->Resize(770,56);
   MassAnalysis_v1->AddFrame(fGroupFrameDC, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
   fGroupFrameDC->MoveResize(8,88,770,56);

   fMainFrame1817->AddFrame(MassAnalysis_v1, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY,1,1,1,1));
   fMainFrame1933->AddFrame(fMainFrame1817, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY));
   fMainFrame1817->MoveResize(0,0,876,539);

   fMainFrame2609->AddFrame(fMainFrame1933, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY));
   fMainFrame1933->MoveResize(0,0,781,452);


   fMainFrame2609->SetMWMHints(kMWMDecorAll,
                        kMWMFuncAll,
                        kMWMInputModeless);
   fMainFrame2609->SetWindowName("MassAnalyser");
   fMainFrame2609->MapSubwindows();

   fMainFrame2609->Resize(fMainFrame2609->GetDefaultSize());
   fMainFrame2609->MapWindow();
   fMainFrame2609->Resize(876,539);
}  


MassAnaGUI :: ~MassAnaGUI(){

	MassAnalysis_v1->Cleanup();
	fMainFrame1817->Cleanup();
	fMainFrame1933->Cleanup();
	fMainFrame2609->Cleanup();
}


void MassAnaGUI ::ActionSetDriftON(){

	fTextButtonOriginal->SetState(kButtonDown);
	fTextButtonCorrected->SetState(kButtonUp);
   ::abort_corr_loop=false;
	std::cout<<"\e[1;33m"<<"Enable drift correct"<<"\e[0m"<<std::endl;

}


void MassAnaGUI ::ActionSetDriftOFF(){

	fTextButtonOriginal->SetState(kButtonUp);
	fTextButtonCorrected->SetState(kButtonDown);
   ::abort_corr_loop=true;
	std::cout<<"\e[1;33m"<<"Stop dirft correct"<<"\e[0m"<<std::endl;

}

void MassAnaGUI:: CloseWindow(){ // close the main window return to root CINT
   ActionSetDriftOFF();
   fMainFrame2609->DeleteWindow();
}

void MassAnaGUI:: ActivateDecoder(){
	string LstFilePATH = fTextEntryPATH->GetText();
	string LstFilename = fTextEntryFilename->GetText();
	gROOT->ProcessLine( Form("ParserMCS6A_2(\"%s\",\"%s\")",LstFilePATH.c_str(),LstFilename.c_str()) );

}

void MassAnaGUI::ActivateDrifCorrect(){
   ActionSetDriftON();
	string LstFilePATH = fTextEntryPATH->GetText();
	string LstFilename = fTextEntryFilename->GetText();

   double T0 = fNumberEntryT0->GetNumber();
   double RefCentoTof = fNumberEntryRefCento->GetNumber();
   int Nbins = fNumberEntryNbins->GetNumber();
   int EveAkill = fNumberEntryEveAkill->GetNumber();
   double RefSigma = fNumberEntryRefSigma->GetNumber();
   double Histo_halfwidth = fNumberEntryHisto_halfwidth->GetNumber();
   fHProgressBar1546->Reset();
   std::cout<<"\e[1;33m"<<"T0 is set to "<<T0<<" [ns]"<<"\e[0m"<<std::endl;
   //const char * command = Form("DriftCorrect(\"%s\",\"%s\",%f,%i,%i,%f,%f,%f,ui->fHProgressBar1546)",LstFilePATH.c_str(),LstFilename.c_str(),RefCentoTof,Nbins,EveAkill,RefSigma,Histo_halfwidth,T0); 
   //gROOT->ProcessLine(command);
   //fHProgressBar1546->Reset();
   //if(t!=NULL)delete t;
  thread t(&MassAnaGUI ::runthread,this,LstFilePATH,LstFilename,RefCentoTof,Nbins,EveAkill,RefSigma,Histo_halfwidth,T0,fHProgressBar1546);
  // thread t(::DriftCorrect,LstFilePATH,LstFilename,RefCentoTof,Nbins,EveAkill,RefSigma,Histo_halfwidth,T0,fHProgressBar1546);
   t.detach();
   gROOT->ProcessLine("\n");
   return;

/*
   if(cc!=NULL) delete cc;
   gROOT->SetBatch();
   cc = new TCanvas();gROOT->SetBatch(kFALSE);
   exc(LstFilePATH,LstFilename,RefCentoTof,Nbins,EveAkill,RefSigma,Histo_halfwidth,T0,fHProgressBar1546);
   runcorrect=true;
   cc->AddExec("exc","exc()");
   //gROOT->SetBatch(kFALSE);*/


}
//string PATH="../",string filename = "mcs_39Kvs143X2plus@500_174927", double centro=17e6, int nbins = 200, int event_akill =600, double ref_sigma=10,double Half_hiswidth = 150,double _inT0=130,TGProgressBar *bar


void MassAnaGUI ::runthread(string _PAHT,string _filename,double _centro, int _nbins, int _event_akill, double _ref_sigma,double _Half_hiswidth,double __intTo,TGProgressBar* _bar){
   ::DriftCorrect(_PAHT, _filename, _centro, _nbins, _event_akill, _ref_sigma, _Half_hiswidth, __intTo, _bar);
   fHProgressBar1546->RaiseWindow();
   //sleep(1);
   gROOT->ProcessLine("cout<<\"done\"<<endl;");
   //fHProgressBar1546->Reset();
}

void MassAnaGUI ::exc(string _PAHT,string _filename,double _centro, int _nbins, int _event_akill, double _ref_sigma,double _Half_hiswidth,double __intTo,TGProgressBar* _bar){
   static string _PAHT_;
   static string _filename_;
   static double _centro_;
   static int _nbins_;
   static int _event_akill_;
   static double _ref_sigma_;
   static double _Half_hiswidth_;
   static double __intTo_;
   static TGProgressBar* _bar_;

   if(_PAHT!=""){
      _PAHT_ = _PAHT;
      _filename_ = _filename;
      _centro_ = _centro;
      _nbins_ = _nbins;
      _event_akill_ = _event_akill;
      _ref_sigma_ = _ref_sigma;
      _Half_hiswidth_ = _Half_hiswidth;
      __intTo_ = __intTo;
      _bar_ = _bar;
   }


   if(runcorrect){
      ::DriftCorrect(_PAHT_, _filename_, _centro_, _nbins_, _event_akill_, _ref_sigma_, _Half_hiswidth_, __intTo_, _bar_);
      runcorrect=false;
   }
   else{;}
}


MassAnaGUI * ui = NULL;

MassAnaGUI * MassAnaGUI2(){
	gROOT->ProcessLine(".L ParserMCS6A_2.C+");
//	gROOT->ProcessLine(".L DriftCorrect2.C+");
	ui = new MassAnaGUI();
	return ui;

}
