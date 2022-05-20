#include "TROOT.h"
#include "TRint.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TSpectrum.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TGraph.h"

//#include "TProfile.h"
//#include "TAxis.h"

#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include "TGProgressBar.h"
#include "funcJohnson.h"

//#define _LASER_


#if defined (MAKECINT)
#pragma link C++ class vector<Long_t>+;   // in order to identify vector<Long_t> 
#endif


using namespace std;

bool abort_corr_loop=false;

double* FeaturePar(TH1D* _inhisto,Long64_t _min_evt, Long64_t _max_evt);
double fit_result(TTree *t1, Long64_t min_evt, Long64_t max_evt, int nbins, double lowx, double highx,Double_t Sigma);
TH1D *h;
TF1 *func=NULL;
funcJohnson* john=NULL;
//TCanvas *c1 = new TCanvas("c1","c1",1000,800);
//TCanvas *c_p1_p2 = new TCanvas("c_p1_p2","t0 correction",1200,1000);
TGraph *g_mean_p1 = new TGraph();
TGraph *g_mean_p2 = new TGraph();
TGraph *g_sigma_p1 = new TGraph();
TGraph *g_sigma_p2 = new TGraph();

double event_akill_global=0;
bool isnewfile=true;

Long64_t stop_point = -1;  // stop_point marker set to -1
int N_byfit=0;      // Number of slices by fit
int N_byappro=0;    // Number of slices by approximating

Long_t empty_slice_gsweeps_L=0; // gsweeps of empty slice
Long_t empty_slice_gsweeps_R=0; // gsweeps of empty slice

int DriftCorrect(string PATH="../",string filename = "mcs_39Kvs143X2plus@500_174927", double centro=17e6, int nbins = 200, int event_akill =600, double ref_sigma=10,double Half_hiswidth = 150,double _inT0=130,TGProgressBar *bar=NULL){

      string inputfile = PATH + "rootfiles/" + filename + ".root";
      TFile *fin = new TFile(inputfile.c_str(),"READ");
      if(!fin->IsOpen()) {cout<<"Can not open file: "<<inputfile<<endl; cout<<"Break!!"<<endl; return 0;}
      TTree *intree = (TTree *)fin->Get("tree");
      TTree* tree_copy = intree;

      //%%%%%% initialize %%%%%%%%%%%%%%%
      isnewfile=true;
 
	event_akill_global = event_akill;

      if(john != NULL) delete john;
      john = new funcJohnson();

      func = john->GetTF1();

      N_byfit=0;
      N_byappro=0;
      //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      cout<<"\e[1;33m reading data from "<<inputfile<<"\e[0m"<<endl;


  Long64_t nentries = intree->GetEntriesFast();
  cout<<"nentries = "<<nentries<<endl;
  //Long64_t stop_point = -1;        // stop_point marker set to -1
  stop_point=-1;
  int counter = 1;    // display purpose 
  int SlideWidth = event_akill;  // number of events to load in each time
  int num_of_SlideWidth = (int) (nentries / SlideWidth);
  cout<<"num_of_SlideWidth = "<<num_of_SlideWidth<<endl;
  double binlowedge1 = centro - Half_hiswidth;    double binlowedge2 = centro - Half_hiswidth;// total histowidth should be 300 usually
  double binupedge1 = centro + Half_hiswidth;     double binupedge2 = centro + Half_hiswidth;
  int    nbins1    = nbins;           int    nbins2    = nbins;
  double sigma_ref_measured1 = ref_sigma;   double sigma_ref_measured2 = ref_sigma; //;6
  int accumu_num = 0;
  double reference1 =0;              double reference2 =0;           // fit first slidewidth data as reference
  double n_slide_central1 = 0;       double n_slide_central2 = 0;
  int    node = 0 ;
                                           //(1000,0,794000,1000,9.006e6,9.007e6) ;  1600,0,793001,100,9.0064e6,9.0067e6                             
                                            //(1000,0,1771000,1000,9.019e6,9.021e6)
                                             //(1000,0,181000,1000,9.0423e6,9.0426e6);

  bool fitTrue = false;
  bool fitTrue3 = false;
//  double offset_stop =0;
//  vector <Double_t> *offset = new vector <Double_t>();

  const int t0_num = 1;         // check effect of 11 kinds of daq delay t0 to drift correction
  // !!!!!!!!!!!! modified to allow setting a single t0
  // in the following:  
  //double t0 = index*t0_delta + t0_delta; index=[0,t0_num) => index=0 only; t0 = t0_delta finally; 
  // by setting t0_delta = setting t0
  ///////////////////////////////////
  double t0_delta = _inT0; //130;    // step size of t0 = 20 ns for multi t0 test situation, t0 delay is about 130 ns usually
  double offset_stop[t0_num] ={0};
  vector <Double_t> **offset = new vector <Double_t> *[t0_num];

	for(int index =0;index<t0_num;index++){offset[index] = new vector <Double_t>[1];}

  if(bar!=NULL) bar->Reset();
 
  for(Long64_t jentry=0; jentry<nentries; jentry++){       //nentries
    if(abort_corr_loop) break;

    if(jentry == (int)counter*nentries/20){  // divide into 20 piese , each piese is equivalent to 5%
      cout<< '\r' << "finish running : "<< counter*5 << "%" <<"\t ratio by fit= "<<N_byfit/(N_byfit+2.0*N_byappro)*100.<<" %";
      cout<<"\t ratio by approximate= "<<2*N_byappro/(N_byfit+2.0*N_byappro)*100.<<" %"<<flush; 
      if(bar != NULL){bar->SetPosition((Float_t)counter*5.0); bar->RaiseWindow();} 
      counter++;
     }

//     cout<<endl;
 //    cout<<"1_nevt = "<<num<<endl; //return 1;
       
     if(jentry == 1){
           cout<<"continue or not: 'n' for break!!!!"<<endl;
             char go;
             cin>>go;
            if( go =='n'){cerr<<"\e[1;33m"<<"program break!!!!!"<<"\e[0m"<<endl; return 0;}
            cout<<"Pending......"<<endl;
      }

    // correct drift time
    if((jentry < SlideWidth) && fitTrue == false){//%%%%%%%%%%% first piece
    cout<<"do 1"<<endl;
    
   // cout<<"do 1 num1 , nevt =  "<<num<< " , "<<nevt<<endl;     
  
      reference1 = fit_result(tree_copy,0,((Long64_t) SlideWidth),nbins1,binlowedge1,binupedge1,sigma_ref_measured1);  //(TTree *t1, Long64_t min_evt, Long64_t max_evt, 
              // cout<<"do 1 num2, nevt = "<<num << " , "<<nevt<<endl;
               cout<<"referen 1 = "<<reference1<<endl;                                         //int nbins, double lowx, double highx);
      reference2 = fit_result(tree_copy,0,((Long64_t) SlideWidth),nbins2,binlowedge2,binupedge2,sigma_ref_measured2);
               cout<<"referen 2 = "<<reference2<<endl; 
            //offset->push_back(1);  // shift==0; ration ==1
         for(int index =0;index<t0_num;index++){offset[index]->push_back(1);}
       fitTrue = true;
                 
    }
    else if( (node<(int)(jentry/SlideWidth)) && (int)(jentry/SlideWidth)< num_of_SlideWidth){ //%%%%%%%%% middle piece
   // cout<<"do 2"<<endl;

          if(stop_point < 0){
              n_slide_central1 = fit_result(tree_copy,jentry,(jentry + (Long64_t) SlideWidth),nbins1,binlowedge1,binupedge1,sigma_ref_measured1);// node here delays by 1;
              n_slide_central2 = fit_result(tree_copy,jentry,(jentry + (Long64_t) SlideWidth),nbins2,binlowedge2,binupedge2,sigma_ref_measured2);
          }else{
              n_slide_central1 = fit_result(tree_copy,stop_point,(jentry + (Long64_t) SlideWidth),nbins1,binlowedge1,binupedge1,sigma_ref_measured1); 
              n_slide_central2 = fit_result(tree_copy,stop_point,(jentry + (Long64_t) SlideWidth),nbins2,binlowedge2,binupedge2,sigma_ref_measured2);
          }


          if(n_slide_central1>0 && n_slide_central2>0){   // good fit
                           accumu_num++;
                 for(int akill=0;akill<accumu_num;akill++){
                    // offset_stop = ((n_slide_central1 - reference1)+(n_slide_central2 - reference2))/2;     //update offset record
				//	offset_stop = ((n_slide_central1 / reference1)+(n_slide_central2 / reference2))/2;
                       	//	offset->push_back( offset_stop ) ;  // update offset
              				for(int index=0;index<t0_num;index++){
              					double t0 = index*t0_delta + t0_delta;   // + t0_delta at the end just want t0 to start from 130 ns 
              					offset_stop[index] = ((n_slide_central1-t0)/(reference1-t0) + (n_slide_central2-t0)/(reference2-t0) )/2;
              					//offset_stop[index] = ((n_slide_central1 - reference1)+(n_slide_central2 - reference2))/2; // by difference method
              					//offset_stop[index] = ((reference1 / n_slide_central1)+(reference2 / n_slide_central2))/2;   // by ratio method
                                     		offset[index]->push_back( offset_stop[index] ) ;  // update offset
              				}

                  }

                        accumu_num =0;          // reset to 0;
                        stop_point = -1;       // reset stop_point

          }
          else if(n_slide_central1==-20 || n_slide_central2==-20){ // empty slice!!!!!
                accumu_num++;
                for(int akill=0;akill<accumu_num;akill++){
                          for(int index=0;index<t0_num;index++){
                               // double t0 = index*t0_delta + t0_delta;   // + t0_delta at the end just want t0 to start from 130 ns 
                               // offset_stop[index] = ((n_slide_central1-t0)/(reference1-t0) + (n_slide_central2-t0)/(reference2-t0) )/2;
                              //  offset_stop[index] = offset[index]->at(offset[index]->size()-1);
                                offset[index]->push_back( offset_stop[index] ) ;  
                          }

                }
                cout<<"\e[1;32m"<<"warning: empty slice founded!!!  "<<empty_slice_gsweeps_L<<" ~ "<<empty_slice_gsweeps_R<<"\e[0m"<<endl;

                accumu_num =0;          // reset to 0;
                stop_point = -1;       // reset stop_point
          }
          else{ // bad fit to accumulate more slideWidth
                       if(accumu_num ==0){ stop_point = jentry;}  // set stop_point
                  accumu_num ++;
                /* if(accumu_num>20) {
                      cerr<<"warning !!!!!! find no peak!!! check point, from nevt= "<<stop_point<<"accumulate num = "<<accumu_num<<endl; 
                           return 0;
                     }*/
          }

       //     cout<<"offset = , node = "<<offset_stop[0]<<"  ,  "<<node+1<<".  nevt range = "<<(node+1)*SlideWidth<<" ~ "<<(node+2)*SlideWidth<<endl; 
       //        cout<<endl;
            //cout<<"jentry= "<<jentry<<endl;

            node++;  //code above only execute once at the beginning of each piece; point to next piece       
              

    }
    else if((int)(jentry/SlideWidth) == num_of_SlideWidth && fitTrue3 == false){ //%%%%%%%%%%% last piese
    //   cout<<"do final"<<endl;

          if(stop_point < 0){
             n_slide_central1 = fit_result(tree_copy,jentry,nentries,nbins1,binlowedge1,binupedge1,sigma_ref_measured1);
             n_slide_central2 = fit_result(tree_copy,jentry,nentries,nbins2,binlowedge2,binupedge2,sigma_ref_measured2);
          }else{
              n_slide_central1 = fit_result(tree_copy,stop_point,nentries,nbins1,binlowedge1,binupedge1,sigma_ref_measured1); 
              n_slide_central2 = fit_result(tree_copy,stop_point,nentries,nbins2,binlowedge2,binupedge2,sigma_ref_measured2);
          }

             // %%%%%%%%%%%%   good fit or bad fit %%%%%%%%%%%%%%
          if(n_slide_central1>0 && n_slide_central2>0){   // good fit
                           accumu_num++;
                 for(int akill=0;akill<accumu_num;akill++){
                   //  offset_stop = ((n_slide_central1 - reference1) + (n_slide_central2 - reference2))/2;
			//	offset_stop = ((n_slide_central1 / reference1)+(n_slide_central2 / reference2))/2;
                   //    offset->push_back( offset_stop ) ;  // update offset

              				for(int index=0;index<t0_num;index++){
              					double t0 = index*t0_delta + t0_delta;  // + t0_delta at the end just want t0 to start from 130 ns
              					offset_stop[index] =  ((n_slide_central1-t0)/(reference1-t0) + (n_slide_central2-t0)/(reference2-t0) )/2;
              					//offset_stop[index] = ((n_slide_central1 - reference1)+(n_slide_central2 - reference2))/2;  // by difference method
              					//offset_stop[index] = ((reference1 / n_slide_central1)+(reference2 / n_slide_central2 ))/2;   // by ratio method
                                     		offset[index]->push_back( offset_stop[index] ) ;  // update offset
              				}

                  }
                        accumu_num =0;          // reset to 0;
                        stop_point = -1;       // reset stop_point

          }
         else if(n_slide_central1==-20 || n_slide_central2==-20){ // empty slice!!!!!
                accumu_num++;
                for(int akill=0;akill<accumu_num;akill++){
                          for(int index=0;index<t0_num;index++){
                               // double t0 = index*t0_delta + t0_delta;   // + t0_delta at the end just want t0 to start from 130 ns 
                               // offset_stop[index] = ((n_slide_central1-t0)/(reference1-t0) + (n_slide_central2-t0)/(reference2-t0) )/2;
                               // offset_stop[index] = offset[index]->at(offset[index]->size()-1);
                                offset[index]->push_back( offset_stop[index] ) ;  // update offset with offset from last slice
                          }

                }
                cout<<"\e[1;32m"<<"warning: empty slice founded!!!  "<<empty_slice_gsweeps_L<<" ~ "<<empty_slice_gsweeps_R<<"\e[0m"<<endl;

                accumu_num =0;          // reset to 0;
                stop_point = -1;       // reset stop_point
          }
          else{ // bad fit 

                           accumu_num++;
                 for(int akill=0;akill<accumu_num;akill++){

                     //  offset->push_back( offset_stop ) ;  // use old offset to approximate
                     //  offset->push_back( offset_stop[index] ) ; 
            				for(int index=0;index<t0_num;index++){
            					offset[index]->push_back( offset_stop[index] ) ;  // update offset with offset from last slice
            				}
                  }
          }


           // cout<<"end of fit : offset = "<<offset_stop[0]<<endl;    
          // cout<<"jentry= "<<jentry<<endl; 

              fitTrue3 = true;  // make sure last piece fit only once
        
    }else{;}
  }// end of for offset extraction


  ////%%%%%%%%%%%%%%%%%%%%%%% end of offset extraction %%%%%%%%%%%%%%%%%%%%%%%%

    if(abort_corr_loop){
       cout<<"external interrupt"<<endl;
       return 0;
    }

     cout<<endl;
     cout<<"nentries = "<<nentries<<endl;
     cout<<"offset num = "<< offset[0]->size()<<endl;
     cout<<endl;
     cout<<endl;
   
      if(bar != NULL) bar->SetPosition((Float_t)counter*5.0); 


      Long64_t nevt = 0;
 //     Long64_t num = 0;
      int nhits[10];
//      int glo = 0;
      vector <Int_t> *channel = new vector <Int_t>();
      vector <Int_t> *edge = new vector <Int_t>();
      vector <Int_t> *value = new vector <Int_t>();
      vector <Double_t> *time = new vector <Double_t>();
      vector <Double_t> *timec = new vector <Double_t>();
      vector <Int_t> *sweeps = new vector <Int_t>();
      vector <Long_t> * sweeps_global = new vector <Long_t>;
      vector <Int_t> *glo = new vector <Int_t>();
      vector <Int_t> *tag = new vector <Int_t>();
      vector <Int_t> *lost = new vector <Int_t>();
      #ifdef _LASER_
      vector <Double_t> *time_gcDrift = new vector<Double_t>();
      #endif
      double offset_newtree = 0;
      
      intree->SetBranchAddress("nevt", &nevt); 
//      intree->SetBranchAddress("nevt", &num); 
      intree->SetBranchAddress("nhits", nhits);
      intree->SetBranchAddress("channel", &channel);
      intree->SetBranchAddress("edge", &edge);
      intree->SetBranchAddress("value", &value);
      intree->SetBranchAddress("time", &time);
      intree->SetBranchAddress("sweeps", &sweeps);
      intree->SetBranchAddress("sweeps_global", &sweeps_global);
      intree->SetBranchAddress("glo", &glo);
      intree->SetBranchAddress("tag", &tag);
      intree->SetBranchAddress("lost", &lost);
      #ifdef _LASER_
      intree->SetBranchAddress("time_gcDrift", &time_gcDrift);
      #endif
  // correct drift time

  string outputfile = PATH + "rootfiles/" + filename +"_dcorrect.root";
  TFile *fout = new TFile(outputfile.c_str(),"RECREATE");
//  TTree *tree = new TTree("tree","analyzed tree");
  TTree **tree = new TTree *[t0_num];
  TH2D *h2p1 = new TH2D("h2p1","h2p1",1000,0,794000,200,9.0063e6,9.0068e6);
  TH2D *h2p2 = new TH2D("h2p2","h2p2",1000,0,794000,100,9.0134e6,9.0139e6);
  TH1D *h1p1 = new TH1D("h1p1","h1p1",200,9.0063e6,9.0068e6);
  TH1D *h1p2 = new TH1D("h1p2","h1p2",100,9.0134e6,9.0139e6);


for(int index_t0=0;index_t0<t0_num;index_t0++){
	double t0 = index_t0 * t0_delta + t0_delta;
  tree[index_t0]= new TTree(Form("tree%d",index_t0),"analyzed tree");
  tree[index_t0]->Branch("nevt",&nevt);
//  tree->Branch("nevt1",&num);
  tree[index_t0]->Branch("nhits",nhits,"nhits[10]/I");
  tree[index_t0]->Branch("channel",&channel);
  tree[index_t0]->Branch("edge",&edge);
  tree[index_t0]->Branch("value",&value); // in unit of 100ps
  tree[index_t0]->Branch("time",&time); // in unit of ns
  tree[index_t0]->Branch("timec",&timec); // in unit of ns
  tree[index_t0]->Branch("sweeps",&sweeps);
  tree[index_t0]->Branch("sweeps_global",&sweeps_global);
  tree[index_t0]->Branch("glo",&glo);
  tree[index_t0]->Branch("tag",&tag);
  tree[index_t0]->Branch("lost",&lost);
  tree[index_t0]->Branch("offset",&offset_newtree);
  #ifdef _LASER_
  tree[index_t0]->Branch("time_gcDrift",&time_gcDrift);
  #endif

 
       counter =0;

  for(Long64_t jentry=0; jentry<nentries; jentry++){
    if(jentry == (int)counter*nentries/20){  // divide into 20 piese , each piese is equivalent to 5%
      cout<< '\r' << "finish running : "<< counter*5 << "%" <<flush;
      counter++;
     }
   
       node = (int)(jentry/SlideWidth);

        intree->GetEntry(jentry);
         
     //   offset_newtree = offset->at(node);


       // offset_newtree = offset->at(node);
        offset_newtree = offset[index_t0]->at(node);
     // calculate timec
    double itime=0;
    for(int j=0; j<((int)time->size()); j++){
       itime = time->at(j);
      //timec->push_back(itime-offset_newtree);  // shift method correct
      if(itime==0){timec->push_back(0);}
      else timec->push_back((itime-t0) / offset_newtree + t0);    // ratio method correct
      // cout<<"offset = "<<offset<<", timec = "<<itime-offset<<endl;
    }// calculate timec



    tree[index_t0]->Fill();
    
    channel->clear();
    edge->clear();

    value->clear();
    time->clear();
    timec->clear();
    sweeps->clear();
    sweeps_global->clear();
    glo->clear();
    tag ->clear();
    lost->clear(); 
    #ifdef _LASER_
    time_gcDrift->clear();
    #endif
    for(int i=0; i<10; i++) nhits[i] = 0;


  }// end of for getEvent loop



  cout<<endl;

  cout<<"\e[1;33m writing tree to "<<outputfile<<"\e[0m"<<endl;
  fout->cd();
  tree[index_t0]->Write();
  

  fout->Close();
  fin->Close();

  
   /*for(int index_t0=0;index_t0<t0_num;index_t0++){
     
       delete tree[index_t0];
    }*/

    delete[] tree;

    cout<<endl;
    cout<<"nentries = "<<nentries<<endl;
    cout<<"num_of_SlideWidth = "<<num_of_SlideWidth<<endl;
    cout<<"ratio by fit= "<<N_byfit/(N_byfit+2.0*N_byappro)*100.<<"%";
    cout<<"\t ratio by approximate= "<<N_byappro*2./(N_byfit+2.0*N_byappro)*100.<<"%"<<endl;
    cout<<"Maximum sweeps = "<<empty_slice_gsweeps_R<<endl;
    if(bar != NULL)bar->SetPosition(100);

/*
   if(index_t0==0){ c_p1_p2->Divide(2,2);}
   //if(index_t0 !=0) {delete fun_g1;delete fun_g2;}
   c_p1_p2->cd(1);
   tree[index_t0]->Draw("timec:nevt>>h2p1","","colz");
   c_p1_p2->cd(2);
   tree[index_t0]->Draw("timec:nevt>>h2p2","","colz");
   c_p1_p2->cd(3);
   tree[index_t0]->Draw("timec>>h1p1","","colz");
   fun_g1 = new TF1("fun_g1","gaus",0,1e9);
   h1p1->Fit(fun_g1,"","",9006.4e3,9006.6e3);
   c_p1_p2->cd(4);
   tree[index_t0]->Draw("timec>>h1p2","","colz");
   fun_g2 = new TF1("fun_g2","gaus",0,1e9);
   h1p2->Fit(fun_g2,"","",9013.4e3,9013.8e3);

   g_mean_p1->SetPoint(index_t0,index_t0*t0_delta,fun_g1->GetParameter(1));
   g_sigma_p1->SetPoint(index_t0,index_t0*t0_delta,fun_g1->GetParameter(2));   
   g_mean_p2->SetPoint(index_t0,index_t0*t0_delta,fun_g2->GetParameter(1));
   g_sigma_p2->SetPoint(index_t0,index_t0*t0_delta,fun_g2->GetParameter(2));

   c_p1_p2->Update();
   if(index_t0==0)  c_p1_p2->WaitPrimitive();
   delete fun_g1;
   delete fun_g2;
*/

} // end of for loop of t0
/*
   c1->Divide(2,2);
   c1->cd(1);
   g_mean_p1->SetTitle("p1 mean;t0 [ns];mean");
   g_mean_p1->Draw("APL*");
   c1->cd(3);
   g_sigma_p1->SetTitle("p1 sigma;t0 [ns];sigma");
   g_sigma_p1->Draw("APL"); 
   c1->cd(2);  
   g_mean_p2->SetTitle("p2 mean;t0 [ns];mean");
   g_mean_p2->Draw("APL*");
   c1->cd(4);
   g_sigma_p2->SetTitle("p2 sigma;t0 [ns];sigma");
   g_sigma_p2->Draw("APL");
*/

  return 1;
}


// dynamic hiso range and dynamic sigam 

double fit_result(TTree *t1, Long64_t min_evt, Long64_t max_evt, int nbins, double lowx, double highx,Double_t Sigma){
    //gROOT->SetBatch(true);
       double n_slide_central = 0;
       static double active_lowx = lowx;
       static double active_highx = highx;
       static bool appro_lasttime = false; // use approximate last time
  
       if(isnewfile){
        active_lowx = lowx;
        active_highx = highx;
        appro_lasttime = false;
         isnewfile=false;
       }
        double Half_hiswidth = (highx-lowx)*0.5;  

    // for search peak and temp result of fitting;//&&&&&&&&&&&  not use now; cover by sigma_StdDev

       // active_lowx = h->GetXaxis()->GetBinCenter(h->GetMaximumBin())-Half_hiswidth;
       // active_highx = n_slide_central+Half_hiswidth;
             
       h = new TH1D("h","h",nbins,active_lowx,active_highx);
        
 
         
       vector <Double_t> *tmp_time = new vector <Double_t>();
       vector <Long_t> * tmp_sweeps_global = new vector <Long_t>;
       t1->SetBranchAddress("time",&tmp_time);
       t1->SetBranchAddress("sweeps_global",&tmp_sweeps_global);
       tmp_time->clear();
       tmp_sweeps_global->clear();

       vector<Long_t> temprecord_sweeps;  // record global sweeps of all hits 


        for(Long64_t index = min_evt;index< max_evt;index++){

               t1->GetEvent(index);
               int nhits = (int) (tmp_time->size());

              for(int n_hit=0;n_hit < nhits;n_hit++){
                  h->Fill( (Double_t) tmp_time->at(n_hit) );
                  temprecord_sweeps.push_back( tmp_sweeps_global->at(n_hit) );
              }
                  
             tmp_time->clear();
             tmp_sweeps_global->clear();


        }     

        empty_slice_gsweeps_L = *min_element(temprecord_sweeps.begin(),temprecord_sweeps.end());
        empty_slice_gsweeps_R = *max_element(temprecord_sweeps.begin(),temprecord_sweeps.end());

        temprecord_sweeps.clear();

         // h->Draw();
//c1->Modified();



        double* getfeature = FeaturePar(h,min_evt,max_evt);  // [0] Amp; [1] Peak center; [2] sigma; [3] bin width

        if(getfeature[6]<(max_evt-min_evt)*0.5*0.1){// empty slice; x0.5 ==> only half of them from tag1
            delete h;
            delete tmp_time;
            delete tmp_sweeps_global;
            return -20;  // value represents empty slice
        }

        //%%%%%%%%%%%%%%%%%%%%   renew histogram with new feature paras  %%%%%%%%%%%%%%%%%%
        /*
                    int bin_N = TMath::Nint((active_highx-active_lowx)/getfeature[3]);
                    if(bin_N<10){
                      cout<<"binN="<<bin_N<<endl;
                      cout<<"exp binwidth = "<<getfeature[3]<<endl;
                      cout<<"exp std ="<<getfeature[2]<<endl;
                      c1->cd();
                      h->Draw();
                    c1->Update();
                    c1->WaitPrimitive();
                    return -10;
                    }
          */

        int Nbins_cal = TMath::Nint((active_highx - active_lowx) / getfeature[3]);

        if (Nbins_cal > 10) {  // use new Nbins define histogram

              delete h;
              h = new TH1D("h","h",Nbins_cal,active_lowx,active_highx);

             for(Long64_t index = min_evt;index< max_evt;index++){

                     t1->GetEvent(index);
                     int nhits = (int) (tmp_time->size());

                 for(int n_hit=0;n_hit < nhits;n_hit++){
                        h->Fill( (Double_t) tmp_time->at(n_hit) );

                    }
                        
                   tmp_time->clear();

              }  

              getfeature = FeaturePar(h,min_evt,max_evt);
        }
        else{
                      cout<<"binN="<<Nbins_cal<<endl;
                      cout<<"exp binwidth = "<<getfeature[3]<<endl;
                      cout<<"exp std ="<<getfeature[2]<<endl;
                      printf("exp fwhm_l = %.4f exp fwhm_r = %.4f\n",getfeature[4],getfeature[5]);
                 /*     c1->cd();
                      h->Draw();
                      c1->Update();
                      c1->WaitPrimitive(); */
        }


        double lambda_old = func->GetParameter(2);
        double gamma_old = func->GetParameter(3);
        double delta_old = func->GetParameter(4);

        func->SetParameters(getfeature[0],getfeature[1],getfeature[2]);



        double fit_rangeL = (getfeature[1] - (getfeature[2]*2.36)*0.5*2.5 );   //0.5* FWHM * 1.5
        fit_rangeL = (fit_rangeL>active_lowx)? fit_rangeL : active_lowx;

        double fit_rangeR = (getfeature[1] + (getfeature[2]*2.36)*0.5*3 );
        fit_rangeR = (fit_rangeR < active_highx)? fit_rangeR : active_highx;




        for(int i=0;i<15;i++){
            h->Fit(func,"LMQN","",fit_rangeL,fit_rangeR);
        }
/*
c1->cd();
h->Draw();
h->Fit(func,"LM","",fit_rangeL,fit_rangeR);
c1->Modified();
c1->Update();
char tempc;
while(1){
  cout<<"go?"<<endl;
  cin>>tempc;
  if(tempc!='n' && tempc!='y') continue;
  else{
    if(tempc=='n') return -10;
    if(tempc=='y')break;
  }
}
*/
        n_slide_central = func->GetMaximumX(fit_rangeL,fit_rangeR,1e-13,500);

        bool condition1 = (lambda_old==0)? true : (TMath::Abs((func->GetParameter(2)-lambda_old)/lambda_old) < 4);  // first fit: lambda_old==0
        bool condition2 = (gamma_old==0)? true : (TMath::Abs((func->GetParameter(3)-gamma_old)/gamma_old) < 9);
        bool condition3 = (delta_old==-1)? true : (TMath::Abs((func->GetParameter(4)-delta_old)/delta_old) < 4);
        bool condition4 = TMath::Abs(n_slide_central - getfeature[1]) <getfeature[2];
        bool condition5 = TMath::Abs((func->GetMaximum(fit_rangeL,fit_rangeR) -getfeature[0])/TMath::Sqrt(getfeature[0])) < 12;

        if(condition1 && condition1 && condition3 && condition4 && condition5){ // fit successful
                active_lowx = n_slide_central-Half_hiswidth;
                active_highx = n_slide_central+Half_hiswidth;
                appro_lasttime = false;
                if(stop_point>0)N_byfit+=2;
                else N_byfit++;

              delete h;
              delete tmp_time;
              return n_slide_central;
        }
        else{ // fit fail

            if(stop_point>0 || t1->GetEntriesFast()==max_evt){ // already fail last slice, ==> approximate by mean or last slice

                  h->GetXaxis()->SetRangeUser(fit_rangeL+(getfeature[2]*2.36)*0.5*1.5,fit_rangeR-(getfeature[2]*2.36)*0.5*1.5);
                  n_slide_central = h->GetMean();

                  active_lowx = n_slide_central-Half_hiswidth;
                  active_highx = n_slide_central+Half_hiswidth;
                  if(appro_lasttime){func->SetParameters(0,0,0,0,-1);}  // reset to initial pars if fitting of continuous slices are both failed
                  else{ func->SetParameters(0,0,lambda_old,gamma_old,delta_old);}
                  appro_lasttime = true;
                  N_byappro++;

                  delete h;
                  delete tmp_time;
                  return n_slide_central;

            }
            else{ // fail current lice ==> accumulate to next slice

                  if(appro_lasttime){func->SetParameters(0,0,0,0,-1);}  // reset to initial pars
                  else{ func->SetParameters(0,0,lambda_old,gamma_old,delta_old);}

                  delete h;
                  delete tmp_time;
                  return -10;
            }

        }
}




double* FeaturePar(TH1D* _inhisto,Long64_t _min_evt, Long64_t _max_evt){
    static double par_return[7];  // Peak height , peak center, sigma, good bin size, FWHM_L, FWHM_R, integral counts
    TH1D* hc = (TH1D*) _inhisto->Clone();
    par_return[6] = hc->Integral(1,hc->GetNbinsX());
    if(par_return[6]<(_max_evt-_min_evt)*0.5*0.1){ // empty slice
       delete hc;
       return par_return;
    }

    hc->Smooth(4);

    int Bin_max = hc->GetMaximumBin();
    par_return[0] = hc->GetBinContent(Bin_max);  // Amp
    par_return[1] = hc->GetBinCenter(Bin_max);   // peak center
    double height_half = par_return[0]*0.5;

    double bin_i_y1;
    double bin_i_y2;
    double candidate_x;
    double FWHM_L=0,FWHM_R=0;

    for(int i=1;i<hc->GetNbinsX();i++){
      bin_i_y1 = hc->GetBinContent(i);
      bin_i_y2 = hc->GetBinContent(i+1);
      candidate_x = (height_half - bin_i_y1) / (bin_i_y2 - bin_i_y1) * hc->GetBinWidth(1) + hc->GetBinCenter(i);
      if(i<Bin_max){
          if(bin_i_y1<height_half && bin_i_y2 >height_half){  FWHM_L = candidate_x; }
      }
      else{
          if(bin_i_y1>height_half && bin_i_y2 <height_half){    FWHM_R = candidate_x; break;  }
      }
    }

    double FWHM = FWHM_R - FWHM_L;

    hc->GetXaxis()->SetRangeUser(FWHM_L-1.5*FWHM,FWHM_R+1.5*FWHM);

    par_return[2] = hc->GetStdDev();

    _inhisto->GetXaxis()->SetRangeUser(FWHM_L-1.5*FWHM,FWHM_R+1.5*FWHM);//FWHM_L-0.25*FWHM,FWHM_R+0.25*FWHM

    double maincounts = _inhisto->Integral();

    //scott rule for bin size:  bin size= 3.49*sigma*(counts)^(-1/3)

    par_return[3] = 3.49 * par_return[2] *TMath::Power(maincounts,-1./3.);

    _inhisto->GetXaxis()->UnZoom();

par_return[4]=FWHM_L;
par_return[5]=FWHM_R;

    delete hc;

    return par_return;

}




#ifndef __CINT__
int main(int argc, char *argv[]){
  if(argc>1){
    DriftCorrect( string(argv[1]) , string(argv[2]),atof(argv[3]),atoi(argv[4]),atoi(argv[5]),atof(argv[6]),atof(argv[7]),atof(argv[8]));
  }
  return 0;
}
#endif







