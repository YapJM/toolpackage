#ifdef _PREVIEWER_
#ifndef _EXEC_H_
#define _EXEC_H_


#include "TSystem.h"
#include "Buttons.h"
#include "TObject.h"
#include "TVirtualX.h"
#include <stack>

TCanvas* c_handle;
bool stack2clear[2] = {1,1};
// fit setting
string AFitOption= "LMEQ";
int AFitSetFunc=2;
char AFit_x_r = 'x';
bool FixFitRange=false;
double fitL_width=0;  // fix fit width for normal fit
double fitR_width=0;

double unbinned_fitL=0;  // fit left edge tof value
double unbinned_fitR=0; // fit right edge tof value


void SetRange(int index,TCanvas* c_get,TH1D* h,double xmin, double xmax);
void SetLogy(TPad* p,bool flag);
void MakerAdjY(TH1D* h_in, TOFMarker* marker,double h_left,double h_right);

//********** for sampling fitting ********************
bool Sampling(int padindex, double _rangeL, double _rangeR);

void Exec(int index=1){
    if (!gPad) {
      Error("exec", "gPad is null, you are not supposed to run this macro");
      return;
    }  

    static double xmouse[2] = {0};
    static int imouse = 0;
    static TH1D* h_history[2] = {};
    static stack<double> lowlimit[2];
    static stack<double> highlimit[2];

    // new histogram setting, clear stack;
    if(index==1){
        if(h_history[index-1] != h_xF){
          h_history[index-1] = h_xF;
          stack2clear[index-1] = true;
        }
        else stack2clear[index-1] = false;

    }
    else if(index ==2){
        if(h_history[index-1] != h_zoom_x){
          h_history[index-1] = h_zoom_x;
          stack2clear[index-1] = true;
        }
        else stack2clear[index-1] = false;
    }
    else{
        Error("exec", "index =1 or 2; exec program only limit to h_xF and h_zoom_x!!");
        return;
    }


    if(stack2clear[index-1]){
      for(unsigned int i=0;i<lowlimit[index-1].size();i++){
        lowlimit[index-1].pop();
        highlimit[index-1].pop();
      }
      stack2clear[index-1]=false;
    }



    c_handle = (TCanvas*)gPad->GetCanvas();
    c_handle->FeedbackMode(kTRUE);


    //TObject *Sele = gPad->GetSelected();
   //if(!select) return;
   //if (!select->InheritsFrom(TH1::Class())) {gPad->SetUniqueID(0); return;}
/*
   TH1D * h_handle = (TH1D*)Sele;
   int index=0;
   const char* getname = h_handle->GetName();

   if(strcmp(getname,"h_xF") ==0){index=0;}
   else if(strcmp(getname,"h_zoom_x") ==0){index=1;}
   else if(strcmp(getname,"h_refF") ==0){index=2;}
   else return;*/

    
  TPad* Pd = (TPad*)c_handle->cd(index);

  //cout<<"\r"<<"index = "<<index<<flush;




  //TPad* Pd = (TPad*) gPad->GetPad(index);

  int pxold = Pd->GetUniqueID();
  int px =  Pd->GetEventX();
  int py =  Pd->GetEventY();

  float uxmin = Pd->GetUxmin();
  float uxmax = Pd->GetUxmax();
  int pxmin = Pd->XtoAbsPixel(uxmin);
  int pxmax = Pd->XtoAbsPixel(uxmax);
  float uymin = Pd->GetUymin();
  float uymax = Pd->GetUymax();
  int pymin = Pd->YtoAbsPixel(uymin);
  int pymax = Pd->YtoAbsPixel(uymax);

  if(px<pxmin) px = pxmin+1;
  if(px>pxmax) px = pxmax-1;
  
  if(pxold) gVirtualX->DrawLine(pxold,pymin,pxold,pymax);
  gVirtualX->DrawLine(px,pymin,px,pymax);
  Pd->SetUniqueID(px);
  Double_t upx =  Pd->AbsPixeltoX(px);
  Double_t x =  Pd->PadtoX(upx);

  EEventType event = static_cast<EEventType> ( Pd->GetEvent());
  if(event == kButton1Down){ // select range
    //cout<<"\r x: "<<Form("%.2f",upx)<<flush;
    xmouse[imouse] = upx;  imouse = (imouse+1)%2;
   /* TObject *select = gPad->GetSelected();
   TH1D * h_handle = (TH1D*)select;
   int index=0;
   const char* getname = h_handle->GetName();

   if(strcmp(getname,"h_xF") ==0){index=0;}
   else if(strcmp(getname,"h_zoom_x") ==0){index=1;}
   else if(strcmp(getname,"h_refF") ==0){index=2;}
   else return;*/
    
    gVirtualX->DrawLine(px,uymin,px,uymax);

   // printf("px=%d , py=%d, upx=%f, x=%f\n",px,py,upx,x);
  }
  else if(event == kButton1Double){
    printf("tof position = %.3f\n\n",x);
    gXposition = x;   //SetROI use

    if(m_ref==0){cout<<"No ref mass setting, no mass calculation"<<endl;}
    else if(tof_ref_cento <1){cout<<"No ref tof setting, no mass calculation"<<endl;}
    else{
      double tem_tof = tof_x_cento[0];
      double tem_tof_err = tof_x_cento_err[0];
      tof_x_cento[0] = x;
      tof_x_cento_err[0]=0.7;
      mass_calculator(1);
      tof_x_cento[0] =tem_tof;
      tof_x_cento_err[0] = tem_tof_err;
      
    }
  }
  else if(event == kKeyPress){ // Key action
      int press = Pd->GetEventX(); //cout<<"press: "<<press<<endl;
      if(press==' '){ // set new range
        double Tmin = TMath::Min(xmouse[0],xmouse[1]);
        double Tmax = TMath::Max(xmouse[0],xmouse[1]);
        histo_zoom_in_x(0,100,Tmin,Tmax);
        printf("used: histo_zoom_in_x(%d,%d,%.1f,%.1f)\n",0,100,Tmin,Tmax);
        c_handle->cd(2)->Modified();
        c_handle->cd(2)->Update();      
        c_handle->Modified();  c_handle->Update();
        
      }else if(press=='z'){ // zoom
        lowlimit[index-1].push(TMath::Min(xmouse[0],xmouse[1]));
        highlimit[index-1].push(TMath::Max(xmouse[0],xmouse[1]));
        SetRange(index,c_handle,h_history[index-1],lowlimit[index-1].top(),highlimit[index-1].top()); 
        if(index==1){ //"h_xF"
            MakerAdjY(h_history[index-1],marker_tof,lowlimit[index-1].top(),highlimit[index-1].top());
            ROIadjY(h_history[index-1],lowlimit[index-1].top(),highlimit[index-1].top());
            c_handle->cd(1)->Modified();
            c_handle->cd(1)->Update();
            c_handle->Modified();  c_handle->Update();
        }     
      }else if(press=='x'){ // unzoom
        if(lowlimit[index-1].size() !=0){
              lowlimit[index-1].pop();
              highlimit[index-1].pop();
        }
        if(lowlimit[index-1].size() ==0){ 
              h_history[index-1]->GetXaxis()->UnZoom(); Pd->Modified();Pd->Update();
              double histo_L = EJE0-300;
              double histo_H = histo_L + sptrFW;
              MakerAdjY(h_history[index-1],marker_tof,histo_L,histo_H);
              ROIadjY(h_history[index-1],histo_L,histo_H);
              c_handle->cd(1)->Modified();
              c_handle->cd(1)->Modified();
              c_handle->cd(1)->Update();
              c_handle->Modified();  c_handle->Update();
        }
        else{
             SetRange(index,c_handle,h_history[index-1],lowlimit[index-1].top(),highlimit[index-1].top()); 
              if(index==1){ //"h_xF"
                  MakerAdjY(h_history[index-1],marker_tof,lowlimit[index-1].top(),highlimit[index-1].top());
                  ROIadjY(h_history[index-1],lowlimit[index-1].top(),highlimit[index-1].top());
                  c_handle->cd(1)->Modified();
                  c_handle->cd(1)->Modified();
                  c_handle->cd(1)->Update();
                  c_handle->Modified();  c_handle->Update();
              } 
        }
      }else if(press=='l'){ // logy
        SetLogy(Pd,!Pd->GetLogy());
        c_handle->Modified(); c_handle->Update();
      }else if(press=='r'){// add new ROI
        if(ROI_initial){
                  ROI_INDEX++;
                  SetROI(gXposition,ROI_WIDTH,(ROI_INDEX-1)%40+1,0,true,0);
                  c_handle->Modified(); c_handle->Update();
                  //printf("%.4f,%f,%d,%d",gXposition,ROI_WIDTH,ROI_INDEX,(ROI_INDEX-1)%20+1); // test purpose
        }
        else{
                  SetROI(gXposition,ROI_WIDTH,ROI_INDEX,0,true,0);
                  c_handle->Modified(); c_handle->Update();
                  //printf("%.4f,%f,%d,%d",gXposition,ROI_WIDTH,ROI_INDEX,(ROI_INDEX-1)%20+1); // test purpose
                  ROI_initial=true;
        }
      }else if(press=='R'){// correct current ROI
                  SetROI(gXposition,ROI_WIDTH,(ROI_INDEX-1)%40+1,0,true,0);
                  c_handle->Modified(); c_handle->Update();
      }else if(press=='s'){
                cout<<"Sampling at pad(4)............"<<endl;
                 if(Sampling(4,-1,-1))cout<<"Successful sampling"<<endl; // sample at tag1 , have to draw arrow to define sampling range

      }else if(press=='S'){
                cout<<"Sampling at pad(2)............"<<endl;
                if(Sampling(2,TMath::Min(xmouse[0],xmouse[1]),TMath::Max(xmouse[0],xmouse[1])))cout<<"Successful sampling"<<endl; // sample at tag1 , have to draw arrow to define sampling range

      }else if(press=='f'){// fitquickly
              c_handle->cd(2);
              if(index==2){
                  double Max_high =0;  // maximum count in selected range
                  double Max_tof=0;    // tof of maximum bin
                  double Tmin = TMath::Min(xmouse[0],xmouse[1]);
                  double Tmax = TMath::Max(xmouse[0],xmouse[1]);
                  int Bin_L = h_zoom_x->FindBin(Tmin);
                  int Bin_R = h_zoom_x->FindBin(Tmax);

                  for(int i=Bin_L;i<=Bin_R;i++){ // find the position and count of maximum bin in the range
                    if( h_zoom_x->GetBinContent(i) > Max_high ){Max_high = h_zoom_x->GetBinContent(i); Max_tof=h_zoom_x->GetBinCenter(i);}
                  }

                  tem_func->SetParameter(0,Max_high);
                  tem_func->SetParameter(1,Max_tof);

                  if(FixFitRange && fitL_width!=0 && fitR_width!=0){ // fix fit range
                      Tmin = Max_tof-fitL_width;
                      Tmax = Max_tof+fitR_width;
                  }
                  else{ // free fix fit range
                      fitL_width = Max_tof - Tmin;  
                      fitR_width = Tmax - Max_tof;
                  }


                  string tem_func_name = tem_func->GetName();
                  cout<<"tem_func_anme = "<<tem_func_name<<endl;

                  //*************** fit with sampling function ******************
                  if(tem_func_name == "fsample"){
                      if(fs->GetOldNPeaks() != fs->NumOfPeaks){fs->Makefitfunc();}

                       if(NumOfPeaks==1){ // one peak
                          fs->SetPars(1,Max_high,Max_tof);
                       }
                       else{// more than one peak

                          for(int i=0;i<NumOfPeaks;i++){
                            cout<<"draw an arrow from top of "<<"\e[1;33m"<<"Peak_"<<i+1<<"\e[0m"<<" to FWHM"<<endl;
                            cout<<"pending......"<<endl;
                            get_para_by_draw(2);
                            fs->SetPars(i+1,tem_high,tem_cento);
                            
                          }

                       }

                      fs->Fit(h_zoom_x,Tmin,Tmax,AFitOption);
                     // tem_func = fs->Getfitfunc();
                      c_handle->Modified(); c_handle->Update();
                      printChi(h_zoom_x,'x',Tmin,Tmax);

                      for(int i=0;i<NumOfPeaks;i++){
                           tof_x_cento[i] = fs->GetTofCenter(i+1); // peakindex from 1
                           tof_x_cento_err[i] = fs->GetTofCenterErr(i+1);
                           printf("cento_%d: %.4f(%.4f)\n",i+1,tof_x_cento[i],tof_x_cento_err[i]);
                           unbinned_fitR = tof_x_cento[i] + fs->range_R; // update fit right edge for unbinned fit
                      }

                          unbinned_fitL = tof_x_cento[0] - fs->range_L;
                          unbinned_fitR = tof_x_cento[NumOfPeaks-1] + fs->range_R; // update fit right edge for unbinned fit

                      if(fs->NumOfPeaks>1){  fs->Draw_subline(c1,2);  }

                  }// end of sampling fit
                  else{//******************** fit with gaus_exp function *********************
                      fitquickly(AFit_x_r,AFitSetFunc,AFitOption.c_str(),true,Tmin,Tmax,50);
                      unbinned_fitL = Tmin;
                      unbinned_fitR = Tmax;
                      c_handle->Modified(); c_handle->Update();
                      printf("using: fitquickly('%c',%d,\"%s\",true,%.4f,%.4f,50)\n",AFit_x_r,AFitSetFunc,AFitOption.c_str(),Tmin,Tmax);  
                  }// end of gaus_exp fit                
              }
              else cout<<"Only canvas 2 has fitting function"<<endl;

      }
     else if(press=='-'){
             double Tmin = TMath::Min(xmouse[0],xmouse[1]);
             double Tmax = TMath::Max(xmouse[0],xmouse[1]);
          	 printf("x1=%.4f;\t x2=%.4f\n",Tmin,Tmax);
          	 printf("distance:%.4f\n",Tmax-Tmin);
      }
      else if(press =='c'){
          printf("Searching for beta - TOF coincidence .......\n\n");
          ShowCoinCondition();
          ShowRIHalflife();
          if(Current_beta_file != ToLoad_beta_file ){  // beta file read to vector
            if(ToLoad_beta_file=="+++"){
              cout<<"warning: no beta file to be load; please load file first"<<endl;
              cout<<"\e[1;33m"<<"beta filename with .lst"<<"\e[0m"<<endl;
              cout<<"using:  LoadNewBetaFile("<<endl;
              return;
            }
            else{
              Read_Beta_lst(FilePath+"../","LST/"+ToLoad_beta_file);
            }
          }
          if(Compare_Beta_Beta()){
                ShowBeta_Beta_Coin();
          }  // which will create tbeta_beta coincident tree, can be an option

          if(Match_hit.size()==0){cout<<"\e[1;23m"<<"No coincidence in store"<<"\e[0m"<<endl; return;}

          double Tmin = TMath::Min(xmouse[0],xmouse[1]);
          double Tmax = TMath::Max(xmouse[0],xmouse[1]);
          FindBetaTof_Coin(Tmin,Tmax, RIHalflive);
      }
      else if(press =='m'){ // mass calculate
        for(int index_m=0;index_m<NumOfPeaks && index_m<10;index_m++){
          cout<<"Unknow mass "<<index_m+1<<":"<<endl;
          mass_calculator(tof_x_cento[index_m],tof_x_cento_err[index_m],laps_x,q_x,tof_ref_cento,tof_ref_cento_err,laps_ref,m_ref,err_ref,bref,q_ref);
          cout<<endl;
        }

      }
      else if(press==','){
          fs->FreeRange = !(fs->FreeRange);
          if(fs->FreeRange){ cout<<"Set to FreeRange== true for sampling fitting, fitting range FREE now!!!"<<endl; FixFitRange =false;}
          else{ cout<<"Set to FreeRange== false"<<endl;
                  cout<<"fitting range fix with  same left and right tail ratio as sampling peak"<<endl;
                  FixFitRange=true;
          }
      }
      else if(press=='u'){
        UnbinnedFit(tem_func,active_tree_name, h_zoom_x, unbinned_fitL,unbinned_fitR);
      }

  }
 
  

  funcS::NumOfPeaks = NumOfPeaks; // keep update to global NumOfPeaks

}//end of Exec


void SetRange(int index,TCanvas* c_get,TH1D* h,double xmin, double xmax){
  h->GetXaxis()->SetRangeUser(xmin,xmax);
  c_get->cd(index)->Modified();
  c_get->cd(index)->Update();      
  c_get->Modified();  c_get->Update();
  return;
}

void SetLogy(TPad* p,bool flag){
  p->SetLogy(flag);
  p->Modified();
  p->Update();      
  return;
}

void MakerAdjY(TH1D* h_in, TOFMarker* marker, double h_left,double h_right){
  string histo_name = h_in->GetName();
  if(histo_name == "h_xF"){
      double Y2Set = h_in->GetBinContent( h_in->GetMaximumBin() );
      if(Y2Set ==0) Y2Set=1.01;
      else Y2Set *=1.01;
      double lineX=0;
      for(int index=0;index<40;index++){
          lineX = marker[index].GetX();
          if(lineX>= h_left && lineX<=h_right){
            marker[index].SetY(Y2Set);
          }
      }
  }
}

bool Sampling(int padindex, double _rangeL, double _rangeR){
    if(fs==NULL){cout<<"sample function is no exist, Abort!!!"<<endl; return false;}
  /*  int padindex=2;
    while(1){
      cout<<"Which Pad for sampling: 2 or 4"<<endl;
      cin>>padindex;
      if(padindex==2 || padindex==4) break;
    }*/

    TH1D* getHisto=NULL;
    if(padindex==2){getHisto=h_zoom_x;}
    else if(padindex==4){
          getHisto=h_zoom_ref;
          cout<<"draw a line for sampling range"<<endl;
          cout<<"pending..."<<endl;
          get_para_by_draw(padindex);
          _rangeL = fitrangeL;
          _rangeR = fitrangeR;
    }
    else{cout<<"Error: padindex!!!"<<endl; return false;}

    fs->Sampling(getHisto,_rangeL,_rangeR);
    fs->RecreateHisto();

    if(fs->SmoothHisto()){
        fs->Makefitfunc();
        tem_func=fs->Getfitfunc(); // get handle of fs
        tof_ref_cento = fs->GetsPeakCenter();
        tof_ref_cento_err = fs->GetsPeakCenter_err();
    }
    else{cout<<"Error in Smooth histogram!!! Abort!!!"<<endl; return false;}

    if(fs->Draw(c1,padindex)){c1->cd(2);}
    else{cout<<"Error, faile to draw sampling line on histogram!!!"<<endl; return false;}

    return true;

}


#endif  // ifndef _Exec_h_

#endif// ifdef _PREVIEWER_
