#include<iostream>
#include<stdlib.h>
#include"TSystem.h"
#include"TTree.h"
#include"TFile.h"
#include<fstream>
#include <vector>
#include <stdio.h>
#include"TMath.h"
#include"TTree.h"
#include"TCanvas.h"
#include"TGraph.h"
#include"TH1.h"
#include"TH2.h"

#include"RefineBeta.h"

#if defined (MAKECINT)
#pragma link C++ class vector<Long_t>+;   // in order to identify vector<Long_t> 
#endif

using namespace std;

// Graph of beta-beta energy distribution and coincidence time difference distribution
TFile* f_beta_tree=NULL;
TTree* tbeta_beta=NULL;
TCanvas* c_beta_beta=NULL;
TH1D* histo_beta_deltaT=NULL; // time difference of coincidence
TH2D* histo_beta_EkeV=NULL;   // E distribution of beta coincidence
TH2D* histo_beta_EkeV_withcut=NULL; // E distribution with cut
TH2D* histo_beta_EkeV_CutAndRej=NULL; // E distribution with cut and reject
TH1D* histo_h1E=NULL;  // Silicon 1 ADC histo
TH1D* histo_h2E=NULL;  // Silicon2 ADC histo



// beta energy calibration
//**********************************
double beta_slope[2]={0};
double beta_intercept[2]={0};

void GetBeta_E_Calibrate_para(){ // [0]=>Silicon1; [1]Silicon2
	printf("Silicon1: beta_slope[0]= %.2f ; beta_intercept[0]=%.2f\n",beta_slope[0],beta_intercept[0]);
	printf("Silicon2: beta_slope[1]= %.2f ; beta_intercept[1]=%.2f\n",beta_slope[1],beta_intercept[1]);
}

void SetBeta_E_Calibrate_para(double _slope1=1,double _intercept1=0,double _slope2=1,double _intercept2=0){ // [0]=>Silicon1; [1]Silicon2
	if(beta_slope[0]!=_slope1) beta_slope[0] =_slope1;
	if(beta_slope[1]!=_slope2) beta_slope[1] =_slope2;
	if(beta_intercept[0]!=_intercept1) beta_intercept[0] =_intercept1;
	if(beta_intercept[1]!=_intercept2) beta_intercept[1] =_intercept2;
	GetBeta_E_Calibrate_para();
}

double ADC2keV_1(int _adc){return beta_slope[0]*_adc + beta_intercept[0];}
double ADC2keV_2(int _adc){return beta_slope[1]*_adc + beta_intercept[1];}
//**************************************


//***************** Beta-beta and betaTof coincident condition *************
Long64_t gate_time=500;   // time window 500 ns
int gate_adc_low=0;         // ADC low limit
int gate_adc_hi=1000;      // ADC High limit

void ShowCoinCondition(){
	printf("coincident condition:  gate time\t gate_adc_low \t gate_adc_hi\n");
	cout<<"\e[1;37m"<<gate_time<<"\t"<<gate_adc_low<<"\t"<<gate_adc_hi<<endl;
	printf("To change condition:     ");
	cout<<"\e[1;37m"<<"SetCoinCondition("<<"\e[0m"<<endl;
	cout<<endl;
}
void SetCoinCondition(Long64_t _gate_time=500, int _gate_adc_low=0, int _gate_adc_hi=1000){
	gate_time= _gate_time;   // time window 500 ns
 	gate_adc_low= _gate_adc_low;         // ADC low limit
 	gate_adc_hi= _gate_adc_hi; 
 	ShowCoinCondition();
}


Long64_t ModuleDelay = 50; // 50[ns] signl to ch4 has a 50 ns delay from synchronizer output !!!!!!!!!!!!!!!!!!!!!
Long64_t RIHalflive = 150 * 1000000;  // radioactive ion half live in ns
void ShowRIHalflife(){
	cout<<"\e[1;37m"<<"RI halflife:  "<<RIHalflive<<" [ns] !!!!"<<"\e[0m"<<endl;
	printf("Set halflife by:  SetRIHalflife()\n");
	cout<<endl;
}
void SetRIHalflife(double InMilliSecond = 150){
		RIHalflive = (Long64_t) InMilliSecond * 1000000;  // to [ns]
		ShowRIHalflife();
}
double Lifetime2Halflife(double Inlivetime_ms=150){
		double  output= Inlivetime_ms*TMath::Log(2);
		cout<<"Halflife T1/2 = "<< output<<" [ms] !!!"<<endl;
		return output;    // in [ms]
}
double Halflife2Lifetime(double Inhalflife_ms=150){
		double output= Inhalflife_ms/TMath::Log(2);
		cout<<"life time tao = "<<output<<" [ms] !!!"<<endl;
		return output;	// in [ms]
}
void SetRIHalflife_byLifeTime(double Inlivetime_ms=150){
	SetRIHalflife(Lifetime2Halflife(Inlivetime_ms));
}
//**************************************************



int syn_CH = 4; // beta-tof synchronizer channel

//vector <Long64_t> *nevt_v = new vector <Long64_t>();
vector<Long64_t> *sweeps_gclock_v = new vector<Long64_t>();
vector <Long64_t> *time_v = new vector <Long64_t>();
vector <Long64_t> *time_relative_v = new vector <Long64_t>();
vector<int> *adc_v = new vector<int>();
vector<int> *channel_v = new vector<int>();

Long64_t nevt_bta=0;  // how many lines or events in a file

string Path_beta_file="./";
string Current_beta_file="---";
string ToLoad_beta_file="+++";

void LoadNewBetaFile(string filename){
	ToLoad_beta_file = filename;
	cout<<"Beta file to be load: "<<filename.c_str()<<endl;
}

void Read_Beta_lst(string PATH,string filename){

	FILE* fin =NULL;
	//string fpath_name = PATH + "LST/" + filename +".txt";
	string fpath_name = PATH + filename;
	fin = fopen(fpath_name.data(),"r");
	if(fin==NULL){
		cout<<"fail to open list file of beta!! break!!"<<endl;
		return;
	}

	nevt_bta=0;
	Long64_t sweeps_gclock=0;
	Long64_t time =0;
	Long64_t current_sweeps_time_gclock = 0;
	int adc=0;
	int channel=0;

/*	//vector <Long64_t> *nevt_v = new vector <Long64_t>();
	vector <Long64_t> *time_v = new vector <Long64_t>();
	vector<int> *adc_v = new vector<int>();
	vector<int> *channel_v = new vector<int>();
*/

	sweeps_gclock_v->clear(); sweeps_gclock_v->shrink_to_fit();
	time_v->clear(); time_v->shrink_to_fit();
	time_relative_v->clear(); time_relative_v->shrink_to_fit();
	adc_v->clear();  adc_v->shrink_to_fit();
	channel_v->clear();  channel_v->shrink_to_fit();
	bool firstloop=false; // turn to true after finish first read

	while(!feof(fin)){
		if(firstloop){
			nevt_bta++;
			time_v->push_back(time);
			channel_v->push_back(channel);
			if(channel==syn_CH){
				 sweeps_gclock++;  // in mrtof.lst file: sweeps start from 1; NOT 0
				 current_sweeps_time_gclock = time;
			}
			sweeps_gclock_v->push_back(sweeps_gclock);
			time_relative_v->push_back(time - current_sweeps_time_gclock);
			adc_v->push_back(adc);
			if(nevt_bta%1000==0)cout<<'\r'<<"nevt_bta= "<<nevt_bta<<flush;
			//if(nevt_bta==10)break;
		}

		fscanf(fin,"%lld,%d,%d",&time,&channel,&adc);
		firstloop=true;

	}
	cout<<endl;
	fclose(fin);

/*
	for(Long64_t index=0;index<100;index++){
		sweeps_gclock = sweeps_gclock_v->at(index);
		time=time_v->at(index);
		adc=adc_v->at(index);
		channel = channel_v->at(index);
		printf("%lld,%lld,%d,%d\n",sweeps_gclock,time,channel,adc);
	}
*/

	/*if(histo_h1E!=NULL) delete histo_h1E;
		histo_h1E = new TH1D("histo_h1E"," Si1 and Si2 ADC ",512,0,1024); histo_h1E->SetLineColor(kBlack);
	if(histo_h2E!=NULL) delete histo_h2E;
		histo_h2E = new TH1D("histo_h2E"," Si1 and Si2 ADC ",512,0,1024); histo_h2E->SetLineColor(kRed);
	for(Long64_t index=0;index<nevt_bta;index++){
		if(channel_v->at(index)==1) histo_h1E->Fill(adc_v->at(index));
		if(channel_v->at(index)==2) histo_h2E->Fill(adc_v->at(index));
	}*/

	if(filename.find("LST/")!=string::npos) filename.erase(filename.begin(),filename.begin()+4);
	Current_beta_file = filename;

	Path_beta_file = PATH;

}


/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  define in RefineBeta.h %%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
struct match_event{
	Long64_t B1_sweeps_gclock;
	Long64_t B1_time;
	Long64_t B1_time_relative;
	int B1_adc;
	Long64_t B2_sweeps_gclock;
	Long64_t B2_time;
	Long64_t B2_time_relative;
	int B2_adc;
	int delta_t;
	void clear(){
		B1_sweeps_gclock=-1;
		B1_time=-1;
		B1_time_relative=-1;
		B1_adc=-1;
		B2_sweeps_gclock=-1;
		B2_time=-1;
		B2_time_relative=-1;
		B2_adc=-1;
		delta_t=-1;
	}
};

struct hit_info{
	bool empty;   // container empty or not empty=> true; full=>false
	Long64_t sweeps_gclock;
	Long64_t time;
	Long64_t time_relative;
	int channel;
	int adc;
};
*/
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


vector<match_event>Match_hit; // store matching pairs

void KeepNewMatching(match_event _NewMatch){ // keep only new matching pairs; skip same ones

	for(unsigned long index=0;index< Match_hit.size();index++){
		Long64_t time1 = Match_hit[index].B1_time;
		Long64_t time2 = Match_hit[index].B2_time;
		if(time1==_NewMatch.B1_time && time2==_NewMatch.B2_time) return; // skip; do not save repeated data
	}

	Match_hit.push_back(_NewMatch); // store new one;
	

}

//void Compare_Beta_Beta(Long64_t gate_time=500,int gate_adc_low=0, int gate_adc_hi=1000){
bool Compare_Beta_Beta(){

	static string file_old;
	static Long64_t gate_time_old=-1;   // for remember the condition of last run, skip "beta Compare" if setting no change
	static int gate_adc_low_old=-1;         // for remember the condition of last run
	static int gate_adc_hi_old=-1;

	if(gate_time_old==gate_time && gate_adc_low_old==gate_adc_low && gate_adc_hi_old==gate_adc_hi && file_old==Current_beta_file){
		cout<<"\e[1;33m"<<"Beta-Beta coincident condition NO change! break"<<"\e[0m"<<endl;
		 return false;
	}

	gate_time_old = gate_time ;
	gate_adc_low_old = gate_adc_low;
	gate_adc_hi_old = gate_adc_hi;
	file_old=Current_beta_file;


	Match_hit.clear();	Match_hit.shrink_to_fit();
	match_event OneMatch;  // for one pair matching hits
	OneMatch.clear();

	hit_info candidate_tem;
	candidate_tem.empty=true;
	Long64_t progress=0; // show progress

	for(Long64_t line_index=0;line_index<nevt_bta;line_index++){
		progress++;
		if(progress%1000==0)cout<<"progress = "<<progress<<"\r"<<flush;
		//int progress_ratio =TMath::Nint((double)(progress/nevt_bta)*100);
		//if(progress_ratio%5==0)cout<<"\r"<<"progress= "<<progress_ratio<<"%"<<flush;

		candidate_tem.time = time_v->at(line_index);
		candidate_tem.time_relative = time_relative_v->at(line_index);
		candidate_tem.adc = adc_v->at(line_index);
		candidate_tem.channel = channel_v->at(line_index);
		candidate_tem.sweeps_gclock = sweeps_gclock_v->at(line_index);

		if(candidate_tem.adc >=gate_adc_low && candidate_tem.adc<=gate_adc_hi && candidate_tem.channel !=syn_CH && candidate_tem.sweeps_gclock!=0){ 
		// energy gate event for first waiting candidate
		// eliminate the event of sweeps_gclock ==0 ==> MCS of TOF not start yet.
			for(Long64_t forward_i=line_index+1;forward_i<nevt_bta;forward_i++){ // seeking coincident forward(towards end of file)
				Long64_t tem_sweeps_gclock= sweeps_gclock_v->at(forward_i);
				Long64_t tem_time = time_v->at(forward_i);
				Long64_t delta_time = TMath::Abs(tem_time-candidate_tem.time);
				Long64_t tem_time_relative = time_relative_v->at(forward_i);
				int tem_adc = adc_v->at(forward_i);
				int tem_channel = channel_v->at(forward_i);
				if(tem_channel == candidate_tem.channel||tem_channel==syn_CH){// sweeps trig or same silicon
					if(delta_time>gate_time)break;
					else continue; 
				}
				else{
					if(delta_time<=gate_time && tem_adc>=gate_adc_low && tem_adc<=gate_adc_hi){// matching
						if(candidate_tem.channel==1){// load data to match vector
							OneMatch.B1_sweeps_gclock = candidate_tem.sweeps_gclock;
							OneMatch.B1_time = candidate_tem.time;
							OneMatch.B1_time_relative = candidate_tem.time_relative;
							OneMatch.B1_adc = candidate_tem.adc;
							OneMatch.B2_sweeps_gclock = tem_sweeps_gclock;
							OneMatch.B2_time = tem_time;
							OneMatch.B2_time_relative = tem_time_relative;
							OneMatch.B2_adc = tem_adc;
							OneMatch.delta_t=(int)(OneMatch.B2_time-OneMatch.B1_time);//delta_time;
						}
						else{
							OneMatch.B1_sweeps_gclock = tem_sweeps_gclock;
							OneMatch.B1_time = tem_time;
							OneMatch.B1_time_relative = tem_time_relative;
							OneMatch.B1_adc = tem_adc;
							OneMatch.B2_sweeps_gclock = candidate_tem.sweeps_gclock;
							OneMatch.B2_time = candidate_tem.time;
							OneMatch.B2_time_relative = candidate_tem.time_relative;
							OneMatch.B2_adc = candidate_tem.adc;
							OneMatch.delta_t=(int)(OneMatch.B2_time-OneMatch.B1_time);//delta_time;
						}
						KeepNewMatching(OneMatch);
						OneMatch.clear();
					}
					else{
						if(delta_time>gate_time)break;
						continue;
					}
				}

			}// loop forward


/*			for(int backward_i=line_index-1;backward_i>=0;backward_i--){ // seeking coincident backward(towards beginning of file)
				Long64_t tem_sweeps_gclock= sweeps_gclock_v->at(backward_i);
				Long64_t tem_time = time_v->at(backward_i);
				Long64_t delta_time = TMath::Abs(tem_time-candidate_tem.time);
				int tem_adc = adc_v->at(backward_i);
				int tem_channel = channel_v->at(backward_i);
				if(tem_channel == candidate_tem.channel||tem_channel==syn_CH){// sweeps trig or same silicon
					if(delta_time>gate_time)break;
					else continue; 
				}
				else{
					if(delta_time<=gate_time && tem_adc>=gate_adc_low && tem_adc<=gate_adc_hi){// matching
						if(candidate_tem.channel==1){// load data to match vector
							OneMatch.B1_sweeps_gclock = candidate_tem.sweeps_gclock;
							OneMatch.B1_time = candidate_tem.time;
							OneMatch.B1_adc = candidate_tem.adc;
							OneMatch.B2_sweeps_gclock = tem_sweeps_gclock;
							OneMatch.B2_time = tem_time;
							OneMatch.B2_adc = tem_adc;
							OneMatch.delta_t=delta_time;
						}
						else{
							OneMatch.B1_sweeps_gclock = tem_sweeps_gclock;
							OneMatch.B1_time = tem_time;
							OneMatch.B1_adc = tem_adc;
							OneMatch.B2_sweeps_gclock = candidate_tem.sweeps_gclock;
							OneMatch.B2_time = candidate_tem.time;
							OneMatch.B2_adc = candidate_tem.adc;
							OneMatch.delta_t=delta_time;
						}
						KeepNewMatching(OneMatch);
						OneMatch.clear();
					}
					else{
						if(delta_time>gate_time)break;
						continue;
					}
				}

			}// loop backward   */


		}
		else{continue;}


	}// end of for loop line_index of file


	/*for(unsigned long result_i=0;result_i<Match_hit.size();result_i++){
		Long64_t B1_t = Match_hit[result_i].B1_time;
		int B1_adc = Match_hit[result_i].B1_adc;
		Long64_t B2_t = Match_hit[result_i].B2_time;
		int B2_adc = Match_hit[result_i].B2_adc;
		printf("%lld,%d \t %lld,%d \t %d\n",B1_t,B1_adc,B2_t,B2_adc,Match_hit[result_i].delta_t);
	}*/
	cout<<endl;
	cout<<"total event before refine= "<<Match_hit.size()<<endl;
	cout<<endl;

	Match_hit.shrink_to_fit();

	vector<match_event> Match_hit_refine;

	RefineMatchHit(Match_hit,Match_hit_refine);

	Match_hit.clear();Match_hit.shrink_to_fit();

	Match_hit = Match_hit_refine;

	Match_hit.shrink_to_fit();

/*
cout<<endl;
	for(unsigned long result_i=0;result_i<Match_hit.size();result_i++){
		Long64_t B1_t = Match_hit[result_i].B1_time;
		int B1_adc = Match_hit[result_i].B1_adc;
		Long64_t B2_t = Match_hit[result_i].B2_time;
		int B2_adc = Match_hit[result_i].B2_adc;
		printf("%lld,%d \t %lld,%d \t %d\n",B1_t,B1_adc,B2_t,B2_adc,Match_hit[result_i].delta_t);
	}*/

	cout<<"total event after refine= "<<Match_hit.size()<<endl;
	cout<<endl;

	ShowCoinCondition();

	return true;

}


bool Beta_Graph_Initialize(){
	//if(c_beta_beta!=NULL) delete c_beta_beta;
	if(histo_beta_deltaT!=NULL) delete histo_beta_deltaT;
	if(histo_beta_EkeV!=NULL) delete histo_beta_EkeV;
	if(histo_beta_EkeV_withcut!=NULL) delete histo_beta_EkeV_withcut;
	if(histo_beta_EkeV_CutAndRej!=NULL) delete histo_beta_EkeV_CutAndRej;
	if(histo_h1E!=NULL) delete histo_h1E;
	if(histo_h2E!=NULL) delete histo_h2E;
	if(tbeta_beta!=NULL) delete tbeta_beta;

	if(c_beta_beta==NULL){
		c_beta_beta = new TCanvas("c_beta_beta","c_beta_beta",1200,800);
		c_beta_beta->Divide(2,2);
	}
	
	c_beta_beta->cd(1);

	if(f_beta_tree!=NULL){
		if(f_beta_tree->IsOpen()){f_beta_tree->Close(); delete f_beta_tree;}
		else{delete f_beta_tree;}
	}

	//************ create file to store beta-beta tree, Avoid crash!!!!  ****************
	string tem_path = Path_beta_file +"../rootfiles/";
	bool PathOK = gSystem->AccessPathName(tem_path.c_str());
	if(!PathOK){ // path exist
		tem_path += "tbeta_beta_tree.root";
		f_beta_tree = new TFile(tem_path.c_str(),"RECREATE");
	}
	else{
		tem_path = Path_beta_file + "tbeta_beta_tree.root";
		f_beta_tree = new TFile(tem_path.c_str(),"RECREATE");
	}

	if(f_beta_tree->IsOpen()){f_beta_tree->cd();}
	else{cout<<"\e[1;33m"<<"Fail to create file to store beta-beta tree!!!! Abort!!!"<<"\e[0m"<<endl; return false;}

	tbeta_beta = new TTree("tbeta_beta","tbeta_beta");

	histo_beta_deltaT = new TH1D("histo_beta_deltaT","Coincidence time difference(Si2-Si1)",gate_time/20,-gate_time,gate_time);//gate_time/40,0,gate_time

	if(beta_slope[0]==0 && beta_slope[1]==0){
			SetBeta_E_Calibrate_para();
	}


	histo_beta_EkeV = new TH2D("histo_beta_EkeV","",200,0.,ADC2keV_1(1024),200,0.,ADC2keV_2(1024));
	//title "coincident Energy Si1 VS Si2" not set for easy to use cut
	histo_beta_EkeV->GetXaxis()->SetTitle("Si_1 [keV]"); histo_beta_EkeV->GetXaxis()->CenterTitle();
	histo_beta_EkeV->GetYaxis()->SetTitle("Si_2 [keV]"); histo_beta_EkeV->GetYaxis()->CenterTitle();

	histo_beta_EkeV_withcut = new TH2D("histo_beta_EkeV_withcut","EkeV_2:EkeV_1",200,0.,ADC2keV_1(1024),200,0.,ADC2keV_2(1024));  // show E-E 2D histo after cut
	histo_beta_EkeV_withcut->GetXaxis()->SetTitle("Si_1 [keV]"); histo_beta_EkeV_withcut->GetXaxis()->CenterTitle();
	histo_beta_EkeV_withcut->GetYaxis()->SetTitle("Si_2 [keV]"); histo_beta_EkeV_withcut->GetYaxis()->CenterTitle();

	histo_beta_EkeV_CutAndRej = new TH2D("histo_beta_EkeV_CutAndRej","EkeV_2:EkeV_1",200,0.,ADC2keV_1(1024),200,0.,ADC2keV_2(1024)); 
	histo_beta_EkeV_CutAndRej->GetXaxis()->SetTitle("Si_1 [keV]"); histo_beta_EkeV_CutAndRej->GetXaxis()->CenterTitle();
	histo_beta_EkeV_CutAndRej->GetYaxis()->SetTitle("Si_2 [keV]"); histo_beta_EkeV_CutAndRej->GetYaxis()->CenterTitle();

	//if(histo_h1E!=NULL) delete histo_h1E;
		histo_h1E = new TH1D("histo_h1E"," Si1 and Si2 ADC ",512,0,1024); histo_h1E->SetLineColor(kBlack);
	//if(histo_h2E!=NULL) delete histo_h2E;
		histo_h2E = new TH1D("histo_h2E"," Si1 and Si2 ADC ",512,0,1024); histo_h2E->SetLineColor(kRed);
	for(Long64_t index=0;index<nevt_bta;index++){
		if(channel_v->at(index)==1) histo_h1E->Fill(adc_v->at(index));
		if(channel_v->at(index)==2) histo_h2E->Fill(adc_v->at(index));
	}

	return true;

}

bool ShowBeta_Beta_Coin(){

	if(!Beta_Graph_Initialize()){
		cout<<"\e[1;33m"<<"Fail to initialize the Beta_Beta_coin Graph beacaue of failure of creating file to save tbeta_beta tree!!"<<"\e[0m"<<endl;
		return false;
	}
	c_beta_beta->cd(1);
	histo_h1E->Draw();
	histo_h2E->Draw("same");
	c_beta_beta->cd(1)->SetLogy();;

	if(Match_hit.size()==0){ cout<<"\e[1;23m"<<"No coincidence in store"<<"\e[0m"<<endl; return false;}


	Long64_t time1=0; Long64_t time1_relative=0; Long64_t sweeps_gclock_1=0; int adc1=0; double EkeV_1=0;
	Long64_t time2=0; Long64_t time2_relative=0; Long64_t sweeps_gclock_2=0; int adc2=0; double EkeV_2=0;
	int delta_t=0;


	tbeta_beta->Branch("sweeps_gclock_1",&sweeps_gclock_1);
	tbeta_beta->Branch("time1",&time1);
	tbeta_beta->Branch("time1_relative",&time1_relative);
	tbeta_beta->Branch("adc1",&adc1);
	tbeta_beta->Branch("EkeV_1",&EkeV_1);
	tbeta_beta->Branch("sweeps_gclock_2",&sweeps_gclock_2);
	tbeta_beta->Branch("time2",&time2);
	tbeta_beta->Branch("time2_relative",&time2_relative);
	tbeta_beta->Branch("adc2",&adc2);
	tbeta_beta->Branch("EkeV_2",&EkeV_2);
	tbeta_beta->Branch("delta_t",&delta_t);

	for(unsigned long index=0;index< Match_hit.size();index++){
			// histo_beta_deltaT->Fill(Match_hit[index].delta_t);	
		sweeps_gclock_1 = Match_hit[index].B1_sweeps_gclock;
		time1 = Match_hit[index].B1_time;
		time1_relative = Match_hit[index].B1_time_relative;
		adc1 = Match_hit[index].B1_adc;
		EkeV_1 = ADC2keV_1(adc1);
		sweeps_gclock_2 = Match_hit[index].B2_sweeps_gclock;
		time2 = Match_hit[index].B2_time;
		time2_relative = Match_hit[index].B2_time_relative;
		adc2 = Match_hit[index].B2_adc;
		EkeV_2 = ADC2keV_2(adc2);
		delta_t = Match_hit[index].delta_t;
		tbeta_beta->Fill();
	}

	c_beta_beta->cd(3);
	tbeta_beta->Draw("delta_t>>histo_beta_deltaT");

	c_beta_beta->cd(2);
	tbeta_beta->Draw("EkeV_2:EkeV_1>>histo_beta_EkeV","","colz");
	histo_beta_EkeV->SetTitle("EkeV_2:EkeV_1");

	cout<<"beta-beta tree is created:  tbeta_beta"<<endl;
	cout<<endl;

	histo_h1E->Write();
	histo_h2E->Write();
	tbeta_beta->Write();
	f_beta_tree->Write();
	tbeta_beta->ResetBranchAddresses();

	return true;

}


TCanvas* c_eject=NULL;
TGraph* Eject0[2] = {NULL,NULL};  // trigger leading edge ==>0; lagging edge ==>1
TGraph* Eject0_o[2] = {NULL,NULL};
TGraph* Eject1[2] = {NULL,NULL};
TGraph* Eject1_o[2] = {NULL,NULL};

void ShowEjectionBeta(double ejection0, double ejection1){  // ejection moment in ns

	ejection0 -= 2000; // mirror open 2us before daq delay
	ejection1 -=2000;

	if(c_eject==NULL){
		c_eject = new TCanvas("c_eject","Beta-Beta coincidence around ejection moment", 1200,800);
		c_eject->Divide(3,2);
	}

	if(Eject0[0]!=NULL){
	 	delete Eject0[0]; delete Eject0[1];
	 	delete Eject0_o[0]; delete Eject0_o[1];
	 	delete Eject1[0]; delete Eject1[1];
	 	delete Eject1_o[0]; delete Eject1_o[1];
	}

	double ejeW = 500 * 1000; // width of ejection trigger ==> 500 us to ns

	Eject0_o[0] = new TGraph(); Eject0_o[1] = new TGraph();  // point of eject
	Eject1_o[0] = new TGraph(); Eject1_o[1] = new TGraph();

	Eject0[0] = new TGraph(); Eject0[1] = new TGraph(); // points of beta-beta
	Eject1[0] = new TGraph(); Eject1[1] = new TGraph();

	for(int index=0;index<2;index++){
		Eject0_o[index]->SetPoint(0,ejection0+index*ejeW,ejection0+index*ejeW);
		if(index==0)Eject0_o[index]->SetPoint(1,ejection0+index*ejeW+3000,ejection0+index*ejeW+3000);
		else Eject0_o[index]->SetPoint(1,ejection0+index*ejeW+4000,ejection0+index*ejeW+4000);
		Eject0_o[index]->SetMarkerStyle(20);
		Eject0_o[index]->SetMarkerSize(1.1);

		if(index==0){
			Eject0_o[index]->SetTitle("time2_relative:time1_relative");
			Eject0_o[index]->GetXaxis()->SetTitle("Si_1 eje0 rising [ns]"); Eject0_o[index]->GetXaxis()->CenterTitle(); 
			Eject0_o[index]->GetXaxis()->SetLimits(ejection0+index*ejeW-1000,ejection0+index*ejeW+6000);
			Eject0_o[index]->GetYaxis()->SetTitle("Si_2 eje0 rising [ns]"); Eject0_o[index]->GetYaxis()->CenterTitle();
			Eject0_o[index]->SetMaximum(ejection0+index*ejeW+6000);
			Eject0_o[index]->SetMinimum(ejection0+index*ejeW-1000);

			Eject0[index]->SetTitle("time2_relative:time1_relative");
			Eject0[index]->GetXaxis()->SetTitle("Si_1 eje0 rising [ns]"); Eject0[index]->GetXaxis()->CenterTitle();
			Eject0[index]->GetXaxis()->SetLimits(ejection0+index*ejeW-1000,ejection0+index*ejeW+6000);
			Eject0[index]->GetYaxis()->SetTitle("Si_2 eje0 rising [ns]"); Eject0[index]->GetYaxis()->CenterTitle();
			Eject0[index]->SetMaximum(ejection0+index*ejeW+6000);
			Eject0[index]->SetMinimum(ejection0+index*ejeW-1000);
			Eject0[index]->SetMarkerStyle(31);
			Eject0[index]->SetMarkerSize(1.1);
			Eject0[index]->SetMarkerColor(kRed);
		}
		else{
			Eject0_o[index]->SetTitle("time2_relative:time1_relative");
			Eject0_o[index]->GetXaxis()->SetTitle("Si_1 eje0 dropping [ns]"); Eject0_o[index]->GetXaxis()->CenterTitle();
			Eject0_o[index]->GetXaxis()->SetLimits(ejection0+index*ejeW-1000,ejection0+index*ejeW+6000);
			Eject0_o[index]->GetYaxis()->SetTitle("Si_2 eje0 dropping [ns]"); Eject0_o[index]->GetYaxis()->CenterTitle();
			Eject0_o[index]->SetMaximum(ejection0+index*ejeW+6000);
			Eject0_o[index]->SetMinimum(ejection0+index*ejeW-1000);

			Eject0[index]->SetTitle("time2_relative:time1_relative");
			Eject0[index]->GetXaxis()->SetTitle("Si_1 eje0 dropping [ns]"); Eject0[index]->GetXaxis()->CenterTitle();
			Eject0[index]->GetXaxis()->SetLimits(ejection0+index*ejeW-1000,ejection0+index*ejeW+6000);
			Eject0[index]->GetYaxis()->SetTitle("Si_2 eje0 dropping [ns]"); Eject0[index]->GetYaxis()->CenterTitle();
			Eject0[index]->SetMaximum(ejection0+index*ejeW+6000);
			Eject0[index]->SetMinimum(ejection0+index*ejeW-1000);
			Eject0[index]->SetMarkerStyle(31);
			Eject0[index]->SetMarkerSize(1.1);
			Eject0[index]->SetMarkerColor(kRed);
		}

		Eject1_o[index]->SetPoint(0,ejection1+index*ejeW,ejection1+index*ejeW);
		if(index==0)Eject1_o[index]->SetPoint(1,ejection1+index*ejeW+3000,ejection1+index*ejeW+3000);
		else Eject1_o[index]->SetPoint(1,ejection1+index*ejeW+4000,ejection1+index*ejeW+4000);
		Eject1_o[index]->SetMarkerStyle(20);
		Eject1_o[index]->SetMarkerSize(1.1);

		if(index==0){	
			Eject1_o[index]->SetTitle("time2_relative:time1_relative"); 
			Eject1_o[index]->GetXaxis()->SetTitle("Si_1 eje1 rising [ns]"); Eject1_o[index]->GetXaxis()->CenterTitle(); 
			Eject1_o[index]->GetXaxis()->SetLimits(ejection1+index*ejeW-1000,ejection1+index*ejeW+6000);
			Eject1_o[index]->GetYaxis()->SetTitle("Si_2 eje1 rising [ns]"); Eject1_o[index]->GetYaxis()->CenterTitle();
			Eject1_o[index]->SetMaximum(ejection1+index*ejeW+6000);
			Eject1_o[index]->SetMinimum(ejection1+index*ejeW-1000);

			Eject1[index]->SetTitle("time2_relative:time1_relative");
			Eject1[index]->GetXaxis()->SetTitle("Si_1 eje1 rising [ns]"); Eject1[index]->GetXaxis()->CenterTitle();
			Eject1[index]->GetXaxis()->SetLimits(ejection1+index*ejeW-1000,ejection1+index*ejeW+6000);
			Eject1[index]->GetYaxis()->SetTitle("Si_2 eje1 rising [ns]"); Eject1[index]->GetYaxis()->CenterTitle();
			Eject1[index]->SetMaximum(ejection1+index*ejeW+6000);
			Eject1[index]->SetMinimum(ejection1+index*ejeW-1000);
			Eject1[index]->SetMarkerStyle(31);
			Eject1[index]->SetMarkerSize(1.1);
			Eject1[index]->SetMarkerColor(kRed);
		}
		else{
			Eject1_o[index]->SetTitle("time2_relative:time1_relative");
			Eject1_o[index]->GetXaxis()->SetTitle("Si_1 eje1 dropping [ns]"); Eject1_o[index]->GetXaxis()->CenterTitle();
			Eject1_o[index]->GetXaxis()->SetLimits(ejection1+index*ejeW-1000,ejection1+index*ejeW+6000);
			Eject1_o[index]->GetYaxis()->SetTitle("Si_2 eje1 dropping [ns]"); Eject1_o[index]->GetYaxis()->CenterTitle();
			Eject1_o[index]->SetMaximum(ejection1+index*ejeW+6000);
			Eject1_o[index]->SetMinimum(ejection1+index*ejeW-1000);

			Eject1[index]->SetTitle("time2_relative:time1_relative");
			Eject1[index]->GetXaxis()->SetTitle("Si_1 eje1 dropping [ns]"); Eject1[index]->GetXaxis()->CenterTitle();
			Eject1[index]->GetXaxis()->SetLimits(ejection1+index*ejeW-1000,ejection1+index*ejeW+6000);
			Eject1[index]->GetYaxis()->SetTitle("Si_2 eje1 dropping [ns]"); Eject1[index]->GetYaxis()->CenterTitle();
			Eject1[index]->SetMaximum(ejection1+index*ejeW+6000);
			Eject1[index]->SetMinimum(ejection1+index*ejeW-1000);
			Eject1[index]->SetMarkerStyle(31);
			Eject1[index]->SetMarkerSize(1.1);
			Eject1[index]->SetMarkerColor(kRed);
		}

	}// end of for format setting


	static vector <double> ejection0_si1_v[2];
	static vector <double> ejection0_si2_v[2];
	static vector <double> ejection1_si1_v[2];
	static vector <double> ejection1_si2_v[2];
	for(int iv=0;iv<2;iv++){
		ejection0_si1_v[iv].clear();
		ejection0_si2_v[iv].clear();
		ejection1_si1_v[iv].clear();
		ejection1_si2_v[iv].clear();
	}

	Long64_t si1_si2_avg=0;

	for(unsigned long index=0;index< Match_hit.size();index++){
			si1_si2_avg = (Match_hit[index].B1_time_relative + Match_hit[index].B2_time_relative)/2;
		for(int j=0;j<2;j++){
			if(si1_si2_avg > ejection0+j*ejeW-1000 && si1_si2_avg < ejection0+j*ejeW+6000){
				ejection0_si1_v[j].push_back( (double) Match_hit[index].B1_time_relative );
				ejection0_si2_v[j].push_back( (double) Match_hit[index].B2_time_relative );
				break;
			}
			if(si1_si2_avg > ejection1+j*ejeW-1000 && si1_si2_avg < ejection1+j*ejeW+6000){
				ejection1_si1_v[j].push_back( (double) Match_hit[index].B1_time_relative);
				ejection1_si2_v[j].push_back( (double) Match_hit[index].B2_time_relative);
				break;
			}
		}

	}

cout<<"eje0 up= "<<ejection0_si1_v[0].size()<<endl;
cout<<"eje0 down= "<<ejection0_si1_v[1].size()<<endl;
cout<<"eje1 up= "<<ejection1_si1_v[0].size()<<endl;
cout<<"eje1 down= "<<ejection1_si1_v[1].size()<<endl;

	for(int j=0;j<2;j++){
		for(unsigned long index=0;index< ejection0_si1_v[j].size();index++){
			Eject0[j]->SetPoint(index,ejection0_si1_v[j][index],ejection0_si2_v[j][index]);
		}
		for(unsigned long index=0;index< ejection1_si1_v[j].size();index++){
			Eject1[j]->SetPoint(index,ejection1_si1_v[j][index],ejection1_si2_v[j][index]);
		}
	}

	c_eject->cd(1)->SetMargin(0.120,0.052,0.1,0.1);
	c_eject->cd(1)->Modified();
	Eject0_o[0]->Draw("AP");
	if(Eject0[0]->GetN())Eject0[0]->Draw("P");


	c_eject->cd(2)->SetMargin(0.120,0.052,0.1,0.1);
	c_eject->cd(2)->Modified();
	Eject0_o[1]->Draw("AP");
	if(Eject0[1]->GetN())Eject0[1]->Draw("P");


	c_eject->cd(3)->SetMargin(0.120,0.052,0.1,0.1);
	c_eject->cd(3)->Modified();
	Eject1_o[0]->Draw("AP");
	if(Eject1[0]->GetN())Eject1[0]->Draw("P");


	c_eject->cd(4)->SetMargin(0.120,0.052,0.1,0.1);
	c_eject->cd(4)->Modified();
	Eject1_o[1]->Draw("AP");
	if(Eject1[1]->GetN())Eject1[1]->Draw("P");


}


