//*************************************************************************
//* Offline decoder
//*
//*************************************************************************
#include <vector>
#include <TRint.h>
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>

#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

//define _LASER_   // switch to old verson(event by sweeps) by disable this line
#define LaserFequence 0.774

#if defined (MAKECINT)
//#pragma extra include "vector";
#pragma link C++ class std::vector<long>+;
#endif 

using namespace std;
/*
int sweepsMax=0; // for using in Online.C to known the maximum sweeps number allowed;
long gsweeps_tof=0;
int gsweeps_old_tof=-1;
bool HasData=false;*/

#ifndef _LASER_
unsigned long int ParserMCS6A_2(string PATH = "../",string filename = "mcs_39Kvs143X2plus@502_160142"){
  const int kMaxBufLen = 0x800000; // 8MB , get line of char in file
  static const int kstart = 6;     // the maximum of channel number, used in a "if()" condition
  unsigned long long int kChannel=0;    // = 0x00000007;
  int kShiftChannel         = 0;
  unsigned long long int kEdge=0;        //= 0x00000008;
  int kShiftEdge            =0;              //3;
  unsigned long long int kTimeData=0;    //= 0xfffffffff0;  //0xfffffff0;
  int kShiftData            = 0;              //4;
  unsigned long long int kSweeps=0;      //= 0x7f0000000000;  //0xffff00000000;
  int kShiftSweeps          = 0;            //40; //32;
  unsigned long long int kTagBit=0;      //= 0xffff000000000000;  //0x7fff000000000000;
  int kShiftTagBit          = 0;             //48;
  unsigned long long int kLost=0;        //= 0x800000000000; //0x8000000000000000;
  int kShiftLost            = 0;             //47;  //63;

  unsigned long long int sweepsclearbit=0;
  unsigned long long int sweepsfull=0;   // value = value of sweeps full + 1;  => 0b 1111 + 1 = 0b 10000;

  Long64_t sweeps_delta_s=0;
  Long64_t sweeps_delta_l=0;


  string inputfile = PATH + "LST/" + filename + ".lst";//"LST/"
  ifstream fin;
  fin.open(inputfile.data(),ios::in);
  if(!fin){
    cerr<<"cannot open LST input file: "<<inputfile<<" from decorder !!! "<<endl;
    return 0;
  }
  

  char *buffer = new char[kMaxBufLen];

  int Nbitshift=0;


  while(!fin.eof()){
      fin.getline(buffer,kMaxBufLen);
      if(fin.eof()){
        cerr<<"cannot find data!!!"<<endl;
	  fin.close();
        return 0;
      }

      //&&&&&&&&&&&&&&&7 read data format description  &&&&&&&&&&&&&&&&&&&&&&
      if(strncmp(buffer,";bit",4)==0){
            unsigned int bitL=0, bitH=0;

            string sbitL;
            sbitL=buffer[4]; 
            if(isdigit(buffer[5])) sbitL = sbitL + buffer[5];
            bitL = atoi(sbitL.c_str());

            string _buffer(buffer);

            string::size_type position =  _buffer.find("..");

            if(position<10){ // using more than one bit
                string sbitH;
                for(int i=(int)position +2;i<11;i++){
                    if(isdigit(buffer[i])) sbitH +=buffer[i];
                    if(buffer[i]==':') break;
                }

                bitH = atoi(sbitH.c_str());
                if(bitH==0){cout<<"error at reading bits format for decorder, break!!"<<endl; fin.close(); return 0;}


            }
            else{ bitH = bitL;} // using only one bit


            unsigned long long int GetBitRange=0;
            unsigned long long int one=1;

            for(unsigned int i=bitL;i<=bitH;i++){GetBitRange+= (one<<i);}  // generate bits to use & for saving corresponding bits



            if(_buffer.find("channel") != _buffer.npos){
                //printf("channel:%d .. %d; %llx ;kshift=%d \n ",bitL,bitH,GetBitRange,Nbitshift);
                kChannel = GetBitRange;
                kShiftChannel = Nbitshift;
            }
            else if(_buffer.find("edge") != _buffer.npos){
                //printf("kEdge:%d .. %d; %llx ;kshift=%d \n ",bitL,bitH,GetBitRange,Nbitshift);  
                kEdge = GetBitRange;
                kShiftEdge = Nbitshift;         
            }
            else if(_buffer.find("timedata") != _buffer.npos){
                 //printf("timedata:%d .. %d; %llx ;kshift=%d \n ",bitL,bitH,GetBitRange,Nbitshift);
                 kTimeData = GetBitRange;
                 kShiftData = Nbitshift;        
            }
            else if(_buffer.find("sweeps") != _buffer.npos){
                
                for(unsigned int i=0;i<(bitH-bitL)+1;i++){sweepsclearbit+= (one<<i);}

                sweepsfull = sweepsclearbit + 1;   // value = value of sweeps full + 1;  => 0b 1111 + 1 = 0b 10000;
		//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
		 //   sweepsMax = sweepsfull;  // only use in Online.C ==> to know the maximum number allowed sweeps number 
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                sweepsclearbit=~sweepsclearbit;
                //printf("sweeps:%d .. %d; %llx ;kshift=%d ;sweepclear=%llx ;sweepsfull value= %llx \n ",bitL,bitH,GetBitRange,Nbitshift,sweepsclearbit,sweepsfull);  
                kSweeps= GetBitRange;
                kShiftSweeps = Nbitshift;

            }
            else if(_buffer.find("tag") != _buffer.npos){
                 //printf("tag:%d .. %d; %llx ;kshift=%d \n ",bitL,bitH,GetBitRange,Nbitshift);    
                 kTagBit = GetBitRange;
                 kShiftTagBit = Nbitshift; 
            }
            else if(_buffer.find("data_lost") != _buffer.npos){
                //printf("data_lost:%d .. %d; %llx ;kshift=%d \n ",bitL,bitH,GetBitRange,Nbitshift); 
                kLost= GetBitRange;
                kShiftLost = Nbitshift;         
            }
            else{cout<<"Error; find unknow property of bits. Break!!!!"<<endl; fin.close(); return 0;}

            Nbitshift+= bitH-bitL +1;

      }


      if(strncmp(buffer,"[DATA]",6)==0){  // strncmp: compare buffer with "[DATA]" for 6 character, if equal will be 0; buffer> data -> value>0
                                                      // buffer< data  -> value<0; based on ACSII sequence
        cout<<"\e[1;33m reading data from "<<inputfile<<"\e[0m"<<endl;
        break;
      }
  }//end of file head reading for data format

  
  Long64_t nevt = 0;
  int nhits[10];
/*
  int& isweep_old= gsweeps_old_tof;// = -1; original
  int iglo =(int) gsweeps_tof; // =0 original
  long& isweep_global = gsweeps_tof; //=0 original*/
  int isweep_old=-1;// gsweeps_old_tof;// = -1; original
  int iglo =0;// (int) gsweeps_tof; // =0 original
  long isweep_global =0;// gsweeps_tof; //=0 original
  vector <Int_t> *channel = new vector <Int_t>();
  vector <Int_t> *edge = new vector <Int_t>();
  vector <long> *value = new vector <long>();
  vector <Double_t> *time = new vector <Double_t>();
  vector <Int_t> *sweeps = new vector <Int_t>();
  vector <long> * sweeps_global = new vector <long>;
  vector <Int_t> *tag = new vector <Int_t>();
  vector <Int_t> *lost = new vector <Int_t>();
vector <Int_t> *glo = new vector <Int_t>();

  string outputfile = PATH + "rootfiles/" + filename + ".root";
  TFile *fout = new TFile(outputfile.c_str(),"RECREATE");
  if(!fout->IsOpen()) {cout<<"Can not Recreate file: "<<outputfile<<" from decorder."<<endl; cout<<"Break!!"<<endl; return 0;}
  TTree *tree = new TTree("tree","rawdata tree");
  tree->Branch("nevt",&nevt,"nevt/L");
  tree->Branch("nhits",nhits,"nhits[10]/I");
  tree->Branch("channel",&channel);
  tree->Branch("edge",&edge);
  tree->Branch("value",&value); // in unit of 100ps
  tree->Branch("time",&time); // in unit of ns
  tree->Branch("sweeps",&sweeps);
  tree->Branch("sweeps_global",&sweeps_global);
tree->Branch("glo",&glo);
  tree->Branch("tag",&tag);
  tree->Branch("lost",&lost);

  unsigned long long int evtdata;
  string datastring;
  char datachar[100];
  istringstream is;
  while(!fin.eof()){

    if(nevt%1000==0) cout<<'\r'<<"read "<<nevt<<" events"<<flush;
    fin>>datachar; 
    datastring = string(datachar);
    is.clear();
    is.str(datastring);
    if(datastring.size()!=16){
	if(nevt==0 && strncmp(datachar,"[DATA]",6)==0) return 0; // in case there is no data after [DATA]
	else{	cout<<"read unknown data"<<endl; continue;}
     }
    is>>hex>>evtdata; 
   // fin >> hex >> evtdata;     // value in evtdata will be ox.... but anyway cout<< evtdata is in decimal; data is always binary in computer
    int ich = (evtdata&kChannel) >> kShiftChannel;  // extract 0...2 bits: channel:1~6 channel
    int iedge = (evtdata&kEdge) >> kShiftEdge;   // extract 4th bit; right shift 3 bit
    long ivalue = (evtdata&kTimeData) >> kShiftData;
    double itime = ivalue*0.1;   // resolution is 0.1 nanosecond
    int isweeps = (evtdata&kSweeps) >> kShiftSweeps;
    int itag = (evtdata&kTagBit) >> kShiftTagBit;
    int ilost = (evtdata&kLost) >> kShiftLost;

    // using tag0
    itag = itag&0x1;   // when itag==1 -> itag=1; when itag==0; -> itag =0; extract the first bit of tag value
    
    if(evtdata==0){
      continue;
    }else if(ich<1 || ich>kstart ){
      cerr<<"read unknow data: "<<evtdata<<endl;
      continue;
    }else if(isweeps != isweep_old){       //(ich == kstart && ivalue == 0){
      // new event
/***********************************************************************
        if(isweeps > isweep_old || isweeps>(isweep_old-10)){  
//cout<<"global sweep = "<<isweep_global<<" , isweeps "<<isweeps<<endl;
          //  isweep_global = (isweep_global&0xFFFFFFFFFFFF0000) + isweeps; // clear lower 4 digit   // 0xFFFFFFFFFFFF0000
          //      iglo = (iglo&0xFFFF0000) + isweeps;

          //  isweep_global = (isweep_global&0xFFFFFFFFFFFFFF80) + isweeps; // clear lower 7 digit (binary)  // 0xFFFFFFFFFFFF0000
          //      iglo = (iglo&0xFFFFFF80) + isweeps;

            isweep_global = (isweep_global&sweepsclearbit) + isweeps; // clear lower 7 digit (binary)  // 0xFFFFFFFFFFFF0000
                iglo = (iglo&sweepsclearbit) + isweeps;

        }else{

             // isweep_global = (isweep_global&0xFFFFFFFFFFFF0000) + isweeps + 0x10000;
//iglo = (iglo&0xFFFF0000) + isweeps + 0x10000;

              //isweep_global = (isweep_global&0xFFFFFFFFFFFFFF80) + isweeps + 0x80;
//iglo = (iglo&0xFFFFFF80) + isweeps + 0x80;

              isweep_global = (isweep_global&sweepsclearbit) + isweeps + (int)sweepsfull;
iglo = (iglo&sweepsclearbit) + isweeps + (int)sweepsfull;
          }
*********************************************************************************/

	if(isweeps > isweep_old){
		sweeps_delta_l = isweeps - isweep_old;   // if same around as isweep_old
		sweeps_delta_s = sweepsfull - isweeps + isweep_old; //if one around before isweep_old; round -1
		if(sweeps_delta_l > sweeps_delta_s){//last round;
			isweep_global = (isweep_global&sweepsclearbit) + isweeps - sweepsfull;
			iglo = (iglo&sweepsclearbit) + isweeps - sweepsfull;
		}
		else{ //same around as isweep_old;
			    isweep_global = (isweep_global&sweepsclearbit) + isweeps; // clear lower 7 digit (binary)  // 0xFFFFFFFFFFFF0000
                iglo = (iglo&sweepsclearbit) + isweeps;
		}
	}
	else{
		sweeps_delta_s = isweep_old - isweeps;  // if same around as isweep_old
		sweeps_delta_l = sweepsfull - isweep_old + isweeps;  //if one around after isweep_old; round+1
		if(sweeps_delta_s>sweeps_delta_l){ //next round
				isweep_global = (isweep_global&sweepsclearbit) + isweeps + sweepsfull;
				iglo = (iglo&sweepsclearbit) + isweeps + sweepsfull;
		}
		else{//same around as isweep_old
				isweep_global = (isweep_global&sweepsclearbit) + isweeps; // clear lower 7 digit (binary)  // 0xFFFFFFFFFFFF0000
		        iglo = (iglo&sweepsclearbit) + isweeps;
		}
	}


      isweep_old = isweeps;  //update isweep_old to new one,  new event: nevt++
//cout<<"global sweep = "<<isweep_global<<" , isweeps "<<isweeps<<endl;
//return ;
                    
      if(nevt>0) tree->Fill(); // do not fill empty data at the beginning

      // clean container
      channel->clear();
      edge->clear();
      value->clear();
      time->clear();
      sweeps->clear();
      sweeps_global->clear();
glo->clear();
      tag->clear();
      lost->clear();
      for(int i=0; i<10; i++) nhits[i] = 0;
      // load data
      nhits[ich]++;
      channel->push_back(ich);
      edge->push_back(iedge);
      value->push_back(ivalue);
      time->push_back(itime);
      sweeps->push_back(isweeps);
      sweeps_global->push_back(isweep_global);
glo->push_back(iglo);
      tag->push_back(itag);
      lost->push_back(ilost);
      nevt++;

    }else{
      if(!fin.eof()){  // prevent to fill repeated data  // multi-hits
      // real time data
      nhits[ich]++;
      channel->push_back(ich);
      edge->push_back(iedge);
      value->push_back(ivalue);
      time->push_back(itime);
      sweeps->push_back(isweeps);
      sweeps_global->push_back(isweep_global);
glo->push_back(iglo);
      tag->push_back(itag);
      lost->push_back(ilost);
       }
    }
  }
  cout<<'\r'<<flush<<"read "<<nevt<<" events"<<endl;
  cout<<"\e[1;33m writing tree to "<<outputfile<<"\e[0m"<<endl; // show color
  fout->cd();
  tree->Fill(); // fill last data
  tree->Write();
  
  delete [] buffer;

   fout->Write();
  //fout->Close();
delete fout;
  fin.clear();
  fin.seekg(0,ios::end);
  unsigned long int FileEnd = fin.tellg();

  fin.close();
   is.clear();
cout<<"check global sweeps: "<<isweep_global<<endl;
  return FileEnd;
}

#endif

/*
#ifndef __CINT__
int main(int argc, char *argv[]){
  if(argc>1){
    ParserMCS6A_2( string(argv[1]) );
  }
  return 0;
}
#endif
*/

#ifdef _LASER_

unsigned long int ParserMCS6A_2(string PATH = "../",string filename = "mcs_39Kvs143X2plus@502_160142"){
  const int kMaxBufLen = 0x800000; // 8MB , get line of char in file
  static const int kstart = 6;     // the maximum of channel number, used in a "if()" condition
  unsigned long long int kChannel=0;    // = 0x00000007;
  int kShiftChannel         = 0;
  unsigned long long int kEdge=0;        //= 0x00000008;
  int kShiftEdge            =0;              //3;
  unsigned long long int kTimeData=0;    //= 0xfffffffff0;  //0xfffffff0;
  int kShiftData            = 0;              //4;
  unsigned long long int kSweeps=0;      //= 0x7f0000000000;  //0xffff00000000;
  int kShiftSweeps          = 0;            //40; //32;
  unsigned long long int kTagBit=0;      //= 0xffff000000000000;  //0x7fff000000000000;
  int kShiftTagBit          = 0;             //48;
  unsigned long long int kLost=0;        //= 0x800000000000; //0x8000000000000000;
  int kShiftLost            = 0;             //47;  //63;

  unsigned long long int sweepsclearbit=0;
  unsigned long long int sweepsfull=0;   // value = value of sweeps full + 1;  => 0b 1111 + 1 = 0b 10000;


  string inputfile = PATH + "LST/" + filename + ".lst";//"LST/"
  ifstream fin;
  fin.open(inputfile.c_str());
  if(!fin){
    cerr<<"cannot open LST input file: "<<inputfile<<" from decorder !!! "<<endl;
    return 0;
  }


  char *buffer = new char[kMaxBufLen];

  int Nbitshift=0;


  while(!fin.eof()){
      fin.getline(buffer,kMaxBufLen);
      if(fin.eof()){
        cerr<<"cannot find data!!!"<<endl;
	  fin.close();
        return 0;
      }

      //&&&&&&&&&&&&&&&7 read data format description  &&&&&&&&&&&&&&&&&&&&&&
      if(strncmp(buffer,";bit",4)==0){
            unsigned int bitL=0, bitH=0;

            string sbitL;
            sbitL=buffer[4]; 
            if(isdigit(buffer[5])) sbitL = sbitL + buffer[5];
            bitL = atoi(sbitL.c_str());

            string _buffer(buffer);

            string::size_type position =  _buffer.find("..");

            if(position<10){ // using more than one bit
                string sbitH;
                for(int i=(int)position +2;i<11;i++){
                    if(isdigit(buffer[i])) sbitH +=buffer[i];
                    if(buffer[i]==':') break;
                }

                bitH = atoi(sbitH.c_str());
                if(bitH==0){cout<<"error at reading bits format for decorder, break!!"<<endl; fin.close(); return 0;}


            }
            else{ bitH = bitL;} // using only one bit


            unsigned long long int GetBitRange=0;
            unsigned long long int one=1;

            for(unsigned int i=bitL;i<=bitH;i++){GetBitRange+= (one<<i);}  // generate bits to use & for saving corresponding bits



            if(_buffer.find("channel") != _buffer.npos){
                //printf("channel:%d .. %d; %llx ;kshift=%d \n ",bitL,bitH,GetBitRange,Nbitshift);
                kChannel = GetBitRange;
                kShiftChannel = Nbitshift;
            }
            else if(_buffer.find("edge") != _buffer.npos){
                //printf("kEdge:%d .. %d; %llx ;kshift=%d \n ",bitL,bitH,GetBitRange,Nbitshift);  
                kEdge = GetBitRange;
                kShiftEdge = Nbitshift;         
            }
            else if(_buffer.find("timedata") != _buffer.npos){
                 //printf("timedata:%d .. %d; %llx ;kshift=%d \n ",bitL,bitH,GetBitRange,Nbitshift);
                 kTimeData = GetBitRange;
                 kShiftData = Nbitshift;        
            }
            else if(_buffer.find("sweeps") != _buffer.npos){
                
                for(unsigned int i=0;i<(bitH-bitL)+1;i++){sweepsclearbit+= (one<<i);}

                sweepsfull = sweepsclearbit + 1;   // value = value of sweeps full + 1;  => 0b 1111 + 1 = 0b 10000;
		//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
		  //  sweepsMax = sweepsfull;  // only use in Online.C ==> to know the maximum number allowed sweeps number 
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                sweepsclearbit=~sweepsclearbit;
                //printf("sweeps:%d .. %d; %llx ;kshift=%d ;sweepclear=%llx ;sweepsfull value= %llx \n ",bitL,bitH,GetBitRange,Nbitshift,sweepsclearbit,sweepsfull);  
                kSweeps= GetBitRange;
                kShiftSweeps = Nbitshift;

            }
            else if(_buffer.find("tag") != _buffer.npos){
                 //printf("tag:%d .. %d; %llx ;kshift=%d \n ",bitL,bitH,GetBitRange,Nbitshift);    
                 kTagBit = GetBitRange;
                 kShiftTagBit = Nbitshift; 
            }
            else if(_buffer.find("data_lost") != _buffer.npos){
                //printf("data_lost:%d .. %d; %llx ;kshift=%d \n ",bitL,bitH,GetBitRange,Nbitshift); 
                kLost= GetBitRange;
                kShiftLost = Nbitshift;         
            }
            else{cout<<"Error; find unknow property of bits. Break!!!!"<<endl; fin.close(); return 0;}

            Nbitshift+= bitH-bitL +1;

      }


      if(strncmp(buffer,"[DATA]",6)==0){  // strncmp: compare buffer with "[DATA]" for 6 character, if equal will be 0; buffer> data -> value>0
                                                      // buffer< data  -> value<0; based on ACSII sequence
        cout<<"\e[1;33m reading data from "<<inputfile<<"\e[0m"<<endl;
        break;
      }
  }// end of file head reading for data format
  


  Long64_t nevt = 0;
  int nhits[10];
  int isweep_old=-1;
int iglo = 0;
  long isweep_global = 0;
 int isweep_laserIn=-1;
 double time_laserIn=-1;
 double itime_gcDrift=0;

 double freqHZ= LaserFequence; // HZ laser frequence
 double cycle=(1/freqHZ)*1e9; // change to ns

  vector <Int_t> *channel = new vector <Int_t>();
  vector <Int_t> *edge = new vector <Int_t>();
  vector <long> *value = new vector <long>();
  vector <Double_t> *time = new vector <Double_t>();
  vector <Int_t> *sweeps = new vector <Int_t>();
  vector <long> * sweeps_global = new vector <long>;
  vector <Int_t> *tag = new vector <Int_t>();
  vector <Int_t> *lost = new vector <Int_t>();
vector <Int_t> *glo = new vector <Int_t>();
  vector <Double_t> *time_gcDrift = new vector<Double_t>();

  string outputfile = PATH + "rootfiles/" + filename + ".root";
  TFile *fout = new TFile(outputfile.c_str(),"RECREATE");
  if(!fout->IsOpen()) {cout<<"Can not Recreate file: "<<outputfile<<" from decorder."<<endl; cout<<"Break!!"<<endl; return 0;}
  TTree *tree = new TTree("tree","rawdata tree");
  tree->Branch("nevt",&nevt,"nevt/L");
  tree->Branch("nhits",nhits,"nhits[10]/I");
  tree->Branch("channel",&channel);
  tree->Branch("edge",&edge);
  tree->Branch("value",&value); // in unit of 100ps
  tree->Branch("time",&time); // in unit of ns
  tree->Branch("time_gcDrift",&time_gcDrift);  // drift time from gas cell
  tree->Branch("sweeps",&sweeps);
  tree->Branch("sweeps_global",&sweeps_global);
tree->Branch("glo",&glo);
  tree->Branch("tag",&tag);
  tree->Branch("lost",&lost);

  unsigned long long int evtdata;
  while(!fin.eof()){

    if(nevt%1000==0) cout<<'\r'<<"read "<<nevt<<" events"<<flush;
    
    fin >> hex >> evtdata;     // value in evtdata will be ox.... but anyway cout<< evtdata is in decimal; data is always binary in computer
    int ich = (evtdata&kChannel) >> kShiftChannel;  // extract 0...2 bits: channel:1~6 channel
    int iedge = (evtdata&kEdge) >> kShiftEdge;   // extract 4th bit; right shift 3 bit
    long ivalue = (evtdata&kTimeData) >> kShiftData;
    double itime = ivalue*0.1;   // resolution is 0.1 nanosecond
    int isweeps = (evtdata&kSweeps) >> kShiftSweeps;
    int itag = (evtdata&kTagBit) >> kShiftTagBit;
    int ilost = (evtdata&kLost) >> kShiftLost;

    // using tag0
    itag = itag&0x1;   // when itag==1 -> itag=1; when itag==0; -> itag =0; extract the first bit of tag value
    
    if(evtdata==0){
      continue;
    }else if(ich<1 || ich>kstart){
      cerr<<"read unknow data: "<<evtdata<<endl;
      continue;
    }else{ 
      // new event per ion
        if(isweeps > isweep_old || isweeps>(isweep_old-10)){  // because sometimes sweeps drop back a bit
//cout<<"global sweep = "<<isweep_global<<" , isweeps "<<isweeps<<endl;
          //  isweep_global = (isweep_global&0xFFFFFFFFFFFF0000) + isweeps; // clear lower 4 digit   // 0xFFFFFFFFFFFF0000
          //      iglo = (iglo&0xFFFF0000) + isweeps;

          //  isweep_global = (isweep_global&0xFFFFFFFFFFFFFF80) + isweeps; // clear lower 7 digit (binary)  // 0xFFFFFFFFFFFF0000
          //      iglo = (iglo&0xFFFFFF80) + isweeps;

            isweep_global = (isweep_global&sweepsclearbit) + isweeps; // clear lower 7 digit (binary)  // 0xFFFFFFFFFFFF0000
                iglo = (iglo&sweepsclearbit) + isweeps;

        }else if(isweeps<=(isweep_old-10)){

             // isweep_global = (isweep_global&0xFFFFFFFFFFFF0000) + isweeps + 0x10000;
//iglo = (iglo&0xFFFF0000) + isweeps + 0x10000;

              //isweep_global = (isweep_global&0xFFFFFFFFFFFFFF80) + isweeps + 0x80;
//iglo = (iglo&0xFFFFFF80) + isweeps + 0x80;

              isweep_global = (isweep_global&sweepsclearbit) + isweeps + (int)sweepsfull;
iglo = (iglo&sweepsclearbit) + isweeps + (int)sweepsfull;
        }

      isweep_old = isweeps;  //update isweep_old to new one,  new event: nevt++
//cout<<"global sweep = "<<isweep_global<<" , isweeps "<<isweeps<<endl;
//return ;
        
     // &&&&&&&& laser event  &&&&&&&&&&&&               
      if(ich==2){
          if(isweep_laserIn==-1){ // no laser in history; the first laser in
            isweep_laserIn = isweeps;
            time_laserIn = itime;
          }
          else if(isweeps-isweep_laserIn > 20){ // real laser in signal; avoid tripple trigger
            isweep_laserIn = isweeps;
            time_laserIn = itime;
          }

      }

      //&&&&&&&&&&&&&&&& calculate transportation time &&&&&&&&&&&&&&&&&&&

      if(isweep_laserIn==-1) itime_gcDrift=0;
      else{
        itime_gcDrift = (isweeps -isweep_laserIn)*25e6+itime-time_laserIn;  // 25 ms = 25e6 ns per sweep
        if(itime_gcDrift>cycle){  // maybe some laser trigger missing due to out of histo range of MCA
            itime_gcDrift = itime_gcDrift-(cycle * (int)(itime_gcDrift/cycle));  // use theoretical laser position before real laserin for update
        }
      }


      if(nevt>0) tree->Fill(); // do not fill empty data at the beginning
      // clean container
      channel->clear();
      edge->clear();
      value->clear();
      time->clear();
      sweeps->clear();
      sweeps_global->clear();
glo->clear();
      tag->clear();
      lost->clear();
      time_gcDrift->clear();
      for(int i=0; i<10; i++) nhits[i] = 0;
      // load data
      nhits[ich]++;
      channel->push_back(ich);
      edge->push_back(iedge);
      value->push_back(ivalue);
      time->push_back(itime);
      sweeps->push_back(isweeps);
      sweeps_global->push_back(isweep_global);
glo->push_back(iglo);
      tag->push_back(itag);
      lost->push_back(ilost);
      time_gcDrift->push_back(itime_gcDrift);
      nevt++;

    }
  }// end of while loop

  cout<<'\r'<<flush<<"read "<<nevt<<" events"<<endl;
  cout<<"\e[1;33m writing tree to "<<outputfile<<"\e[0m"<<endl; // show color
  fout->cd();
  //tree->Fill(); // fill last data
  tree->Write();
  
  delete [] buffer;

  fout->Write();
 // fout->Close();
delete fout;
  fin.clear()
  fin.seekg(0,ios::end);
  unsigned long int FileEnd = fin.tellg();
  fin.close();
  return  FileEnd;
}
#endif

