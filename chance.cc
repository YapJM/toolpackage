#include<iostream>
#include<stdlib.h>
#include<vector>
#include"TMath.h"
#include"TRandom3.h"
#include<ctime>

using namespace std;


void FakeCoincidence(int numtof=3,int numbeta=2,double Time_measure = 20,double coingate=1, int Maxloop=5){
	
		TRandom3* ranTOF = new TRandom3((unsigned int)time(0));
		TRandom3* ranBeta= new TRandom3(7192+(unsigned int)time(0)/209);

		/*	double Time_measure = 20; // second;
			double coingate= 1;  // second;

			int numtof=3;  // num of tof events
			int numbeta=2;  // num of beta events*/

		double* tof_time = new double[numtof];
		int* sequence_tof_list = new int[numtof];
		double* tof_time_s = new double[numtof];

		double* beta_time = new double[numbeta];
		int* sequence_beta_list = new int[numbeta];
		double* beta_time_s  = new double[numbeta];

		vector<double> tof_list;
		vector<double> beta_list;

		int numcase = TMath::Min(numtof,numbeta)+1;

		int* result_counts = new int[numcase];

		for(int index = 0; index<numcase;index++){
			result_counts[index]=0;
		}


		for(int nloop=0;nloop< Maxloop; nloop++){
				tof_list.clear();
				beta_list.clear();

				if(nloop%1000==0) cout<<"processing......"<<(double)nloop/Maxloop*100<<"%"<<"\r"<<flush;

				for(int index_tof=0; index_tof<numtof; index_tof++){

					tof_time[index_tof] =ranTOF->Uniform(Time_measure);

				}

				for(int index_beta=0;index_beta<numbeta;index_beta++){
					beta_time[index_beta] = ranBeta->Uniform(Time_measure);

				}


				TMath::Sort(numtof,tof_time,sequence_tof_list,kFALSE);
				TMath::Sort(numbeta,beta_time,sequence_beta_list,kFALSE);


				for(int index_tof=0; index_tof<numtof; index_tof++){

					tof_time_s[index_tof] = tof_time[sequence_tof_list[index_tof]];
					//cout<<tof_time_s[index_tof]<<"\t";

				}
				//cout<<endl;

				for(int index_beta=0; index_beta<numbeta; index_beta++){

					beta_time_s[index_beta] = beta_time[sequence_beta_list[index_beta]];
					//cout<<beta_time_s[index_beta]<<"\t";

				}
				//cout<<endl;


					for(int index_beta=0;index_beta<numbeta;index_beta++){

						for(int index_tof=0;index_tof<numtof;index_tof++){
							if(tof_time_s[index_tof]<beta_time_s[index_beta] && index_tof<numtof-1) continue;
							else if(tof_time_s[index_tof]>=beta_time_s[index_beta]){
										if( index_tof>=1){index_tof--;}
										else{break;}
							}
														
								
							
							if(beta_time_s[index_beta]-tof_time_s[index_tof] < coingate){ //candidate
								bool newtof = true, newbeta =true;
								for(int ilist=0;ilist<(int)tof_list.size();ilist++){
									if(tof_time_s[index_tof] == tof_list[ilist]){ newtof=false;break;} // repeat tof event
								}
								for(int jlist=0;jlist<(int)beta_list.size();jlist++){
									if(beta_time_s[index_beta] == beta_list[jlist]){newbeta=false;break;}
								}
								if(newtof && newbeta){
									tof_list.push_back(tof_time_s[index_tof]);
									beta_list.push_back(beta_time_s[index_beta]);
								}

							}
								
								break; // to next beta
							
						}// tof loop

					}// beta loop


					result_counts[tof_list.size()]++;
		}// end of nloop

			cout<<endl;
			cout<<"Num of coincidence \t Happens in loop \t Probability"<<endl;
			double mean=0;
			for(int index=0;index<numcase;index++){
				printf("%d \t %d \t %.2f%%\n",index,result_counts[index],(double)result_counts[index]/Maxloop*100);
				mean+=index*(double)result_counts[index]/Maxloop;
			}
			cout<<"mean = "<<mean<<endl;

	delete ranTOF;
	delete ranBeta;
	delete[] tof_time;
	delete[] sequence_tof_list;
	delete[] tof_time_s;
	delete[] beta_time;
	delete[] sequence_beta_list;
	delete[] beta_time_s;
	delete[] result_counts;




}