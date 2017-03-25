#include<iostream>
#include<vector>
#include<sstream>
#include<fstream>
#include<cmath>
#include<string>
#include<stdlib.h>

using namespace std;

int main(int argc,char** argv){

	if(argc==1) {
		cout<<"-------------------------------------"<<endl;
		cout<<"GetShuffle"<<endl;
                cout<<"G. Diana."<<endl;
		cout<<"-------------------------------------"<<endl;
                cout<<"\tUsage:"<<endl;
		cout<<"\tCCap <Genotype> <group> <folder name> <label> <gridsize>"<<endl;
		exit(0);
	}


	ifstream infile[6];
	stringstream inpdf_filename;
	inpdf_filename<<"/home/diana/workspace/Analysis/Information/"
		      <<"inpdf_"<<argv[1]<<".dat";
	ifstream inpdf_file(inpdf_filename.str().c_str());
	const int gridsize=atoi(argv[5]);
	const int group=atoi(argv[2]);
	int i,l,j,k,l1;
	int counter=0;

        stringstream filename;

	string foodlist[6] = {"1", "2e+07", "6.3e+07", "6.3e+08", "2e+09", "1.1e+10"};

	for(i=0;i<6;i++){
		filename.str("/home/diana/workspace/Analysis/R_projects/");
		filename.seekp(0, ios_base::end);
		filename<<argv[3]<<"/"<<argv[4]<<'_'<<argv[1]<<"_"<<foodlist[i]<<"_GS"<<gridsize<<"_group"<<group<<".dat";
		infile[i].open(filename.str().c_str());
	}
	
	ofstream outfile_info("info.dat");

	double joint[6][gridsize][gridsize][gridsize];
	double jointS[6][gridsize][gridsize][gridsize];

	for(i=0;i<6;i++){
		for(j=0;j<gridsize;j++){
			for(k=0;k<gridsize;k++){
				for(l=0;l<gridsize;l++){
					infile[i]>>joint[i][j][k][l];
				}
			}
		}
	}

	double input_pdf[6];
	for(k=0;k<6;k++){
		inpdf_file>>input_pdf[k];
	}
	inpdf_file.close();

	double pg=0;
	double pcg[6];
	double logsum[6];
	double Z=0;
	double Info=0;
	double Info_vec[6];
	int xyID=atoi(argv[2]);
	double PmargX[6][gridsize];
	double PmargY[6][gridsize];
	double PmargZ[6][gridsize];
	double PmargOverFood=0;


	
// NORMALIZATION and Initialization
	for(k=0;k<6;k++){
		Z=0;
		for(i=0;i<gridsize;i++){
		for(j=0;j<gridsize;j++){
		for(l=0;l<gridsize;l++){
			Z+=joint[k][i][j][l];
		}
		}
		}

		for(i=0;i<gridsize;i++){
		for(j=0;j<gridsize;j++){
		for(l=0;l<gridsize;l++){
			joint[k][i][j][l]/=Z;
		}
		}
		}
	}

	for(i=0;i<6;i++){
		Info_vec[i]=0;
	}


	for(i=0;i<gridsize;i++){
		for(k=0;k<6;k++){
			PmargX[k][i]=0;
			PmargY[k][i]=0;
			PmargZ[k][i]=0;
		}
	}

	for(k=0;k<6;k++){
	for(i=0;i<gridsize;i++){
	for(j=0;j<gridsize;j++){
	for(l=0;l<gridsize;l++){
		PmargX[k][i]+=joint[k][i][j][l];
		PmargY[k][i]+=joint[k][j][i][l];
		PmargZ[k][i]+=joint[k][j][l][i];
	}}}}
    

	for(i=0;i<gridsize;i++){
	for(j=0;j<gridsize;j++){
	for(l=0;l<gridsize;l++){
		for(k=0;k<6;k++){
			jointS[k][i][j][l]=PmargX[k][i]*PmargY[k][j]*PmargZ[k][l];
		}
	}}}
    
	for(i=0;i<gridsize;i++){
	for(j=0;j<gridsize;j++){
	for(l=0;l<gridsize;l++){
	    PmargOverFood=0;
		for(k=0;k<6;k++) PmargOverFood+=input_pdf[k]*jointS[k][i][j][l];
		for(k=0;k<6;k++){
			if(jointS[k][i][j][l]>0 && PmargOverFood>0) {
				Info+=jointS[k][i][j][l]*input_pdf[k]*log(jointS[k][i][j][l]/PmargOverFood)/log(2.);
				Info_vec[k]+=jointS[k][i][j][l]*log(jointS[k][i][j][l]/PmargOverFood)/log(2.);
			}
		}
	}
	}
	}

	for(i=0;i<6;i++) cout<<Info_vec[i]<<' ';
	cout<<Info<<endl;
}






					



