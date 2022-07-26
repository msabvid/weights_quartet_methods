/*Postgres:
  referencia: http://www.linux.com/community/blogs/133-general-linux/562240-postgresql-c-tutorial
  compilar: g++ readFasta.cpp -o readFasta -I/usr/local/include/ -lpqxx -lpq
  */




//#include "declaration.h"
#include <stdlib.h>
#include <stdio.h>
#include <algorithm>
#include <armadillo>
#include <string>
#include <fstream>
#include <iostream>
#include <map>
#include <sys/types.h>
#include <sys/stat.h>
#include <pwd.h>


//per connectar-me a postgres
//#include <pqxx/pqxx>

//#include <direct.h>

using namespace std;    //using namespace std;
//using namespace pqxx; //postgres
//using namespace arma; // for 'armadillo';

string nuc[4]={"A","C","G","T"};
int N=4;
//vector <string> observations=createpatterns(4);

struct Alignment {
    unsigned int num_taxa;
    unsigned int seq_len;
    vector<string> taxa; // the taxa names
	map<string,string> seqs;  // the sequences
  //    unsigned int num_states;
//    map<char,int> table;
};


//reads from fasta file fname and saves the alignment into an alignment struct
Alignment readFASTA (string fname) {     
    Alignment align;
    string line;
    // vector<string> species;
	fstream file;
	string speciesname;
	string dna="";
	string dnaaux;
	int numchar;
	int i;
	vector<int> al_length;
	//opening file
	file.open(fname.c_str(), fstream::in); //open(pointer to filename, type i.e. in/out)
	if (file == NULL) {
	  cout <<"cannot open file " << fname.c_str() << "\n";
		exit(1);
	}
	// get the first line from the FASTA file
	getline(file, line, '\n');//reads the first line
	//cout << line; //test
	// exit(1);
	if (line[0]!='>'){
          cout <<"error: not a FASTA file \n";
          exit(0);
    }
	numchar=line.length()-1;
	speciesname=line.substr(1,numchar);
    if (speciesname.length()==4) {
        speciesname = speciesname.substr(0,3)+"0"+speciesname.substr(3,4);
    }
	align.taxa.push_back(speciesname);
	// cout << align.taxa.front();
	if (align.seqs[speciesname].length() > 0) {
		cerr << "Warning - found 2+ sequences with same name " << speciesname << endl;        
    }
    while(! file.eof()){
		getline(file, line, '\n');

		if (line[0]!='>'){
			dna += line;
			//  cout << line;
		}else{ 
			//cout << dna <<endl;
			align.seqs[speciesname]=dna;
			al_length.push_back(dna.length());
			numchar=line.length()-1;
			speciesname=line.substr(1,numchar);
            if (speciesname.length()==4) {
                speciesname = speciesname.substr(0,3)+"0"+speciesname.substr(3,4);
            }
			align.taxa.push_back(speciesname);
			dna.clear();

			if (align.seqs[speciesname].length() > 0) {
				cerr << "Warning - found 2+ sequences with same name " << speciesname << endl;        
			}
		}
    }
    align.seqs[speciesname]=dna;
    al_length.push_back(dna.length());
    file.close();
    align.num_taxa=align.taxa.size();
    if(al_length.size() != align.num_taxa){
        cout << "error in getting the alignment!!!";
    }
    for (i=0; i< al_length.size()-1; i++){
        if (al_length[i] != al_length[i+1]){
           cout << "sequence " <<align.taxa[i+1] <<" doesn't have same length \n";
            exit(0);
          }
    }
    align.seq_len=al_length.back();
    sort(align.taxa.begin(), align.taxa.end());
    return align;
}



string convertInt(int number) {
   std::stringstream ss;//create a stringstream
   ss << number;//add number to the stream
   return ss.str();//return a string with the contents of the stream
}


void writeFASTA (Alignment align, string fname) {
	
	string outputName;
	string q="";
	ofstream outputFile;
	string taxonName;
	const char *c;
	for (int i=0; i<align.num_taxa-3; i++) {
		for (int j=i+1; j<align.num_taxa-2; j++) {
			for (int k=j+1; k<align.num_taxa-1; k++) {
				for (int l=k+1;l<align.num_taxa; l++) {
					q = convertInt(i+1) + "_" + convertInt(j+1) + "_" + convertInt(k+1) + "_" + convertInt(l+1);
                    outputName = fname + q + ".fa";
                    //cout << outputName << "\n";
					//ofstream outputFile;
					outputFile.open(outputName.c_str());
                    if (outputFile==NULL) {
                        cout << "cannot open file " << outputName.c_str() << "\n";
                        exit(1);
                    }
					if (align.taxa[i].size()>=3 && align.taxa[i].substr(0,3)=="seq") {
						c = align.taxa[i].substr(3).c_str();
						taxonName = convertInt(atoi(c)+1);
					} else {
						taxonName = align.taxa[i];
					}
					//outputFile << ">" << align.taxa[i] << "\n";
                    //outputFile << ">" << taxonName << "\n";
                    outputFile << ">1" << "\n";
					outputFile << align.seqs[align.taxa[i]] << "\n";
					if (align.taxa[j].size()>=3 && align.taxa[j].substr(0,3)=="seq") {
						c = align.taxa[j].substr(3).c_str();
						taxonName = convertInt(atoi(c)+1);
					} else {
						taxonName = align.taxa[j];
					}
					//outputFile << ">" << align.taxa[j] << "\n";
                    //outputFile << ">" << taxonName << "\n";
                    outputFile << ">2" << "\n";
					outputFile << align.seqs[align.taxa[j]] << "\n";
					if (align.taxa[k].size()>=3 && align.taxa[k].substr(0,3)=="seq") {
						c = align.taxa[k].substr(3).c_str();
						taxonName = convertInt(atoi(c)+1);
					} else {
						taxonName = align.taxa[k];
					}
					//outputFile << ">" << align.taxa[k] << "\n";
                    //outputFile << ">" << taxonName << "\n";
                    outputFile << ">3" << "\n";
					outputFile << align.seqs[align.taxa[k]] << "\n";
					if (align.taxa[l].size()>=3 && align.taxa[l].substr(0,3)=="seq") {
						c = align.taxa[l].substr(3).c_str();
						taxonName = convertInt(atoi(c)+1);
					} else {
						taxonName = align.taxa[l];
					}					
					//outputFile << ">" << align.taxa[l] << "\n";
                    //outputFile << ">" << taxonName << "\n";
                    outputFile << ">4" << "\n";
					outputFile << align.seqs[align.taxa[l]] << "\n";
					outputFile.close();
				}
			}
		}
	}
}


int main(int argc, char* argv[]) {


    string fname;
    string pathin;
    Alignment align;
    unsigned pos;
    string pathout;

    //if (!R.empty()) {
    //    for (result::const_iterator c = R.begin(); c != R.end(); ++c) {

    //fname = c[1].as(string());
    fname = argv[2];

    //pathin = c[0].as(string());
    pathin = argv[1];

    pathin = pathin + "/" + fname;

    align = readFASTA(pathin);
    pos = fname.find(".fa");
    fname = fname.substr(0,pos);

    //cmdmkdir = "mkdir /home/marc/MAMME-TFM/data/AlignmentsQuartets/" + fname;
    //system(cmdmkdir.c_str());

    int myuid;

    passwd *mypasswd;

    myuid = getuid();

    mypasswd = getpwuid(myuid);
    pathout = mypasswd->pw_dir;
    pathout = pathout + "/PhyloReconstruction";
    int status = mkdir(pathout.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    pathout = pathout + "/"+fname;
    status = mkdir(pathout.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    pathout = pathout + "/InputData";
    status = mkdir(pathout.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    pathout = pathout + "/";
    writeFASTA(align, pathout);

	return(0);
}

