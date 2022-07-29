/*
  compilar: g++ NJ.cpp functions.o graphic.o -o NJ -O1 -larmadillo
  */




#include "declaration.h"
#include <stdlib.h>
#include <stdio.h>
#include <algorithm>
#include <armadillo>
#include <string>
#include <fstream>
#include <iostream>
#include <time.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <pwd.h>


//per connectar-me a postgres
//#include <pqxx/pqxx>

//#include <direct.h>

using namespace std;    //using namespace std;
//using namespace pqxx; //postgres
using namespace arma; // for 'armadillo';

//string nuc[4]={"A","C","G","T"};
//int N=4;
////vector <string> observations=createpatterns(4);

//struct Alignment {
//    unsigned int num_taxa;
//    unsigned int seq_len;
//    vector<string> taxa; // the taxa names
//	map<string,string> seqs;  // the sequences
//  //    unsigned int num_states;
////    map<char,int> table;
//};


//reads from fasta file fname and saves the alignment into an alignment struct
Alignment MARCreadFASTA (string fname) {
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
    return align;
}








int main(int argc, char* argv[]) {

    //postgres---------------
    string filename;
    //  int n_partitions;

    string AlignmentOriginal;
    Alignment align;
    int n_taxa;

    ofstream outfile;
    unsigned pos;

    string pathout;
    string commandmkdir;

    string outputname;
    string quartet;

    float weight_01,weight_02,weight_03;
    float dist_01, dist_02, dist_03, sum;


    //fname = c[1].as(string());
    string fname = argv[2];

    //pathin = c[0].as(string());
    string pathin = argv[1];

    pathin = pathin + "/" + fname;

    align = MARCreadFASTA(pathin);
    n_taxa = align.num_taxa;

    pos = fname.find(".fa");
    fname = fname.substr(0,pos);

    string path_weight;

    int myuid;

    passwd *mypasswd;

    myuid = getuid();

    mypasswd = getpwuid(myuid);
    pathout = mypasswd->pw_dir;
    pathin = pathout + "/PhyloReconstruction/" + fname + "/InputData/";
    pathout = pathout + "/PhyloReconstruction/"+fname+"/Weights";
    int status = mkdir(pathout.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    pathout = pathout + "/NJ";
    status = mkdir(pathout.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);


    pathout = pathout + "/";

    for (int i1=0; i1<n_taxa-3; ++i1) {
        for (int i2=i1+1; i2<n_taxa-2; ++i2) {
            for (int i3=i2+1; i3<n_taxa-1; ++i3) {
                for (int i4=i3+1; i4<n_taxa; ++i4) {
                    quartet = convertInt(i1+1) + "_" + convertInt(i2+1) + "_" + convertInt(i3+1) + "_" + convertInt(i4+1);
                    outputname = pathout + quartet + ".txt";
                    outfile.open(outputname.c_str()); // cft means "corrected frobenius tail"
                    filename = pathin + quartet + ".fa";
                    //cout << filename << "\n";

                    // NEIGHBOUR-JOINING : PL_DISTANCE
                    align=readFASTA(filename);
                    mat dist= zeros<mat>(4,4); // distance matrix
                    dist=distance_PL(align);

                    dist_01=(1./4)*(dist(0,2)+dist(0,3)+dist(1,2)+dist(1,3))-(1./2)*(dist(0,1)+dist(2,3));
                    weight_01=exp(dist_01);
                    dist_02=(1./4)*(dist(0,1)+dist(0,3)+dist(1,2)+dist(2,3))-(1./2)*(dist(0,2)+dist(1,3));
                    weight_02=exp(dist_02);
                    dist_03=(1./4)*(dist(0,1)+dist(0,2)+dist(1,3)+dist(2,3))-(1./2)*(dist(0,3)+dist(1,2));
                    weight_03=exp(dist_03);
                    sum=weight_01+weight_02+weight_03;
                    outfile << weight_01/sum << endl;
                    outfile << weight_02/sum << endl;
                    outfile << weight_03/sum << endl;
                    outfile.close();
                }
            }
        }
    }


	return(0);
}

