/*compilar: g++ ML.cpp functions.o graphic.o -o ML -O1 -larmadillo
  executar: ./ML
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


// global variables
//string nuc[4]={"A","C","G","T"}; already in functions.cpp




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


bool MLconvergence(string pathin)
{
    ifstream infile(pathin.c_str());
    string w;
    int count=0;
    while (infile >> w) {
        count = count+1;
    }
    cout << count << endl;
    if (count==3) {
        return(true);
    } else {
        return(false);
    }
}


int main(int argc, char* argv[]) {

    cout << "holan" << endl;

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
    pathout = pathout + "/ML";
    status = mkdir(pathout.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);


    pathout = pathout + "/";

    for (int i1=0; i1<n_taxa-3; ++i1) {
        for (int i2=i1+1; i2<n_taxa-2; ++i2) {
            for (int i3=i2+1; i3<n_taxa-1; ++i3) {
                for (int i4=i3+1; i4<n_taxa; ++i4) {
                    quartet = convertInt(i1+1) + "_" + convertInt(i2+1) + "_" + convertInt(i3+1) + "_" + convertInt(i4+1);
                    filename = pathin + quartet + ".fa";
                    cout << filename << endl;
                    commandmkdir = "./paml4species_killing.pl " + filename;
                    system(commandmkdir.c_str());
                    path_weight = pathout + "/" + quartet + ".txt";
                    if (!MLconvergence(path_weight)){
                        return(0);
                    } else {
                        cout << "YIIIHAAAAA" << endl;
                    }
                }
            }
        }
    }

    return(0);

}
