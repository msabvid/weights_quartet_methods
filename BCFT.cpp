/* dist_prog_wg
// This program generates 4 linear graphics for a given value of length and par_a (introduced as arguments)
// 2 graphics are for the distance to symmetric positive definite matrices: one with the percetange of success, and the other with the actual values of this distance
// 2 graphics are for the Frob. norm: one with the percentatge of sucess, and the other with the actual values of the tail of singular values. 

versió del 05/VII/2013
 *
Copyright (C) 2013 Casanellas and Fernandez-Sanchez.
This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
*/


/*compilar: g++ -std=c++0x BCFT.cpp functions.o graphic.o -o BCFT -O1 -larmadillo
  executar: ./BCFT 1
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

//#include <pqxx/pqxx>


using namespace std;	
using namespace arma; // for 'armadillo';
//using namespace pqxx; //postgres


/*
struct Alignment {
    unsigned int num_taxa;
    unsigned int seq_len;
    vector<string> taxa; // the taxa names
    map<string,string> seqs;  // the sequences
  //    unsigned int num_states;
//    map<char,int> table;
};*/


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



/*
string convertInt(int number) {
   std::stringstream ss;//create a stringstream
   ss << number;//add number to the stream
   return ss.str();//return a string with the contents of the stream
}*/





// global variables
//string nuc[4]={"A","C","G","T"}; already in functions.cpp


int main(int argc, char* argv[])
{
    //string filename=argv[2];
    //unsigned pos = filename.find(".");
    string filename;
    int n_partitions=atoi(argv[1]);


    //declaracio vars per despres postgres--------------
    //string AlignmentOriginal;
    Alignment align;
    int n_taxa;

    ofstream outfile;
    unsigned pos;

    string pathout;

    string outputname;
    string quartet;

    float coord1, coord2, coord3, sum;
    float tail1, tail2, tail3;
    int weight_12, weight_13, weight_14;

//    //postgres
//    connection C("dbname=tfm user=marc password=mochuelo hostaddr=127.0.0.1 port=5432");
//    string tableName("original_alignments");
//    if (C.is_open()) {

//        cout << "We are connected to " << C.dbname() << endl;

//    } else {

//        cout << "We are not connected!" << endl;
//    }

//    nontransaction N(C);

//    result R(N.exec("select absolute_path, fasta_name from "+tableName + " where absolute_path = '" + argv[2] +"'"));
//    if (!R.empty()) {
//        for (result::const_iterator c = R.begin(); c != R.end(); ++c) {


    //fname = c[1].as(string());
    string fname = argv[3];

    //pathin = c[0].as(string());
    string pathin = argv[2];

    pathin = pathin + "/" + fname;

    align = MARCreadFASTA(pathin);
    n_taxa = align.num_taxa;

    pos = fname.find(".fa");
    fname = fname.substr(0,pos);


    int myuid;

    passwd *mypasswd;

    myuid = getuid();

    mypasswd = getpwuid(myuid);
    pathout = mypasswd->pw_dir;
    pathin = pathout + "/PhyloReconstruction/" + fname + "/InputData/";
    pathout = pathout + "/PhyloReconstruction/"+fname+"/Weights";
    int status = mkdir(pathout.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    pathout = pathout + "/BCFT";
    status = mkdir(pathout.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);



//    pathin = "/home/marc/MAMME-TFM/data/AlignmentsQuartets/" + AlignmentOriginal + "/";
//    pathout = "/home/marc/MAMME-TFM/data/Weights/" + AlignmentOriginal;


//    commandmkdir = "mkdir " + pathout;
//    system(commandmkdir.c_str());
//    pathout = pathout + "/binary_cfrobtail";
//    commandmkdir = "mkdir " + pathout;
//    system(commandmkdir.c_str());

    pathout = pathout + "/";
    weight_12=0;
    weight_13=0;
    weight_14=0;

    for (int i1=0; i1<n_taxa-3; ++i1) {
     for (int i2=i1+1; i2<n_taxa-2; ++i2) {
         for (int i3=i2+1; i3<n_taxa-1; ++i3) {
             for (int i4=i3+1; i4<n_taxa; ++i4) {
                 quartet = convertInt(i1+1) + "_" + convertInt(i2+1) + "_" + convertInt(i3+1) + "_" + convertInt(i4+1);
                 outputname = pathout + quartet + ".txt";
                 outfile.open(outputname.c_str()); // cft means "corrected frobenius tail"
                 filename = pathin + quartet + ".fa";
                 //cout << filename << "\n";
                 weight_12=0;
                 weight_13=0;
                 weight_14=0;
                 map<string,float> tensor=my_readFASTA(filename); // obtain the tensor from the fasta file


                 // CORRECTED FROBENIUS TAIL
                 tail1=corrected(tensor,n_partitions,{0,1},{2,3});
                 tail2=corrected(tensor,n_partitions,{0,2},{1,3});
                 tail3=corrected(tensor,n_partitions,{0,3},{1,2});
                 coord1=tail2*tail3;
                 coord2=tail1*tail3;
                 coord3=tail1*tail2;
                 if (coord1 >=max(coord2,coord3)) {weight_12=1;}
                 if (coord2 >=max(coord1,coord3)) {weight_13=1;}
                 if (coord3 >=max(coord1,coord2)) {weight_14=1;}
                 outfile << weight_12 << endl;
                 outfile << weight_13 << endl;
                 outfile << weight_14 << endl;
                 outfile.close();
             }
         }
     }
    }


//        }
//    }





   return(0);
}



