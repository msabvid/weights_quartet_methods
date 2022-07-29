/*
  Compilar: g++ -std=c++0x mainPhyloReconstruction.cpp -o mainPhyloReconstruction
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

using namespace std;



int main(int argc, char* argv[]) {
//    argv[1] = path del fasta original
//    argv[2] = nom del fasta original
//    argv[3] = arbre original tree: datacc, datacd o datadd


    string pathin = argv[1];
    string fname = argv[2];
    string original_tree = argv[3];

    string command;

    cout << "READFASTA" << endl;
    command = "./readFasta " + pathin + " " + fname;
    system(command.c_str());

    cout << "BCFT" << endl;
    command = "./BCFT 1 " + pathin + " " + fname;
    system(command.c_str());

    cout << "CFT" << endl;
    command = "./CFT 1 " + pathin + " " + fname;
    system(command.c_str());

    cout << "NJ" << endl;
    command = "./NJ " + pathin + " " + fname;
    system(command.c_str());

//    command = "./ML " + pathin + " " + fname;
//    system(command.c_str());


    command = "./QP BCFT " + pathin + " " + fname + " " + original_tree;
    system(command.c_str());
    command = "./QP CFT " + pathin + " " + fname + " " + original_tree;
    system(command.c_str());
    command = "./QP NJ " + pathin + " " + fname + " " + original_tree;
    system(command.c_str());
//    command = "./QP ML " + pathin + " " + fname + " " + original_tree;
//    system(command.c_str());


    command = "./WO BCFT " + pathin + " " + fname + " " + original_tree;
    system(command.c_str());
    command = "./WO CFT " + pathin + " " + fname + " " + original_tree;
    system(command.c_str());
    command = "./WO NJ " + pathin + " " + fname + " " + original_tree;
    system(command.c_str());
//    command = "./WO ML " + pathin + " " + fname + " " + original_tree;
//    system(command.c_str());


    command = "./Willson BCFT " + pathin + " " + fname + " " + original_tree;
    system(command.c_str());
    command = "./Willson CFT " + pathin + " " + fname + " " + original_tree;
    system(command.c_str());
    command = "./Willson NJ " + pathin + " " + fname + " " + original_tree;
    system(command.c_str());
//    command = "./Willson ML " + pathin + " " + fname + " " + original_tree;
//    system(command.c_str());


    return(0);
}
