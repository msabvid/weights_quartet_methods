/******Weight Optimization******/

/*Compilar:
  g++ -std=c++0x WO.cpp -o WO -I/usr/include/python2.7/ -lpython2.7
  */

#include<Python.h>
#include<iostream>
#include<fstream>
#include<string>
#include<vector>
#include<random>
#include<set>
#include<algorithm>
#include<sstream>
#include<fstream>
#include<ctime>
#include<map>
//#include <pqxx/pqxx>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <pwd.h>



using namespace std;
//using namespace pqxx; //postgres



struct node {
	int key;
	node *father;
	node *left;
	node *right;
	double weight;
};


bool ComparebyWeight(node *a, node *b) {
	return(a->weight > b->weight);			
}



class Tree {

	
	public:
		node *root;
		int nLeaves;
		vector<node*> leaves;
		
		Tree(int n) {
			root = new node;
			root->father = NULL;
			root->weight = -1;
			root->key = 0;
			root->left = NULL;
			root->right = NULL;
			nLeaves = n;
			for (int i=0; i<nLeaves; i++) {
				leaves.push_back(new node);
				leaves[i] = NULL;
			}
		}
		
		//~Tree();
		
		node* splitEdge(node *intNode, int key){
			node *newNode = new node;
			newNode->father = intNode->father;
			newNode->weight = 0;
			newNode->key = key;
			newNode->left = intNode;
			newNode->right = NULL;
			if (intNode->father->left == intNode) {
				intNode->father->left = newNode;
			} else {
				intNode->father->right = newNode;
			}
			intNode->father = newNode;
			return(newNode);
		}
		
		void insertNode(node *intNode, int key, string side, bool IsLeaf){
			if (side == "right") {
				intNode->right = new node;
				intNode->right->father = intNode;
				intNode->right->weight = 0;
				intNode->right->key = key;
				intNode->right->left = NULL;
				intNode->right->right = NULL;
				if (IsLeaf) {
					leaves[key-1] = intNode->right;			
				}
			} else {
				intNode->left = new node;
				intNode->left->father = intNode;
				intNode->left->weight = 0;
				intNode->left->key = key;
				intNode->left->left = NULL;
				intNode->left->right = NULL;
				if (IsLeaf) {
					leaves[key-1] = intNode->left;			
				}				
			}
		}
		
		void FirstQuartet(vector<int> quartet) {
			insertNode(root, nLeaves +1, "left", false);
			insertNode(root, nLeaves +2, "right", false);
			insertNode(root->left, quartet[0], "left", true);
			insertNode(root->left, quartet[1], "right", true);
			insertNode(root->right, quartet[2], "left", true);
			insertNode(root->right, quartet[3], "right", true);
		}
		
		vector<node*> GetPath(int key1, int key2){
			vector<node*> up1;
			vector<node*> up2;
			vector<node*> output;
			node *node1 = leaves[key1-1];
			node *node2 = leaves[key2-1];
			node *parser;
			parser = node1;
			while (parser!=NULL) {
				up1.push_back(parser);
				parser = parser->father;
			}
			parser = node2;
			//up2.push_back(parser);
			while (parser!=NULL) {
				up2.push_back(parser);
				parser = parser->father;
			}
			
			bool found = false;
			int pos_up1=0;
			int pos_up2=0;
			while (!found) {
				while (!found && pos_up2<up2.size()) {
					if (up1[pos_up1]->key == up2[pos_up2]->key) {
						found = true;
					} else {
						pos_up2 = pos_up2+1;
					}
				}
				if (!found && (pos_up1+1)<up1.size()) {
					pos_up1 = pos_up1+1;
					pos_up2 = 0;
				}
			}
			for (int i=0; i<pos_up1; i++) {output.push_back(up1[i]);}
			for (int i=0; i<pos_up2; i++) {output.push_back(up2[i]);}
			output.push_back(up2[pos_up2]);
			return(output);
		}
		
		
		node* FindMiddleNode(int key1, int key2, int key3) {
			vector<node*> pathkey12 = GetPath(key1, key2);
			vector<node*> pathkey13 = GetPath(key1, key3);
			vector<node*> pathkey23 = GetPath(key2, key3);
			vector<node*> intersection1; 
			
			//set_intersection(pathkey12.begin(), pathkey12.end(), pathkey13.begin(), pathkey13.end(), back_inserter(intersection1));
			int j;
			bool intersected;
			for (int i=0; i<pathkey12.size(); i++) {
				j = 0;
				intersected = false;
				while (j<pathkey13.size() && !intersected) {
					if (pathkey12[i]->key == pathkey13[j]->key) {
						intersected = true;
						intersection1.push_back(pathkey12[i]);
					} else {j = j+1;}
				}  
			}
			
			vector<node*> intersection2;
			//set_intersection(pathkey23.begin(), pathkey23.end(), intersection1.begin(), intersection1.end(), back_inserter(intersection2));
			for (int i=0; i<pathkey23.size(); i++) {
				j = 0;
				intersected = false;
				while (j<intersection1.size() && !intersected) {
					if (pathkey23[i]->key == intersection1[j]->key) {
						intersected = true;
						intersection2.push_back(pathkey23[i]);
					} else {j = j+1;}
				}  
			}
			if (intersection2.size()>1) {
				cout << "Algo malament està passant!!\n";
			}
			return(intersection2[0]);
		}
		
		
		void InsertWeightsSubtree(node *start, node *u, bool uIsLeaf, double w) {
			if (start->father!=NULL) {
				start->weight += w;
			}
			if (start!=u || (start==u && !uIsLeaf)) {
				if (start->left!=NULL) {
					InsertWeightsSubtree(start->left, u, uIsLeaf, w);
				}
				if (start->right!=NULL) {
					InsertWeightsSubtree(start->right, u, uIsLeaf, w);
				}
			}
		}
		
		void InsertWeights(int key, node *u, double w) {
			node *parser = leaves[key-1];
			while (parser->father!=u && parser->father!=NULL) {
				parser = parser->father;
			}
			if (parser->father==u) {
				InsertWeightsSubtree(parser, u, false, w);
			} else {
				InsertWeightsSubtree(root, u, true, w);
			}
		}
		
		
		void InsertWeightsSplit(int key1, int key2, int key3, node *u, double w1, double w2, double w3) {
			//node *u = FindMiddleNode(key1, key2, key3);
			InsertWeights(key1, u, w1);
			InsertWeights(key2, u, w2);
			InsertWeights(key3, u, w3);
		}
		
		
	
		vector<node*> HighestWeight(node *start) {
			vector<node*> output;
			vector<node*> outputLeft;
			vector<node*> outputRight;
			if (start->left==NULL) {
				output.push_back(start);
				return(output);
			} else if (start->father!=NULL) {
				outputLeft = HighestWeight(start->left);
				outputRight = HighestWeight(start->right);
				output.push_back(start);
				output.insert(output.end(), outputLeft.begin(), outputLeft.end());
				output.insert(output.end(), outputRight.begin(), outputRight.end());
				sort(output.begin(), output.end(), ComparebyWeight);
				output.erase(output.begin()+2, output.end());
				return(output);
			} else {
				outputLeft = HighestWeight(start->left);
				outputRight = HighestWeight(start->right);		
				output = outputLeft;
				output.insert(output.end(), outputRight.begin(), outputRight.end());
				sort(output.begin(), output.end(), ComparebyWeight);
				output.erase(output.begin()+2, output.end());
				return(output);
			}
		}
		
		
		double safety(node *node1, node *node2) {
			double output = (node1->weight - node2->weight)/(node1->weight + node2->weight);
			return(output);
		}
		
				
		void InitEdges(node *start) {
			if (start->left!=NULL) {
				start->left->weight = 0;
				InitEdges(start->left);
			}
			if (start->right!=NULL) {
				start->right->weight = 0;
				InitEdges(start->right);
			}
		}

		string convertInt(int number) {
		   std::stringstream ss;//create a stringstream
		   ss << number;//add number to the stream
		   return ss.str();//return a string with the contents of the stream
		}

		
		
		string resultQP(node *start) {
			string output;
			if (start->left != NULL && start->right !=NULL) {
				output = "(";
				output.append(resultQP(start->left));
				output.append(", ");
				output.append(resultQP(start->right));
				output.append(")");
				return(output);
			} else {
				output = convertInt(start->key);
				return(output);
			}
		}

        string resultWOLong(node *start, bool first) {
            string output;
            if (start->left != NULL && start->right !=NULL) {
                output = "(";
                output.append(resultWOLong(start->left, false));
                output.append(",");
                output.append(resultWOLong(start->right, false));
                output.append(")");
                if (!first) {
                    output.append(":1");
                }
                return(output);
            } else {
                output = "t"+convertInt(start->key) + ":1";
                return(output);
            }
        }

        string getresultWO() {
            return(resultWOLong(root, true));
        }
		
		void printQP() {
			cout << resultQP(root) << "\n";
		}
		
		string allKeys(node *start) {
			string output;
			output = convertInt(start->key);
			if (start->left !=NULL) {
				output.append(allKeys(start->left));
			} 
			if (start->right != NULL) {
				output.append(allKeys(start->right));
			}
			return(output);
		}
		
	
		
		void PrintTree() {
			cout << allKeys(root) << "\n";
		}
					
			
};


string convertInt(int number) {
   std::stringstream ss;//create a stringstream
   ss << number;//add number to the stream
   return ss.str();//return a string with the contents of the stream
}

vector<double> QuartetReader(vector<int> quartet, string pathin)
{
    sort(quartet.begin(), quartet.end());
    vector<double> output;
    string strFile="";
    for (int i=0; i<quartet.size(); i++){
        strFile.append(convertInt(quartet[i]));
        if (i<quartet.size()-1) {strFile.append("_");}
    }
    strFile.append(".txt");
    strFile = pathin + strFile;
    ifstream infile(strFile.c_str());
    double w;
    while (infile >> w) {
        output.push_back(w);
    }
    return(output);
}





struct Alignment {
    unsigned int num_taxa;
    unsigned int seq_len;
    vector<string> taxa; // the taxa names
    map<string,string> seqs;  // the sequences
  //    unsigned int num_states;
//    map<char,int> table;
};


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


string WOmethod(string method, string originalPath, string fileName) {

    //string method = argv[1];
    //string fastaName = argv[2];

    string pathin = originalPath;
    Alignment align = MARCreadFASTA(pathin);
    int n = align.num_taxa;

    //cout << "n=" << n << endl;

    unsigned pos = fileName.find(".fa");
    string fastaName = fileName.substr(0,pos);

    int myuid;

    passwd *mypasswd;

    myuid = getuid();

    mypasswd = getpwuid(myuid);
    pathin = mypasswd->pw_dir;

    pathin = pathin + "/PhyloReconstruction/"+fastaName+"/Weights/"+method + "/";

    //srand(time(NULL));
    //int n=5;
	Tree PhylogeneticTree(n);
	int RandomIndex;
	vector<int> R;
	vector<int> S;
	vector<int> quartet;
	double aux;
	for (int i=1; i<=n; i++) {
		R.push_back(i);
	}

	//-------First w4tree-------//
	for (int i=1; i<=4; ++i) {
		RandomIndex = rand() % R.size();
		S.push_back(*(R.begin()+RandomIndex));
		R.erase(R.begin()+RandomIndex);
	}
	sort(S.begin(), S.end());
	quartet = S;
    vector<double> weights = QuartetReader(quartet, pathin);
	if (weights[1]>=weights[0] && weights[1]>=weights[2]) {
		aux = quartet[2];
		quartet[2] = quartet[1];
		quartet[1] = aux;
	} else if (weights[2]>=weights[0] && weights[2]>=weights[1]) {
		aux = quartet[3];
		quartet[3] = quartet[1];
		quartet[1] = aux;
	}
	PhylogeneticTree.FirstQuartet(quartet);
	//---------------------------//
	
	int count = n+3;
	//node *lowest;
	vector<node*> Highest;
	node *SafestNode;
	node *NewIntNode;
	node *u;
	double MaxSafety;
	int index_MaxSafety;
	while(!R.empty()) {
		for (int l=0; l<R.size(); l++) {
			PhylogeneticTree.InitEdges(PhylogeneticTree.root);
			for (int i=0; i<S.size()-2; i++) {
				for (int j=i+1; j<S.size()-1; j++) {
					for (int k=j+1; k<S.size(); k++) {
						weights.clear();
						quartet.clear();
						quartet.push_back(S[i]);
						quartet.push_back(S[j]);
						quartet.push_back(S[k]);
						//Middle node!!!!--------------
						u = PhylogeneticTree.FindMiddleNode(quartet[0], quartet[1], quartet[2]);
						//-----------------------------
						quartet.push_back(R[l]);
						sort(quartet.begin(), quartet.end());
                        weights = QuartetReader(quartet, pathin);
						if (quartet[0] == R[l]) {
							PhylogeneticTree.InsertWeightsSplit(quartet[1], quartet[2], quartet[3], u, weights[0], weights[1], weights[2]);
						} else if (quartet[1] == R[l]) {
							PhylogeneticTree.InsertWeightsSplit(quartet[0], quartet[2], quartet[3], u, weights[0], weights[2], weights[1]);
						} else if (quartet[2] == R[l]) {
							PhylogeneticTree.InsertWeightsSplit(quartet[0], quartet[1], quartet[3], u, weights[1], weights[2], weights[0]);
						} else {
							PhylogeneticTree.InsertWeightsSplit(quartet[0], quartet[1], quartet[2], u, weights[2], weights[1], weights[0]);
						}
					}
				}
			}
			Highest = PhylogeneticTree.HighestWeight(PhylogeneticTree.root);
			if (l==0) {
				MaxSafety = PhylogeneticTree.safety(Highest[0], Highest[1]);
				index_MaxSafety = 0;
				SafestNode = Highest[0];
			} else if (PhylogeneticTree.safety(Highest[0], Highest[1]) > MaxSafety){
				MaxSafety = PhylogeneticTree.safety(Highest[0], Highest[1]);
				index_MaxSafety = l;
				SafestNode = Highest[0];
			}
		}
		NewIntNode = PhylogeneticTree.splitEdge(SafestNode, count); 
		PhylogeneticTree.insertNode(NewIntNode, R[index_MaxSafety], "right", true);
		count = count+1;
		S.push_back(R[index_MaxSafety]);
		R.erase(R.begin() + index_MaxSafety);
		sort(S.begin(), S.end());
	}
    //cout << "The END" << endl;
    return(PhylogeneticTree.getresultWO());
}


long getDistance(string WOtree, string originalTree) {
    string ccTree = "(((((t1:1,t2:0.1):0.1,t3:0.1):0.1,(t4:0.1,t5:0.9):0.1):0.1,t6:0.8):0.1,(t7:0.8,((t8:0.9,t9:0.1):0.1,(t10:0.1,(t11:0.1,t12:0.9):0.1):0.1):0.1):0.1)";
    string cdTree = "(((((t1:0.9,t2:0.1):0.1,t3:0.1):0.1,(t4:0.1,t5:0.9):0.1):0.1,t6:0.8):0.1,(t7:0.8,(t8:0.1,((t9:0.9,t10:0.1):0.1,(t11:0.1,t12:0.9):0.1):0.1):0.1):0.1)";
    string ddTree = "(((((t1:0.9,t2:0.1):0.1,(t3:0.1,t4:0.9):0.1):0.1,t5:0.1):0.1,t6:0.8):0.1,(t7:0.8,(t8:0.1,((t9:0.9,t10:0.1):0.1,(t11:0.1,t12:0.9):0.1):0.1):0.1):0.1)";

    Py_Initialize();

    PyRun_SimpleString("import dendropy");
    string aux;
    if (originalTree=="datacc") {
        aux = "s1 = \"" +ccTree+"\"";
        PyRun_SimpleString(aux.c_str());
    } else if (originalTree=="datacd") {
        aux = "s1 = \"" +cdTree+"\"";
        PyRun_SimpleString(aux.c_str());
    } else {
        aux = "s1 = \"" +ddTree+"\"";
        PyRun_SimpleString(aux.c_str());
    }
    //PyRun_SimpleString("s2 = \"((t5:2.161175,t6:0.161175):0.392293,((t4:0.104381,(t2:0.075411,t1:0.075411):1):0.065840,t3:0.170221):0.383247)\"");
    aux = "s2 = \"" +WOtree+"\"";
    PyRun_SimpleString(aux.c_str());
    PyRun_SimpleString("tree1 = dendropy.Tree.get_from_string(s1, 'newick')");
    PyRun_SimpleString("tree2 = dendropy.Tree.get_from_string(s2, 'newick')");
    PyRun_SimpleString("result = tree1.symmetric_difference(tree2)");
    PyObject * module = PyImport_AddModule("__main__"); // borrowed reference
    //PyObject * module1 = PyImport_AddModule("dendropy");

    assert(module);                                     // __main__ should always exist
    PyObject * dictionary = PyModule_GetDict(module);   // borrowed reference
    assert(dictionary);                                 // __main__ should have a dictionary
    PyObject * result
    = PyDict_GetItemString(dictionary, "result");     // borrowed reference

    assert(result);                                     // just added result
    assert(PyInt_Check(result));                        // result should be an integer
    long result_value = PyInt_AS_LONG(result);          // already checked that it is an int

    std::cout << result_value << std::endl;

    Py_Finalize();
    return(result_value);
}






int main(int argc, char* argv[]) {
    cout << "WWWWWOOOOOOOO" << endl;

    //argv[1] = weighted_method
    //argv[2] = path
    //argv[3] = aliniament/filename
    //argv[4] = original tree
    //argv[5] = llargada de l'alineament
    string weighted_method = argv[1];
    string path = argv[2];
    string filename = argv[3];
    path = path+ "/" + filename;
    string original_tree = argv[4];
    //string alignment_length = argv[5];
    //string resultWOmethod;
    int dist;
    srand(time(NULL));





    Py_Initialize();

    PyRun_SimpleString("import dendropy");

    PyRun_SimpleString("trees = dendropy.TreeList()");
    string pyCommand;
    string resultWOmethod;
    for (int i=0; i<100; ++i) {
        //cout << "i=" << i << endl;
        resultWOmethod = WOmethod(weighted_method, path, filename);
        //cout << resultWOmethod << endl;
        pyCommand = "s  = \"" +resultWOmethod+"\"";
        PyRun_SimpleString(pyCommand.c_str());
        pyCommand = "trees.append(dendropy.Tree.get_from_string(s, 'newick'))";
        PyRun_SimpleString(pyCommand.c_str());
        //InsertQuery = "insert into temp values ('"+resultWOmethod+"')";
        //TEMP.exec(InsertQuery);
    }
    //TEMP.commit();
    pyCommand = "con_tree = trees.consensus(min_freq=0.5)";
    PyRun_SimpleString(pyCommand.c_str());
    pyCommand = "str_con_tree = con_tree.as_string('newick')";
    PyRun_SimpleString(pyCommand.c_str());

    PyObject * module = PyImport_AddModule("__main__"); // borrowed reference
    assert(module);                                     // __main__ should always exist
    PyObject * dictionary = PyModule_GetDict(module);   // borrowed reference
    assert(dictionary);                                 // __main__ should have a dictionary
    PyObject * str_con_tree = PyDict_GetItemString(dictionary, "str_con_tree");     // borrowed reference
    assert(str_con_tree);                                     // just added result
    //assert(PyInt_Check(con_tree));                        // result should be an integer
    string con_tree = PyString_AS_STRING(str_con_tree);          // already checked that it is an int
    unsigned pos1 = con_tree.find("(");
    con_tree = con_tree.substr(pos1);
    con_tree = con_tree.replace(con_tree.end()-1, con_tree.end(), "");
    //cout << "con_tree = " << con_tree << endl;
    Py_Finalize();

    string strDist = convertInt(getDistance(con_tree, original_tree));

    int myuid;

    passwd *mypasswd;

    myuid = getuid();

    mypasswd = getpwuid(myuid);
    string pathout = mypasswd->pw_dir;
    unsigned pos = filename.find(".fa");
    filename = filename.substr(0,pos);
    pathout = pathout + "/PhyloReconstruction/"+filename+"/Result";
    int status = mkdir(pathout.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    pathout = pathout +"/"+weighted_method;
    status = mkdir(pathout.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

    string outputName;
    ofstream outputFile;
    pathout = pathout + "/WO.txt";
    outputFile.open(pathout.c_str());
    outputFile << strDist << ";" << con_tree << "\n";
    outputFile.close();



    return(0);
}






