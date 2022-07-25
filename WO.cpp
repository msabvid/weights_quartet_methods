/******Weight Optimization******/





#include<iostream>
#include<fstream>
#include<string>
#include<vector>
#include<random>
#include<set>
#include<algorithm>
#include <sstream>
#include<fstream>
#include<ctime>

using namespace std;


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

vector<double> QuartetReader(vector<int> quartet)
{
	sort(quartet.begin(), quartet.end());
	vector<double> output;
	string strFile="";
	for (int i=0; i<quartet.size(); i++){
		strFile.append(convertInt(quartet[i]));
		if (i<quartet.size()-1) {strFile.append("_");}
	}
	strFile.append(".txt");
	ifstream infile(strFile);
	double w;
	while (infile >> w) {
		output.push_back(w);
	}
	return(output);
}


int main() {
	srand(time(0));
	int n=5;
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
	vector<double> weights = QuartetReader(quartet);
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
						weights = QuartetReader(quartet);
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
	PhylogeneticTree.printQP();
	return(0);
}
