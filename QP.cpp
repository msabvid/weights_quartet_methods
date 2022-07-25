
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
	//node *left;
	//node *rigth;
	double weight;
};


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
			return(output);
		}
		
		
	
		node* HighestWeight(node *start) {
			double valueLeft;
			double valueRight;
			if (start->father==NULL) { 
				//root
				if (HighestWeight(start->left)>=HighestWeight(start->right)) {
					return(HighestWeight(start->left));
				} else {
					return(HighestWeight(start->right));
				}
			} else {	
				if (start->left!=NULL && start-> right!=NULL) {
					if (start->weight>=HighestWeight(start->left)->weight && start->weight>=HighestWeight(start->right)->weight) {
						return(start);
					} else if (HighestWeight(start->left)->weight >= start->weight && HighestWeight(start->left)->weight >= HighestWeight(start->right)->weight) {
						return HighestWeight(start->left);
					} else if (HighestWeight(start->right)->weight >= start->weight && HighestWeight(start->right)->weight >= HighestWeight(start->left)->weight) {
						return HighestWeight(start->right);
					}
				} else if (start->left!=NULL) {
					if (start->weight >= HighestWeight(start->left)->weight) {
						return start;
					} else {
						return HighestWeight(start->left);
					}
				} else if (start->right!=NULL) {
					if (start->weight >= HighestWeight(start->right)->weight) {
						return start;
					} else {
						return HighestWeight(start->right);
					}
				} else {
					return(start);
				}
			}
		}
		
				
		node* LowestWeight(node *start) {
			double valueLeft;
			double valueRight;
			if (start->father==NULL) { 
				//root
				if (LowestWeight(start->left)<=LowestWeight(start->right)) {
					return LowestWeight(start->left);
				} else {
					return LowestWeight(start->right);
				}
			} else {	
				if (start->left!=NULL && start-> right!=NULL) {
					if (start->weight<=LowestWeight(start->left)->weight && start->weight<=LowestWeight(start->right)->weight) {
						return start;
					} else if (LowestWeight(start->left)->weight<=start->weight && LowestWeight(start->left)->weight<=LowestWeight(start->right)->weight) {
						return LowestWeight(start->left);
					} else if (LowestWeight(start->right)->weight<=start->weight && LowestWeight(start->right)->weight<=LowestWeight(start->left)->weight) {
						return LowestWeight(start->right);
					}
				} else if (start->left!=NULL) {
					if (start->weight<=LowestWeight(start->left)->weight) {
						return start;
					} else {
						return LowestWeight(start->left);
					}
				} else if (start->right!=NULL) {
					if (start->weight<=LowestWeight(start->right)->weight) {
						return start;
					} else {
						return LowestWeight(start->right);
					}
				} else {
					return(start);
				}
			}
		}

		
		void InsertWeights(int key1, int key2, double w) {
			vector <node*> path = GetPath(key1, key2);
			for (int i=0; i<path.size(); i++) {
				path[i]->weight += w;
			}
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
	node *lowest;
	node *NewIntNode;
	while(!R.empty()) {
		PhylogeneticTree.InitEdges(PhylogeneticTree.root);
		RandomIndex = rand() % R.size();
		for (int i=0; i<S.size()-2; i++) {
			for (int j=i+1; j<S.size()-1; j++) {
				for (int k=j+1; k<S.size(); k++) {
					weights.clear();
					quartet.clear();
					quartet.push_back(S[i]);
					quartet.push_back(S[j]);
					quartet.push_back(S[k]);
					quartet.push_back(R[RandomIndex]);
					sort(quartet.begin(), quartet.end());
					weights = QuartetReader(quartet);
					if (quartet[0] == R[RandomIndex]) {
						PhylogeneticTree.InsertWeights(quartet[2], quartet[3], weights[0]);
						PhylogeneticTree.InsertWeights(quartet[1], quartet[3], weights[1]);
						PhylogeneticTree.InsertWeights(quartet[1], quartet[2], weights[2]);
					} else if (quartet[1] == R[RandomIndex]) {
						PhylogeneticTree.InsertWeights(quartet[2], quartet[3], weights[0]);
						PhylogeneticTree.InsertWeights(quartet[0], quartet[2], weights[1]);
						PhylogeneticTree.InsertWeights(quartet[0], quartet[3], weights[2]);
					} else if (quartet[2] == R[RandomIndex]) {
						PhylogeneticTree.InsertWeights(quartet[0], quartet[1], weights[0]);
						PhylogeneticTree.InsertWeights(quartet[1], quartet[3], weights[1]);
						PhylogeneticTree.InsertWeights(quartet[0], quartet[3], weights[2]);
					} else {
						PhylogeneticTree.InsertWeights(quartet[0], quartet[1], weights[0]);
						PhylogeneticTree.InsertWeights(quartet[0], quartet[2], weights[1]);
						PhylogeneticTree.InsertWeights(quartet[1], quartet[2], weights[2]);
					}
				}
			}
		}
		lowest = PhylogeneticTree.LowestWeight(PhylogeneticTree.root);
		//cout << lowest->key << "\n";
		NewIntNode = PhylogeneticTree.splitEdge(lowest, count); 
		PhylogeneticTree.insertNode(NewIntNode, R[RandomIndex], "right", true);
		//PhylogeneticTree.insertNode(PhylogeneticTree.splitEdge(PhylogeneticTree.LowestWeight(PhylogeneticTree.root), count), R[RandomIndex], "right", true);
		//cout << "NewIntNode " << NewIntNode->key << "; New Leaf " << NewIntNode->right->key << "\n";
		count = count+1;
		S.push_back(R[RandomIndex]);
		R.erase(R.begin()+RandomIndex);
		sort(S.begin(), S.end());
	}
	PhylogeneticTree.printQP();
	return(0);
}
