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


struct w4Tree {
	vector<int> quartet;
	vector<double> weights;
};



class Tree {

	
	public:
		node *root;
		int nLeaves;
		vector<node*> leaves;
		vector<w4Tree> J;
		
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
		
		void QuartetReader(vector<int> quartet){
			w4Tree q;
			sort(quartet.begin(), quartet.end());
			q.quartet = quartet;
			//vector<double> output;
			string strFile="";
			for (int i=0; i<quartet.size(); i++){
				strFile.append(convertInt(quartet[i]));
				if (i<quartet.size()-1) {strFile.append("_");}
			}
			strFile.append(".txt");
			ifstream infile(strFile);
			double w;
			while (infile >> w) {
				q.weights.push_back(w);
			}
			J.push_back(q);
		}
		
		
		void EraseLeaf(node *position, int key) {
			//if (position==position->father->left) {
				//position->father->left = position->left;
			//} else {
				//position->father->right = position->left;
			//}
			//delete position->right;
			//delete position;
			if (position->father == position->father->father->left) {
				position->father = position->father->father;
				position->father->left = position;
			} else {
				position->father = position->father->father;
				position->father->right = position;
			}
			leaves[key-1]=NULL;
		}
		
		void InsertLeaf(node *position, int key) {
			node *NewIntNode;
			NewIntNode = splitEdge(position, -1); 
			insertNode(NewIntNode, key, "right", true);
			//NewIntNode->right = new node;
			//NewIntNode->right->father = intNode;
			//NewIntNode->right->weight = 0;
			//NewIntNode->right->key = key;
			//NewIntNode->right->left = NULL;
			//NewIntNode->right->right = NULL;
		}
		
		
		double LocalInconsistency(node *position, int key) {

			InsertLeaf(position, key);

			double max;
			bool first = true;
			double e;
			for (int i=0; i<J.size(); i++) {
				if (find(J[i].quartet.begin(),J[i].quartet.end(), key)!=J[i].quartet.end()) {
					e = excess(J[i]);
					if (first) {
						max = e;
						first = false;
					} else if (max<e) {
						max = e;
					}
				}
			}
			EraseLeaf(position, key);

			return(max);
		}
		
		double excess(w4Tree w4T){
			double output;
			vector<int> q = InducedTree(w4T.quartet);
			int ind = IndexSplit(q);
			vector<double> orderedWeights = w4T.weights;
			sort(orderedWeights.begin(), orderedWeights.end());
			if (w4T.weights[ind] = orderedWeights[0]) {
				output = orderedWeights[0]-orderedWeights[1];
			} else {
				output = w4T.weights[ind] - orderedWeights[0];
			}
			return(output);			
		}
		
		
		int IndexSplit(vector<int> q) {
			vector<int> q_aux = q;
			sort(q_aux.begin(), q_aux.end());
			int *min;
			if (q[0] == q_aux[0]) {
				min = &q[0];
			} else if (q[1] == q_aux[0]) {
				min = &q[1];
			} else if (q[2] == q_aux[0]) {
				min = &q[2];
			} else {
				min = &q[3];
			}
			if (q_aux[0]==q[0] || q_aux[0]==q[2]) {
				if (q_aux[1]==*(min+1)) {
					return 0;
				} else if (q_aux[2]==*(min+1)) {
					return 1;
				} else {
					return 2;
				} 
			} else {
				if (q_aux[1]==*(min-1)) {
					return 0;
				} else if (q_aux[2]==*(min-1)) {
					return 1;
				} else {
					return 2;
				} 
			}
		}
		
		
		node* IntersectionPaths(vector<node*> a, vector<node*> b){
			node *output = NULL;
			bool found = false;
			int i=0;
			int j=0;
			while (!found && i<a.size()) {
				while(!found && j<b.size()) {
					if (a[i]->key==b[j]->key) {
						found=true;
						output = a[i];
					} else {
						j=j+1;
					}
				}
				i=i+1;
			} 
			return(output);						
		}
		 
		
		vector<int> InducedTree(vector<int> q) {
			vector<node*> Path12 = GetPath(q[0],q[1]);
			vector<node*> Path34 = GetPath(q[2],q[3]);
			vector<node*> Path13 = GetPath(q[0],q[2]);
			vector<node*> Path24 = GetPath(q[1],q[3]);
			vector<int> output;
			if (IntersectionPaths(Path12, Path34)==NULL) {
				return(q);
			} else if (IntersectionPaths(Path13, Path24)==NULL) {
				output.push_back(q[0]);
				output.push_back(q[2]);
				output.push_back(q[1]);
				output.push_back(q[3]);
				return(output);
			} else {
				output.push_back(q[0]);
				output.push_back(q[3]);
				output.push_back(q[1]);
				output.push_back(q[2]);
				return(output);
			}
		}
		
		
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
		
				
		void InsertWeights(int key, node *start) {
			if (start->father!=NULL) {
				start->weight = LocalInconsistency(start, key);
				//cout << "hola\n";
				//exit(0);
			}
			if (start->left!=NULL) {
				//cout << "hola\n";
				InsertWeights(key, start->left);
				InsertWeights(key, start->right);
			}

		}
		
		

		
		vector<node*> LowestWeight(node *start) {
			vector<node*> output;
			vector<node*> outputLeft;
			vector<node*> outputRight;
			if (start->left==NULL) {
				output.push_back(start);
				return(output);
			} else if (start->father!=NULL) {
				outputLeft = LowestWeight(start->left);
				outputRight = LowestWeight(start->right);
				output.push_back(start);
				output.insert(output.end(), outputLeft.begin(), outputLeft.end());
				output.insert(output.end(), outputRight.begin(), outputRight.end());
				sort(output.begin(), output.end(), ComparebyWeight);
				//output.erase(output.begin()+2, output.end());
				output.erase(output.begin(), output.end()-2);
				return(output);
			} else {
				outputLeft = LowestWeight(start->left);
				outputRight = LowestWeight(start->right);		
				output = outputLeft;
				output.insert(output.end(), outputRight.begin(), outputRight.end());
				sort(output.begin(), output.end(), ComparebyWeight);
				//output.erase(output.begin()+2, output.end());
				output.erase(output.begin(), output.end()-2);
				return(output);
			}
		}
		
		
		double safety(node *node1, node *node2) {
			double output = (node1->weight - node2->weight);
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

//vector<double> QuartetReader(vector<int> quartet)
//{
	//sort(quartet.begin(), quartet.end());
	//vector<double> output;
	//string strFile="";
	//for (int i=0; i<quartet.size(); i++){
		//strFile.append(convertInt(quartet[i]));
		//if (i<quartet.size()-1) {strFile.append("_");}
	//}
	//strFile.append(".txt");
	//ifstream infile(strFile);
	//double w;
	//while (infile >> w) {
		//output.push_back(w);
	//}
	//return(output);
//}


int main() {
	srand(time(0));
	int n=5;
	Tree PhylogeneticTree(n);
	int RandomIndex;
	vector<int> R;
	vector<int> S;
	vector<w4Tree> J;
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
	PhylogeneticTree.QuartetReader(quartet);
	vector<double> weights = PhylogeneticTree.J[0].weights;
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
	PhylogeneticTree.printQP();
	
	//---------------------------//
	//exit(0);
	int count = n+3;
	vector<node*> Lowest;
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
						quartet.push_back(R[l]);
						sort(quartet.begin(), quartet.end());
						PhylogeneticTree.QuartetReader(quartet);
					}
				}
			}

			PhylogeneticTree.InsertWeights(R[l], PhylogeneticTree.root);
			
			Lowest = PhylogeneticTree.LowestWeight(PhylogeneticTree.root);
						

			if (l==0) {
				MaxSafety = PhylogeneticTree.safety(Lowest[0], Lowest[1]);
				index_MaxSafety = 0;
				SafestNode = Lowest[1];
			} else if (PhylogeneticTree.safety(Lowest[0], Lowest[1]) > MaxSafety){
				MaxSafety = PhylogeneticTree.safety(Lowest[0], Lowest[1]);
				index_MaxSafety = l;
				SafestNode = Lowest[1];
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
