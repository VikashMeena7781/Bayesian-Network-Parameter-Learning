#include <iostream>
#include <string>
#include <vector>
#include <list>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include<bits/stdc++.h>


// Format checker just assumes you have Alarm.bif and Solved_Alarm.bif (your file) in current directory
using namespace std;

// Our graph consists of a list of nodes where each node is represented as follows:
class Graph_Node{

private:
	string Node_Name;  // Variable name 
	vector<int> Children; // Children of a particular node - these are index of nodes in graph.
	vector<string> Parents; // Parents of a particular node- note these are names of parents
	int nvalues;  // Number of categories a variable represented by this node can take
	vector<string> values; // Categories of possible values
	vector<float> CPT; // conditional probability table as a 1-d array . Look for BIF format to understand its meaning

public:
	// Constructor- a node is initialised with its name and its categories
    Graph_Node(string name,int n,vector<string> vals)
	{
		Node_Name=name;
	
		nvalues=n;
		values=vals;
		

	}
	string get_name()
	{
		return Node_Name;
	}
	vector<int> get_children()
	{
		return Children;
	}
	vector<string> get_Parents()
	{
		return Parents;
	}
	vector<float> get_CPT()
	{
		return CPT;
	}
	int get_nvalues()
	{
		return nvalues;
	}
	vector<string> get_values()
	{
		return values;
	}
	void set_CPT(vector<float> new_CPT)
	{
		CPT.clear();
		CPT=new_CPT;
	}
    void set_Parents(vector<string> Parent_Nodes)
    {
        Parents.clear();
        Parents=Parent_Nodes;
    }
    // add another node in a graph as a child of this node
    int add_child(int new_child_index )
    {
        for(int i=0;i<Children.size();i++)
        {
            if(Children[i]==new_child_index)
                return 0;
        }
        Children.push_back(new_child_index);
        return 1;
    }



};
 // The whole network represted as a list of nodes
class network{
	list <Graph_Node> Pres_Graph;

public:
	int addNode(Graph_Node node)
	{
		Pres_Graph.push_back(node);
		return 0;
	}
    
    
	int netSize()
	{
		return Pres_Graph.size();
	}

	list<Graph_Node>& getPresGraph() {
        return Pres_Graph;
    }
    // get the index of node with a given name
    int get_index(string val_name)
    {
        list<Graph_Node>::iterator listIt;
        int count=0;
        for(listIt=Pres_Graph.begin();listIt!=Pres_Graph.end();listIt++)
        {
            if(listIt->get_name().compare(val_name)==0)
                return count;
            count++;
        }
        return -1;
    }
// get the node at nth index
    list<Graph_Node>::iterator get_nth_node(int n)
    {
       	list<Graph_Node>::iterator listIt;
        int count=0;
        for(listIt=Pres_Graph.begin();listIt!=Pres_Graph.end();listIt++)
        {
            if(count==n)
                return listIt;
            count++;
        }
        return listIt; 
    }
    //get the iterator of a node with a given name
    list<Graph_Node>::iterator search_node(string val_name)
    {
        list<Graph_Node>::iterator listIt;
        for(listIt=Pres_Graph.begin();listIt!=Pres_Graph.end();listIt++)
        {
            if(listIt->get_name().compare(val_name)==0)
                return listIt;
        }
        cout<<"node not found\n";
        return listIt;
    }
};



network read_network()
{
	network Alarm;
	string line;
	int find=0;
  	ifstream myfile("alarm.bif"); 
  	string temp;
  	string name;
  	vector<string> values;
  	
    if (myfile.is_open())
    {
    	while (!myfile.eof())
    	{
    		stringstream ss;
      		getline (myfile,line);
      		
      		
      		ss.str(line);
     		ss>>temp;
     		
     		
     		if(temp.compare("variable")==0)
     		{
                    
     				ss>>name;
     				getline (myfile,line);
                   
     				stringstream ss2;
     				ss2.str(line);
     				for(int i=0;i<4;i++)
     				{
     					ss2>>temp;
     				}
     				values.clear();
     				while(temp.compare("};")!=0)
     				{
     					values.push_back(temp);
     					
     					ss2>>temp;
    				}
     				Graph_Node new_node(name,values.size(),values);
     				int pos=Alarm.addNode(new_node);

     				
     		}
     		else if(temp.compare("probability")==0)
     		{
                    
     				ss>>temp;
     				ss>>temp;
     				
                    list<Graph_Node>::iterator listIt;
                    list<Graph_Node>::iterator listIt1;
     				listIt=Alarm.search_node(temp);
                    int index=Alarm.get_index(temp);
                    ss>>temp;
                    values.clear();
     				while(temp.compare(")")!=0)
     				{
                        listIt1=Alarm.search_node(temp);
                        listIt1->add_child(index);
     					values.push_back(temp);
     					
     					ss>>temp;

    				}
                    listIt->set_Parents(values);
    				getline (myfile,line);
     				stringstream ss2;
                    
     				ss2.str(line);
     				ss2>> temp;
                    
     				ss2>> temp;
                    
     				vector<float> curr_CPT;
                    string::size_type sz;
     				while(temp.compare(";")!=0)
     				{
                        
     					curr_CPT.push_back(atof(temp.c_str()));
     					
     					ss2>>temp;
                       
                        

    				}
                    
                    listIt->set_CPT(curr_CPT);


     		}
            else
            {
                
            }
    		
    	}
    	if(find==1)
    	myfile.close();
  	}else{
		cout << "Unable to open file";
	}
  	return Alarm;
}


vector<vector<float>> finaL_prob()
{
	network Alarm;
	string line;
	int find=0;
  	ifstream myfile("gold_alarm.bif"); 
  	string temp;
  	string name;
  	vector<string> values;
	vector<vector<float>> final_prob;

  	
    if (myfile.is_open())
    {
    	while (!myfile.eof())
    	{
    		stringstream ss;
      		getline (myfile,line);
      		
      		
      		ss.str(line);
     		ss>>temp;
     		
     		if(temp.compare("probability")==0)
     		{
                    
     				ss>>temp;
     				ss>>temp;
     				
                    // list<Graph_Node>::iterator listIt;
                    // list<Graph_Node>::iterator listIt1;
     				// listIt=Alarm.search_node(temp);
                    // int index=Alarm.get_index(temp);
                    ss>>temp;
                    values.clear();
     				while(temp.compare(")")!=0)
     				{
                        // listIt1=Alarm.search_node(temp);
                        // listIt1->add_child(index);
     					// values.push_back(temp);
     					
     					ss>>temp;

    				}
                    // listIt->set_Parents(values);
    				getline (myfile,line);
     				stringstream ss2;
                    
     				ss2.str(line);
     				ss2>> temp;
                    
     				ss2>> temp;
                    
     				vector<float> curr_CPT;
                    string::size_type sz;
     				while(temp.compare(";")!=0)
     				{
                        
     					curr_CPT.push_back(atof(temp.c_str()));
     					
     					ss2>>temp;
    				}
                    final_prob.push_back(curr_CPT);

     		}
            else
            {
                
            }
    		
    	}
    	// if(find==1)
    	// myfile.close();
  	}else{
		cout << "Unable to open file";
	}
  	return final_prob;
}


float get_cpt(int col,Graph_Node &node,vector<int> coeff,network &Alarm){
	vector<string> parents = node.get_Parents();
	int temp=1;
	int row=0;
	// cout<<coeff.size()<<" "<<parents.size()<<endl;
	if(coeff.size()!=parents.size()){
		cout<<"Error in CPT size and Parents size\n";
		cout<<"Node: "<<node.get_name()<<endl;
		cout<<"Node_parents"<<endl;
		for(string x : parents){
			cout<<x<<" ";
		}
		cout<<endl;
		exit(-1);
	}
	for(int i=parents.size()-1;i>=0;i--){
		int values= Alarm.search_node(parents[i])->get_nvalues();
		// temp*=values;
		row+=coeff[i]*temp;
		temp*=values;
	}
	if(node.get_nvalues()*row+col>=node.get_CPT().size()){
		cout<<"Index out of range in get_cpt func!";
		exit(0);
	}
	if(node.get_CPT()[node.get_nvalues()*row+col]==0){
		cout<<"Node_Name: "<<node.get_name()<<endl;
		cout<<"cpt: ";
		for(auto x : node.get_CPT()){
			cout<<x<<" ";
		}	
		cout<<"\n";
		exit(0);
	}
	return node.get_CPT()[node.get_nvalues()*row+col];
}

int update_cpt_indv(int col, Graph_Node &node,vector<int> coeff,network &Alarm, float new_cpt){
	vector<string> parents = node.get_Parents();
	int temp=1;
	int row=0;
	for(int i=parents.size()-1;i>=0;i--){
		int values= Alarm.search_node(parents[i])->get_nvalues();
		// temp*=values;
		row+=(coeff[i]*temp);
		temp*=values;
	}
	if(node.get_nvalues()*row+col>=node.get_CPT().size()){
		cout<<"Index out of range in update_cpt_indv func!"<<endl;
		cout<<"Index:"<<" "<<node.get_nvalues()*row+col<<" Size: "<<node.get_CPT().size()<<endl;
		cout<<col<<" "<<row<<" "<<node.get_nvalues()<<endl;
		for(int x : coeff){
			cout<<x<<endl;
		}
		for(string x : parents){
			cout<<Alarm.search_node(x)->get_nvalues()<<endl;
		}
		exit(0);
	}else{
		node.get_CPT()[node.get_nvalues()*row+col]=new_cpt;
		// cout<<node.get_CPT()[node.get_nvalues()*row+col]<<endl;
		return node.get_nvalues()*row+col;
	}
	// node.get_CPT()[node.get_nvalues()*row+col]=new_cpt;
}

pair<int,int> Numofsample(vector<string> &values,vector<vector<string>> &Datafile,string temp,vector<pair<int,string>> &computed_val,int col,vector<int> indxes){
	int count=0;
	int count2=0;
	int row=0;
	// cout<<temp<<endl;
	if(values.size()==0){
		for(auto &umap : Datafile){
			if(umap[col]==temp){
				count2++;
			}else{
				if(computed_val[row].first==col && computed_val[row].second == temp){
					count2++;
				}
			}
			count++;
			row++;
		}
		return {count2,count};
	}
	for(auto &umap : Datafile){
		bool search=true;
		int indx=0;
		for(string x : values){
			if(!(umap[indxes[indx]]==x || (computed_val[row].first== indxes[indx] && (computed_val[row].second==x)))){
				search=false;
				break;
			}
			indx++;
		}
		if(search){
			count++;
			if(umap[col]==temp){
				count2++;
			}
			else{
				if(computed_val[row].first==col && computed_val[row].second == temp){
					count2++;
				}
			}
		}
		row++; 
	}
	return {count2,count};
}


struct Entry {
    string value;
    int index;
};

void generatePermutationsWithIndices(vector<vector<string>>& temp, vector<vector<Entry>>& result, vector<Entry> current, int depth) {
    if (depth == temp.size()) {
        result.push_back(current);
        return;
    }
    for (int i = 0; i < temp[depth].size(); i++) {
        Entry entry;
        entry.value = temp[depth][i];
        entry.index = i;
        current.push_back(entry);
        generatePermutationsWithIndices(temp, result, current, depth + 1);
        current.pop_back();
    }
}

void update_cpt(network &Alarm,vector<vector<string>> &Datafile,vector<pair<int,string>> &computed_val){
	auto &nodes = Alarm.getPresGraph();
	int col_of_node=0;
	for (auto node = nodes.begin(); node!= nodes.end(); ++node) {
		auto parents = node->get_Parents();
		// int col=0;
		vector<vector<Entry>> result;
		vector<vector<string>> temp;
		vector<int> indexs;
		for(auto parent : parents){
			auto a = Alarm.search_node(parent)->get_values();
			int indx = Alarm.get_index(parent);
			indexs.push_back(indx);
			vector<string> tmp;
			for(string x : a){
				tmp.push_back(x);
			}
			temp.push_back(tmp);
		}
		generatePermutationsWithIndices(temp, result, vector<Entry>(), 0);
		int col=0;
		vector<float> updated_cpt(node->get_CPT().size());
		for(string x : node->get_values()){
			for(const auto& dataset : result){
				vector<string> values;
				vector<int> coeff;
				for(const Entry& entry : dataset){
					// cout<<entry.index<<" "<<entry.value<<endl;
					values.push_back(entry.value);
					coeff.push_back(entry.index);
				}
				auto pare = Numofsample(values,Datafile,x,computed_val,col_of_node,indexs);
				float prob;
				// if(pare.first==0){
				// 	prob=0;
				// 	// cout<<"Node Name :"<<node->get_name()<<" Value: "<<x<<endl;
				// 	// cout<<"Node col: "<<col_of_node<<endl;
				// 	// cout<<"Parents Value: "<<endl;
				// 	// int i=0;
				// 	// for(string x : values){
				// 	// 	cout<<x<<" "<<Alarm.get_index(Alarm.get_nth_node(indexs[i])->get_name())<<endl;
				// 	// 	i++;
				// 	// }
				// 	// cout<<pare.second<<endl;
				// 	// exit(0);
				// }else{
					// if(pare.second==0){
						// cout<<"Node Name :"<<node->get_name()<<"Value: "<<x<<endl;
						// cout<<"Parents Value: ";
						// cout<<"Node col: "<<col_of_node<<endl;
						// for(string x : values){
						// 	cout<<x<<" ";
						// }
					// 	exit(0);
					// }
				prob = static_cast<float>(pare.first+0.1) / static_cast<float>(pare.second+0.1*node->get_nvalues());
				// }
				// cout<<pare.first<<" "<<pare.second<<" "<<prob<<endl;
				// float rounded_prob = round(prob * 1e6) / 1e6;
				int indx=update_cpt_indv(col,*node,coeff,Alarm,prob);
				if(indx>=updated_cpt.size()){
					cout<<"Index out of range : "<<node->get_name();
					exit(-1);
				}
				updated_cpt[indx]=prob;
				// if(pare.first==pare.second){
			}
			col++;
		}
		node->set_CPT(updated_cpt);
		col_of_node++;
	}
}


// void evaluate_q_mark(vector<string> &mp1,network &Alarm, int row, vector<pair<int,string>> &computed_val){
// 	auto const &node = Alarm.get_nth_node(computed_val[row].first);
// 	float maxi=-1;
// 	string best_val;
// 	int col=0;
// 	for(string val : node->get_values()){
// 		vector<int> coeff;
// 		for(auto &parent : node->get_Parents()){
// 			auto tmp = Alarm.search_node(parent)->get_values();
// 			int parent_indx = Alarm.get_index(parent);
// 			int cof=0;
// 			for(string x : tmp){
// 				if(mp1[parent_indx]==x){
// 					coeff.push_back(cof);
// 					break;
// 				}
// 				cof++;
// 			}
// 		}
// 		float prob = get_cpt(col,*node,coeff,Alarm);
// 		col++;
// 		for(auto &parent : node->get_Parents()){
// 			auto tmp = Alarm.search_node(parent);
// 			int indx = Alarm.get_index(parent);
// 			int temp_col=0;
// 			for(string x : tmp->get_values()){
// 				if(x==mp1[indx]){
// 					break;
// 				}
// 				temp_col++;
// 			}
// 			vector<int> temp_coeff;
// 			for(auto &paren : tmp->get_Parents()){
// 				auto paren_tmp = Alarm.search_node(paren);
// 				int paren_indx = Alarm.get_index(paren);
// 				int cof=0;
// 				for(string y : paren_tmp->get_values()){
// 					if(mp1[paren_indx]==y){
// 						temp_coeff.push_back(cof);
// 						break;
// 					}
// 					cof++;
// 				}   
// 			}
// 			prob*=get_cpt(temp_col,*tmp,temp_coeff,Alarm);
// 		}
// 		if(prob>=maxi){
// 			maxi=prob;
// 			best_val=val;
// 		}
// 		if(prob==0) break;
// 	}
// 	// cout<<"Node_name: "<<node->get_name()<<endl;
// 	// cout<<"Prev_Val: "<<computed_val[row].second<<" "<<"New_val: "<<best_val<<" Prob: "<<maxi<<endl;
// 	computed_val[row].second=best_val;
// 	// cout<<best_val<<" "<<maxi<<endl;
// 	// exit(0);
// }

void evaluate_q_mark(vector<string> &mp1,network &Alarm, int row, vector<pair<int,string>> &computed_val){
	auto const &node = Alarm.get_nth_node(computed_val[row].first);
	float maxi=-1;
	string best_val;
	int col=0;
	for(auto val : node->get_values()){
		vector<int> coeff;
		for(auto &parent : node->get_Parents()){
			auto tmp = Alarm.search_node(parent)->get_values();
			int parent_indx = Alarm.get_index(parent);
			int cof=0;
			for(string x : tmp){
				if(mp1[parent_indx]==x){
					coeff.push_back(cof);
					break;
				}
				cof++;
			}
		}
		float prob = static_cast<float>(get_cpt(col,*node,coeff,Alarm));
		// cout<<prob<<endl;
		// cout<<"Node_name: "<<node->get_name()<<endl;
		// cout<<"Node_val: "<<val<<endl;
		// cout<<"Node_cpt: ";
		// for(auto x : node->get_CPT()){
		// 	cout<<x<<" ";
		// }
		// cout<<endl;
		float temp=1;
		for(int children : node->get_children()){
			auto child = Alarm.get_nth_node(children);
			int temp_col=0;
			for(string x : child->get_values()){
				if(x==mp1[children]){
					break;
				}
				temp_col++;
			}
			vector<int> coeff_child;
			for(auto &child_parent : child->get_Parents()){
				auto tmp = Alarm.search_node(child_parent)->get_values();
				int child_parent_indx = Alarm.get_index(child_parent);
				int cof2=0;
				if(child_parent==node->get_name()){
					coeff_child.push_back(col);
				}else{
					for(string x : tmp){
						if(mp1[child_parent_indx]==x){
							coeff_child.push_back(cof2);
							break;
						}
					cof2++;
					}
				}
				
			}
			temp *= static_cast<float>(get_cpt(temp_col,*child,coeff_child,Alarm));
			if(static_cast<float>(get_cpt(temp_col,*child,coeff_child,Alarm))==0 || temp==0){
				cout<<"Node: "<<node->get_name()<<endl;
				cout<<"Current_cpt: "<<get_cpt(temp_col,*child,coeff_child,Alarm)<<endl;
				cout<<"Cpt: ";
				for(auto x : node->get_CPT()) {
					cout << x << ' ';
				}
				cout<<endl;
				exit(0);
			}
		}
		prob*=temp;
		if(prob>=maxi){
			maxi=prob;
			best_val=val;
		}
		col++;
		// cout<<val<<" Prob: "<<prob<<endl;
		// if(prob==0)
	}
	// cout<<"Node_Name: "<<node->get_name()<<endl;
	// cout<<"Old_val: "<<computed_val[row].second<<" New_val: "<<best_val<<" Best_prob: "<<maxi<<endl;
	computed_val[row].second=best_val;
}

void compute_value_quesMark(int temp,vector<string> &mp1,network &Alarm,int row,vector<pair<int,string>> &computed_val){
	auto const &node = Alarm.get_nth_node(temp);
	vector<int> coeff;
	for(auto &parent : node->get_Parents()){
		auto tmp = Alarm.search_node(parent)->get_values();
		int parent_indx = Alarm.get_index(parent);
		int cof=0;
		for(string x : tmp){
			if(mp1[parent_indx]==x){
				coeff.push_back(cof);
				break;
			}
			cof++;
		}
	}
	int maxi=-1;
	int col=0;
	string value;
	// cout<<"Node_Name: "<<node->get_name()<<endl;
	// cout<<"Last value:"<<computed_val[row].second<<endl;
	// cout<<"Possible_values: "<<endl;
	for(auto val : node->get_values()){
		int cpt = get_cpt(col,*node,coeff,Alarm);
		if(cpt>maxi){
			value=val;
		}
		cout<<val<<" "<<cpt<<endl;
		col++;
		maxi=max(maxi,cpt);
	}
	// cout<<"Best Value: "<<value<<endl;
	if(row>=computed_val.size()){
		cout<<"Row Index Out Of Range in compute_value_qmark fun"<<endl;
		exit(-1);
	}
	computed_val[row].first=temp;
	computed_val[row].second=value;
	// cout<<computed_val[row].second<<endl;
	
}

pair<vector<pair<int,string>>,vector<vector<string>>> ReadDataset(network &Alarm){
	vector<pair<int,string>> computed_val;
	vector<vector<string>> Data;
	ifstream myfile("records.dat");
	if (!myfile.is_open()) {
		cout << "Unable to open file";
	}else{
		string line;
		while (getline(myfile, line)) {
			istringstream iss(line); // Create a string stream for the current line
			string word;
			int index=0;
			int temp=-1;
			// unordered_map<string,bool> mp1;
			vector<string> tmp;
			// unordered_map<int,string> umap;
			while (iss >> word) {
				// string x = "?";
				if(word.size()==3){
					// cout<<word<<endl;
					temp=index;
					// exit(0);
				}
				// umap[index]=word;
				// cout<<word<<" ";
				index++;
				// mp1[word]=true;
				// exit(0);
				tmp.push_back(word);
			}
			// cout<<endl;
			// exit(0);
			if(temp!=-1){
				// auto const &node = Alarm.get_nth_node(temp);
				// vector<int> coeff;
				// for(auto &parent : node->get_Parents()){
				// 	auto tmp1 = Alarm.search_node(parent)->get_values();
				// 	int parent_indx = Alarm.get_index(parent);
				// 	int cof=0;
				// 	// exit(0);
				// 	for(string x : tmp1){
				// 		if(tmp[parent_indx]==x){
				// 			coeff.push_back(cof);
				// 			break;
				// 		}
				// 		cof++;
				// 	}
				// 	// exit(0);
				// }
				// // exit(0);
				// if(coeff.size()!=node->get_Parents().size()){
				// 	cout<<"Error in CPT size and Parents size in ReadDatabase func..\n";
				// 	cout<<coeff.size()<<" "<<node->get_Parents().size()<<endl;
				// 	cout<<"Node: "<<node->get_name()<<endl;
				// 	cout<<"Node_parents"<<endl;
				// 	for(string x : node->get_Parents()){
				// 		int parent_indx = Alarm.get_index(x);
				// 		auto tmp1 = Alarm.search_node(x)->get_values();
				// 		// exit(0);
				// 		// cout<<tmp1.size();
				// 		for(string x : tmp1){
				// 			cout<<x<<endl;
				// 			if(parent_indx>=tmp.size()) cout<<"Error";
				// 			exit(0);
				// 			if(tmp[parent_indx]==x){
				// 				cout<<x<<" "<<endl;
				// 			}else{
				// 				cout<<"Original: "<<x<<" "<<"Required "<<tmp1[parent_indx]<<endl;
				// 			}
				// 			exit(0);
				// 		}
						
				// 	}
				// 	// cout<<endl;
				// 	//
				// 	exit(0);
				// }
				// // exit(0);
				// // cout<<coeff.size()<<" "<<node->get_Parents().size()<<endl;
				// int maxi=-1;
				// int col=0;
				// string value;
				// for(auto &val : node->get_values()){
				// 	// exit(0);
				// 	int cpt = get_cpt(col,*node,coeff,Alarm);
				// 	// exit(0);
				// 	if(cpt>=maxi){
				// 		value=val;
				// 	}
				// 	col++;
				// 	maxi=max(maxi,cpt);
				// }
				auto const &node = Alarm.get_nth_node(temp);
				computed_val.push_back({temp,node->get_values()[0]});
			}else{
				computed_val.push_back({-1,"null"});
			}
			Data.push_back(tmp);
		}
		

	}
	myfile.close();

	return {computed_val,Data};
}


void write_network_to_bif(network Alarm)
{
    ofstream outfile("Solved_Alarm.bif");
    if (!outfile.is_open())
    {
        cout << "Unable to open the output file";
        return;
    }

    outfile << "// Bayesian Network in the Interchange Format" << endl;

    auto nodes = Alarm.getPresGraph();
    for (auto node = nodes.begin(); node != nodes.end(); ++node)
    {
        outfile << "variable \"" << node->get_name() << "\" {" << endl;
        outfile << "\ttype discrete[" << node->get_nvalues() << "] { ";
        for (int i = 0; i < node->get_nvalues(); ++i)
        {
            outfile << "\"" << node->get_values()[i] << "\"";
            if (i < node->get_nvalues() - 1)
            {
                outfile << " ";
            }
        }
        outfile << " };" << endl;
        outfile << "\tproperty \"position = (x, y)\" ;" << endl;
        outfile << "}" << endl;
    }

    for (auto node = nodes.begin(); node != nodes.end(); ++node)
    {
        outfile << "probability (";
        outfile << " \"" << node->get_name() << "\"";
        auto parents = node->get_Parents();
        for (auto it = parents.begin(); it != parents.end(); ++it)
        {
            outfile << " \"" << *it << "\"";
        }
        outfile << " ) {" << endl;
        outfile << "\ttable ";
        for (int i = 0; i < node->get_CPT().size(); ++i)
        {
            outfile << node->get_CPT()[i];
            if (i < node->get_CPT().size() - 1)
            {
                outfile << " ";
            }
        }
        outfile << ";" << endl;
        outfile << "}" << endl;
    }
    outfile.close();
}


int main()
{	
	chrono::seconds runTime(120);
    auto startTime = chrono::high_resolution_clock::now();
	network Alarm;
	Alarm=read_network();
	auto pare = ReadDataset(Alarm);
	vector<pair<int,string>> computed_val;
	vector<vector<string>> mps;
	computed_val=pare.first;
	mps=pare.second;
	int itr=0;
	vector<vector<float>> nex = finaL_prob();
	while(itr<3){
		auto currentTime = chrono::high_resolution_clock::now();
		auto elapsedTime = chrono::duration_cast<chrono::seconds>(currentTime - startTime);
		if (elapsedTime >= runTime) {
			break;
		}
		int row=0;
		for(auto &mp : mps){
			int col=computed_val[row].first;
			if(col!=-1){
				evaluate_q_mark(mp,Alarm,row,computed_val);
			}
			row++;
		}
		update_cpt(Alarm,mps,computed_val);
		itr++;
	}
	// write_network_to_bif(Alarm);


	vector<vector<float>> final_prob;
	 for(int i = 0; i < nex.size();i++){
		auto x = nex[i];
		vector<vector<float>> new1;
		int mul = 1;
		// cout<<i<<endl;
		// cout<<node_at_index[i]<<endl;
		vector<string> parents = Alarm.get_nth_node(i)->get_Parents();		
		if (parents.size()==0){final_prob.push_back(x);continue;}
		for(auto y : parents){
			mul *= Alarm.search_node(y)->get_nvalues();
		}
		// cout<<mul<<endl;
		// cout<<"new1: "<<endl;
		for(int r = 0;r < x.size();r+=mul){
			// cout<<"temp: "<<endl;
			vector<float> temp;
			for(int l = r;l < r + mul;l++){temp.push_back(x[l]);}
			// for(auto z : temp){cout<<z<<" ";}
			new1.push_back(temp);
		}
		// cout<<endl;
		// cout<<"hi2";
		vector<float> temp2;
		for(int r = 0;r < mul;r++){
			for(int l = 0;l < new1.size();l++){
				temp2.push_back(new1[l][r]);
			}
		}
		final_prob.push_back(temp2);
	}
	auto &nodes = Alarm.getPresGraph();
	float score=0;
	int i=0;
	for(auto node=nodes.begin();node!=nodes.end();node++){
		// cout<<"Node Name: "<<node->get_name()<<endl;
		// cout<<"Cpt: ";
		int j=0;
		for(auto x : node->get_CPT()){
			// cout<<x<<" ";
			score+=abs(x-final_prob[i][j]);
			j++;
		}
		// cout<<endl;
		i++;
	}
	cout<<"Total Iterations: "<<itr<<endl;
	cout<<"Score: "<<score<<endl;
	auto currentTime = chrono::high_resolution_clock::now();
	auto elapsedTime = chrono::duration_cast<chrono::seconds>(currentTime - startTime);
	cout<<"Time: "<<elapsedTime.count()<<endl;

}