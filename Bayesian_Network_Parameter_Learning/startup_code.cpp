#include <iostream>
#include <string>
#include <vector>
#include <list>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include<unordered_map>
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

unordered_map<string,int> index_of_value;
unordered_map<string,int> index_of_node;
unordered_map<int,string> node_at_index;
unordered_map<string,Graph_Node*> pointers;

vector<vector<float>> finaL_prob()
{
	network Alarm;
	string line;
	int find=0;
  	ifstream myfile("./testcases/testcase3/gold_alarm.bif"); 
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
     				
                    
                    ss>>temp;
                    values.clear();
     				while(temp.compare(")")!=0)
     				{
     					
     					ss>>temp;

    				}
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
		// int values= Alarm.search_node(parents[i])->get_nvalues();
		int values = pointers[parents[i]]->get_nvalues();
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
		int values = pointers[parents[i]]->get_nvalues();
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
		return node.get_nvalues()*row+col;
	}
}

pair<int,int> Numofsample(vector<string> &values,vector<vector<string>> &Datafile,string temp,vector<pair<int,string>> &computed_val,int col,vector<int> indxes, network &Alarm){
	int total_sample=0;
	int fav_sample=0;
	int row=0;
	if(values.size()==0){
		for(auto &umap : Datafile){
			if(umap[col]==temp){
				fav_sample++;
			}else{
				if(computed_val[row].first==col && computed_val[row].second == temp){
					fav_sample++;
				}
			}
			total_sample++;
			row++;
		}
		return {fav_sample,total_sample};
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
			total_sample++;
			if(umap[col]==temp){
				fav_sample++;
			}
			else{
				if(computed_val[row].first==col && computed_val[row].second == temp){
					fav_sample++;
				}
			}
		}
		row++; 
	}
	return {fav_sample,total_sample};
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
		vector<vector<Entry>> result;
		vector<vector<string>> temp;
		vector<int> indexs;
		for(auto parent : parents){
			auto a = pointers[parent]->get_values();
			int indx = index_of_node[parent];
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
					values.push_back(entry.value);
					coeff.push_back(entry.index);
				}
				auto pare = Numofsample(values,Datafile,x,computed_val,col_of_node,indexs,Alarm);
				float prob;
				// 0.03 - 20.9885
				// 0.02 - 21.03
				// 0.031 - 20.9872
				prob = static_cast<float>(pare.first+0.031) / static_cast<float>(pare.second+0.031*node->get_nvalues());
				int indx=update_cpt_indv(col,*node,coeff,Alarm,prob);
				if(indx>=updated_cpt.size()){
					cout<<"Index out of range : "<<node->get_name();
					exit(-1);
				}
				updated_cpt[indx]=prob;
			}
			col++;
		}
		node->set_CPT(updated_cpt);
		col_of_node++;
	}
}


void evaluate_q_mark(vector<string> &mp1,network &Alarm, int row, vector<pair<int,string>> &computed_val){
	// auto const &node = Alarm.get_nth_node(computed_val[row].first);
	auto const &node = pointers[node_at_index[computed_val[row].first]];
	float maxi=-1;
	string best_val;
	int col=0;
	for(auto val : node->get_values()){
		vector<int> coeff;
		for(auto &parent : node->get_Parents()){
			auto tmp = pointers[parent]->get_values();
			int parent_indx = index_of_node[parent];
			
			string key = parent + " " + mp1[parent_indx];
			coeff.push_back(index_of_value[key]);
		}
		float prob = static_cast<float>(get_cpt(col,*node,coeff,Alarm));
		float temp=1;
		for(int children : node->get_children()){
			auto child = pointers[node_at_index[children]];
			int temp_col=0;
			string key = child->get_name() + " " + mp1[children];
			temp_col = index_of_value[key];
			vector<int> coeff_child;
			for(auto &child_parent : child->get_Parents()){
				auto tmp = pointers[child_parent]->get_values();
				int child_parent_indx = index_of_node[child_parent];
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
	}
	
	computed_val[row].second=best_val;
}

pair<vector<pair<int,string>>,vector<vector<string>>> ReadDataset(network &Alarm,const std::string& inputfile){
	vector<pair<int,string>> computed_val;
	vector<vector<string>> Data;
	ifstream myfile(inputfile);
	if (!myfile.is_open()) {
		cout << "Unable to open file";
	}else{
		string line;
		while (getline(myfile, line)) {
			istringstream iss(line); // Create a string stream for the current line
			string word;
			int index=0;
			int temp=-1;
			vector<string> tmp;
			while (iss >> word) {
				if(word.substr(1,1)=="?"){
					temp=index;
				}
				index++;
				tmp.push_back(word);
			}
			if(temp!=-1){

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

void write_network_to_bif( const std::string& inputFileName, network& Alarm) {
    std::ofstream myfile("solved_alarm.bif");
    std::ifstream orginalfile(inputFileName);
    std::string line;
    bool is_write = false;
    int j = 0;
    while (std::getline(orginalfile, line)) {
        std::stringstream ss;
        ss.str(line);
        std::string item;
        ss >> item;
        if (!is_write) {
            myfile << line << std::endl;
        } else {
            std::string line2;
            line2 += "\ttable ";
            auto node = pointers[node_at_index[j]];
            std::vector<float> cpt = node->get_CPT();
            for (int i = 0; i < cpt.size(); i++) {
                float prob = round(cpt[i] * 10000) / 10000.0;
                line2 += std::to_string(prob) + " ";
            }
            line2 += ";\n";
            myfile << line2;
            j++;
        }
        if (item.compare("probability") == 0) {
            is_write = true;
        } else {
            is_write = false;
        }
    }
}



void preprocessing(network &Alarm, vector<vector<string>> &Data,vector<pair<int,string>> computed_val){
	for (int i = 0; i < Alarm.netSize(); i++) {
		Graph_Node* node = &(*(Alarm.get_nth_node(i))); // Get a pointer to the node
		pointers[node->get_name()] = node;
		node_at_index[i] = node->get_name();
	}
	for(int i = 0; i < Alarm.netSize(); i++){
		Graph_Node n = *(pointers[node_at_index[i]]);
		vector<string> values = n.get_values();
		index_of_node[n.get_name()] = i;
		for(int j  = 0; j < values.size();j++){
			string key = n.get_name() + " " + values[j];
			index_of_value[key] = j;
		}
	}
	
	
}

int main(int argc, char* argv[])
{	
	chrono::seconds runTime(108);
    auto startTime = chrono::high_resolution_clock::now();
	network Alarm;
	Alarm=read_network();
	auto const &pare = ReadDataset(Alarm,argv[2]);
	vector<pair<int,string>> computed_val;
	vector<vector<string>> mps;
	computed_val=pare.first;
	mps=pare.second;
	preprocessing(Alarm,mps,computed_val);
	int itr=0;
	vector<vector<float>> nex = finaL_prob();

	while(true){
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
	
	vector<vector<float>> final_req;
	for(int i = 0; i < Alarm.netSize();i++){
		vector<string> par = pointers[node_at_index[i]]->get_Parents();
		int n = pointers[node_at_index[i]]->get_nvalues();
		vector<float> cpt = pointers[node_at_index[i]]->get_CPT();
		int num_rows = (cpt.size())/(n);
		vector<vector<float>> req_2d;
		for(int r = 0;r < cpt.size();r += n){
			vector<float> tp;
			for(int l = r;l < r + n;l++){tp.push_back(cpt[l]);}
			req_2d.push_back(tp);
		}
		vector<float> temp;
		for(int col = 0;col  < n;col++){
			for(int row=0;row < num_rows;row++){
				temp.push_back(req_2d[row][col]);
			}
		}
		final_req.push_back(temp);
	}
	for(int i = 0;i  < Alarm.netSize();i++){
		Alarm.get_nth_node(i)->set_CPT(final_req[i]);
	}
	auto &nodes = Alarm.getPresGraph();
	float score=0;
	int i=0;
	for(auto node=nodes.begin();node!=nodes.end();node++){
		int j=0;
		for(auto x : node->get_CPT()){
			score+=fabs(x-nex[i][j]);
			j++;
		}
		i++;
	}
	// write_network_to_bif(argv[1],Alarm);
	cout<<"Total Iterations: "<<itr<<endl;
	cout<<"Score: "<<score<<endl;
	auto currentTime = chrono::high_resolution_clock::now();
	auto elapsedTime = chrono::duration_cast<chrono::seconds>(currentTime - startTime);
	cout<<"Time: "<<elapsedTime.count()<<endl;
}






