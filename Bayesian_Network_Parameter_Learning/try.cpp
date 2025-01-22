#include <iostream>
#include <string>
#include <vector>
#include <list>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include<unordered_map>
#include <thread>
#include <random>
#include <iomanip>
using namespace std;
// int number_of_different_q_marks = 0;

// Format checker just assumes you have Alarm.bif and Solved_Alarm.bif (your file) in current directory

float roundToDecimalPlaces(float value, int decimalPlaces) {
    // double multiplier = std::pow(10.0, decimalPlaces);
    return round(value * 10000) / 10000.0;
}
// Our graph consists of a list of nodes where each node is represented as follows:
class Graph_Node{

private:
	string Node_Name;  // Variable name 
	vector<int> Children; // Children of a particular node - these are index of nodes in graph.
	vector<string> Parents; // Parents of a particular node- note these are names of parents
	int nvalues;  // Number of categories a variable represented by this node can take
	vector<string> values; // Categories of possible values
	vector<float> CPT; // conditional probability table as a 1-d array . Look for BIF format to understand its meaning4
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
unordered_map<string,Graph_Node*> pointers;
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
    	while (! myfile.eof() )
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

int find_index(string val,string par,network& Alarm,unordered_map<string,int> &index_of_value, unordered_map<string,int> &index_of_node, unordered_map<int,string> &node_at_index){
	string key = par + " " + val;
	if (index_of_value.find(key) == index_of_value.end()){cout<<"ERROR"<<endl;return -1;}
	return index_of_value[key];

}

int find_row(Graph_Node node,vector<string>& values,network& Alarm,unordered_map<string,int> &index_of_value, unordered_map<string,int> &index_of_node, unordered_map<int,string> &node_at_index){
	int mul = 1;
	int n = values.size();
	vector<string> par = node.get_Parents();
	int ans = 0;
	for(int i = 0; i < n; i++){
		string val = values[n-i-1];
		ans += find_index(val,par[n-i-1],Alarm,index_of_value, index_of_node, node_at_index)*mul;
		int idx = index_of_node[par[n-i-1]];
		mul *= (pointers[node_at_index[idx]])->get_nvalues();
	}
	return ans;
}

int find_main_idx(Graph_Node node,string value_node,vector<string> values,network& Alarm,unordered_map<string,int> &index_of_value, unordered_map<string,int> &index_of_node, unordered_map<int,string> &node_at_index){
	int j = find_index(value_node,node.get_name(),Alarm,index_of_value, index_of_node, node_at_index);
	return ((node.get_nvalues())*find_row(node,values,Alarm,index_of_value, index_of_node, node_at_index)  + j);
}

float find_CPT_probability(network& Alarm,string node_name,string value,vector<string>& parent_values,unordered_map<string,int> &index_of_value, unordered_map<string,int> &index_of_node, unordered_map<int,string> &node_at_index){
	int id = find_main_idx(*pointers[node_name],value,parent_values,Alarm,index_of_value, index_of_node, node_at_index);
	return roundToDecimalPlaces(pointers[node_name]->get_CPT()[id],4);
}

float question_mark_helper(network& Alarm,vector<string> &vec,unordered_map<string,int> &index_of_value, unordered_map<string,int> &index_of_node, unordered_map<int,string> &node_at_index,string node_name){
	float ans = 1;
	int sz = Alarm.netSize();
	vector<int> children = pointers[node_name]->get_children();
	children.push_back(index_of_node[node_name]);
	for(int i = 0; i < children.size();i ++){
		vector<string> values;
		Graph_Node req = *(pointers[node_at_index[children[i]]]);
		vector<string> par = req.get_Parents();
		if (par.size() == 0){
			string name = req.get_name();
			string value_of_node = vec[index_of_node[name]];
			int idx = find_index(value_of_node,name,Alarm,index_of_value, index_of_node, node_at_index);
			vector<float> CPT = req.get_CPT();
			ans *= CPT[idx];
			continue;
		}
		for(int j = 0; j < par.size();j++){
			values.push_back(vec[index_of_node[par[j]]]);
		}
		ans *= find_CPT_probability(Alarm,req.get_name(),vec[index_of_node[req.get_name()]],values,index_of_value, index_of_node, node_at_index);
		ans = roundToDecimalPlaces(ans,4);
	}
	ans = roundToDecimalPlaces(ans,4);
	return ans;
}

string return_best_val(string node_name,vector<string> &vec,network& Alarm,unordered_map<string,int> &index_of_value, unordered_map<string,int> &index_of_node, unordered_map<int,string> &node_at_index){
	Graph_Node node = *(pointers[node_name]);
	vector<string> vals = node.get_values();
	string best = "";
	float probability = -1;
	for(int i = 0; i < vals.size();i++){
		vec[index_of_node[node_name]] = vals[i];
		float temp = question_mark_helper(Alarm,vec,index_of_value, index_of_node, node_at_index,node_name);
		if (temp > probability){
			probability = temp;
			best = vals[i];
		}
	}
	return best;
}
unordered_map<int,pair<pair<bool,int>,vector<string>>> ReadData(network& Alarm,const string file){
	ifstream myfile(file);
		string line;
		int line_n = 0;
		unordered_map<int,pair<pair<bool,int>,vector<string>>> mapper;
		while (getline(myfile, line)) {
			istringstream iss(line); // Create a string stream for the current line
			string word;
			vector<string> lst;
			bool is_q_mark = false;
			int i = 0;
			int id_req = -1;
			while (iss >> word) {
				int len = word.size();
				lst.push_back(word);
				word = word.substr(1,len-2);
				if (word == "?"){
                    // cout<<"in"<<endl;
                    is_q_mark = true;id_req = i;
				
				}
				
				i++;
			}
			// cout<<endl;
			mapper[line_n] = {{is_q_mark,id_req},lst}; 

			line_n++;
		}
		myfile.close();
		return mapper;
	// }
	
}

unordered_map<int,string> find_best_values(unordered_map<int,pair<pair<bool,int>,vector<string>>> &dataset,network& Alarm,unordered_map<string,int> &index_of_value, unordered_map<string,int> &index_of_node, unordered_map<int,string> &node_at_index){
	int sz = dataset.size();
	unordered_map<int,string> ans;
	for(int i = 0; i < sz;i++){
		if ((dataset[i].first).first == true){
			Graph_Node node = *(pointers[node_at_index[(dataset[i].first).second]]);
			string node_name = node.get_name();
			ans[i] = return_best_val(node_name,(dataset[i].second),Alarm,index_of_value, index_of_node, node_at_index);
		}
		else{ans[i]="";}
		
	}

	return ans;
}

float find_prob_from_data(network& Alarm,Graph_Node node,string node_val,vector<string> &parents_vals,unordered_map<int,pair<pair<bool,int>,vector<string>>> &dataset,unordered_map<int,string> &ans,unordered_map<string,int> &index_of_value, unordered_map<string,int> &index_of_node, unordered_map<int,string> &node_at_index){
	vector<string> parents = node.get_Parents();
	int idx_node = index_of_node[node.get_name()];
	vector<int> indices = {};
	for(int i = 0;i < parents.size();i++){indices.push_back(index_of_node[parents[i]]);}
	int sz  = dataset.size();
	float total_count = 0.02;
	float favourable_count = 0.02;
	for(int i = 0; i < sz;i++){
		vector<string> req = (dataset[i].second);
		bool to_include = true;
		for(int j = 0;j  < indices.size();j++){
			int id = indices[j];
			if (req[id] == parents_vals[j] || (req[id] == "?" && ans[i] == parents_vals[j])){continue;}
			else{to_include=false;break;}
		}
		if (to_include){total_count++;
			if (req[idx_node] == node_val || (req[idx_node] == "?" && ans[i] == node_val)){favourable_count++;}
		}
	}
	if (total_count==0){return 0;}
	return roundToDecimalPlaces(float(float(favourable_count)/float(total_count)),4);
}	

vector<int> find_coeff(network& Alarm,vector<string> &Parents,int row_idx){
	int maxi  = 1;
	vector<int> req;
	for(int i = 1; i < Parents.size();i++){
		Graph_Node node = *(pointers[Parents[i]]);
		maxi *= node.get_nvalues();
		req.push_back(node.get_nvalues());
	}
	vector<int> coeff;
	int i = 0;
	req.push_back(1);
	while(coeff.size() != Parents.size()){
		int cf = row_idx/maxi;
		coeff.push_back(cf);
		row_idx %= maxi;
		maxi /= req[i];
		i++;
	}
	return coeff;
}

void update_cpt(string node_name,network& Alarm,unordered_map<int,pair<pair<bool,int>,vector<string>>> &dataset,unordered_map<int,string> &ans,unordered_map<string,int> &index_of_value, unordered_map<string,int> &index_of_node, unordered_map<int,string> &node_at_index){
	Graph_Node node = *(pointers[node_name]);
	vector<string> values = node.get_values();
	vector<string> parents = node.get_Parents();
	vector<float> CPT = node.get_CPT();
	for(int index = 0;index < CPT.size();index++){
		int n = node.get_nvalues();
		int j = (index)%n;
		int row_idx = (index - j)/n;
		vector<int> coeff = find_coeff(Alarm,parents,row_idx);
		vector<string> req;
		for(int i = 0; i < coeff.size();i++){
			req.push_back((pointers[parents[i]]->get_values())[coeff[i]]);
		}
		CPT[index] = find_prob_from_data(Alarm,node,values[j],req,dataset,ans,index_of_value, index_of_node, node_at_index);
	}
	(pointers[node_name])->set_CPT(CPT);
}

void update(network& Alarm,unordered_map<int,pair<pair<bool,int>,vector<string>>> &dataset,unordered_map<int,string> &ans,unordered_map<string,int> &index_of_value, unordered_map<string,int> &index_of_node, unordered_map<int,string> &node_at_index){
	int sz = Alarm.netSize();
	for(int i = 0; i < sz; i++){
		update_cpt(node_at_index[i],Alarm,dataset,ans,index_of_value, index_of_node, node_at_index);
	}
}

float compute_init_probab(string node_name,string value_of_node,vector<int> coeff,unordered_map<int,pair<pair<bool,int>,vector<string>>> &dataset,unordered_map<string,int> &index_of_value, unordered_map<string,int> &index_of_node, unordered_map<int,string> &node_at_index,vector<string> &best_poss){
	vector<string> par = pointers[node_name]->get_Parents();
	vector<string> values_par;
	for(int i = 0; i < par.size();i++){
		values_par.push_back(pointers[par[i]]->get_values()[coeff[i]]);
	}
	float total_count = 0.01;
	float favourable_count = 0.01;
	for(int i = 0; i < dataset.size();i++){
		bool include = true;
		vector<string> req = dataset[i].second;
		for(int j = 0; j < par.size();j++){
			int id = index_of_node[par[j]];
			if (req[id] == values_par[j] || req[id].substr(1,req[id].size()-2) == "?" && values_par[j] == best_poss[id]){continue;}
			else{include=false;break;}
		}
		if (include){
			total_count++;
			if (req[index_of_node[node_name]] == value_of_node || req[index_of_node[node_name]].substr(1,req[index_of_node[node_name]].size()-2) == "?" && value_of_node == best_poss[index_of_node[node_name]]){favourable_count++;}
		}
	}
	return float(favourable_count/total_count);
}
vector<string> find_init_best_values(unordered_map<int,pair<pair<bool,int>,vector<string>>> &dataset,network& Alarm,unordered_map<string,int> &index_of_value, unordered_map<string,int> &index_of_node, unordered_map<int,string> &node_at_index){
	vector<string> req;
	for(int i = 0; i < Alarm.netSize();i++){
		int n = pointers[node_at_index[i]]->get_nvalues();
		vector<int> vec(n,0);
		for(int j = 0;j  < dataset.size();j++){
			string val = dataset[i].second[i];
			if (val.substr(1,val.size()-2) != "?"){
				string key = node_at_index[i] + " " + val;
				vec[index_of_value[key]]++;
			}
		}
		int ct = -1;
		int id = -1;
		for(int k = 0; k  < vec.size();k++){
			if (vec[k] > ct){id = k;ct = vec[k];}
		}
		req.push_back(pointers[node_at_index[i]]->get_values()[id]);
	}
	return req;
}
void initialise(network& Alarm,unordered_map<string,int> &index_of_value, unordered_map<string,int> &index_of_node, unordered_map<int,string> &node_at_index,unordered_map<int,pair<pair<bool,int>,vector<string>>> &dataset,vector<string> &best_pred){
	int sz = Alarm.netSize();
	random_device rd;  // Used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); // Mersenne Twister 19937 generator
    std::uniform_real_distribution<double> distribution(0.0, 1.0); //
	for(int i = 0;i < sz;i++){
		vector<float> CPT = pointers[node_at_index[i]]->get_CPT();
		for(int j = 0; j < CPT.size();j++){
			vector<string> par = pointers[node_at_index[i]]->get_Parents();
				CPT[j] = 0;}
		// }
		(pointers[node_at_index[i]])->set_CPT(CPT);
	}
}


string write_infile(Graph_Node &node){
    string line2;
    line2+="\ttable ";
    // auto node = pointers[node_at_index[j]];
    vector<float> cpt = node.get_CPT();
    for(int i=0;i<cpt.size();i++){
        float prob = round(cpt[i]*10000)/10000.0;
        // cout<<prob<<endl;
        line2+=to_string(prob)+ " ";
    }
    line2+=";\n";
    return line2;
}

void network_to_bif_format(network &Alarm,unordered_map<int,string> &node_at_index,const string input){
    ofstream myfile("solved_alarm.bif");
	ifstream orginalfile(input);
	string line;
	bool should_write=false;
	int j=0;
	while(getline(orginalfile,line)){
		stringstream ss;
		ss.str(line);
		string item;
		ss>>item;
		if(!should_write){
			myfile<<line<<endl;
		}else{
			auto node = pointers[node_at_index[j]];
			myfile << write_infile(*node);
			j++;
		}
		if(item.compare("probability")==0){
			should_write=true;
		}else{
			should_write=false;
		}
	}
    


}
int main(int argc, char* argv[])
{
    chrono::seconds runTime(97);
    auto startTime = chrono::high_resolution_clock::now();
	network Alarm;
	Alarm=read_network();
	unordered_map<string,int> index_of_value;
	unordered_map<string,int> index_of_node;
	unordered_map<int,string> node_at_index;
    // exit(0);
	unordered_map<int,pair<pair<bool,int>,vector<string>>> dataset = ReadData(Alarm,argv[2]);
    // exit(0);
	for (int i = 0; i < Alarm.netSize(); i++) {
        Graph_Node* node = &(*(Alarm.get_nth_node(i))); 
        pointers[node->get_name()] = node;
        node_at_index[i] = node->get_name();
	}
    // exit(0);
	vector<string> init_best = find_init_best_values(dataset,Alarm,index_of_value, index_of_node,node_at_index);
    // exit(0);
	initialise(Alarm,index_of_value, index_of_node, node_at_index,dataset,init_best);
	int sz = Alarm.netSize();
	
	for(int i = 0; i < sz; i++){
		Graph_Node n = *(pointers[node_at_index[i]]);
		vector<string> values = n.get_values();
		index_of_node[n.get_name()] = i;
		for(int j  = 0; j < values.size();j++){
			string key = n.get_name() + " " + values[j];
			index_of_value[key] = j;
		}
	}
	
    int itr=0;
	while(true){
		unordered_map<int,string> bests = {};
        auto currentTime = chrono::high_resolution_clock::now();
		auto elapsedTime = chrono::duration_cast<chrono::seconds>(currentTime - startTime);
		if (elapsedTime >= runTime) {
			break;
		}
		bests = find_best_values(dataset,Alarm,index_of_value, index_of_node, node_at_index);
		update(Alarm,dataset,bests,index_of_value, index_of_node, node_at_index);
        itr++;
	}

	vector<vector<float>> req;
	for(int i = 0; i < sz; i++){
		// cout<<i<<endl;
		vector<float> temp;
		Graph_Node n = *(pointers[node_at_index[i]]);
		vector<float> CPT = n.get_CPT();
		for(auto x : CPT){temp.push_back(x);}
		// cout<<endl;
		req.push_back(temp);
	}
	vector<vector<float>> final_req;
	for(int i = 0; i < sz;i++){
		vector<string> par = pointers[node_at_index[i]]->get_Parents();
		// int mul =1;
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
	for(int i = 0; i < sz;i++){
		Alarm.get_nth_node(i)->set_CPT(final_req[i]);
	}
    network_to_bif_format(Alarm,node_at_index,argv[1]);
	// vector<vector<float>> nex = finaL_prob();

    auto currentTime = chrono::high_resolution_clock::now();
	auto elapsedTime = chrono::duration_cast<chrono::seconds>(currentTime - startTime);
	// cout<<"Score: "<<answer<<endl;
    cout<<"Total iteration: "<<itr<<endl;
    cout<<"Total_time: "<<elapsedTime.count()<<endl;


}