
#include "Circut_Simulator.hpp"
#include <iostream>
#include <string>
#include <vector>
#include <sstream>
#include <Eigen>
#include <bits/stdc++.h>
#include <numeric>

using namespace std;
using namespace Eigen;


bool isComponent(string x)
{
	if(isalpha(x[0])){
		return true;}
	return false;}

VectorXd matrixSolve(MatrixXd m,VectorXd v)
{
	MatrixXd inverse = m.inverse();
	return inverse *v;}


ostream& operator<<(ostream& os, const Component& c)
{
    os << c.type << ':' << c.name << ':' << c.A.label << ':' << c.B.label << ':' << c.value << endl;
    return os;
}
ostream& operator<<(ostream& os, const Node& c)
{

    os << "Number: " << c.number  << " label: "<< c.label << " supernode: "<< c.super << endl;
    return os;
}

bool operator == (Node const &n1, Node const &n2)
{
     return(n1.label == n2.label);
}

vector<Node> findNodes(vector<Component> list)
{
	vector<Node> Nodes;
	for (Component part : list){
		if(find (Nodes.begin(), Nodes.end(), part.A) == Nodes.end()){
			if (part.A.label=="0")
				Nodes.insert(Nodes.begin(),part.A);
			else
				Nodes.push_back(part.A);}
		if(find (Nodes.begin(), Nodes.end(), part.B) == Nodes.end()){
			if (part.B.label=="0")
				Nodes.insert(Nodes.begin(),part.B);
			else
				Nodes.push_back(part.B);}

	}
	return Nodes;
}
vector<Component> patchSupernodes(vector<Component> list)
{
    vector<Component> out = list;
    for(int i=0;i<out.size();i++){
            if(out[i].type == 'V' || out[i].type == 'v'){
                int topnode = 99;
                int botnode = 99;
                //cout << out[i].A.number << " b is: " << out[i].B.number << endl;
                if(out[i].A.super>out[i].B.super)
                {topnode = out[i].A.super;
                botnode =out[i].B.super;}
                else {topnode = out[i].B.super;
                botnode =out[i].A.super;}
                 cout << "topnode is :" << topnode << " botnode is: " << botnode << endl;
                 
                for(int x=0;x<out.size();x++){
                    //cout << "loop: " << x << endl;
                    if(out[x].A.super == botnode){
                    out[x].A.super = topnode;
                    cout << "set " << out[x].name << out[x].A.number << "to " <<  out[x].A.super << endl;
                    }
                    else
                    {//out[x].A.super = out[x].A.super;
                    cout <<  out[x].name << "has " << out[x].A.number << " not equal " << botnode << endl;       
                    }
                    
                    if(out[x].B.super == botnode){
                    out[x].B.super = topnode;
                    cout << "set " << out[x].name << out[x].B.number << " to " <<  out[x].B.super << endl;
                    }
                    else
                    {//out[x].B.super = out[x].A.number;
                    cout << out[x].B.number << " not equal " << botnode << endl;
                    }

                }
                }
                
                
            }
    return out;
}

vector<Component> patchComponents(vector<Component> list)
{
	vector<Component> out = list;
	vector<Node> nodes = findNodes(out);
        int i;
        if(nodes[0].label=="0")
            i = 0;
        else
            i = 1;
        
	for(i;i<nodes.size();i++){
		for(int j = 0;j<out.size();j++){
			if ((out[j].A.label)==nodes[i].label){
				out[j].A.number=i;
                                out[j].A.super=i;
				//cout << "set " << out[j].A.label <<  " " << out[j].A.number << endl;
			}
			if ((out[j].B.label)==nodes[i].label){
                                out[j].B.number=i;
                                out[j].B.super=i;
                                //cout << "set " << out[j].B.label <<  " " << out[j].B.number << endl;
						}
		}
	}    
        //return out;
	return patchSupernodes(out);
}

vector<Component> readInput()
{
	string x;
	vector<string> strings;
	vector<Component> components;
	while(getline(cin, x))
	{
		strings.push_back(x);
	}

	for(auto line : strings){
		if(isComponent(line)){
			string properties[4];
			stringstream ss(line);
			int count=0;
			 while (ss >> properties[count]){
					count++; }
			string name;
			if(isalnum((properties[0])[1])){name = properties[0];}
			else{name =(properties[0]).substr (3,(properties[0].length())-1);}

			Component c1((properties[0])[0],name,properties[1],properties[2],stof(properties[3]));
			components.push_back(c1);
		}
		}


	return components;

}
int main()
{
    
    vector<Component> test = readInput();
    vector<Component> out = patchComponents(test);
    vector<Node> nlist = findNodes(out);
    
    for (auto x : nlist)
    { cout << "Name is :" << x.label <<  " number is :" << x.number << endl; 
    }
    for(auto x : out)
    {
            cout << x.name << endl << "A :"<< x.A << "B :" << x.B << endl;
    }


    int noden = compute_noden(findNodes(out));
    pair<MatrixXd, vector<float>> knowns = conductance_current (out, noden);
    cout << knowns.first << endl;
    VectorXd v(noden);
    for(int i =0;i<noden;i++){v(i) = knowns.second[i];}
    cout << v << endl << endl;
 
    cout << matrixSolve(knowns.first,v);
        
}
