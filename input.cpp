
#include "Circut_Simulator.hpp"
#include <iostream>
#include <string>
#include <vector>
#include <sstream>
#include <Eigen>
#include <bits/stdc++.h>
#include <numeric>
#include <cmath>

using namespace std;
using namespace Eigen;


bool isComponent(string x)
{
	if(isalpha(x[0])){
		return true;}
	return false;}
bool isCmd(string x)
{
	if(x[0] == '.' && x[1]!='e'){
		return true;}
	return false;}

VectorXd matrixSolve(MatrixXd m,VectorXd v)
{
	MatrixXd inverse = m.inverse();
	return inverse *v;}


ostream& operator<<(ostream& os, const Component& c)
{
    if(c.isSignal)
    {os << c.type << ':' << c.name << ':' << c.A.label << ':' << c.B.label << ':' << c.DCOff << ':' << c.amplitude << ':' << c.frequency << endl;}
    else
    {os << c.type << ':' << c.name << ':' << c.A.label << ':' << c.B.label << ':' << c.value << endl;}
    return os;
}
ostream& operator<<(ostream& os, const Node& c)
{

    os << "Number: " << c.number  << " label: "<< c.label << " supernode: "<< c.super << endl;
    return os;
}
ostream& operator<<(ostream& os, const Simulation& c)
{

    os << "Type: " << c.type  << " stop: "<< c.stop << " timestep: "<< c.step << endl;
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
            if(out[i].type == 'V' || out[i].type == 'v' || out[i].type == 'C' || out[i].type == 'c' || out[i].type == 'I' || out[i].type == 'i'){
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

bool isData(char c){
    if(isdigit(c) || c == '.')
        return false;
    return true;}

bool isSci(char c){
    if(c == 'p' || c == 'P' ||c == 'N' || c == 'u' || c == 'U' || c == 'm' || c == 'K' || c == 'k' || c == 'M' || c == 'G')
        return false;
    return true;}

string removeChar(string s, char style)
{
    string x = s;
    if(style == 'D'){   
        x.erase(std::remove_if(x.begin(), x.end(), isData), x.end());}
    
    if(style == 'S'){
    x.erase(std::remove_if(x.begin(), x.end(), isSci), x.end());}
    
    return x;
}

float procData(string x)
{

    float num = stof(removeChar(x,'D'));
    if(removeChar(x,'S') != "")
    {
        string sci = removeChar(x,'S');
        if(tolower(sci[0]) == tolower('k'))
            return num * 1000;
        else if(tolower(sci[0]) == tolower('p'))
            return num * pow(10,-12);
        else if(tolower(sci[0]) == tolower('n'))
            return num * pow(10,-9);
        else if(tolower(sci[0]) == tolower('u'))
            return num * pow(10,-6);
        else if(sci[0] == 'm')
            return num * pow(10,-3);
        else if(sci[0] == 'M')
            return num * pow(10,6);
        else if(sci[0] == tolower('G'))
            return num * pow(10,-3);
        else{
            return num;
        }
               
    }
    return num;
}

void writeFile(Node n)
{
    cout << "time" << '\t' << n.label << endl;
}

void writeFile(float time, float voltage)
{
    cout << time << '\t' << voltage << endl;
}


pair<vector<Component>, Simulation> readInput()
{
	string x;
	vector<string> strings;
	vector<Component> components;
        Simulation sim;
	while(getline(cin, x))
	{
		strings.push_back(x);
	}

	for(auto line : strings){
		if(isComponent(line)){
			vector<string> properties;
			stringstream ss(line);
			int count=0;
			 while (ss >> x){
                             properties.push_back(x);
					}
			string name;
			if(isalnum((properties[0])[1])){name = properties[0];}
			else{name =(properties[0]).substr (3,(properties[0].length())-1);}
                        if(properties.size()<5){
			Component c1((properties[0])[0],name,properties[1],properties[2],procData(properties[3]));
                        components.push_back(c1);}
                        else{
			Component c1((properties[0])[0],name,properties[1],properties[2],procData(properties[3]),procData(properties[4]),procData(properties[5]));
                        components.push_back(c1);}
			
                        
		}
                if(isCmd(line))
                {
                    vector<string> properties;
			stringstream ss(line);
			int count=0;
			 while (ss >> x){
                             properties.push_back(x);
					}
                        string type =(properties[0]).substr (1,(properties[0].length())-1);
                        if(type == "op")
                        {sim.type=type;}
                        if(type == "tran"){
                        sim.type=(properties[0]).substr (1,(properties[0].length())-1);
                        sim.stop=(procData(properties[2]));
                        sim.step=(procData(properties[4]));
                        }
                }
	}
        


	return make_pair(components, sim);

}
int main()
{
    
    pair<vector<Component>, Simulation> testm = readInput();
    vector<Component> out = patchComponents(testm.first);
    vector<Node> nlist = findNodes(out);

    VectorXd component_currents = VectorXd::Zero (out.size());

    /*/
    cout << testm.second;
    for(auto x : out)
    {
        cout << x;
    }
    
    for (auto x : nlist)
    { cout << "Name is :" << x.label <<  " number is :" << x.number << endl; 
    }
    for(auto x : out)
    {
            cout << x.name << endl << "A :"<< x.A << "B :" << x.B << endl;
    }
    /*/
   
    int noden = compute_noden(nlist);

    cout << endl;

    pair<MatrixXd, VectorXd> knowns = conductance_current (out, noden, 0, 0, component_currents);

    cout << endl;

    //cout << knowns.first << endl;
    
    /*/ not quit sure what this is
    VectorXd v(noden);
    for(int i =0;i<noden;i++){v(i) = knowns.second[i];}
    cout << v << endl << endl;
    /*/

    test(noden, knowns.first, knowns.second);

    cout << endl;

    VectorXd pastnodes;

    pastnodes = matrixSolve(knowns.first, knowns.second);

    cout << pastnodes << endl;

    cout << endl;

    component_currents = comp_currents (out, nlist, pastnodes, 0);

    cout << component_currents << endl;    
    
    cout << "f ";

    /*/
    VectorXd transrhs;
    for(int i = 1; i <= 10; i++)
    {
        cout << "test" << endl;
        transrhs = VectorUpdate (out, noden, i, pastnodes, component_currents, 1);
        knowns = conductance_current (out, noden, i, 1, component_currents);
        pastnodes << matrixSolve(knowns.first, knowns.second);
        cout << endl;
        cout << pastnodes;
    }
    /*/
    
}
