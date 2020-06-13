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


bool isComponent(const string & x)
{
	if(isalpha(x[0])){
		return true;}
	return false;}
bool isCmd(const string & x)
{
	if(x[0] == '.' && x[1]!='e'){
		return true;}
	return false;}

VectorXd matrixSolve(const MatrixXd & m, const VectorXd & v)
{
	MatrixXd inverse = m.inverse();
	return inverse *v;}


ostream& operator<<(ostream& os, const Node& c)
{

    os << "Number: " << c.number  << " label: "<< c.label << " supernode: "<< c.super << " reactiveSuper: " << c.reactiveSuper;
    return os;
}

ostream& operator<<(ostream& os, const Component& c)
{
    if(c.isSignal)
    {os << c.type << ':' << c.name << ':' << c.A.label << ':' << c.B.label << ':' << c.DCOff << ':' << c.amplitude << ':' << c.frequency << endl;}
    else
    {os << c.type << ':' << c.name << ':' << c.A.label << ':' << c.B.label << ':' << c.value << endl;}
    return os;
}
ostream& operator<<(ostream& os, const Simulation& c)
{

    os << "Type: " << c.type  << " stop: "<< c.stop << " timestep: "<< c.step << endl;
    return os;
}
string nodeName(const int & i, const vector<Component> & out)
{
    for(auto x : out){
        if(x.A.number == i)
            return x.A.label;
        else if(x.B.number == i)
            return x.B.label;
    }
    return "not found node number";
}


bool operator == (Node const &n1, Node const &n2)
{
     return(n1.label == n2.label);
}

/*/
vector<Node> findNodes(vector<Component> list)
{
	vector<Node> Nodes;
    vector<Node> children;
	for (Component part : list){
		if(find (Nodes.begin(), Nodes.end(), part.A) == Nodes.end()){
			if (part.A.label=="0"){
				Nodes.insert(Nodes.begin(),part.A);}
            else
            {
                Nodes.push_back(part.A);
            }
        }
		if(find (Nodes.begin(), Nodes.end(), part.B) == Nodes.end()){
			if (part.B.label=="0"){
				Nodes.insert(Nodes.begin(),part.B);}
            else if (part.B.poser && find(children.begin(), children.end(), part.B) == children.end())
            {
                children.push_back(part.B);
            }
			else
				{Nodes.push_back(part.B);}
        }
	}
    for(auto x : children)
    {
        if(find (Nodes.begin(), Nodes.end(), x) != Nodes.end()){
                Nodes.erase(std::remove(Nodes.begin(), Nodes.end(), x), Nodes.end());
                Nodes.insert(Nodes.begin()+1,x);
                cerr << "inserted " << x.label <<endl;
            }
        else{
             Nodes.insert(Nodes.begin()+1,x);
             cerr << "inserted " << x.label << endl;
        }
    }
	return Nodes;
}
/*/
vector<Node> findNodes(vector<Component> list)
{
	vector<Node> Nodes;
    vector<Node> children;
	for (Component part : list){
		if(find (Nodes.begin(), Nodes.end(), part.A) == Nodes.end()){
			if (part.A.label=="0"){
				Nodes.insert(Nodes.begin(),part.A);}
            else
            {
                Nodes.push_back(part.A);
            }
        }
		if(find (Nodes.begin(), Nodes.end(), part.B) == Nodes.end()){
			if (part.B.label=="0"){
				Nodes.insert(Nodes.begin(),part.B);}
            else if (part.B.poser && find(children.begin(), children.end(), part.B) == children.end())
            {
                children.push_back(part.B);
            }
			else
				{Nodes.push_back(part.B);}
        }
	}
    for(auto x : children)
    {
        if(find (Nodes.begin(), Nodes.end(), x) != Nodes.end()){
                Nodes.erase(std::remove(Nodes.begin(), Nodes.end(), x), Nodes.end());
                Nodes.insert(Nodes.begin()+1,x);
                cerr << "inserted " << x.label <<endl;
            }
        else{
             Nodes.insert(Nodes.begin()+1,x);
             cerr << "inserted " << x.label << endl;
        }
    }
	return Nodes;
}

vector<Component> patchSupernodes(vector<Component> list)
{
    vector<Component> out = list;
    for(int i=0;i<out.size();i++){
            if(out[i].type == 'V' || out[i].type == 'v' || out[i].type == 'C' || out[i].type == 'c'|| out[i].type == 'D' || out[i].type == 'd'){
                int topnode = 99;
                int botnode = 99;
                //cout << out[i].A.number << " b is: " << out[i].B.number << endl;
                if(out[i].A.super>out[i].B.super)
                {topnode = out[i].A.super;
                botnode =out[i].B.super;}
                else {topnode = out[i].B.super;
                botnode =out[i].A.super;}
                //cerr << "topnode is :" << topnode << " botnode is: " << botnode << endl;

                for(int x=0;x<out.size();x++){
                    //cout << "loop: " << x << endl;
                    if(out[x].A.super == botnode){
                    out[x].A.super = topnode;
                    out[x].A.reactiveSuper=(out[i].type == 'C' || out[i].type == 'c');
                    //cout << "set " << out[x].name << out[x].A.number << "to " <<  out[x].A.super << endl;
                    }
                    else
                    {//out[x].A.super = out[x].A.super;
                    //cout <<  out[x].name << "has " << out[x].A.number << " not equal " << botnode << endl;
                    }

                    if(out[x].B.super == botnode){
                    out[x].B.super = topnode;
                    out[x].B.reactiveSuper=(out[i].type == 'C' || out[i].type == 'c');
                    //cout << "set " << out[x].name << out[x].B.number << " to " <<  out[x].B.super << endl;
                    }
                    else
                    {//out[x].B.super = out[x].A.number;
                    //cout << out[x].B.number << " not equal " << botnode << endl;
                    }

                }
                }


            }
    return out;
}

vector<Component> patchSupernodeInductor(vector<Component> list)
{
    vector<Component> out = list;
    for(int i=0;i<out.size();i++){
    out[i].A.super = out[i].A.number;
    out[i].B.super = out[i].B.number;
    }
    for(int i=0;i<out.size();i++){
            if(out[i].type == 'V' || out[i].type == 'v' || out[i].type == 'C' || out[i].type == 'c'|| out[i].type == 'D' || out[i].type == 'd'){
                int topnode = 99;
                int botnode = 99;
                //cout << out[i].A.number << " b is: " << out[i].B.number << endl;
                if(out[i].A.super>out[i].B.super)
                {topnode = out[i].A.super;
                botnode =out[i].B.super;}
                else {topnode = out[i].B.super;
                botnode =out[i].A.super;}
                //cerr << "topnode is :" << topnode << " botnode is: " << botnode << endl;

                for(int x=0;x<out.size();x++){
                    //cout << "loop: " << x << endl;
                    if(out[x].A.super == botnode){
                    out[x].A.super = topnode;
                    out[x].A.reactiveSuper=(out[i].type == 'C' || out[i].type == 'c');
                    //cout << "set " << out[x].name << out[x].A.number << "to " <<  out[x].A.super << endl;
                    }
                    else
                    {//out[x].A.super = out[x].A.super;
                    //cout <<  out[x].name << "has " << out[x].A.number << " not equal " << botnode << endl;
                    }

                    if(out[x].B.super == botnode){
                    out[x].B.super = topnode;
                    out[x].B.reactiveSuper=(out[i].type == 'C' || out[i].type == 'c');
                    //cout << "set " << out[x].name << out[x].B.number << " to " <<  out[x].B.super << endl;
                    }
                    else
                    {//out[x].B.super = out[x].A.number;
                    //cout << out[x].B.number << " not equal " << botnode << endl;
                    }

                }
                }


            }
    return out;
}

vector<Component> patchComponents(const vector<Component> & list)
{
	vector<Component> out = list;
	vector<Node> nodes = findNodes(out);
    for (auto x : nodes)
    {
        cerr << x.label <<endl;
    }
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

bool isData(const char & c){
    if(isdigit(c) || c == '.')
        return false;
    return true;}

bool isSci(const char & c){
    if(c == 'p' || c == 'P' ||c == 'N' || c == 'u' || c == 'U' || c == 'm' || c == 'K' || c == 'k' || c == 'M' || c == 'G' || c == 'µ')
        return false;
    return true;}

string removeChar(const string & s, const char & style)
{
    string x = s;
    if(style == 'D'){
        x.erase(std::remove_if(x.begin(), x.end(), isData), x.end());}

    if(style == 'S'){
    x.erase(std::remove_if(x.begin(), x.end(), isSci), x.end());}

    return x;
}

float procData(const string & x)
{

    float num = stof(removeChar(x,'D'));
    if(removeChar(x,'S') != "")
    {
        string sci = removeChar(x,'S');
        if(tolower(sci[0]) == tolower('k'))
            num = num * 1000;
        else if(tolower(sci[0]) == tolower('p'))
            num = num * pow(10,-12);
        else if(tolower(sci[0]) == tolower('n'))
            num = num * pow(10,-9);
        else if((tolower(sci[0]) == tolower('u')) || sci[0] == 'µ')
            num = num * pow(10,-6);
        else if(sci[0] == 'm')
            num = num * pow(10,-3);
        else if(sci[0] == 'M')
            num = num * pow(10,6);
        else if(sci[0] == tolower('G'))
            num =num * pow(10,-3);
        else{
            return -1;
        }

    }
    size_t found = x.find('-');
    if (found != string::npos)
        {return num * -1;}
    else
        return num;
}

void writeOP(const vector<Node>& nlist, const vector<Component>& out,const VectorXd& pastnodes, const VectorXd& component_currents)
{

    bool fakes[nlist.size()+out.size()-1] = {false};
    int count = 0;
    cout << nlist[1].label;
    for(int i = 2;i<nlist.size();i++) {
        if(!nlist[i].poser){
        cout << '\t' << nlist[i].label;
        }
        else{fakes[count]=1;}
        count ++;
    }
    for(int i = 0;i<out.size();i++) {
        if(!out[i].poser){
        cout << '\t' << out[i].name;
        }
        else{fakes[count]=1;}
        count ++;
    }
    cout << endl;
    count = 0;
    for(int i = 0;i<pastnodes.size();i++) {
        if(!fakes[count]){
        cout << '\t' << pastnodes(i);
        }
        count ++;
    }
     for(int i = 0;i<component_currents.size();i++) {
        if(!fakes[count]){
        cout << '\t' << component_currents(i);
        }
        count++;
     }
    cout << endl;

}
void writeOPReadable(const vector<Node>& nlist, const vector<Component>& out,const VectorXd& pastnodes, const VectorXd& component_currents)
{
    for(int i = 0;i<pastnodes.size();i++) {
        if(!nlist[i+1].poser){
        cout << "V(" << nlist[i+1].label << ") : ";
        cout << pastnodes(i) << endl;}
    }
    for(int i = 0;i<out.size();i++) {
        if(!out[i].poser){
        cout << "I(" << out[i].name << ") : ";
        cout << component_currents(i) << endl;
        }
    }
    cout << endl;
}
void writeOPZero(const VectorXd& pastnodes, const VectorXd& component_currents)
{
    cout << 0;
    for(int i = 0;i<pastnodes.size();i++) {
        cout << '\t' << pastnodes(i);
    }
     for(int i = 0;i<component_currents.size();i++) {
        cout << '\t' << component_currents(i);
    }
    cout << endl;
}

void writeTranHeaders(const vector<Node>& nlist, const vector<Component>& out,const VectorXd& pastnodes, const VectorXd& component_currents)
{
    bool fakes[nlist.size()+out.size()-1] = {false};
    int count = 0;
    cout << "time";
    for(int i = 1;i<nlist.size();i++) {
        if(!nlist[i].poser){
        cout << '\t' << nlist[i].label;
        }
        else{fakes[count]=1;}
        count ++;
    }
    for(int i = 0;i<out.size();i++) {
        if(!out[i].poser){
        cout << '\t' << out[i].name;
        }
        else{fakes[count]=1;}
        count ++;
    }
    cout << endl;
    cout << 0;
    count = 0;
    for(int i = 0;i<pastnodes.size();i++) {
        if(!fakes[count]){
        cout << '\t' << pastnodes(i);
        }
        count ++;
    }
     for(int i = 0;i<component_currents.size();i++) {
        if(!fakes[count]){
        cout << '\t' << component_currents(i);
        }
        count++;
     }
    cout << endl;

}

void writeTran(const vector<Node>& nlist, const vector<Component>& out,const VectorXd& pastnodes,const VectorXd& component_currents, float time)
{
     bool fakes[nlist.size()+out.size()-1] = {false};
    int count = 0;
    for(int i = 1;i<nlist.size();i++) {
        if(!nlist[i].poser){;}
        else{fakes[count]=1;}
        count ++;
    }
    for(int i = 0;i<out.size();i++) {
        if(!out[i].poser){;}
        else{fakes[count]=1;}
        count ++;
    }
    cout << time;
    count = 0;
    for(int i = 0;i<pastnodes.size();i++) {
        if(!fakes[count]){
        cout << '\t' << pastnodes(i);
        }
        count ++;
    }
     for(int i = 0;i<component_currents.size();i++) {
        if(!fakes[count]){
        cout << '\t' << component_currents(i);
        }
        count++;
     }
    cout << endl;
}


pair<vector<Component>, Simulation> readInput()
{
	string x;
	vector<string> strings;
	vector<Component> components;
    int fakeNode = 1;
    int fakeRes = 1;
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
                        if((properties[0])[0]=='D'||(properties[0])[0]=='d'){
                            Diode d1('D',name,properties[1],properties[2],properties[3]);
                            components.push_back(d1);
                        }
                        else if((properties[0])[0]=='C'||(properties[0])[0]=='c'){
                            Component c1('C',name,properties[1],"fakeNode"+to_string(fakeNode),procData(properties[3]));
                            Component c2('R',"fakeRes"+to_string(fakeRes),"fakeNode"+to_string(fakeNode),properties[2],(float)(0.0001/procData(properties[3])));
                            c2.poser=true;
                            c1.B.poser=true;
                            c2.A.poser=true;
                            components.push_back(c1);
                            components.push_back(c2);
                            fakeNode++;
                            fakeRes++;
                        }
                        else if((properties[3]).find('{') != std::string::npos){
                            Component v1(toupper((properties[0])[0]),name,properties[1],properties[2],(properties[3])[1]);
                            components.push_back(v1);
                        }
                        else if((properties[3]).find('{') != std::string::npos){
                            Component v1(toupper((properties[0])[0]),name,properties[1],properties[2],(properties[3])[1]);
                            components.push_back(v1);
                        }

                        else if(properties.size()<5){
			                Component c1(toupper((properties[0])[0]),name,properties[1],properties[2],procData(properties[3]));
                            components.push_back(c1);
                        }
                        else{

                            string DCoff = properties[3].substr(4,properties[3].length());
			                Component c1(toupper((properties[0])[0]),name,properties[1],properties[2],procData(DCoff),procData(properties[4]),procData(properties[5]));
                            components.push_back(c1);
                        }

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
                        if(type == "step")
                        {
                            sim.isStep = true;
                            Param p1 = {properties[2][0],procData(properties[3]),procData(properties[4]),procData(properties[5])};
                            sim.steps.push_back(p1);
                            cerr << "Pusshed back";
                        }

                        if(type == "tran"){
                            if(properties.size()>3){
                                sim.type=type;
                                sim.stop=(procData(properties[2]));
                                sim.step=(procData(properties[4]));
                            }
                            else{
                                sim.type=type;
                                sim.stop=(procData(properties[1]));
                                sim.step = 99999999;
                            }

                        }

                }
    }
	return make_pair(components, sim);

}

vector<Component> reorderVoltages(const vector<Component> &in)
{
    vector<Component> res;
    vector<Component> vol;
    vector<Component> duo;
    for(auto x : in)
    {
        if(x.type=='V' || x.type=='v')
            vol.push_back(x);
        else
            res.push_back(x);
    }
    reverse(vol.begin(),vol.end());
    duo.reserve(res.size() + vol.size() ); // preallocate memory
    duo.insert( duo.end(), res.begin(), res.end() );
    duo.insert( duo.end(), vol.begin(), vol.end() );
    return duo;
}

void stepOP(vector<Component> &list, const Simulation & sim, const int & i, const vector<Node> & nlist)
{
    int j = 0;
    int counter =0;
    Param values = sim.steps[i];
    while(counter<= values.stop){
;
        for(int i =0; i<list.size();i++)
        {
            if(list[i].isVar){
                if(list[i].var == values.var){
                    list[i].value = (values.start + j*values.interval);
                }
            }
        }
        int noden = compute_noden(nlist);
        pair<VectorXd, VectorXd> knowns = no_prior_change (list, nlist, noden);
        writeOPReadable(nlist, list, knowns.first, knowns.second);
        cout << endl << endl;
        j++;
        counter = values.start+(values.interval*j);
    }
}
void stepTran(vector<Component> &list, const Simulation & sim, const int & i, const vector<Node> & nlist, const float & duration, const float & interval)
{
    int j = 0;
    int counter =0;
    Param values = sim.steps[i];
    while(j< (ceil((values.stop-values.start)/values.interval))+1){
        cout << "Step Information : " << values.var << "=" << values.start+(values.interval*j);
        cout << " Run : " << j+1<< "/" << (ceil((values.stop-values.start)/values.interval))+1 << endl;
        for(int i =0; i<list.size();i++)
        {
            if(list[i].isVar){
                if(list[i].var == values.var){
                    list[i].value = (values.start + j*values.interval);
                }
            }
        }
        int noden = compute_noden(nlist);
        pair<VectorXd, VectorXd> knowns = no_prior_change (list, nlist, noden);
        transient (list, nlist, noden, duration, interval, knowns.first, knowns.second);
        cout << endl << endl;
        j++;
        counter = values.start+(values.interval*j);
    }
}
int main()
{
    pair<vector<Component>, Simulation> testm = readInput();
    Simulation sim =testm.second;
    vector<Component> out = patchComponents(testm.first);
    vector<Node> nlist = findNodes(out);
    for(auto x : nlist)
    {
        //cerr << x << endl;
    }
    //cerr << endl;
    for(auto x : out)
    {
        cerr << x;
        cerr << "A is "<< x.A << " superlabel: " << nodeName(x.A.super,out) << endl;
        cerr << "B is "<< x.B << " superlabel: " << nodeName(x.B.super,out) << endl;
    }





    //main running part
    if(sim.isStep){
        if(sim.type=="op"){
            for(int i = 0; i<sim.steps.size();i++){
                stepOP(out,sim,i,nlist);
                }
        }
        else if (sim.type=="tran"){
            float duration = sim.stop;
            float interval = 0.000001;
            for(int i = 0; i<sim.steps.size();i++){
                stepTran(out,sim,0,nlist,duration,interval);
            }


        }
    }

    else if(sim.type=="op"){
        int noden = compute_noden(nlist);
        pair<VectorXd, VectorXd> knowns = no_prior_change (out, nlist, noden);
        writeOPReadable(nlist, out, knowns.first, knowns.second);

    }
    else if (sim.type=="tran"){

        int noden = compute_noden(nlist);
        pair<VectorXd, VectorXd> knowns = no_prior_change (out, nlist, noden);
        float duration = sim.stop;
        float interval = 0.000001;
        transient (out, nlist, noden, duration, interval, knowns.first, knowns.second);
    }
    else
    {
        cout << "Invalid simulation command";
    }
}
