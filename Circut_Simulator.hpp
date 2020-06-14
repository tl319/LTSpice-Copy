/*
 * lmao gotem
 */

/*
 * File:   Circut_Simulator.hpp
 * Author: Jason Tang, David Cormier, Thomas Van Issum
 *
 * Created on May 25, 2020, 12:16 PM
 */

#ifndef CIRCUT_SIMULATOR_HPP
#define CIRCUT_SIMULATOR_HPP

#include <iostream>
#include <string>
#include <vector>
#include <sstream>
#include <Eigen>
#include <bits/stdc++.h>
#include <cmath>
#include <chrono>
using namespace std;
using namespace Eigen;

struct Node
{
	int number;
	string label;
	int super;
	bool reactiveSuper = false;
	bool poser = false;
	Node(){number = -1; label = "N/A"; super = -99;}
	Node(string l){number = -1; super = -1; label=l;}
	Node(int n, string l, int s)
	{number = n; super = s; label=l;
	}


};

struct Param
{
	char var;
	float start;
	float stop;
	float interval;

	Param(char x, float st, float ed, float inter)
	{var =x;
	start= st;
	stop = ed;
	interval = inter;
	}

};

struct Simulation
{
	string type;
	float stop = -1;
    float step = -1;
	bool isStep =false;

	Simulation(){;}
	Simulation(string t, float stp, float sep)
	{type = t;
	stop = stp;
	sep = sep;
	}

	//only used in  step
	vector<Param> steps;
};

struct Component
{
	char type;
	string name;
	Node A;
	Node B;
	float value;
	bool poser = false;
	char var;
	//no zero when dealing with a time dependent signal
	bool isSignal = false;
	float DCOff=0;
	float amplitude=0;
	float frequency=0;

	// used for step param
	bool isVar = false;
		Component(){type = '0';name ="N/A";A.label="N/A";B.label="N/A";value=0;}
		Component(char t,string n, string nA,string nB,float v){
		type = t;
		name = n;
		A.label = nA;
		B.label = nB;
		value = v;
	}

        Component(char t,string n, Node nA,Node nB,float v){
		type = t;
		name = n;
		A = nA;
		B = nB;
		value = v;
	}
        Component(char t,string n, string nA,string nB,float DC, float a, float f){
		type = t;
		name = n;
		A.label = nA;
		B.label = nB;
		DCOff=DC;
		amplitude=a;
		frequency=f;
		isSignal = true;
	}
	Component(char t,string n, string nA,string nB,char s){
	isVar =true;
	type = t;
	name = n;
	var = s;
	A.label = nA;
	B.label = nB;
	}

};

struct Diode: Component{

	float is;
		Diode(char t,string n, string nA,string nB,string s){
		type = t;
		name = n;
		A.label = nA;
		B.label = nB;
		if(s=="1N4148"){
			is = 2.52*pow(2.52,-9);
		}
	}


};

struct Variable: Component{

	char var;

	Variable(char t,string n, string nA,string nB,char s){
	isVar =true;
	type = t;
	name = n;
	var = s;
	A.label = nA;
	B.label = nB;
	}
	char getVar()
	{
		return var;
	}


};

// *** Component Management Functions ***

pair<vector<Component>, Simulation> readInput();
//reads a netlist from stdin and sorts them into a vec of components and a simulation;

vector<Node> findNodes(const vector<Component> & list);
// returns a list of node. Duplicates are not counted and 0 is set to ground if it exits

vector<Component> patchComponents(const vector<Component> & list);
// returns a new vector of components that have their number and super vaules updated. Has to be new, as functions take in a copy for the vector.




//  *** Matrix Functions ***
VectorXd matrixSolve(const MatrixXd & m, const VectorXd & v);
//takes in the vector of currents and conductance matrix, and spits out the node voltage vector.

pair<MatrixXd, VectorXd> conductance_current (const vector<Component> & comps, const int & noden);
//return complete conductance matrix and current vector

void writeFile(const Node & n);
void writeFile(const float & time, const float & voltage);
//writes out in the format specified by LTspice/matlab
//input of node just writes the header *** NEEDED ***
//use the function after to write time and voltage.

VectorXd VectorUpdate (const vector<Component> & comps, const int & noden, const float & time, const VectorXd & pastnodes, const VectorXd & component_currents,
const float & interval, const vector<int> & c_vs_row);
//updates the rhs ("current") vector during transient simulations

pair<MatrixXd, vector<int>> MatrixUpdate (vector<Component> & comps, const int & noden, const float & interval);
//update matrix to reflect the change in behaviour of reactive components

//is the node vector actually used by the functions within? (and cannot be rpelaced by the label int)
pair<VectorXd, VectorXd> no_prior_change (const vector<Component> & comps, const vector<Node> & nodes, const int & noden);
//return voltage and current vectors for operating point or first point of transient analysis

void transient (vector<Component> & comps, const vector<Node> & nodes, const int & noden, const float & duration,
const float & interval, const VectorXd & pastnodes, const VectorXd & pastcurrents);




// *** Quality of Life functions ***
ostream& operator<<(ostream& os, const Component& c);
ostream& operator<<(ostream& os, const Node& c);
// allows you to cout both components and nodes.

float procData(const string & x);
//used to extract data ie. remove braces get K

void testBetter(const int & noden, const MatrixXd & conducts, const VectorXd & currents, const vector<Node> & nodes);
//print out conductance matrix and current vector
void writeTran(const vector<Node>& nlist, const vector<Component>& out,const VectorXd& pastnodes,const VectorXd& component_currents, float time);

void writeTranHeaders(const vector<Node>& nlist, const vector<Component>& out,const VectorXd& pastnodes, const VectorXd& component_currents);
void writeOPZero(const VectorXd& pastnodes, const VectorXd& component_currents);

void writeOP(const vector<Node>& nlist, const vector<Component>& out,const VectorXd& pastnodes, const VectorXd& component_currents);

void writeOPReadable(const vector<Node>& nlist, const vector<Component>& out,const VectorXd& pastnodes, const VectorXd& component_currents);



// *** Functions used by functions. Only use if necessary ***
bool operator == (Node const &n1, Node const &n2);
//used in the node finder function.

vector<Component> patchSupernodes(const vector<Component> & list);
//used in patchNodes to correctly fill out the values of the supernodes after the inital numbering operation

bool isComponent(const string & x);
//returns true if a line from a string (probs from stdin) is a component or not

string removeChar(const string & s, const char & style);
// removes unused charaters. Input 'D' extracts data input 'S' extracts scientific notation
bool isSci(const char & c);
bool isData(const char & c);
//used in above function

string nodeName(const int & i, const vector<Component> & out);
vector<Component> patchSupernodeInductor(const vector<Component> & list);

int compute_noden(const vector<Node> & nodes);
//calculate number of non-ground nodes from node vector

int nA(const Component & c);
int nB(const Component & c);
int SnA(const Component & c);
int SnB(const Component & c);
//read node/supernode A/B of a given component

vector<Component> common_node (const vector<Component> & comps, const Component & C, const Node & A);
//find all other components attached to a component C at a given node A

int component_index (const vector<Component> & comps, const Component & C);
//returns the index of a component in the component vector, useful since this and component current vector are ordered identically

VectorXd recursive_currents (const vector<Component> & comps, const vector<Node> & nlist, const VectorXd & nodev, const VectorXd & prevnodev,
const float & interval, const VectorXd & past_currents, const bool & op);
//compute currents accross "insufficient" (V, C, D) components in series with other such components

float recursive_basecase (const int & i, const Component & C, const vector<Component> & comps, const vector<Node> & nlist, const VectorXd & nodev,
const VectorXd & prevnodev, const float & interval, VectorXd & comp_currents, const Node & used_node, const bool & op);
//used in above

#endif
