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
using namespace std;
using namespace Eigen;

struct Node
{
	int number;
	string label;
	int super;
	bool reactiveSuper = false;
	Node(){number = -1; label = "N/A"; super = -99;}
	Node(string l){number = -1; super = -1; label=l;}
	Node(int n, string l, int s)
	{number = n; super = s; label=l;
	}
        

};
struct Simulation
{
	string type;
	float stop = -1;
    float step = -1;
        
	Simulation(){;}
	Simulation(string t, float stp, float sep)
        {type = t;
        stop = stp;
        sep = sep;
        }
};

struct Component
{
	char type;
	string name;
	Node A;
	Node B;
	float value;
        
        //no zero when dealing with a time dependent signal
        bool isSignal = false;
        float DCOff=0;
        float amplitude=0;
        float frequency=0;

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
       
};

// *** Component Management Functions ***

pair<vector<Component>, Simulation> readInput();
//reads a netlist from stdin and sorts them into a vec of components and a simulation;

vector<Node> findNodes(vector<Component> list);
// returns a list of node. Duplicates are not counted and 0 is set to ground if it exits

vector<Component> patchComponents(vector<Component> list);
// returns a new vector of components that have their number and super vaules updated. Has to be new, as functions take in a copy for the vector.




//  *** Matrix Functions ***
VectorXd matrixSolve(MatrixXd m,VectorXd v);
//takes in the vector of currents and conductance matrix, and spits out the node voltage vector. 

pair<MatrixXd, VectorXd> conductance_current (vector<Component> comps, int noden);
//return complete conductance matrix and current vector

void writeFile(Node n);
void writeFile(float time, float voltage);
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

vector<pair<VectorXd, VectorXd>> transient (vector<Component> & comps, const vector<Node> & nodes, const int & noden, const float & duration, 
const float & interval, const VectorXd & pastnodes, const VectorXd & pastcurrents);




// *** Quality of Life functions ***
ostream& operator<<(ostream& os, const Component& c);
ostream& operator<<(ostream& os, const Node& c);
// allows you to cout both components and nodes.

float procData(string x);
//used to extract data ie. remove braces get K

void test(int noden, MatrixXd conducts, VectorXd currents);
//print out conductance matrix and current vector

void writeTran(const VectorXd& pastnodes,const VectorXd& component_currents, float time);

void writeTranHeaders(const vector<Node>& nlist, const vector<Component>& out);

void writeOP(const vector<Node>& nlist, const vector<Component>& out,const VectorXd& pastnodes, const VectorXd& component_currents);

void writeOPReadable(const vector<Node>& nlist, const vector<Component>& out,const VectorXd& pastnodes, const VectorXd& component_currents);



// *** Functions used by functions. Only use if necessary ***
bool operator == (Node const &n1, Node const &n2);
//used in the node finder function.

vector<Component> patchSupernodes(vector<Component> list);
//used in patchNodes to correctly fill out the values of the supernodes after the inital numbering operation

bool isComponent(string x);
//returns true if a line from a string (probs from stdin) is a component or not

string removeChar(string s, char style);
// removes unused charaters. Input 'D' extracts data input 'S' extracts scientific notation
bool isSci(char c);
bool isData(char c);
//used in above function

int compute_noden(vector<Node>  nodes);
//calculate number of non-ground nodes from node vector

int nA(Component c);
int nB(Component c);
int SnA(Component c);
int SnB(Component c);
//read node/supernode A/B of a given component

vector<Component> common_node (const vector<Component> & comps, const Component & C, const Node & A);
//find all other components attached to a component C at a given node A

int component_index (const vector<Component> & comps, const Component & C);
//returns the index of a component in the component vector, useful since this and component current vector are ordered identically

float vs_current (vector<Component> comps, Component C, vector<bool> & computed, VectorXd currents, Node N);
//used to compute current through voltage sources or capacitors

VectorXd recursive_currents (const vector<Component> & comps, const vector<Node> & nlist, const VectorXd & nodev, const float & interval, 
const VectorXd & past_currents);
//compute currents accross "insufficient" (V, C, D) components in series with other such components

float recursive_basecase (const int & i, const Component & C, const vector<Component> & comps, const vector<Node> & nlist, VectorXd nodev, 
const float & interval, vector<bool> & computed, VectorXd & comp_currents, const Node & used_node);
//used in above

vector<bool> incorrect_assumptions(VectorXd comp_currents, vector<Component> comps);
//check currents accross nonlinear components to determine their mode

pair<MatrixXd, vector<int>> CorrectAssumptions (vector<Component> comps, int noden, vector<bool> incorrect_assumptions);
//rewrite matrix with new assumptions
//obsolete, change matrix to vector

#endif