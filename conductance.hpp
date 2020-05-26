#include <iostream>
#include <string>
#include <vector>
#include <numeric>
#include <Eigen>

using namespace Eigen;
using namespace std;

struct Node
{
	int number;
	string label;
	int super;
	Node(){number = -1; label = "N/A"; super = -99;}
	Node(string l){number = -1; super = -1; label=l;}
};

struct Component
{
	char type;
	string name;
	Node A;
	Node B;
	float value;

	Component(char t,string n, string nA,string nB,float v){
		type = t;
		name = n;
		A.label = nA;
		B.label = nB;
		value = v;
	}       
};

//calculate number of non-ground nodes from node vector
int compute_noden(vector<int> nodes);

//return complete conductance matrix and current vector
pair<MatrixXd, vector<float>> conductance_current (vector<Component> comps, int noden);

//print out conductance matrix and current vector
void test(int noden, MatrixXd conducts, vector<float> currents);

//read node/supernode A/B of a given component
int nA(Component c);
int nB(Component c);
int SnA(Component c);
int SnB(Component c);

//function declarations

void test(int noden, MatrixXd conducts, vector<float> currents)
{
    //print each line of the matrix
    for(int i = 0; i<noden; i++)
    {
        for(int j = 0; j<noden; j++)
        {
            cout << conducts(i, j) << " ";
        }
        cout << endl;
    }

    //print current vector
    for(int i = 0; i<noden; i++)
    {
        cout << currents[i] << endl;
    }

    cout << endl;
}

int nA(Component c)
{
    return (c.A.number);
}

int nB(Component c)
{
    return (c.B.number);
}

int SnA(Component c)
{
    return (c.A.super);
}

int SnB(Component c)
{
    return (c.B.super);
}

int compute_noden(vector<Node> nodes)
{
    //number of non-ground nodes
    int noden = nodes.size();
    //loop redundant if vector is ordered
    for(int i = 0; i<nodes.size(); i++)
    {
        if( nodes[i].number == 0)
        {
            noden--;
        }
    }

    return noden;
}

pair<MatrixXd, vector<float>> conductance_current (vector<Component> comps, int noden)
{

    //conductance matrix
    MatrixXd conducts = MatrixXd::Zero (noden, noden);
    //rhs vector
    vector<float> currents(noden, 0);

    //used to tell where the first empty spot in the node vector is
    int nodes_added = 0;
    //used to tell if an empty spot must be written to
    bool found = 0;

    float conductance;

    //"locked" matrix lines so that if they're written to by a v source they aren't written to again 
    vector<int> locked (noden, 0);

    //cycle through each component
    for(int i = 0; i<comps.size(); i++)
    {

        //dealing with resistors
        if(comps[i].type == 'R')
        {
            conductance = 1/comps[i].value;

            //if row nA has not already been edited to represent a voltage source and nA is not ground, 
            //then add and subtract conductance in columns nA and nB 
            //Row corresponds to supernode (node if no supernode) and column to actual node
            if(locked[ nA(comps[i])-1] == 0)
            {
                if( nA(comps[i]) != 0)
                {
                    conducts ( SnA(comps[i]) -1, nA(comps[i]) -1) += conductance;
                }
                if( nA(comps[i]) != 0 && nB(comps[i]) != 0)
                {
                    conducts ( SnA(comps[i]) -1, nB(comps[i]) -1) -= conductance;
                }
            }

            //Ditto for nB
            if(locked[ nB(comps[i])-1] == 0)
            {
                if( nB(comps[i]) != 0)
                {
                    conducts ( SnB(comps[i]) -1, nB(comps[i]) -1) += conductance;
                }
                if( nA(comps[i]) != 0 && nB(comps[i]) != 0)
                {
                    conducts ( SnB(comps[i]) -1, nA(comps[i]) -1) -= conductance;
                }
            }

        }

        //add case for nA = nB = 0?
        if(comps[i].type == 'V')
        {
            if( nB(comps[i]) == 0)
            {
                locked[nA(comps[i])-1] = 1;
                for(int j = 0; j<noden; j++)
                {
                    if(j == (nA(comps[i]) -1))
                    {
                        conducts (j, j) = 1;
                        currents[j] = comps[i].value; 
                    } else {
                        conducts (nA(comps[i]) -1, j) = 0;
                    }
                }
            }

            if(nA(comps[i]) == 0)
            {
                locked[nB(comps[i])-1] = 1;
                for(int j = 0; j<noden; j++)
                {
                    if(j == (nB(comps[i]) -1))
                    {
                        conducts (j, j) = 1;
                        currents[j] = (-1)*comps[i].value; 
                    } else {
                        conducts (nB(comps[i]) -1, j) = 0;
                    }
                }
            }

            /*/
            if(nA(comps[i]) != 0 && nB(comps[i]) != 0)
            {
                if( nA(comps[i]) > nB(comps[i]) )
                locked[comps[i].nA-1] = 1;
                for(int j = 0; j<noden; j++)
                {
                    if(j == (comps[i].nA -1))
                    {
                        conducts (nA(comps[i]) -1, j) = 1;
                        currents[j] = comps[i].value; 
                    } else {
                        conducts (nA(comps[i]) -1, j) = 0;
                    }
                }
            }
            /*/
        }

            //test(noden, conducts, currents);
    }

    return {conducts, currents};

}
