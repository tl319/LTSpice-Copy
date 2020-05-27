#include "Circut_Simulator.hpp"
#include <iostream>
#include <string>
#include <vector>
#include <numeric>
#include <Eigen>

using namespace Eigen;
using namespace std;


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

        //dealing with voltage sources                  add case for nA = nB = 0?
        if(comps[i].type == 'V')
        {
            //negative terminal to ground
            if( nB(comps[i]) == 0 && nA(comps[i]) != 0)
            {
                //prevent ulterior editing of the matrix row assigned to represent the voltage source
                locked[nA(comps[i])-1] = 1;
                //cycle through all columns to  edit the row corresponding to the positive terminal node
                for(int j = 0; j<noden; j++)
                {
                    //at the column corresponding to the positive terminal, write 1
                    if(j == (nA(comps[i]) -1))
                    {
                        conducts (j, j) = 1;
                        //also write the source voltage to this index in the rhs vector
                        currents[j] = comps[i].value; 
                    } else {
                        //other columns are set to 0
                        conducts (nA(comps[i]) -1, j) = 0;
                    }
                }
            }

            //positive terminal to ground
            if(nA(comps[i]) == 0 && nB(comps[i]) != 0)
            {
                //prevent ulterior editing of the matrix row assigned to represent the voltage source
                locked[nB(comps[i])-1] = 1;
                //cycle through all columns to  edit the row corresponding to the negative terminal node
                for(int j = 0; j<noden; j++)
                {
                    //at the column corresponding to the negative terminal, write 1
                    if(j == (nB(comps[i]) -1))
                    {
                        conducts (j, j) = 1;
                        //also write -1* the source voltage to this index in the rhs vector
                        currents[j] = (-1)*comps[i].value; 
                    } else {
                        //other columns are set to 0
                        conducts (nB(comps[i]) -1, j) = 0;
                    }
                }
            }

            //floating voltage source
            if(nA(comps[i]) != 0 && nB(comps[i]) != 0)
            {
                //assign the row corresponding to the lowest numbered node as that representing the voltage source
                int row;

                if( nA(comps[i]) > nB(comps[i]) )
                {
                    locked[nB(comps[i])-1] = 1;
                    row = nB(comps[i]) -1;
                } else {
                    locked[nA(comps[i])-1] = 1;
                    row = nA(comps[i]) -1;
                }

                //in that row of the rhs vector, write the current source value
                currents[row] = comps[i].value; 

                //write the 1, -1 and 0s in the appropriate columns
                for(int j = 0; j<noden; j++)
                {
                    if(j == (nA(comps[i]) -1))
                    {
                        conducts (row, j) = 1;
                    } else {
                        if(j == (nB(comps[i]) -1))
                        {
                            conducts (row, j) = (-1);
                        } else {
                            conducts (row, j) = 0;
                        }
                    }
                }
            }

            //short circuit
            if(nA(comps[i]) == 0 && nB(comps[i]) == 0)
            {
                cout << "bruh";
                assert(0);
            }
            
        }

        if(comps[i].type == 'I')
        {
            if(nA(comps[i]) != 0)
            {
                currents[nA(comps[i]) -1] += comps[i].value;
            }

            if(nB(comps[i]) != 0)
            {
                currents[nB(comps[i]) -1] -= comps[i].value;
            }
        }

        //test(noden, conducts, currents);
    }

    return {conducts, currents};

}
