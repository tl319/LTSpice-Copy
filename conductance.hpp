#include <iostream>
#include <string>
#include <vector>
#include <numeric>
#include <Eigen>

using namespace Eigen;
using namespace std;

struct Component
{
    char type;
    string name;
    int nA;
    int nB;
    float value;
};

//calculate number of non-ground nodes from node vector
int compute_noden(vector<int> nodes);

//return complete conductance matrix and current vector
pair<MatrixXd, vector<float>> conductance_current (vector<Component> comps, int noden);

//print out conductance matrix and current vector
void test(int noden, MatrixXd conducts, vector<float> currents);

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

int compute_noden(vector<int> nodes)
{
    //number of non-ground nodes
    int noden = nodes.size();
    //loop redundant if vector is ordered
    for(int i = 0; i<nodes.size(); i++)
    {
        if(nodes[i] == 0)
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

            if(locked[comps[i].nA-1] == 0)
            {
                if(comps[i].nA != 0)
                {
                    conducts (comps[i].nA -1, comps[i].nA -1) += conductance;
                }
                if(comps[i].nA != 0 && comps[i].nB != 0)
                {
                    conducts (comps[i].nA -1, comps[i].nB -1) -= conductance;
                }
            }

            if(locked[comps[i].nB-1] == 0)
            {
                if(comps[i].nB != 0)
                {
                    conducts (comps[i].nB -1, comps[i].nB -1) += conductance;
                }
                if(comps[i].nA != 0 && comps[i].nB != 0)
                {
                    conducts (comps[i].nB -1, comps[i].nA -1) -= conductance;
                }
            }

        }

//add case for nA = nB = 0?
        if(comps[i].type == 'V')
        {
            if(comps[i].nB == 0)
            {
                locked[comps[i].nA-1] = 1;
                for(int j = 0; j<noden; j++)
                {
                    if(j == (comps[i].nA -1))
                    {
                        conducts (j, j) = 1;
                        currents[j] = comps[i].value; 
                    } else {
                        conducts (comps[i].nA -1, j) = 0;
                    }
                }
            }

            if(comps[i].nA == 0)
            {
                locked[comps[i].nB-1] = 1;
                for(int j = 0; j<noden; j++)
                {
                    if(j == (comps[i].nB -1))
                    {
                        conducts (j, j) = 1;
                        currents[j] = (-1)*comps[i].value; 
                    } else {
                        conducts (comps[i].nB -1, j) = 0;
                    }
                }
            }

            /*/
            if(comps[i].nA != 0 && comps[i].nB != 0)
            {
                locked[comps[i].nA-1] = 1;
                for(int j = 0; j<noden; j++)
                {
                    if(j == (comps[i].nA -1))
                    {
                        conducts (comps[i].nA -1, j) = 1;
                        currents[j] = comps[i].value; 
                    } else {
                        conducts (comps[i].nA -1, j) = 0;
                    }
                }
            }
            /*/
        }

            //test(noden, conducts, currents);
    }

    return {conducts, currents};

}
