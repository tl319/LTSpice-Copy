#include <iostream>
#include <string>
#include <vector>
#include <numeric>
#include <Eigen>

using namespace Eigen;
using namespace std;

int main()
{
    struct Component
    {
	    char type;
	    string name;
	    int nA;
	    int nB;
	    float value;
    };

    Component R1 {'R', "R1", 1, 4, 1};
    Component R2 {'R', "R2", 1, 2, 2};
    Component R3 {'R', "R3", 2, 3, 3};
    Component R4 {'R', "R4", 3, 4, 4};

    vector<Component> comps {R1, R2, R3, R4};

    //MatrixXd populate_conductance_matrix (vector<Compenent> comps)
    //{
        //number of non-ground nodes
        //only for testing
        int noden = 4;
        //conductance matrix
        MatrixXd conducts = MatrixXd::Zero (noden, noden);
        //lhs vector
        vector<int> nodes {1, 2, 3, 4};
        //rhs vector
        vector<float> currents(noden);

        //used to tell where the first empty spot in the node vector is
        int nodes_added = 0;
        //used to tell if an empty spot must be written to
        bool found = 0;

        float conductance;

        //indices of the two component nodes in the node vector
        int Aval;
        int Bval;

        //cycle through each component
        for(int i = 0; i<comps.size(); i++)
        {
            Aval = 0;
            Bval = 0;

            //dealing with resistors
            if(comps[i].type == 'R')
            {
                conductance = 1/comps[i].value;

                conducts (comps[i].nA -1, comps[i].nA -1) += conductance;
                conducts (comps[i].nA -1, comps[i].nB -1) -= conductance;

                conducts (comps[i].nB -1, comps[i].nB -1) += conductance;
                conducts (comps[i].nB -1, comps[i].nA -1) -= conductance;
            }

            if(comps[i].type == 'V')
            {
                
            }
        }

    //}

    //testing
    for(int i = 0; i<noden; i++)
    {
        for(int j = 0; j<noden; j++)
        {
            cout << conducts(i, j) << " ";
        }
        cout << endl;
    }

}