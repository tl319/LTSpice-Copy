#include "Circut_Simulator.hpp"
#include <iostream>
#include <string>
#include <vector>
#include <numeric>
#include <Eigen>

using namespace Eigen;
using namespace std;

pair<MatrixXd, vector<int>> MatrixUpdate (const vector<Component> & comps, const int & noden);

//function declarations

void testBetter(const int & noden, const MatrixXd & conducts, const VectorXd & currents, const vector<Node> & nodes)
{
    //print each line of the matrix
    cerr << conducts << "      *      " << endl;

    for(auto x : nodes)
    {
        cerr << "("<< x.label << ")        =      " << endl;
    };
    //print current vector
    cerr <<  currents << endl;
}

int nA(const Component & c)
{
    return (c.A.number);
}

int nB(const Component & c)
{
    return (c.B.number);
}

int SnA(const Component & c)
{
    return (c.A.super);
}

int SnB(const Component & c)
{
    return (c.B.super);
}

//compute number of non-ground nodes
int compute_noden(const vector<Node> & nodes)
{
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

pair<MatrixXd, VectorXd> conductance_current(const vector<Component> & comps, const int & noden)
{
    //value to take from sources
    float val;

    //conductance matrix
    MatrixXd conducts = MatrixXd::Zero (noden, noden);
    //rhs vector
    VectorXd currents = VectorXd::Zero (noden);

    float conductance;

    //semiconductor thermal voltage
    float vt = 0.025;

    //"locked" matrix lines so that if they're written to by a v source they aren't written to again
    vector<bool> locked (noden, 0);

    //cycle through each component
    for(int i = 0; i<comps.size(); i++)
    {

        //value is either component value member or signal at t = 0
        if( (comps[i].isSignal == 1) && (comps[i].type == 'V') )
        {
            val = comps[i].DCOff;
            //cerr << "DCoff= " << comps[i].DCOff << " val= " << val << endl;
        } else {
            val = comps[i].value;
        }

        //dealing with resistors
        if(comps[i].type == 'R')
        {
            conductance = 1/val;

            //node numbers that only take the supernode value if the supernode exists in DC operation (isnt from a capacitor)
            int DCsuperA;
            int DCsuperB;
            if(comps[i].A.reactiveSuper == true)
            {
                DCsuperA = nA(comps[i]);
            } else {
                DCsuperA = SnA(comps[i]);
            }

            if(comps[i].B.reactiveSuper == true)
            {
                DCsuperB = nB(comps[i]);
            } else {
                DCsuperB = SnB(comps[i]);
            }

            //if row nA has not already been edited to represent a voltage source and nA is not ground,
            //then add and subtract conductance in columns nA and nB
            //Row corresponds to supernode (node if no supernode) and column to actual node
            if(DCsuperA != 0)
            {
                if(locked[ DCsuperA-1] == 0)
                {
                    if( nA(comps[i]) != 0)
                    {
                        conducts ( DCsuperA -1, nA(comps[i]) -1) += conductance;
                    }
                    if( nA(comps[i]) != 0 && nB(comps[i]) != 0)
                    {
                        conducts ( DCsuperA -1, nB(comps[i]) -1) -= conductance;
                    }
                }
            }

            if( DCsuperB != 0)
            {
                //Ditto for nB
                if(locked[ DCsuperB -1] == 0)
                {
                    if( nB(comps[i]) != 0)
                    {
                        conducts ( DCsuperB  -1, nB(comps[i]) -1) += conductance;
                    }
                    if( nA(comps[i]) != 0 && nB(comps[i]) != 0)
                    {
                        conducts ( DCsuperB -1, nA(comps[i]) -1) -= conductance;
                    }
                }
            }
        }

        //dealing with voltage sources
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
                        currents(j) = val;
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
                        currents(j) = (-1)*val;
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

                if( nA(comps[i]) > nB(comps[i]) && (locked[ nB(comps[i])-1] == 0)|| locked[ nA(comps[i])-1] == 1 )
                {
                    locked[nB(comps[i])-1] = 1;
                    //cerr << "locked " << nB(comps[i])-1 <<endl;
                    row = nB(comps[i]) -1;
                } else {
                    locked[nA(comps[i])-1] = 1;
                    row = nA(comps[i]) -1;
                    //cerr << "locked " << nA(comps[i])-1 <<endl;
                }

                //in that row of the rhs vector, write the current source value
                currents(row) = val;

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
            if(nA(comps[i]) ==  nB(comps[i]))
            {
                cout << "short circuit";
                assert(0);
            }
        }

        //current sources
        if(comps[i].type == 'I')
        {
            if(nA(comps[i]) != 0 && locked[nA(comps[i])-1] == 0)
            {
                currents(nA(comps[i]) -1) -= val;
            }

            if(nB(comps[i]) != 0 && locked[nB(comps[i])-1] == 0)
            {
                currents(nB(comps[i]) -1) += val;
            }
        }

        //capacitors are open circuits at t = 0/.op so are not processed
        //inductors
       if(comps[i].type == 'L')
        {
            //LTSpice treats inductors as 1mOhm resistors in DC coditions
            conductance = 1000;

            if(SnA(comps[i]) != 0)
            {
                if(locked[ SnA(comps[i])-1] == 0)
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
            }

            if(SnB(comps[i]) != 0)
            {
                //Ditto for nB
                if(locked[ SnB(comps[i])-1] == 0)
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
        }

        //cerr << comps[i].name << endl;
        //testBetter(noden, conducts, currents, findNodes(comps));
    }
    //cout << "help" << endl;
    return {conducts, currents};
}

//update matrix to reflect the change in behaviour of reactive components
//the second pair member is used to store the row values representing source and capacitor voltages
pair<MatrixXd, vector<int>> MatrixUpdate (vector<Component> & comps, const int & noden, const float & interval)
{
    //code copied wholesale from above, trim redundancies later

    //value to take from sources
    float val;

    //conductance matrix
    MatrixXd conducts = MatrixXd::Zero (noden, noden);

    float conductance;

    //mechanism to ensure that the 1/0/-1 and voltage values are written to the same row
    //probably need smth similar for current sources (no?)
    vector<int> c_vs_row (comps.size());

    //"locked" matrix lines so that if they're written to by a v source they aren't written to again
    vector<bool> locked (noden, 0);

    //used to "flip" inductors
    Node fakeswitch;

    for(int i = 0; i<comps.size(); i++)
    {
        val = comps[i].value;

        //dealing with voltage sources
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

                if( nA(comps[i]) > nB(comps[i]) && locked[ nB(comps[i])-1] == 0 || locked[ nA(comps[i])-1] == 1 )
                {
                    locked[nB(comps[i])-1] = 1;
                    row = nB(comps[i]) -1;
                } else {
                    locked[nA(comps[i])-1] = 1;
                    row = nA(comps[i]) -1;
                }

                c_vs_row[i] = row;//store the row used to represent the voltage source

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
        }

        //dealing with resistors
        if(comps[i].type == 'R')
        {
            conductance = 1/val;

            //if row nA has not already been edited to represent a voltage source and nA is not ground,
            //then add and subtract conductance in columns nA and nB
            //Row corresponds to supernode (node if no supernode) and column to actual node
            if(SnA(comps[i]) != 0)
            {
                if(locked[ SnA(comps[i])-1] == 0)
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
            }

            if(SnB(comps[i]) != 0)
            {
                //Ditto for nB
                if(locked[ SnB(comps[i])-1] == 0)
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
        }

        if(comps[i].type == 'L')
        {
            //rearrange nodes to correct inductor current sign
            fakeswitch = comps[i].A;
            comps[i].A = comps[i].B;
            comps[i].B = fakeswitch;

            conductance = interval/(comps[i].value);

            //if row nA has not already been edited to represent a voltage source and nA is not ground,
            //then add and subtract conductance in columns nA and nB
            //Row corresponds to supernode (node if no supernode) and column to actual node
            if(SnA(comps[i]) != 0)
            {
                if(locked[ SnA(comps[i])-1] == 0)
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
            }

            if(SnB(comps[i]) != 0)
            {
                //Ditto for nB
                if(locked[ SnB(comps[i])-1] == 0)
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
        }

        if(comps[i].type == 'C')
        {

            //B to ground
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
                    } else {
                        //other columns are set to 0
                        conducts (nA(comps[i]) -1, j) = 0;
                    }
                }
            }

            //A to ground
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

                if( nA(comps[i]) > nB(comps[i]) && locked[ nB(comps[i])-1] == 0 || locked[ nA(comps[i])-1] == 1 )
                {
                    locked[nB(comps[i])-1] = 1;
                    row = nB(comps[i]) -1;
                } else {
                    locked[nA(comps[i])-1] = 1;
                    row = nA(comps[i]) -1;
                }
                //store the row used to represent the capacitor voltage
                c_vs_row[i] = row;

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
        }
        //cerr << "tran " << comps[i].name << endl;
        //cerr << conducts << endl;
    }
    return {conducts, c_vs_row};
}
