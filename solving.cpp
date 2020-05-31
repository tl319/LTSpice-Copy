#include "Circut_Simulator.hpp"

using namespace std;
using namespace Eigen;

VectorXd VectorUpdate (vector<Component> comps, int noden, float time, VectorXd pastnodes, float interval)
{
    //assign the row corresponding to the lowest numbered node as that representing the voltage source
    //for floating voltage sources
    int row;

    float val;

    //rhs vector
    VectorXd currents = VectorXd::Zero (noden);

    for(int i = 0; i<comps.size(); i++)
    {
        if(comps[i].type == 'V' || comps[i].type == 'I')
        {
            val = comps[i].DCOff + (comps[i].amplitude)*sin( (comps[i].frequency)*time );
        } else {
            val = comps[i].value;
        }

        //dealing with voltage sources                  add case for nA = nB = 0?
        if(comps[i].type == 'V')
        {
            //negative terminal to ground
            if( nB(comps[i]) == 0 && nA(comps[i]) != 0)
            {
                //write the source voltage to this index in the rhs vector
                currents(nA(comps[i]) -1) = val; 
            }

            //positive terminal to ground
            if(nA(comps[i]) == 0 && nB(comps[i]) != 0)
            {
                //also write -1* the source voltage to this index in the rhs vector
                currents(nB(comps[i]) -1) = (-1)*val; 
            }

            //floating voltage source
            //to optimise, make the row value saved after initial calculation
            if(nA(comps[i]) != 0 && nB(comps[i]) != 0)
            {
                if( nA(comps[i]) > nB(comps[i]) )
                {
                    //locked[nB(comps[i])-1] = 1;
                    row = nB(comps[i]) -1;
                } else {
                    //locked[nA(comps[i])-1] = 1;
                    row = nA(comps[i]) -1;
                }
                //in that row of the rhs vector, write the current source value
                currents(row) = val; 
            }    
        }

        //dealing with current sources
        if(comps[i].type == 'I')
        {
            if(nA(comps[i]) != 0)
            {
                currents(nA(comps[i]) -1) += val;
            }

            if(nB(comps[i]) != 0)
            {
                currents(nB(comps[i]) -1) -= val;
            }
        }

        if(comps[i].type == 'L')
        {
            if(nA(comps[i]) != 0)
            {
                currents(nA(comps[i]) -1) += (pastnodes(nA(comps[i]) -1) - pastnodes(nB(comps[i]) -1)) * interval/val;
            }

            if(nB(comps[i]) != 0)
            {
                currents(nB(comps[i]) -1) -= (pastnodes(nA(comps[i]) -1) - pastnodes(nB(comps[i]) -1)) * interval/val;
            }
        }
    }
    
    return currents;
}