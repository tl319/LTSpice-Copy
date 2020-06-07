#include "Circut_Simulator.hpp"

using namespace std;
using namespace Eigen;

//!!
//remember to change all arguments to const refs when possible
//!!
//"pre-compute" values that are computed more than once per scope
//!!

VectorXd VectorUpdate (vector<Component> comps, int noden, float time, VectorXd pastnodes, VectorXd component_currents, float interval, vector<int> c_vs_row, vector<bool> incorrect_assumptions)
{
    //assign the row corresponding to the lowest numbered node as that representing the voltage source
    //for floating voltage sources
    int row;

    //rhs vector
    VectorXd currents = VectorXd::Zero (noden);

    //component value
    float val;
    for(int i = 0; i<comps.size(); i++)
    {
        if(comps[i].type == 'V' && comps[i].isSignal == 1)
        {
            val = comps[i].DCOff + (comps[i].amplitude)*sin( (comps[i].frequency)*2*M_PI*time );
        } else {
            val = comps[i].value;
        }

        //dealing with voltage sources              
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
                row = c_vs_row[i];

                //in that row of the rhs vector, write the current source value
                currents(row) = val; 
            }    
        }

        //current sources
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

        //optimise by taking values from current vector
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

        if(comps[i].type == 'C')
        {
            //negative terminal to ground
            if( nB(comps[i]) == 0 && nA(comps[i]) != 0)
            {
                //write the source voltage to this index in the rhs vector
                currents(nA(comps[i]) -1) = (component_currents(i)) * interval/val; 
            }

            //positive terminal to ground
            //not too sure about this one
            if(nA(comps[i]) == 0 && nB(comps[i]) != 0)
            {
                //also write -1* the source voltage to this index in the rhs vector
                currents(nB(comps[i]) -1) = (-1)*(component_currents(i)) * interval/val;
            }

            //floating voltage source
            //to optimise, make the row value saved after initial calculation
            if(nA(comps[i]) != 0 && nB(comps[i]) != 0)
            {
                row = c_vs_row[i];

                //in that row of the rhs vector, write the current source value
                currents(row) = (component_currents(i)) * interval/val; 
            }    
        }

        if(comps[i].type == 'D' && incorrect_assumptions[i] == 0)
        {
                row = c_vs_row[i];

                //in that row of the rhs vector, write 0.7
                currents(row) = 0.7; 
        }

    }    
    return currents;
}

//after calculating voltages and currents assuming all NL components are active, determine which actually are
vector<bool> incorrect_assumptions(VectorXd component_currents, vector<Component> comps)
{
    vector<bool> incorrect_assumptions (comps.size()+1, 0);

    for(int i = 0; i< comps.size(); i++)
    {
        //cout << comps[i].name << " " << comp_currents[i] << endl;
        if(comps[i].type == 'D' && component_currents[i] < 0)
        {
            //cout << "oof" << endl;
            incorrect_assumptions[comps.size()] = 1;
            incorrect_assumptions[i] = 1;
        }
    }

    return incorrect_assumptions;
}

//return voltage and current vector for correct nonlinear modes
pair<VectorXd, VectorXd> adjust_modes(MatrixXd lhs, VectorXd rhs, const vector<Component> & comps, const vector<Node> & nodes) 
{
    vector<bool> oof;
    VectorXd nodev;
    pair<VectorXd, vector<bool>> component_currents;
    VectorXd new_rhs;

    //probably should take as argument
    int noden = compute_noden(nodes);

    pair<MatrixXd, vector<int>> new_mat;

    nodev = matrixSolve(lhs, rhs);

    //what happens with the increment in operating point?
    //component_currents = comp_currents(comps, nodes, nodev, 1);

    component_currents = sufficient_currents (comps, nodes, nodev, 1);

    series_currents (comps, nodes, nodev, 1, component_currents.second, component_currents.first);

    oof = incorrect_assumptions(component_currents.first, comps);

    while(oof[comps.size()] == 1)
    {
        cout << "while ";
        new_mat = CorrectAssumptions (comps, noden, oof);
        new_rhs = VectorUpdate (comps, noden, 1, nodev, component_currents.first, 1, new_mat.second, oof);
        nodev = matrixSolve(new_mat.first, new_rhs);
        //component_currents = comp_currents(comps, nodes, nodev, 1);
        component_currents = sufficient_currents (comps, nodes, nodev, 1);
        series_currents (comps, nodes, nodev, 1, component_currents.second, component_currents.first);
        oof = incorrect_assumptions(component_currents.first, comps);
    }

    return{nodev, component_currents.first};
}

//compute currents accross each "sufficient" (R, L, I) component, ie currents that can be directly computed without othe currents
pair<VectorXd, vector<bool>> sufficient_currents (vector<Component> comps, vector<Node> nlist, VectorXd nodev, float interval)
{
    //register components take care of, to differetiate non calculated values from 0 currents
    vector<bool> computed (comps.size(), 0);

    VectorXd currents = VectorXd::Zero (comps.size());

    //voltages at each node
    float VA;
    float VB;

    for(int i = 0; i<comps.size(); i++)
    {
        if(nA(comps[i]) != 0)
        {
            VA = nodev(nA(comps[i]) - 1);
        } else {
            VA = 0;
        }

        if(nB(comps[i]) != 0)
        {
            VB = nodev(nB(comps[i]) - 1);
        } else {
            VB = 0;
        }

        //straightforward inclusion of known or calculated currents
        if(comps[i].type == 'R')
        {
            currents(i) = ( VA - VB )/comps[i].value;
            computed[i] = 1;
        }
        if(comps[i].type == 'I')
        {
            currents(i) = comps[i].value;
            computed[i] = 1;
        }
        if(comps[i].type == 'L')
        {
            currents(i) = ( VA - VB )*interval/comps[i].value;
            computed[i] = 1;
        }
    }
    return {currents, computed};
}

//used in vs_current to determine other components connected to C at A
vector<Component> common_node (vector<Component> comps, Component C, Node A)
{
    vector<Component> shared_node;

    for(int i = 0; i<comps.size(); i++)
    {
        cout << "for ";
        if( ( nA(comps[i]) == A.number || nB(comps[i]) == A.number ) && comps[i].name != C.name )
        {
            shared_node.push_back(comps[i]);
        }
    }
    cout << endl;
    return shared_node;
}

//returns the index of a component in  the component vector
int component_index (vector<Component> comps, Component C)
{
    for(int i = 0; i<comps.size(); i++)
    {
        if(comps[i].name == C.name)
        {
            return i;
        }
    }
}

//compute currents accross "insufficient" (V, C, D) components in series with sufficient components (ie traversed by a known current)
void series_currents (vector<Component> comps, vector<Node> nlist, VectorXd nodev, float interval, vector<bool> & computed, VectorXd & comp_currents)
{
    vector<Component> same_node;
    bool acceptable = 1;
    float current = 0;
    Node used_node;

    for(int i = 0; i<comps.size(); i++)
    {
        if(comps[i].type == 'V' || comps[i].type == 'C' || comps[i].type == 'D')
        {
            same_node = common_node(comps, comps[i], comps[i].A);
            used_node = comps[i].A;
            for(int j = 0; j<same_node.size(); j++)
            {
                if(same_node[j].type == 'V' || same_node[j].type == 'C' || same_node[j].type == 'D')
                {
                    acceptable = 0;
                    break;
                }
            }

            if(acceptable == 0)
            {
                same_node = common_node(comps, comps[i], comps[i].B);
                used_node = comps[i].B;
                for(int j = 0; j<same_node.size(); j++)
                {
                    if(same_node[j].type == 'V' || same_node[j].type == 'C' || same_node[j].type == 'D')
                    {
                        acceptable = 0;
                        break;
                    }
                }
            }

            if(acceptable == 1)
            {
                for(int j = 0; j<same_node.size(); j++)
                {
                    if(same_node[j].type == 'I')
                    {
                        if(used_node.label == same_node[i].A.label)
                        {
                            current -= same_node[i].value;
                        } else {
                            current += same_node[i].value;
                        }
                    }

                    //do inductors

                    if(same_node[j].type == 'R')
                    {
                        current -= abs( comp_currents( component_index(comps, same_node[j]) ) );
                    }
                    cout << comps[i].name << " " << same_node[j].name << " " << current << endl;
                }

                comp_currents( i ) = current;
                current = 0;
            }
            acceptable = 1;
        }
    } 
}