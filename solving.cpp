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
            cout << "oof" << endl;
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

    recursive_currents (comps, nodes, nodev, 1, component_currents.second, component_currents.first);

    cout << "D I " << component_currents.first(2) << endl;

    oof = incorrect_assumptions(component_currents.first, comps);

    while(oof[comps.size()] == 1)
    {
        new_mat = CorrectAssumptions (comps, noden, oof);
        new_rhs = VectorUpdate (comps, noden, 1, nodev, component_currents.first, 1, new_mat.second, oof);
        nodev = matrixSolve(new_mat.first, new_rhs);
        //component_currents = comp_currents(comps, nodes, nodev, 1);
        component_currents = sufficient_currents (comps, nodes, nodev, 1);
        series_currents (comps, nodes, nodev, 1, component_currents.second, component_currents.first);
        recursive_currents (comps, nodes, nodev, 1, component_currents.second, component_currents.first);
        oof = incorrect_assumptions(component_currents.first, comps);
    }

    return{nodev, component_currents.first};
}

//collapse all 3 main current functions into 1

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
            if(VA > 0)
            {
                currents(i) = ( VA - VB )/comps[i].value;
                computed[i] = 1;
            } else {
                currents(i) = ( VB - VA )/comps[i].value;
                computed[i] = 1;
            }

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
        if( ( nA(comps[i]) == A.number || nB(comps[i]) == A.number ) && comps[i].name != C.name )
        {
            shared_node.push_back(comps[i]);
        }
    }
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


//probably some stuff to trim
void pseudo_basecase(Component comp, vector<Component> comps, vector<Node> nlist, VectorXd nodev, float interval, vector<bool> & computed, VectorXd & comp_currents)
{
    if(comp.type == 'V' || comp.type == 'C' || comp.type == 'D')
    {
        vector<Component> same_node;
        bool acceptable = 1;
        float total_current = 0;
        float current = 0;
        Node used_node;

        same_node = common_node(comps, comp, comp.A);
        used_node = comp.A;
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
            acceptable = 1;
            same_node = common_node(comps, comp, comp.B);
            used_node = comp.B;
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
                    if(used_node.label == same_node[j].A.label)
                    {
                        current = (-1)*same_node[j].value;
                    } else {
                        current = same_node[j].value;
                    }

                    if(used_node.label == comp.A)
                    {
                        total_current += current;
                    } else {
                        total_current -= current;
                    }
                }

                //do inductors

                if(same_node[j].type == 'R')
                {
                    //why does this work? come back later
                    //also fix for some cases
                    total_current -= comp_currents( component_index(comps, same_node[j]) );
                }
            }

            if(comp.type == 'D')
            {
                total_current = total_current*(-1);
            }

            comp_currents( component_index(comps, comp) ) = total_current;
            computed[component_index(comps, comp)] = 1;
            total_current = 0;
        }
        acceptable = 1;
    }
}

//compute currents accross "insufficient" (V, C, D) components in series with sufficient components (ie traversed by a known current)
void series_currents (vector<Component> comps, vector<Node> nlist, VectorXd nodev, float interval, vector<bool> & computed, VectorXd & comp_currents)
{
    for(int i = 0; i<comps.size(); i++)
    {
        pseudo_basecase(comps[i], comps, nlist,nodev, interval, computed, comp_currents);
    } 
}

float recursive_basecase (Component C, vector<Component> comps, vector<Node> nlist, VectorXd nodev, float interval, vector<bool> & computed, VectorXd & comp_currents)
{
    vector<Component> same_node;
    bool acceptable = 1;
    float total_current;
    float current = 0;
    Node used_node;

    if(C.type == 'V' || C.type == 'C' || C.type == 'D')
    {
        same_node = common_node(comps, C, C.A);
        used_node = C.A;

        for(int j = 0; j<same_node.size(); j++)
        {
            total_current = 0;
            if( computed[component_index(comps, same_node[j])] == 1 )
            {
                if(same_node[j].type == 'R')
                {
                    //why does this work? come back later
                    //also fix for some cases
                    total_current -= comp_currents( component_index(comps, same_node[j]) );
                }

                if(same_node[j].type == 'I')
                {
                    if(used_node.label == same_node[j].A.label)
                    {
                        current = (-1)*same_node[j].value;
                    } else {
                        current = same_node[j].value;
                    }

                    if(used_node.label == C.A)
                    {
                        total_current += current;
                    } else {
                        total_current -= current;
                    }
                }

                if(same_node[j].type == 'V' || same_node[j].type == 'C' || same_node[j].type == 'D')
                {
                    if(used_node.label == same_node[j].A.label)
                    {
                        current = (-1)*comp_currents( component_index(comps, same_node[j]) );
                    } else {
                        current = comp_currents( component_index(comps, same_node[j]) );
                    }

                    if(used_node.label == C.A)
                    {
                        total_current += current;
                    } else {
                        total_current -= current;
                    }
                }
            } else {
                if(used_node.label == same_node[j].A.label)
                {
                    current = (-1)*recursive_basecase (same_node[j], comps, nlist, nodev, interval, computed, comp_currents);
                } else {
                    current = recursive_basecase (same_node[j], comps, nlist, nodev, interval, computed, comp_currents);
                }

                if(used_node.label == C.A)
                {
                    total_current += current;
                } else {
                    total_current -= current;
                }
            }
        }
    }
    if(C.type == 'D')
    {
        total_current = total_current*(-1);
    }
    return total_current;
}

//compute currents accross "insufficient" (V, C, D) components in series with other such components
void recursive_currents (vector<Component> comps, vector<Node> nlist, VectorXd nodev, float interval, vector<bool> & computed, VectorXd & comp_currents)
{
    for(int i = 0; i<comps.size(); i++)
    {
        if( computed[i] == 0 )
        {
            comp_currents( i ) = recursive_basecase (comps[i], comps, nlist, nodev, interval, computed, comp_currents);
            computed[component_index(comps, comps[i])] = 1;
        }
    }
}