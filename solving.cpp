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
    VectorXd component_currents;
    VectorXd new_rhs;

    //probably should take as argument
    int noden = compute_noden(nodes);

    pair<MatrixXd, vector<int>> new_mat;

    nodev = matrixSolve(lhs, rhs);

    //what happens with the increment in operating point?
    //component_currents = comp_currents(comps, nodes, nodev, 1);

    component_currents = recursive_currents (comps, nodes, nodev, 1);

    cout << "D I " << component_currents(2) << endl;

    oof = incorrect_assumptions(component_currents, comps);

    while(oof[comps.size()] == 1)
    {
        new_mat = CorrectAssumptions (comps, noden, oof);
        new_rhs = VectorUpdate (comps, noden, 1, nodev, component_currents, 1, new_mat.second, oof);
        nodev = matrixSolve(new_mat.first, new_rhs);
        //component_currents = comp_currents(comps, nodes, nodev, 1);
        component_currents = recursive_currents (comps, nodes, nodev, 1);
        oof = incorrect_assumptions(component_currents, comps);
    }

    return{nodev, component_currents};
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

//compute currents accross "insufficient" (V, C, D) components in series with other such components
VectorXd recursive_currents (vector<Component> comps, vector<Node> nlist, VectorXd nodev, float interval)
{
    //register components take care of, to differetiate non calculated values from 0 currents
    vector<bool> computed (comps.size(), 0);

    VectorXd comp_currents = VectorXd::Zero (comps.size());

    for(int i = 0; i<comps.size(); i++)
    {
        if( computed[i] == 0 )
        {
            comp_currents( i ) = recursive_basecase (i, comps[i], comps, nlist, nodev, interval, computed, comp_currents);
            computed[component_index(comps, comps[i])] = 1;
        }
    }
    return comp_currents;
}

float recursive_basecase (int i, const Component & C, const vector<Component> & comps, const vector<Node> & nlist, VectorXd nodev, const float & interval, vector<bool> & computed, VectorXd & comp_currents)
{
    vector<Component> A_node;
    vector<Component> B_node;
    vector<Component> same_node;
    bool acceptable = 1;
    float total_current;
    float current = 0;
    Node used_node;

    //voltages at each node
    float VA;
    float VB;
    
    if(nA(C) != 0)
    {
        VA = nodev(nA(C) - 1);
    } else {
        VA = 0;
    }

    if(nB(C) != 0)
    {
        VB = nodev(nB(C) - 1);
    } else {
        VB = 0;
    }

    if(C.type == 'I')
    {
        total_current = C.value;
        computed[i] = 1;
    }

    if(C.type == 'R' || C.type == 'L')
    {

        //straightforward inclusion of known or calculated currents
        if(C.type == 'R')
        {
            if(VA > 0)
            {
                total_current = ( VA - VB )/C.value;
                computed[i] = 1;
            } else {
                total_current = ( VB - VA )/C.value;
                computed[i] = 1;
            }

        }

        if(C.type == 'L')
        {
            total_current = ( VA - VB )*interval/C.value;
            computed[i] = 1;
        }
    }

    //probably can replace some common_node calls with i
    if(C.type == 'V' || C.type == 'C' || C.type == 'D')
    {
        total_current = 0;
        //does this actually save time?
        A_node = common_node(comps, C, C.A);
        B_node = common_node(comps, C, C.B);
        if(A_node.size() < B_node.size())
        {
            used_node = C.A;
            same_node = A_node;
        } else {
            used_node = C.B;
            same_node = B_node;
        }

        cout << C.name << " same_node size " << same_node.size() << endl;

        for(int j = 0; j<same_node.size(); j++)
        {
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
                    current = (-1)*recursive_basecase (i, same_node[j], comps, nlist, nodev, interval, computed, comp_currents);
                } else {
                    current = recursive_basecase (i, same_node[j], comps, nlist, nodev, interval, computed, comp_currents);
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