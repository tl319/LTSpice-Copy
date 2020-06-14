#include "Circut_Simulator.hpp"

using namespace std;
using namespace Eigen;
using namespace std::chrono;

//update the current vector to reflect reactive component behaviour and changes in source values
VectorXd VectorUpdate (const vector<Component> & comps, const int & noden, const float & time, const VectorXd & pastnodes, const VectorXd & component_currents,
const float & interval, const vector<int> & c_vs_row)
{
    //assign the row corresponding to the lowest numbered node as that representing the voltage source
    //for floating voltage sources
    int row;

    //rhs vector
    VectorXd currents = VectorXd::Zero (noden);

    //component value
    float val;

    //voltages at each node
    float VA;
    float VB;

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
                currents(nA(comps[i]) -1) -= val;
            }

            if(nB(comps[i]) != 0)
            {
                currents(nB(comps[i]) -1) += val;
            }
        }

        //taking values from current vector causes a phase shift
        if(comps[i].type == 'L')
        {
            if(nA(comps[i]) != 0)
            {
                currents(nA(comps[i]) -1) += component_currents(i);
            }

            if(nB(comps[i]) != 0)
            {
                currents(nB(comps[i]) -1) -= component_currents(i);
            }
        }

        if(comps[i].type == 'C')
        {
            //assign node voltage values
            //nA/nB read the node number
            if(nA(comps[i]) != 0)
            {
                VA = pastnodes(nA(comps[i])-1);
            } else {
                VA = 0;
            }

            if(nB(comps[i]) != 0)
            {
                VB = pastnodes(nB(comps[i])-1);
            } else {
                VB = 0;
            }

            //B to ground
            if( nB(comps[i]) == 0 && nA(comps[i]) != 0)
            {
                //write the source voltage to this index in the rhs vector
                currents(nA(comps[i]) -1) = VA - VB - component_currents(i)*interval/val;
            }

            //A to ground
            //not too sure about this one
            if(nA(comps[i]) == 0 && nB(comps[i]) != 0)
            {
                //also write -1* the source voltage to this index in the rhs vector
                currents(nB(comps[i]) -1) += (-1)*(VA - VB - (component_currents(i))*interval/val);
            }

            //not grounded
            //to optimise, make the row value saved after initial calculation
            if(nA(comps[i]) != 0 && nB(comps[i]) != 0)
            {
                row = c_vs_row[i];

                //in that row of the rhs vector, write the current source value
                currents(row) = VA - VB - component_currents(i)*interval/val;
            }
        }
        //cerr << comps[i].name << endl;
        //cerr << currents << endl;
    }
    return currents;
}

//return voltage and current vectors for operating point or first point of transient analysis
pair<VectorXd, VectorXd> no_prior_change (const vector<Component> & comps, const vector<Node> & nodes, const int & noden)
{
    VectorXd nodev;
    //perhaps avoid creating by separating recursive_currents used for op and trans
    VectorXd component_currents = VectorXd::Zero (comps.size());

    //generate conductance matrix and rhs vector
    pair<MatrixXd, VectorXd> knowns = conductance_current (comps, noden);
    //testBetter(noden, knowns.first, knowns.second, nodes);

    //compute node voltages
    nodev = matrixSolve(knowns.first, knowns.second);
    //generate
    VectorXd prevnodev = VectorXd::Zero(nodev.size());

    component_currents = recursive_currents (comps, nodes, nodev, prevnodev, 0, component_currents, true);

    return{nodev, component_currents};
}

//the function is a void as the data from transient timesteps is not needed after it is printed to the output file
void transient (vector<Component> & comps, const vector<Node> & nodes, const int & noden, const float & duration,
const float & interval, const VectorXd & pastnodes, const VectorXd & pastcurrents)
{
    VectorXd nodev = pastnodes;
    VectorXd prevnodev = VectorXd::Zero(nodev.size());
    VectorXd component_currents = pastcurrents;
    vector<pair<VectorXd, VectorXd>> values;
    VectorXd rhs = VectorXd::Zero (comps.size());

    //DGenerate a single conductance matrix
    pair<MatrixXd, vector<int>> Mat = MatrixUpdate (comps, noden, interval);
    //cerr << Mat.first << endl;
    //cerr << endl;

    writeTranHeaders(nodes, comps,pastnodes,pastcurrents);

    //cycle through timesteps 1 to n of the simulation
    for(float i = interval; i<=duration; i += interval)
    {
        rhs = VectorUpdate (comps, noden, i, nodev, component_currents, interval, Mat.second);
        //store previous node voltages for use in backwards Euler
        prevnodev = nodev;
        //calculate node voltages and component currents, then print to cout
        nodev = matrixSolve(Mat.first, rhs);
        component_currents = recursive_currents (comps, nodes, nodev, prevnodev, interval, component_currents, false);
        writeTran(nodes, comps, nodev, component_currents, i);
    }
}


//used in recursive_basecase to determine other components connected to C at A
vector<Component> common_node (const vector<Component> & comps, const Component & C, const Node & A)
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
int component_index (const vector<Component> & comps, const Component & C)
{
    for(int i = 0; i<comps.size(); i++)
    {
        if(comps[i].name == C.name)
        {
            return i;
        }
    }
}

//compute currents accross each component
VectorXd recursive_currents (const vector<Component> & comps, const vector<Node> & nlist, const VectorXd & nodev, const VectorXd & prevnodev,
const float & interval, const VectorXd & past_currents, const bool & op)
{
    VectorXd comp_currents = VectorXd::Zero (comps.size());

    for(int i = 0; i<comps.size(); i++)
    {
        comp_currents( i ) = recursive_basecase (i, comps[i], comps, nlist, nodev, prevnodev, interval, comp_currents, comps[i].A, op);
    }
    return comp_currents;
}

//edits to computed below should be redundant
float recursive_basecase (const int & i, const Component & C, const vector<Component> & comps, const vector<Node> & nlist, const VectorXd & nodev,
const VectorXd & prevnodev, const float & interval, VectorXd & comp_currents, const Node & used_node, const bool & op)
{
    vector<Component> same_node;
    bool acceptable = 1;
    float total_current;
    float current = 0;
    Node rec_used;

    //node (a/0 or b/1) at which we're adding/subtracting the currents into the component
    bool node_add;

    //voltages at each node
    float VA;
    float VB;

    //voltages at each node of a capacitor model resistor
    float CVA;
    float CVB;

    //capacitor model resistor
    Component capr;

    //in DC operation, LT spice treats inductors as 1mOhm resistors
    //computers are bad at dividing by small numbers so we're directly defining the conductance
    float gl = 1000;

    //straightforward inclusion of known or calculated currents
    if(C.type == 'I')
    {
        total_current = C.value;
    }

    if(C.type == 'R' || C.type == 'L')
    {
        //assign node voltage values
        //nA/nB read the node number
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

        //apply Ohm's law
        if(C.type == 'R')
        {
            total_current = ( VA - VB )/C.value;
        }

        //in DC treat inductors as resistors
        if(C.type == 'L' && op == true)
        {
            total_current = ( VA - VB )*gl;
        }

        //in AC apply add the previous current (Ieq) to the current through geq, in accordance with the backwards Euler model
        if(C.type == 'L' && op == false)
        {
            total_current = comp_currents( component_index( comps, C ) )+( prevnodev(nB(C) - 1) - prevnodev(nA(C) - 1) )*interval/(C.value);
        }
    }

    //the current through a capacitor is that through the resistor in the backwards Euler model
    if(C.type == 'C' && op == false)
    {
        capr = comps[ component_index( comps, C ) + 1 ];
        if(nA(capr) != 0)
        {
            VA = nodev(nA(capr) - 1);
        } else {
            VA = 0;
        }

        if(nB(capr) != 0)
        {
            VB = nodev(nB(capr) - 1);
        } else {
            VB = 0;
        }

        total_current = ( VA - VB )/capr.value;
    }

    //recursively call the function on adjacent nodes
    if(C.type == 'V')
    {
        total_current = 0;

        same_node = common_node(comps, C, used_node);

        for(int j = 0; j<same_node.size(); j++)
        {
                //switch the used_node to the other node of same_node[j]
                if(used_node.number == same_node[j].A.number)
                {
                    rec_used = same_node[j].B;
                } else {
                    rec_used = same_node[j].A;
                }

                //add or subtract the recursively obtained current depending on which node connects the comonents, as this affects the sign given to currents

                if(used_node.number == same_node[j].A.number)
                {
                    current = (-1)*recursive_basecase (i, same_node[j], comps, nlist, nodev, prevnodev, interval, comp_currents, rec_used, op);
                } else {
                    current = recursive_basecase (i, same_node[j], comps, nlist, nodev, prevnodev, interval, comp_currents, rec_used, op);
                }

                if(used_node.number == C.A.number)
                {
                    total_current += current;
                } else {
                    total_current -= current;
                }
        }
    }
    return total_current;
}
