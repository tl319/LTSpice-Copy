#include "Circut_Simulator.hpp"

using namespace std;
using namespace Eigen;

//!!
//remember to change all arguments to const refs when possible
//!!
//"pre-compute" values that are computed more than once per scope
//!!

VectorXd VectorUpdate (vector<Component> comps, int noden, float time, VectorXd pastnodes, VectorXd comp_currents, float interval, vector<int> c_vs_row)
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
                currents(nA(comps[i]) -1) = (comp_currents(i)) * interval/val; 
            }

            //positive terminal to ground
            //not too sure about this one
            if(nA(comps[i]) == 0 && nB(comps[i]) != 0)
            {
                //also write -1* the source voltage to this index in the rhs vector
                currents(nB(comps[i]) -1) = (-1)*(comp_currents(i)) * interval/val;
            }

            //floating voltage source
            //to optimise, make the row value saved after initial calculation
            if(nA(comps[i]) != 0 && nB(comps[i]) != 0)
            {
                row = c_vs_row[i];

                //in that row of the rhs vector, write the current source value
                currents(row) = (comp_currents(i)) * interval/val; 
            }    
        }
    }    
    return currents;
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

//used to compute current through voltage sources or capacitors
//perhaps optimimse by replacing N with an int & by "looking both ways"
float vs_current (vector<Component> comps, Component C, vector<bool> & computed, VectorXd currents, Node N)
{
    vector<Component> shared_node;
    //vector<Component> shared_node_B;
    float total_current = 0;

    Node checked_node;

    //to avoid infinite recursion, consistently check currents into "the other node"
    if( N.number == nA(C) )
    {
        checked_node = C.B;
    } else {
        assert(N.number == nB(C));
        checked_node = C.A;
    }
    shared_node = common_node(comps, C, checked_node);
    
    //sum up the currents entering node A of the component
    for(int j = 0; j<shared_node.size(); j++)
    {
        //check that the current has been computed
        if( computed[component_index( comps, shared_node[j] )] == 1 )
        {
            total_current += currents( component_index( comps, shared_node[j] ) );
        } else {
            total_current += vs_current (comps, shared_node[j], computed, currents, checked_node);
        }
    }
    //checking known currents at node B
    /*/
    if(cv_check == 0)
    {
        cv_check = 1;
        total_current = 0;
        shared_node_B = common_node(comps, C, C.B);
        for(int j = 0; j<shared_node_B.size(); j++)
        {
            if( computed[component_index( comps, shared_node_B[j] )] == 1 )
            {
                total_current += currents( component_index( comps, shared_node_B[j] ) );
            } else {
                cv_check = 0;
            }
        }
    }
    /*/
    //recurrent calls

    computed[component_index(comps, C)] = 1;    
    return total_current;
}

//compute currents accross each component, set interval to 0 if at start of transient or operating point
VectorXd comp_currents (vector<Component> comps, vector<Node> nlist, VectorXd nodev, float interval)
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

    for(int i = 0; i<comps.size(); i++)
    {
        //recursively obtaining currents through components functioning as voltage sources (i.e. other currents into one of their nodes)
        if(comps[i].type == 'C' || comps[i].type == 'V')
        {
            //pass in node A by default, perhaps possible to optimise
            currents(i) = (-1) * vs_current(comps, comps[i], computed, currents, comps[i].A);
        }
    }
    return currents;
}

//after calculating voltages and currents assuming all NL components are active, determine which actually are
vector<bool> incorrect_assumptions(VectorXd comp_currents, vector<components> comps)
{
    vector<bool> incorrect_assumptions (comps.size(), 0);

    for(int i = 0; i< comps.size(); i++)
    {
        if(comps[i].type == 'D' && comp_currents[i] < 0)
        {
            incorrect_assumptions[i] = 1l
        }
    }

    return incorrect_assumptions;
}