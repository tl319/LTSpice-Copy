#include "Circut_Simulator.hpp"

using namespace std;
using namespace Eigen;

//!!
//remember to change all arguments to const refs when possible
//!!
//"pre-compute" values that are computed more than once per scope
//!!

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
            //negative terminal to ground
            if( nB(comps[i]) == 0 && nA(comps[i]) != 0)
            {
                //write the source voltage to this index in the rhs vector
                currents(nA(comps[i]) -1) = nA(comps[i]) - nB(comps[i]) - (component_currents(i))*interval/val;
            }

            //positive terminal to ground
            //not too sure about this one
            if(nA(comps[i]) == 0 && nB(comps[i]) != 0)
            {
                //also write -1* the source voltage to this index in the rhs vector
                currents(nB(comps[i]) -1) += (-1)*(nA(comps[i]) - nB(comps[i]) - (component_currents(i))*interval/val);
            }

            //floating voltage source
            //to optimise, make the row value saved after initial calculation
            if(nA(comps[i]) != 0 && nB(comps[i]) != 0)
            {
                row = c_vs_row[i];

                //in that row of the rhs vector, write the current source value
                currents(row) = (pastnodes(nA(comps[i])-1) - pastnodes(nB(comps[i])-1) - (component_currents(i))*interval/val) ;
                /*/
                cout << "VA= " << nA(comps[i]) << endl;
                cout << "VB= " << nB(comps[i]) << endl;
                cout << "comp current n-1 " << currents(row) << endl;
                cout << "cap voltage= " << currents(row) << endl;
                cout << "cap voltage= " << currents(row) << endl;
                /*/
            }
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
            //cerr << "oof" << endl;
            incorrect_assumptions[comps.size()] = 1;
            incorrect_assumptions[i] = 1;
        }
    }

    return incorrect_assumptions;
}

//return voltage and current vectors for operating point or first point of transient analysis
pair<VectorXd, VectorXd> no_prior_change (const vector<Component> & comps, const vector<Node> & nodes, const int & noden)
{
    VectorXd nodev;
    //perhaps avoid creating by separating recursive_currents used for op and trans
    VectorXd component_currents = VectorXd::Zero (comps.size());

    pair<MatrixXd, VectorXd> knowns = conductance_current (comps, noden);
    //test(noden, knowns.first, knowns.second);



    //cerr << "nodev" << endl;
    nodev = matrixSolve(knowns.first, knowns.second);
    VectorXd prevnodev = VectorXd::Zero(nodev.size());
    //cerr << "comp I" << endl;
    component_currents = recursive_currents (comps, nodes, nodev, prevnodev, 0, component_currents);

    //cerr << "return" << endl;
    return{nodev, component_currents};
}

//might be preferable to cout values directly (using Jason's function) after they're calculated to avoid cycling through the duration twice
//and have a void function if possible
vector<pair<VectorXd, VectorXd>> transient (vector<Component> & comps, const vector<Node> & nodes, const int & noden, const float & duration,
const float & interval, const VectorXd & pastnodes, const VectorXd & pastcurrents)
{
    VectorXd nodev = pastnodes;
    VectorXd component_currents = pastcurrents;
    vector<pair<VectorXd, VectorXd>> values;
    VectorXd rhs = VectorXd::Zero (comps.size());

    VectorXd prevnodev = VectorXd::Zero(nodev.size());
    vector<Component> newcomps = patchSupernodeInductor(comps);
    pair<MatrixXd, vector<int>> Mat = MatrixUpdate (newcomps, noden, interval);
    //cerr << Mat.first << endl;
    writeTranHeaders(nodes, comps);
    writeOPZero(pastnodes, pastcurrents);
    //begin one interval after 0
    //i is time in seconds
    for(float i = interval; i<duration; i += interval)
    {
        rhs = VectorUpdate (newcomps, noden, i, nodev, component_currents, interval, Mat.second);

        //cout << "rhs" << endl;
        //cout << rhs << endl;
        //cout << "done" << endl;

        /*/
        cout << "rhs:" << endl;
        for(int i = 0; i<rhs.size(); i++)
        {
            cout << rhs[i] << endl;
        }
        cout << "f" << endl;
        /*/

        prevnodev = nodev;
        nodev = matrixSolve(Mat.first, rhs);
        component_currents = recursive_currents (newcomps, nodes, nodev, prevnodev, interval, component_currents);

        /*/
        cout << "comp_i" << endl;
        cout << component_currents << "e" << endl;
        /*/

        values.push_back( {nodev, component_currents} );
        writeTran(nodev, component_currents, i);
    }

    return values;
}

//used in vs_current to determine other components connected to C at A
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
const float & interval, const VectorXd & past_currents)
{
    //register components take care of, to differetiate non calculated values from 0 currents
    vector<bool> computed (comps.size(), 0);

    VectorXd comp_currents = VectorXd::Zero (comps.size());

    for(int i = 0; i<comps.size(); i++)
    {
        //if( computed[i] == 0 )
        //{
            if(comps[i].type == 'L')
            {
                comp_currents( i ) = past_currents(i) + recursive_basecase (i, comps[i], comps, nlist, nodev, prevnodev, interval, computed, comp_currents,
                comps[i].A);
            } else {
                //cout << "b4 base" << endl;
                comp_currents( i ) = recursive_basecase (i, comps[i], comps, nlist, nodev, prevnodev, interval, computed, comp_currents, comps[i].A);
                //cout << comp_currents << endl;
                //cout << "after base" << endl;
            }

            //computed[component_index(comps, comps[i])] = 1;
        //}
    }
    return comp_currents;
}

//edits to computed below should be redundant
float recursive_basecase (const int & i, const Component & C, const vector<Component> & comps, const vector<Node> & nlist, const VectorXd & nodev,
const VectorXd & prevnodev, const float & interval, vector<bool> & computed, VectorXd & comp_currents, const Node & used_node)
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

    if(C.type == 'I')
    {
        total_current = C.value;
        computed[ component_index(comps, C) ] = 1;
    }

    if(C.type == 'R' || C.type == 'L')
    {
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

        //straightforward inclusion of known or calculated currents
        if(C.type == 'R')
        {
            //cout << "b4 R" << endl;
            //this must also be fixed
            total_current = ( VA - VB )/C.value;
            computed[ component_index(comps, C) ] = 1;
            /*/
            if(VA > 0)
            {
                total_current = ( VA - VB )/C.value;
                computed[i] = 1;
            } else {
                total_current = ( VB - VA )/C.value;
                computed[i] = 1;
            }
            /*/
            //cout << "after R" << endl;
        }

        if(C.type == 'L')
        {
            total_current = total_current = comp_currents( component_index( comps, C ) )+( prevnodev(nB(C) - 1) - prevnodev(nA(C) - 1) )*interval/(C.value);
            computed[ component_index(comps, C) ] = 1;
        }
    }

    //probably can replace some common_node calls with i
    if(C.type == 'V' || C.type == 'C')
    {
        //cout << "ree" << endl;
        total_current = 0;

        same_node = common_node(comps, C, used_node);

        for(int j = 0; j<same_node.size(); j++)
        {
            //cout << "for" << endl;
            //cout << C.name << " " << used_node.label << endl;
            if( computed[component_index(comps, same_node[j])] == 1 )
            {
                //cout << same_node[j].name << " computed" << endl;
                if(same_node[j].type == 'R')
                {
                    //cout << "rec R" << endl;
                    //why does this work? come back later
                    //also fix for some cases
                    //cout << used_node.label << " " << same_node[j].A.label << endl;
                    if(used_node.number == nA(C))
                    {
                        node_add = 0;
                    } else {
                        node_add = 1;
                    }

                    if(used_node.number == same_node[j].A.number && node_add == 0 || used_node.number == same_node[j].B.number && node_add == 1)
                    {
                        total_current -= comp_currents( component_index(comps, same_node[j]) );
                    } else {
                        total_current += comp_currents( component_index(comps, same_node[j]) );
                    }
                    //cout << "rec R end" << endl;
                }

                if(same_node[j].type == 'I' || same_node[j].type == 'L' )
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

                if(same_node[j].type == 'V' || same_node[j].type == 'C')
                {
                    //cout << "rec V"<< endl;
                    //still some redundant recalculations
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
                    //cout << "rec V end"<< endl;
                }
            } else {
                //cout << same_node[j].name << " compute rec" << endl;
                //cout << "rec else" << endl;

                //switch the used_node to the other node of same_node[j]
                if(used_node.label == same_node[j].A.label)
                {
                    rec_used = same_node[j].B;
                } else {
                    rec_used = same_node[j].A;
                }

                if(used_node.label == same_node[j].A.label)
                {
                    current = (-1)*recursive_basecase (i, same_node[j], comps, nlist, nodev, prevnodev, interval, computed, comp_currents, rec_used);
                } else {
                    current = recursive_basecase (i, same_node[j], comps, nlist, nodev, prevnodev, interval, computed, comp_currents, rec_used);
                }

                if(used_node.label == C.A)
                {
                    total_current += current;
                } else {
                    total_current -= current;
                }
                //cout << "rec else end" << endl;
            }
        }
        //cout << "eer" << endl;
       computed[ component_index(comps, C) ] = 1;
    }

    //cout << C.name << " " << total_current << endl;
    comp_currents(component_index(comps, C)) = total_current;
    return total_current;
}
