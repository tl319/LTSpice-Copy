#include "conductance.hpp"

using namespace Eigen;
using namespace std;

/*/
//move functions to header

//how are variable sources implemented
vector<pair<MatriXxd, vector<float>>> transient(vector<Components> comps, int noden, float length, float interval)
{
    vector<pair<MatriXxd, vector<float>>> v;
    for(float t = 0; t < length; t += interval)
    {
        v.push_back(conductance_current(comps, noden));
    }
} 

vector<float> solve_eqn(pair<MatrixXd, vector<float>> knowns)
{
    //invert matrix

    //multiply currents by inverse
}
/*/

int main()
{
    Node a {1, "a", 1};
    Node b {2, "b", 2};
    Node c {3, "c", 3};
    Node d {4, "d", 4};
    Node e {0, "e", 0};

    Component R1 {'R', "R1", a, b, 2};
    Component V1 {'V', "V1", a, e, 5};
    Component R2 {'R', "R2", c, b, 3};
    Component R4 {'R', "R4", e, d, 4};
    Component I1 {'I', "I1", c, d, 6};

    vector<Component> comps {R1, V1, R2, R4, I1};

    vector<Node> nodes {e, a, b, c, d};

    int noden = compute_noden(nodes);

    pair<MatrixXd, vector<float>> knowns = conductance_current (comps, noden);

    test(noden, knowns.first, knowns.second);
}