#include "conductance.hpp"

using namespace Eigen;
using namespace std;

int main()
{
    Node a {1, "a", 1};
    Node b {2, "b", 3};
    Node c {3, "c", 3};
    Node d {0, "d", 0};

    Component R1 {'R', "R1", a, b, 1};
    Component R2 {'V', "V1", a, d, 2};
    Component R3 {'V', "V2", c, b, 4};
    Component R4 {'R', "R4", c, d, 3};

    vector<Component> comps {R1, R2, R3, R4};

    vector<Node> nodes {d, a, b, c};

    int noden = compute_noden(nodes);

    pair<MatrixXd, vector<float>> knowns = conductance_current (comps, noden);

    test(noden, knowns.first, knowns.second);
}