#include "conductance.hpp"

using namespace Eigen;
using namespace std;

int main()
{
    Component R1 {'R', "R1", 1, 2, 1};
    Component R2 {'V', "V1", 0, 1, 2};
    Component R3 {'R', "R3", 3, 2, 4};
    Component R4 {'R', "R4", 3, 0, 3};

    vector<Component> comps {R1, R2, R3, R4};

    vector<int> nodes {0, 1, 2, 3};

    int noden = compute_noden(nodes);

    pair<MatrixXd, vector<float>> knowns = conductance_current (comps, noden);

    test(noden, knowns.first, knowns.second);
}