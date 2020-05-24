#include <iostream>
#include <string>
#include <vector>
#include <numeric>

using namespace std;

int main()
{
    int noden = 4;

    vector<string> nodes(noden);
    iota(nodes.begin(), nodes.end(), 0);

    for(int i = 0; i<nodes.size(); i++)
    {
        cout << nodes[i] << " ";
    }
}