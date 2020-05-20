//============================================================================
// Name        : HelloWorld.cpp
// Author      :
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>
#include <string>
#include <vector>
#include <sstream>
#include <Eigen>

using namespace std;
using namespace Eigen;
struct Component
{
	char type;
	string name;
	string nA;
	string nB;
	float value;
};

bool isComponent(string x)
{
	if(isalpha(x[0])){
		return true;}
	return false;
}

VectorXd matrixSolve(MatrixXd m,VectorXd v)
{
	MatrixXd inverse = m.inverse();
	return v * inverse;


}


ostream& operator<<(ostream& os, const Component& c)
{
    os << c.type << ':' << c.name << ':' << c.nA << ':' << c.nB << ':' << c.value << endl;
    return os;
}
vector<Component> readInput()
{
	string x;
	vector<string> strings;
	vector<Component> components;
	while(getline(cin, x))
	{
		strings.push_back(x);
	}

	for(auto line : strings){
		if(isComponent(line)){
			string properties[4];
			stringstream ss(line);
			int count=0;
			 while (ss >> properties[count]){
					count++; }
			string name;
			if(isalnum((properties[0])[1])){name = properties[0];}
			else{name =(properties[0]).substr (3,(properties[0].length())-1);}
			Component c1 = {(properties[0])[0],name,properties[1],properties[2],stof(properties[3])};
			components.push_back(c1);

		}
		}


	return components;

}



