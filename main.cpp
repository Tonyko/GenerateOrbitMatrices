#include<iostream>
#include<vector>

#include "libraryZ3.hpp"

using namespace std;

int n = 15;
int k = 6;
int l = 1;
int m = 3;
int p = 3;
int f = 3;
int o = (n - f) / p;

int isoConst = 5;

vector<vector<int>> linEqFixSolution;
vector<vector<int>> linEqOrbitSolution0;
vector<vector<int>> linEqOrbitSolution2;

int main() {
    cout << "Program start " << n << " " << k << " " << l << " " << m << " " << p << " " << f << endl;
    cout << "Number of fixed points and orbits: " << f << " " << o << endl;
    cout << "Size of orbit matrix: " << (f + o) << "x" << (f + o) << endl;

    // Create and solve linear equations for fixed rows
    // x0 + x1 = f
    // y0 + y1 = o
    // x1 + p*y1 = k

    vector<vector<int>> linearSystem {{1, 1, 0, 0}, {0, 0, 1, 1}, {0, 1, 0, p}};
    vector<int> vector {f, o, k};

    // filter: matrix B must have 0 -> x0 >= 1, is graph with zero's on diagonalise

    linEqFixSolution = solveFixedEq(linearSystem, vector, 4).filter(x => x(0) != 0)


    vector<vector<int>> linearSystem {{1,1,0,0,0,0}, {0,0,1,1,1,1}, {0,1,0,1,2,3}, {0,3,0,1,4,9}};
    vector<int> vect{3,4,6,12}; 
    vector<vector<int>> solution = solveLinearEquations(linearSystem, vect, 0, 3);

    return 0;
}
