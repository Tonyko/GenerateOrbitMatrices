#include<vector>
using namespace std;

vector<vector<int>> solve_linear_equations(vector<vector<int>> linEquations, vector<int> vect, int from, int to);
vector<vector<int>> solve_linear_equations_with_constraints_fix(vector<vector<int>> linEquations, vector<int> vect, int from, int to, int row, int f, int fixDegree, int orbitDegree);
vector<vector<int>> solve_linear_equations_with_constraints_orbit(vector<vector<int>> linEquations, vector<int> vect, int from, int to, int value);
