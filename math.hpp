#include <vector>
#include <boost/rational.hpp>

using namespace std;
using namespace boost;

void print_matrix(vector<vector<int>> matrix);
void print_matrix(vector<vector<rational<int>>> matrix);
void print_vector(vector<int> vect);
void print_vector(vector<rational<int>> vect);

vector<vector<int>> addition_matrix(vector<vector<int>> matrixA, vector<vector<int>> matrixB);
vector<vector<int>> multiply_const_matrix(vector<vector<int>> matrix, int c);
vector<vector<int>> create_matrix(int size, int c);
vector<vector<int>> create_indexes(int size);
vector<int> create_vector(int size);
vector<vector<rational<int>>> reducedRowEchelonForm(vector<vector<int>> matrix);
vector<int> checkSysOrbit(vector<vector<rational<int>>> mat, vector<int> solution);
vector<int> checkSysFix(vector<vector<rational<int>>> mat, vector<int> solution);
vector<int> next(vector<int> v, int b);