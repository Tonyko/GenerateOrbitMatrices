#include<vector>
#include<iostream>

#include "math.hpp"

using namespace std;

void print_matrix(vector<vector<int>> matrix) {
    for (int i=0; i<matrix.size(); i++) {
        for (int j=0; j<matrix.size(); j++) {
            cout << matrix[i][j] << " ";
        }
        cout << endl;
    }
}

void print_vector(vector<int> vect) {
    for (int i=0; i<vect.size(); i++) cout << vect[i] << " ";
    cout << endl;
}
