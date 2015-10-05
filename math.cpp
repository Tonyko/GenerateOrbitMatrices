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

vector<vector<int>> addition_matrix(vector<vector<int>> matrixA, vector<vector<int>> matrixB) {
    vector<vector<int>> tmpMatrix = matrixA;
    for (int i=0; i<matrixA.size(); i++) {
        for (int j=0; j<matrixA[i].size(); j++) {
            tmpMatrix[i][j] = matrixA[i][j] + matrixB[i][j];
        }
    }

    return tmpMatrix;
}

vector<vector<int>> multiply_const_matrix(vector<vector<int>> matrix, int c) {
    vector<vector<int>> tmpMatrix = matrix;
    for (int i=0; i<matrix.size(); i++) {
        for (int j=0; j<matrix[i].size(); j++) {
            tmpMatrix[i][j] = matrix[i][j] * c;
        }
    }

    return tmpMatrix;
}