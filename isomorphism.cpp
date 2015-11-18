#include <vector>
#include <iostream>

#include "isomorphism.hpp"

using namespace std;

vector<vector<int>> fix_col_permutation(vector<int> row, vector<vector<int>> indexes) {
    vector<vector<int>> tmpMatrix;
    for (int i=0; i<indexes.size(); i++) {
        vector<int> list1;
        vector<int> list2;
        for (int j=0; j<indexes[i].size(); j++) {
            if (row[indexes[i][j]] == 0) list1.push_back(indexes[i][j]);
            if (row[indexes[i][j]] == 1) list2.push_back(indexes[i][j]);
        }
        if (!list1.empty()) tmpMatrix.push_back(list1);
        if (!list2.empty()) tmpMatrix.push_back(list2);
    }

    return tmpMatrix;
}

vector<vector<int>> orbit_col_permutation(vector<int> row, vector<vector<int>> indexes) {
    vector<vector<int>> tmpMatrix;
    for (int i=0; i<indexes.size(); i++) {
        vector<int> list1;
        vector<int> list2;
        vector<int> list3;
        vector<int> list4;
        for (int j=0; j<indexes[i].size(); j++) {
            if (row[indexes[i][j]] == 0) list1.push_back(indexes[i][j]);
            if (row[indexes[i][j]] == 1) list2.push_back(indexes[i][j]);
            if (row[indexes[i][j]] == 2) list3.push_back(indexes[i][j]);
            if (row[indexes[i][j]] == 3) list4.push_back(indexes[i][j]);
        }
        if (!list1.empty()) tmpMatrix.push_back(list1);
        if (!list2.empty()) tmpMatrix.push_back(list2);
        if (!list3.empty()) tmpMatrix.push_back(list3);
        if (!list4.empty()) tmpMatrix.push_back(list4);
    }

    return tmpMatrix;
}

vector<int> print_cert_matrix(vector<vector<int>> matrix) {
    vector<int> tmpVector;
    for (int i=0; i<matrix.size(); i++) {
        for(int j=i; j<matrix[i].size(); j++) {
            tmpVector.push_back(matrix[i][j]);
        }
    }

    return tmpVector;
}
