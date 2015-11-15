#include <vector>
#include <boost/rational.hpp>
#include <iostream>

#include "math.hpp"

using namespace std;
using namespace boost;

void print_matrix(vector<vector<int>> matrix) {
    for (int i=0; i<matrix.size(); i++) {
        for (int j=0; j<matrix[i].size(); j++) {
            cout << matrix[i][j] << " ";
        }
        cout << endl;
    }
}

void print_matrix(vector<vector<rational<int>>> matrix) {
    for (int i=0; i<matrix.size(); i++) {
        for (int j=0; j<matrix[i].size(); j++) {
            cout << matrix[i][j] << " ";
        }
        cout << endl;
    }
}

void print_vector(vector<int> vect) {
    for (int i=0; i<vect.size(); i++) cout << vect[i] << " ";
    cout << endl;
}

void print_vector(vector<rational<int>> vect) {
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

vector<vector<int>> create_matrix(int size, int c) {
    vector<vector<int>> tmpMatrix;
    for (int i=0; i<size; i++) {
        vector<int> tmpVector;
        for (int j=0; j<size; j++) {
            tmpVector.push_back(c);
        }
        tmpMatrix.push_back(tmpVector);
    }

    return tmpMatrix;
}

vector<int> create_vector(int size) {
    vector<int> tmpVector;
    for(int i=0; i<size; i++) tmpVector.push_back(i);

    return tmpVector;
}

vector<vector<int>> create_indexes(int size) {
    vector<vector<int>> tmpIndexes;
    vector<int> tmpVector;
    for(int i=0; i<size; i++) tmpVector.push_back(i);
    tmpIndexes.push_back(tmpVector);

    return tmpIndexes;
}

vector<int> next(vector<int> v, int b) {
  int index=-1;
  bool t=false;
  for (int i=0; i<v.size(); i++) {
    if (v[i] < b-1) {
      t=true;
      index=i;
      break;
    }
  }
  if (index != -1) {
    fill(v.begin(), v.begin() + index, 0);
    v[index]++;
    return v;
  }
  else {
    vector<int> res;
    for (int i=0; i<v.size(); i++) res.push_back(-1);
    return res;
  }
}

vector<int> checkSysFix(vector<vector<rational<int>>> mat, vector<int> solution) {
    vector<int> res;
    bool test = true;
    for (int i=0; i<mat.size(); i++) {
        if (test == false) break;
        rational<int> sum(0);
        rational<int> vec = mat[i][mat[i].size() - 1];
        for (int j=0; j<mat[i].size()-1; j++) {
            sum += mat[i][j] * rational<int>(solution[j]);
        }
        if ((vec - sum != 0) && (vec - sum != 1)){
            test = false;
        }
        else {
            if (vec - sum == 0) res.push_back(0);
            if (vec - sum == 1) res.push_back(1);
        }
    }

    if (test == false) {
        res.clear();
        res.push_back(9);
    }
    return res;
}

vector<int> checkSysOrbit(vector<vector<rational<int>>> mat, vector<int> solution) {

    vector<int> res;
    bool test = true;
    for (int i=0; i<mat.size(); i++) {
        if (test == false) break;
        rational<int> sum(0);
        rational<int> vec = mat[i][mat[i].size() - 1];
        for (int j=0; j<mat[i].size()-1; j++) {
            sum += mat[i][j] * rational<int>(solution[j]);
        }
        if ((vec - sum != 0) && (vec - sum != 1) && (vec - sum != 2) && (vec - sum != 3)){
            test = false;
        }
        else {
            if (vec - sum == 0) res.push_back(0);
            if (vec - sum == 1) res.push_back(1);
            if (vec - sum == 2) res.push_back(2);
            if (vec - sum == 3) res.push_back(3);
        }
    }
    if (test == false) {
        res.clear();
        res.push_back(9);
    }

    return res;
}


vector<vector<rational<int>>> reducedRowEchelonForm(vector<vector<int>> matrix) {
    int lead = 0;
    int rowCount = matrix.size();
    int colCount = matrix[0].size();

    vector<vector<rational<int>>> newMatrix;
    for (int i=0; i<matrix.size(); i++) {
        vector<rational<int>> tmpVect;
        for (int j=0; j<matrix[i].size(); j++) {
            tmpVect.push_back(rational<int>(matrix[i][j]));
        }
        newMatrix.push_back(tmpVect);
    }  

    for (int r=0; r<rowCount; r++) {
        int i=r;

        while (lead < colCount && newMatrix[i][lead] == 0) {
            i++;
            if (i == rowCount) {
                i = r;
                lead++;
            }
        }

        if (lead < colCount) {
            // swap rows
            vector<rational<int>> tmp;
            tmp = newMatrix[i];
            newMatrix[i] = newMatrix[r];
            newMatrix[r] = tmp;

            rational<int> lv = newMatrix[r][lead];
            for (int index=0; index<newMatrix[r].size(); index++) {
                newMatrix[r][index] /= lv;
            }

            for (int i=0; i<rowCount; i++) {
                if (i != r) {
                    lv = newMatrix[i][lead];
                    for (int index=0; index<newMatrix[i].size(); index++) {
                        newMatrix[i][index] -= lv * newMatrix[r][index];
                    }
                }
            }
            lead++;
        }
    }

    return newMatrix;
}
