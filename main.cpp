#include <vector>
#include <algorithm>
#include <string>

#include <sys/stat.h>
#include <sys/types.h>

#include <iostream>
#include <fstream>
#include <boost/rational.hpp>
#include <tuple>

#include "libraryZ3.hpp"
#include "math.hpp"
#include "test.hpp"
#include "isomorphism.hpp"
#include "main.hpp"

using namespace std;
using namespace boost;

// strongly regular graph parameters
int n;
int k;
int l;
int m;
int p;
int f;
int o;

int isoConst = 5;

vector<vector<int>> linEqFixSolution;
vector<vector<int>> linEqOrbitSolution0;
vector<vector<int>> linEqOrbitSolution2;

vector<vector<vector<int>>> solutions;

// const
vector<vector<int>> iMatrix;
vector<vector<int>> jMatrix;
vector<vector<int>> constMatrix;
vector<int> nMatrix;
vector<int> lambdaSet;
vector<int> muSet;

// from Behbahani thesis
void solve_linear_equations() {
    {
        // Create and solve linear equations for fixed rows using z3 library
        // x0 + x1 = f
        // y0 + y1 = o
        // x1 + p*y1 = k
        vector<vector<int>> linearSystem {{1, 1, 0, 0}, {0, 0, 1, 1}, {0, 1, 0, p}};
        vector<int> vect {f, o, k};
        vector<vector<int>> tmpSolutions = solve_linear_equations(linearSystem, vect, 0, max(f, o));

        // filter: matrix must have 0 on diagonalise -> x0 >= 1
        for (int i=0; i<tmpSolutions.size(); i++) {
            if (tmpSolutions[i][0] != 0) linEqFixSolution.push_back(tmpSolutions[i]);
        }
    }

    {
        // Create and solve linear equations for fixed rows
        // x0 + x1 = f
        // y0 + y1 + y2 + y3 = o
        // x1 + y1 + 2 y2 + 3 y3 = k
        // p x1 + y1 + 4 y2 + 9 y3 = s_rr / p, where s_rr = (k - m) n_r[p] + m n_r^2[p] + (l-m)c_rr n_r[p] and c_rr is 0 1 2
        // there is lemma in theory c_rr could be only even, that 0 or 2 that implies 2 possible linear systems
        // s = (k-m)*p + m*p*p + (l-m)*c_rr*p
        vector<vector<int>> linearSystem = {{1, 1, 0, 0, 0, 0}, {0, 0, 1, 1, 1, 1}, {0, 1, 0, 1, 2, 3}, {0, 3, 0, 1, 4, 9}};
        int s = (k - m) * p + m * p * p + (l - m) * 0 * p;
        vector<int> vect {f, o, k, s / 3};
        vector<vector<int>> tmpSolutions = solve_linear_equations(linearSystem, vect, 0, max(f, o));

        // filter: if c_rr = 0 then y0 >0
        for (int i=0; i<tmpSolutions.size(); i++) {
            if (tmpSolutions[i][2] != 0) linEqOrbitSolution0.push_back(tmpSolutions[i]);
        }
    }

    {
        vector<vector<int>> linearSystem = {{1, 1, 0, 0, 0, 0}, {0, 0, 1, 1, 1, 1}, {0, 1, 0, 1, 2, 3}, {0, 3, 0, 1, 4, 9}};
        int s = (k - m) * p + m * p * p + (l - m) * 2 * p;
        vector<int> vect {f, o, k, s / 3};
        linEqOrbitSolution2 = solve_linear_equations(linearSystem, vect, 0, max(f, o));
    }
}

// create constant matrices
void const_matrices() {
    for (int i=0; i<f+o; i++) {
        vector<int> tmpVector;
        for(int j=0; j<f+o; j++) {
            if (i==j) tmpVector.push_back(k - m);
            else tmpVector.push_back(0);
        }
        iMatrix.push_back(tmpVector);
    }

    for (int i=0; i<f+o; i++) {
        vector<int> tmpVector;
        for(int j=0; j<f+o; j++) {
            if (j<f) tmpVector.push_back(1);
            else tmpVector.push_back(p);
        }
        jMatrix.push_back(tmpVector);
    }    

    constMatrix = addition_matrix(iMatrix, multiply_const_matrix(jMatrix, m));

    nMatrix = jMatrix[0];

    int i = l % p;
    while (i <= n) {
      lambdaSet.push_back(i);
      i += p;
    }

    int ii = m % p;
    while (ii <= n) {
      muSet.push_back(ii);
      ii += p;
    }
}

void generate_orbit_matrix(vector<vector<int>> matrix, int row, vector<vector<int>> indexes) {

    int numberFix0 = 0;
    int number0 = 0;
    int number1 = 0;
    int number2 = 0;
    int number3 = 0;

    vector<int> tmpSol;
    for (int i=0; i<f; i++) {
        if (matrix[row][i] == 0) numberFix0++;
        tmpSol.push_back(matrix[row][i]);
    }
    for (int i=f; i<row; i++) {
        if (matrix[row][i] == 0) number0++;
        if (matrix[row][i] == 1) number1++;
        if (matrix[row][i] == 2) number2++;
        if (matrix[row][i] == 3) number3++;
        tmpSol.push_back(matrix[row][i]);
    }

    // create lin. equations from matrix equation -> Behbahani
    // s_ir = mi n_i n_r + (l - m) c_ir n_r
    // s_ir = sum_k 1 to o {c_ik c_rk n_k}
    // number of lin. equations = row - 1 (s_12, s13)
    vector<vector<int>> linEq;
    for (int i=0; i<row; i++) {
        int value = (m * nMatrix[i] * nMatrix[row]) + ((l -m) * nMatrix[row] * matrix[i][row]);

        vector<int> equation;
        for (int j=0; j<(f+o)-row; j++) equation.push_back(0);

        for (int j=0; j<row; j++) {
            if (tmpSol[j] != 0) {
                value -= tmpSol[j] * matrix[i][j] * nMatrix[j];
            }
        }

        for (int j=row; j<f+o; j++) {
            equation[j - row] = matrix[i][j] * nMatrix[j];
        }

        equation.push_back(value);
        linEq.push_back(equation);
    }

    tuple<bool, vector<vector<rational<int>>>, vector<int>> rat = create_rational_system(linEq, tmpSol, row);
    bool test = get<0>(rat); 
    vector<vector<rational<int>>> mat2 = get<1>(rat);
    vector<int> sys = get<2>(rat);

    if (test == true) return;

    vector<int> v;
    vector<int> end;

    for (int i=0; i<mat2[0].size()-1; i++) {
        v.push_back(0);
        end.push_back(-1);
    }
    if (mat2[0].size() == 1) {
        vector<int> res = checkSysOrbit(mat2, v);
    
        vector<int> sol(v);
        for (int i=0; i<sys.size(); i++) sol.insert(sol.begin() + sys[i], res[i]);

        vector<int> possibleRow;
        for(int j=0; j<tmpSol.size(); j++) possibleRow.push_back(tmpSol[j]);
        for(int j=0; j<sol.size(); j++) possibleRow.push_back(sol[j]);

        int number00 = 0;
        int number11 = 0;
        int number22 = 0;
        int number33 = 0;

        bool testSeq = false;
        for (int j=0; j<indexes.size(); j++) {
            int num = -1;
            if (testSeq == true) break;
            for (int k=0; k<indexes[j].size(); k++) {
                if (indexes[j][k] + f > row) {
                    if (possibleRow[f + indexes[j][k]] >= num) {
                        num = possibleRow[f + indexes[j][k]];
                        if (possibleRow[f + indexes[j][k]] == 0) number00++;
                        if (possibleRow[f + indexes[j][k]] == 1) number11++;
                        if (possibleRow[f + indexes[j][k]] == 2) number22++;
                        if (possibleRow[f + indexes[j][k]] == 3) number33++;
                    }
                    else testSeq = true;
                }
            }
        }
        if (testSeq == false) {
            for (int jj=0; jj<linEqOrbitSolution0.size(); jj++) {
                if (possibleRow[row] == 0 &&
                    ((number2 + number22) + (3*(number3 + number33)) + (f-numberFix0)) == m &&
                    linEqOrbitSolution0[jj][0] == numberFix0 &&
                    linEqOrbitSolution0[jj][1] == f-numberFix0 &&
                    linEqOrbitSolution0[jj][2] == number0 + number00 + 1 &&
                    linEqOrbitSolution0[jj][3] == number1 + number11 &&
                    linEqOrbitSolution0[jj][4] == number2 + number22 &&
                    linEqOrbitSolution0[jj][5] == number3 + number33) {

                    vector<vector<int>> newMatrix(matrix);
                    newMatrix[row] = possibleRow;
                    for (int j=row+1; j<f+o; j++) newMatrix[j][row] = newMatrix[row][j];
                    
                    // test on isomorphism
                    bool testIso = false;
                    vector<int> cert = print_cert_matrix(newMatrix);
                    for (int j=0; j<indexes.size(); j++) {
                        if (testIso == true) break;
                        if (indexes[j].size() != 1) {
                            vector<int> tmp = indexes[j];
                            do {
                                vector<int> listPermutation = create_vector(f + o);
                                for (int k=0; k<indexes[j].size(); k++) listPermutation[f + indexes[j][k]] = f + tmp[k];
                                vector<int> cert2;
                                for (int k=0; k<newMatrix.size(); k++) {
                                    for (int l=k; l<newMatrix.size(); l++) cert2.push_back(newMatrix[listPermutation[k]][listPermutation[l]]);
                                }
                                if (operator<(cert2, cert)) {
                                    testIso = true;
                                    break;
                                }
                            } while(next_permutation(tmp.begin(),tmp.end()));
                        }
                    }
                    if (testIso == false) {
                        if (row + 1 == f + o) {
                            solutions.push_back(newMatrix);
                        }
                        else {
                            vector<int> orbitRow(possibleRow.begin() + f, possibleRow.end());
                            generate_orbit_matrix(newMatrix, row + 1, orbit_col_permutation(orbitRow, indexes));
                        }
                    }
                }
            }

            for (int jj=0; jj<linEqOrbitSolution2.size(); jj++) {
                if (possibleRow[row] == 2 &&
                    ((number2 + number22) + (3*(number3 + number33)) + (f-numberFix0) + 1) == l &&
                    linEqOrbitSolution2[jj][0] == numberFix0 &&
                    linEqOrbitSolution2[jj][1] == f-numberFix0 &&
                    linEqOrbitSolution2[jj][2] == number0 + number00 &&
                    linEqOrbitSolution2[jj][3] == number1 + number11 &&
                    linEqOrbitSolution2[jj][4] == number2 + number22 + 1 &&
                    linEqOrbitSolution2[jj][5] == number3 + number33) {
                    
                    vector<vector<int>> newMatrix(matrix);
                    newMatrix[row] = possibleRow;
                    for (int j=row+1; j<f+o; j++) newMatrix[j][row] = newMatrix[row][j];
                    
                    // test on isomorphism
                    bool testIso = false;
                    vector<int> cert = print_cert_matrix(newMatrix);
                    for (int j=0; j<indexes.size(); j++) {
                        if (testIso == true) break;
                        if (indexes[j].size() != 1) {
                            vector<int> tmp = indexes[j];
                            do {
                                vector<int> listPermutation = create_vector(f + o);
                                for (int k=0; k<indexes[j].size(); k++) listPermutation[f + indexes[j][k]] = f + tmp[k];
                                vector<int> cert2;
                                for (int k=0; k<newMatrix.size(); k++) {
                                    for (int l=k; l<newMatrix.size(); l++) cert2.push_back(newMatrix[listPermutation[k]][listPermutation[l]]);
                                }
                                if (operator<(cert2, cert)) {
                                    testIso = true;
                                    break;
                                }
                            } while(next_permutation(tmp.begin(),tmp.end()));
                        }
                    }

                    if (testIso == false) {
                        if (row + 1 == f + o) {
                            solutions.push_back(newMatrix);
                        }
                        else {
                            vector<int> orbitRow(possibleRow.begin() + f, possibleRow.end());
                            generate_orbit_matrix(newMatrix, row + 1, orbit_col_permutation(orbitRow, indexes));
                        }
                    }
                }
            }
        }
    }
    else {
        while (operator!=(v, end)) {

            vector<int> res = checkSysOrbit(mat2, v);
            if ((res.size() == 1) && (res[0] == 9)) {
                v = next(v, 4);
                continue;
            } 
            
            vector<int> sol(v);

            for (int i=0; i<sys.size(); i++) sol.insert(sol.begin() + sys[i], res[i]);
            
            vector<int> possibleRow;
            for(int j=0; j<tmpSol.size(); j++) possibleRow.push_back(tmpSol[j]);
            for(int j=0; j<sol.size(); j++) possibleRow.push_back(sol[j]);

            int number00 = 0;
            int number11 = 0;
            int number22 = 0;
            int number33 = 0;
            
            bool testSeq = false;
            for (int j=0; j<indexes.size(); j++) {
                int num = -1;
                if (testSeq == true) break;
                for (int k=0; k<indexes[j].size(); k++) {
                    if (indexes[j][k] + f > row) {
                        if (possibleRow[f + indexes[j][k]] >= num) {
                            num = possibleRow[f + indexes[j][k]];
                            if (possibleRow[f + indexes[j][k]] == 0) number00++;
                            if (possibleRow[f + indexes[j][k]] == 1) number11++;
                            if (possibleRow[f + indexes[j][k]] == 2) number22++;
                            if (possibleRow[f + indexes[j][k]] == 3) number33++;
                        }
                        else testSeq = true;
                    }
                }
            }
            
            if (testSeq == true) {
                v = next(v, 4);
                continue;
            }

            for (int jj=0; jj<linEqOrbitSolution0.size(); jj++) {
                if (possibleRow[row] == 0 &&
                    ((number2 + number22) + (3*(number3 + number33)) + (f-numberFix0)) == m &&
                    linEqOrbitSolution0[jj][0] == numberFix0 &&
                    linEqOrbitSolution0[jj][1] == f-numberFix0 &&
                    linEqOrbitSolution0[jj][2] == number0 + number00 + 1 &&
                    linEqOrbitSolution0[jj][3] == number1 + number11 &&
                    linEqOrbitSolution0[jj][4] == number2 + number22 &&
                    linEqOrbitSolution0[jj][5] == number3 + number33) {

                    vector<vector<int>> newMatrix(matrix);
                    newMatrix[row] = possibleRow;
                    for (int j=row+1; j<f+o; j++) newMatrix[j][row] = newMatrix[row][j];
                    
                    // test on isomorphism
                    bool testIso = false;
                    vector<int> cert = print_cert_matrix(newMatrix);
                    for (int j=0; j<indexes.size(); j++) {
                        if (testIso == true) break;
                        if (indexes[j].size() != 1) {
                            vector<int> tmp = indexes[j];
                            do {
                                vector<int> listPermutation = create_vector(f + o);
                                for (int k=0; k<indexes[j].size(); k++) listPermutation[f + indexes[j][k]] = f + tmp[k];
                                vector<int> cert2;
                                for (int k=0; k<newMatrix.size(); k++) {
                                    for (int l=k; l<newMatrix.size(); l++) cert2.push_back(newMatrix[listPermutation[k]][listPermutation[l]]);
                                }
                                if (operator<(cert2, cert)) {
                                    testIso = true;
                                    break;
                                }
                            } while(next_permutation(tmp.begin(),tmp.end()));
                        }
                    }
                    if (testIso == false) {
                        if (row + 1 == f + o) {
                            solutions.push_back(newMatrix);
                        }
                        else {
                            vector<int> orbitRow(possibleRow.begin() + f, possibleRow.end());
                            generate_orbit_matrix(newMatrix, row + 1, orbit_col_permutation(orbitRow, indexes));
                        }
                    }
                }
            }

            for (int jj=0; jj<linEqOrbitSolution2.size(); jj++) {
                if (possibleRow[row] == 2 &&
                    ((number2 + number22) + (3*(number3 + number33)) + (f-numberFix0) + 1) == l &&
                    linEqOrbitSolution2[jj][0] == numberFix0 &&
                    linEqOrbitSolution2[jj][1] == f-numberFix0 &&
                    linEqOrbitSolution2[jj][2] == number0 + number00 &&
                    linEqOrbitSolution2[jj][3] == number1 + number11 &&
                    linEqOrbitSolution2[jj][4] == number2 + number22 + 1 &&
                    linEqOrbitSolution2[jj][5] == number3 + number33) {
                    
                    vector<vector<int>> newMatrix(matrix);
                    newMatrix[row] = possibleRow;
                    for (int j=row+1; j<f+o; j++) newMatrix[j][row] = newMatrix[row][j];
                    
                    // test on isomorphism
                    bool testIso = false;
                    vector<int> cert = print_cert_matrix(newMatrix);
                    for (int j=0; j<indexes.size(); j++) {
                        if (testIso == true) break;
                        if (indexes[j].size() != 1) {
                            vector<int> tmp = indexes[j];
                            do {
                                vector<int> listPermutation = create_vector(f + o);
                                for (int k=0; k<indexes[j].size(); k++) listPermutation[f + indexes[j][k]] = f + tmp[k];
                                vector<int> cert2;
                                for (int k=0; k<newMatrix.size(); k++) {
                                    for (int l=k; l<newMatrix.size(); l++) cert2.push_back(newMatrix[listPermutation[k]][listPermutation[l]]);
                                }
                                if (operator<(cert2, cert)) {
                                    testIso = true;
                                    break;
                                }
                            } while(next_permutation(tmp.begin(),tmp.end()));
                        }
                    }

                    if (testIso == false) {
                        if (row + 1 == f + o) {
                            solutions.push_back(newMatrix);
                        }
                        else {
                            vector<int> orbitRow(possibleRow.begin() + f, possibleRow.end());
                            generate_orbit_matrix(newMatrix, row + 1, orbit_col_permutation(orbitRow, indexes));
                        }
                    }
                }
            }

            v = next(v, 4);
        }
    }
}

void generate_fix_matrices(vector<vector<int>> matrix, int row, vector<vector<int>> fixIndexes, vector<vector<int>> orbitIndexes) {
    if (row == 0) {
        for (int j=0; j<linEqFixSolution.size(); j++) {
            vector<int> newRow;
            for (int i=0; i<(f-linEqFixSolution[j][1]); i++) newRow.push_back(0);
            for (int i=0; i<linEqFixSolution[j][1]; i++) newRow.push_back(1);
            for (int i=0; i<(o-linEqFixSolution[j][3]); i++) newRow.push_back(0);
            for (int i=0; i<(linEqFixSolution[j][3]); i++) newRow.push_back(1);

            vector<vector<int>> newMatrix(matrix);
            newMatrix[row] = newRow;

            // transpose matrix
            for(int i=row+1; i<f; i++) newMatrix[i][row] = newMatrix[row][i];
            for(int i=f; i<f+o; i++) {
                if (newMatrix[row][i] == 1) newMatrix[i][row] = p;
                else newMatrix[i][row] = 0;
            }
            
            if (row + 1 == f) {
                vector<int> orbitRow(newRow.begin() + f, newRow.end());
                generate_orbit_matrix(newMatrix, f, fix_col_permutation(orbitRow, orbitIndexes));
            }
            else {
                vector<int> fixRow(newRow.begin(), newRow.begin() + f);
                vector<int> orbitRow(newRow.begin() + f, newRow.end());
                generate_fix_matrices(newMatrix, row + 1, fix_col_permutation(fixRow, fixIndexes), fix_col_permutation(orbitRow, orbitIndexes));
            }
        }
    }
    else {
        // create lin. equations from matrix equation -> Behbahani
        // s_ir = mi n_i n_r + (l - m) c_ir n_r
        // s_ir = sum_k 1 to o {c_ik c_rk n_k}
        // number of lin. equations = row - 1 (s_12, s13)

        vector<vector<int>> linEq;

        vector<int> tmpSol;
        int tmpSumFix = 0;

        for (int i=0; i<row; i++) {
            tmpSol.push_back(matrix[row][i]);
            tmpSumFix += matrix[row][i];
        }

        for (int i=0; i<row; i++) {
            int value = (m * nMatrix[i] * nMatrix[row]) + ((l -m) * nMatrix[row] * matrix[i][row]);

            vector<int> equation;
            for (int j=0; j<(f+o)-row; j++) equation.push_back(0);

            for (int j=0; j<row; j++) {
                if (tmpSol[j] == 1) value -= matrix[i][j] * nMatrix[j];
            }

            for (int j=row; j<f+o; j++) {
                equation[j - row] = matrix[i][j] * nMatrix[j];
            }

            equation.push_back(value);
            linEq.push_back(equation);
        }

        tuple<bool, vector<vector<rational<int>>>, vector<int>> rat = create_rational_system(linEq, tmpSol, row);

        bool test = get<0>(rat); 
        vector<vector<rational<int>>> mat2 = get<1>(rat);
        vector<int> sys = get<2>(rat);

        if (test == true) return;

        vector<int> v(mat2[0].size() - 1);
        vector<int> end; 
        for (int i=0; i<mat2[0].size()-1; i++) end.push_back(-1);

        while (operator!=(v, end)) {

            vector<int> res = checkSysFix(mat2, v);

            if ((res.size() == 1) && (res[0] == 9)) {
                v = next(v, 2);
                continue;
            }

            vector<int> sol(v);
            for (int i=0; i<sys.size(); i++) sol.insert(sol.begin() + sys[i], res[i]);

            vector<int> possibleRow;
            for(int j=0; j<tmpSol.size(); j++) possibleRow.push_back(tmpSol[j]);
            for(int j=0; j<sol.size(); j++) possibleRow.push_back(sol[j]);

            int sumFix = tmpSumFix;
            int sumOrbit = 0;

            bool testSeq = false;

            // test on decreasing sequence
            for (int j=0; j<fixIndexes.size(); j++) {
                int num = -1;
                if (testSeq == true) continue;
                for (int k=0; k<fixIndexes[j].size(); k++) {
                    if (fixIndexes[j][k] > row) {
                        if (possibleRow[fixIndexes[j][k]] >= num) {
                            num = possibleRow[fixIndexes[j][k]];
                            sumFix += possibleRow[fixIndexes[j][k]];
                        }
                        else testSeq = true;
                    }
                }
            }
            if (testSeq == true) {
                v = next(v, 2);
                continue;
            }

            for (int j=0; j<orbitIndexes.size(); j++) {
                int num = -1;
                if (testSeq == true) continue;
                for (int k=0; k<orbitIndexes[j].size(); k++) {
                    if (f + orbitIndexes[j][k] > row) {
                        if (possibleRow[f + orbitIndexes[j][k]] >= num) {
                            num = possibleRow[f + orbitIndexes[j][k]];
                            sumOrbit += possibleRow[f + orbitIndexes[j][k]];
                        }
                        else testSeq = true;
                    }
                }
            }
            if (testSeq == true) {
                v = next(v, 2);
                continue;
            }

            bool testDegree = false;
            for (int j=0; j<linEqFixSolution.size(); j++) {
                if ((possibleRow[row] == 0) && (sumFix == linEqFixSolution[j][1]) && (sumOrbit == linEqFixSolution[j][3])) testDegree = true;
            }

            if (testDegree == false) {
                v = next(v, 2);
                continue;
            }

            // transpose matrix
            matrix[row] = possibleRow;
            for(int j=row+1; j<f; j++) matrix[j][row] = matrix[row][j];
            for(int j=f; j<f+o; j++) {
                if (matrix[row][j] == 1) matrix[j][row] = p;
                else matrix[j][row] = 0;
            }

            // lambda/ mu test
            bool testLM = false;
            for (int j=0; j<row; j++) {
                int count = 0;

                for (int k=0; k<f; k++) {
                    if ((matrix[row][k] == 1) && (matrix[j][k] == 1)) count++;
                }

                if (matrix[j][row] == 1) {
                    // lambdaSet
                    if(std::find(lambdaSet.begin(), lambdaSet.end(), count) == lambdaSet.end()) {
                        testLM = true;
                    }
                }
                else {
                    // muSet
                    if(std::find(muSet.begin(), muSet.end(), count) == muSet.end()) {
                        testLM = true;
                    }
                }
            }
            if (testLM == true) {
                v = next(v, 2);
                continue;
            }

            // test on isomorphism
            bool testIso = false;
            if (row < isoConst) {
                vector<int> cert = print_cert_matrix(matrix);
                
                for (int j=0; j<fixIndexes.size(); j++) {
                    if (testIso == true) break;
                    if (fixIndexes[j].size() != 1) {
                        vector<int> tmp = fixIndexes[j];
                        do {
                            vector<int> listPermutation = create_vector(f + o);
                            for (int k=0; k<fixIndexes[j].size(); k++) listPermutation[fixIndexes[j][k]] = tmp[k];
                            vector<int> cert2;
                            for (int k=0; k<matrix.size(); k++) {
                                for (int l=k; l<matrix.size(); l++) cert2.push_back(matrix[listPermutation[k]][listPermutation[l]]);
                            }
                            if (operator<(cert2, cert)) {
                                testIso = true;
                                break;
                            }
                        } while(next_permutation(tmp.begin(),tmp.end()));
                    }
                }
                if (testIso == true) {
                    v = next(v, 2);
                    continue;
                }

                for (int j=0; j<orbitIndexes.size(); j++) {
                    if (testIso == true) break;
                    if (orbitIndexes[j].size() != 1) {
                        vector<int> tmp = orbitIndexes[j];
                        do {
                            vector<int> listPermutation = create_vector(f + o);
                            for (int k=0; k<orbitIndexes[j].size(); k++) listPermutation[f + orbitIndexes[j][k]] = f + tmp[k];
                            vector<int> cert2;
                            for (int k=0; k<matrix.size(); k++) {
                                for (int l=k; l<matrix.size(); l++) cert2.push_back(matrix[listPermutation[k]][listPermutation[l]]);
                            }
                            if (operator<(cert2, cert)) {
                                testIso = true;
                                break;
                            }
                        } while(next_permutation(tmp.begin(),tmp.end()));
                    }
                }
                if (testIso == true) {
                    v = next(v, 2);
                    continue;
                }
            }

            if (row + 1 == f) {
                vector<int> orbitRow(possibleRow.begin() + f, possibleRow.end());
                generate_orbit_matrix(matrix, f, fix_col_permutation(orbitRow, orbitIndexes));
            }
            else {
                vector<int> fixRow(possibleRow.begin(), possibleRow.begin() + f);
                vector<int> orbitRow(possibleRow.begin() + f, possibleRow.end());
                generate_fix_matrices(matrix, row + 1, fix_col_permutation(fixRow, fixIndexes), fix_col_permutation(orbitRow, orbitIndexes));
            }           
            
            v = next(v, 2);
        }
    }
}

void save_matrices(bool debug) {

    string name =  to_string(n) + "," + to_string(k) + "," + to_string(l) + "," + to_string(m) + "," + to_string(f) + "," + to_string(p);
    string fold = "matrices/srg(" + name + ")";

    mkdir("matrices", 0777);
    mkdir(fold.c_str(), 0777);

    string fil1 = fold + "/RAW-solutions(" + name + ").txt";
    ofstream file;
    file.open(fil1);

    string fil2 = fold + "/RAW-solutions(GAP)(" + name + ").txt";
    ofstream file2;
    file2.open(fil2);
    file2 << "sol := [";

    vector<vector<vector<int>>> newSolutions;
    for (int i=0; i<solutions.size(); i++) {
        vector<vector<int>> tmpMatrix = solutions[i];
        for (int j=0; j<f; j++) {
            for (int k=f; k<f+o; k++) {
                if (tmpMatrix[j][k] == 1) tmpMatrix[j][k] = 3;
            }
        }
        for (int j=f; j<f+o; j++) {
            for (int k=0; k<f; k++) {
                if (tmpMatrix[j][k] == 3) tmpMatrix[j][k] = 1;
            }
        }

        for (int j=0; j<f+o; j++) {
            for (int k=0; k<f+o; k++) {
                file << tmpMatrix[j][k];
                if (k != f + o - 1) file << " ";
            }
            if (j != f + o - 1) file << "; ";
        }
        file << "\n";
        
        file2 << "[";
        for (int j=0; j<f+o; j++) {
            file2 << "[";
            for (int k=0; k<f+o; k++) {
                file2 << tmpMatrix[j][k];
                if (k != f + o - 1) file2 << ", ";
            }
            file2 << "]";
            if (j != f + o - 1) file2 << ", ";
        }
        file2 << "]";
        if (i != solutions.size() - 1) file2 << ", ";
    }
    file2 << "];";

    file.close();
    file2.close();

    string output = fold + "/solutions(" + name + ").txt"; 
    string output2 = fold + "/number-solutions(" + name + ").txt"; 

    string sage = "sage Test.sage '" + fil2 + "' '" + output + "' '" + output2 + "'";

    if (debug == false) sage += "> tmp";
    system(sage.c_str());
}

int load_number_solutions() {
    string name =  to_string(n) + "," + to_string(k) + "," + to_string(l) + "," + to_string(m) + "," + to_string(f) + "," + to_string(p);
    string fold = "matrices/srg(" + name + ")/number-solutions(" + name + ").txt";
    string line;
    ifstream myfile (fold);

    int number;

    if (myfile.is_open()) {
        getline (myfile,line);
        number = atoi(line.c_str());
        
        myfile.close();
    }

    return number;
}

int generate(int nn, int kk, int ll, int mm, int pp, int ff, bool debug) {
    n = nn;
    k = kk;
    l = ll;
    m = mm;
    p = pp;
    f = ff;
    o = (n - f) / p;

    linEqFixSolution.clear();
    linEqOrbitSolution0.clear();
    linEqOrbitSolution2.clear();

    solutions.clear();

    // const
    iMatrix.clear();
    jMatrix.clear();
    constMatrix.clear();
    nMatrix.clear();
    lambdaSet.clear();
    muSet.clear();

    if (debug == true) {
        cout << "Program start " << n << " " << k << " " << l << " " << m << " " << p << " " << f << endl;
        cout << "Number of fixed points and orbits: " << f << " " << o << endl;
        cout << "Size of orbit matrix: " << (f + o) << "x" << (f + o) << endl;
    }

    solve_linear_equations();

    if (debug == true) {
        cout << "Number of solutions on fixed rows: " << linEqFixSolution.size() << endl;
        cout << "Number of solutions on orbit rows: " << linEqOrbitSolution0.size() << " " << linEqOrbitSolution2.size() << endl;
    }

    const_matrices();

    generate_fix_matrices(create_matrix(f + o, 5), 0, create_indexes(f), create_indexes(o));

    if (debug == true) {
        cout << "Number of solutions: " << solutions.size() << endl;
    }

    save_matrices(debug);

    return load_number_solutions();
}

int main() {
    // generate(15, 6, 1, 3, 3, 3, true);
    // generate(28, 12, 6, 4, 3, 10, true); 
    // generate(35, 16, 6, 8, 3, 2, true);
    // generate(37, 18, 8, 9, 3, 1, true);
    // generate(41, 20, 9, 10, 3, 5, true);
    small_tests();
    return 0;
}
