#include <iostream>
#include <vector>
#include <algorithm>

#include "libraryZ3.hpp"
#include "math.hpp"

using namespace std;

// strongly regular graph parameters
int n = 25;
int k = 12;
int l = 5;
int m = 6;
int p = 3;
int f = 1;
int o = (n - f) / p;

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

    vector<vector<int>> linEquations;
    vector<int> vect;

    for (int i=0; i<row; i++) {
        int value = (m * nMatrix[i] * nMatrix[row]) + ((l -m) * nMatrix[row] * matrix[i][row]);

        vector<int> equation;
        for (int j=0; j<(f+o)-row; j++) equation.push_back(0);

        for (int j=0; j<row; j++) {
            if (tmpSol[j] != 0) value -= tmpSol[j] * matrix[i][j] * nMatrix[j];
        }

        for (int j=row; j<f+o; j++) {
            equation[j - row] = matrix[i][j] * nMatrix[j];
        }

        vect.push_back(value);
        linEquations.push_back(equation);
    }

    vector<vector<int>> sol0 = solve_linear_equations_with_constraints_orbit(linEquations, vect, 0, 3, 0);
    for (int i=0; i<sol0.size(); i++) {
        vector<int> possibleRow;
        for (int j=0; j<tmpSol.size(); j++) possibleRow.push_back(tmpSol[j]);
        for (int j=0; j<sol0[i].size(); j++) possibleRow.push_back(sol0[i][j]);

        int number00 = 0;
        int number11 = 0;
        int number22 = 0;
        int number33 = 0;

        bool testSeq = false;
        for (int j=0; j<indexes.size(); j++) {
            int num = -1;
            if (testSeq == true) continue;
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
        if (testSeq == true) continue;

        for (int jj=0; jj<linEqOrbitSolution0.size(); jj++) {
            if (linEqOrbitSolution0[jj][0] == numberFix0 &&
                linEqOrbitSolution0[jj][1] == f-numberFix0 &&
                linEqOrbitSolution0[jj][2] == number0 + number00 + 1 &&
                linEqOrbitSolution0[jj][3] == number1 + number11 &&
                linEqOrbitSolution0[jj][4] == number2 + number22 &&
                linEqOrbitSolution0[jj][5] == number3 + number33) {
                
                matrix[row] = possibleRow;
                for (int j=row+1; j<f+o; j++) matrix[j][row] = matrix[row][j];
                
                // test on isomorphism
                bool testIso = false;
                vector<int> cert = print_cert_matrix(matrix);
                for (int j=0; j<indexes.size(); j++) {
                    if (testIso == true) break;
                    if (indexes[j].size() != 1) {
                        vector<int> tmp = indexes[j];
                        do {
                            vector<int> listPermutation = create_vector(f + o);
                            for (int k=0; k<indexes[j].size(); k++) listPermutation[f + indexes[j][k]] = f + tmp[k];
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
                if (testIso == true) continue;

                if (row + 1 == f + o) {
                    solutions.push_back(matrix);
                }
                else {
                    vector<int> orbitRow(possibleRow.begin() + f, possibleRow.end());
                    generate_orbit_matrix(matrix, row + 1, orbit_col_permutation(orbitRow, indexes));
                }
            }
        }
    }

    vector<vector<int>> sol2 = solve_linear_equations_with_constraints_orbit(linEquations, vect, 0, 3, 2);
    for (int i=0; i<sol2.size(); i++) {
        vector<int> possibleRow;
        for (int j=0; j<tmpSol.size(); j++) possibleRow.push_back(tmpSol[j]);
        for (int j=0; j<sol2[i].size(); j++) possibleRow.push_back(sol2[i][j]);

        int number00 = 0;
        int number11 = 0;
        int number22 = 0;
        int number33 = 0;

        bool testSeq = false;
        for (int j=0; j<indexes.size(); j++) {
            int num = -1;
            if (testSeq == true) continue;
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
        if (testSeq == true) continue;

        for (int jj=0; jj<linEqOrbitSolution2.size(); jj++) {
            if (linEqOrbitSolution2[jj][0] == numberFix0 &&
                linEqOrbitSolution2[jj][1] == f-numberFix0 &&
                linEqOrbitSolution2[jj][2] == number0 + number00 &&
                linEqOrbitSolution2[jj][3] == number1 + number11 &&
                linEqOrbitSolution2[jj][4] == number2 + number22 + 1 &&
                linEqOrbitSolution2[jj][5] == number3 + number33) {
                
                matrix[row] = possibleRow;
                for (int j=row+1; j<f+o; j++) matrix[j][row] = matrix[row][j];
                
                // test on isomorphism
                bool testIso = false;
                vector<int> cert = print_cert_matrix(matrix);
                for (int j=0; j<indexes.size(); j++) {
                    if (testIso == true) break;
                    if (indexes[j].size() != 1) {
                        vector<int> tmp = indexes[j];
                        do {
                            vector<int> listPermutation = create_vector(f + o);
                            for (int k=0; k<indexes[j].size(); k++) listPermutation[f + indexes[j][k]] = f + tmp[k];
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
                if (testIso == true) continue;

                if (row + 1 == f + o) {
                    solutions.push_back(matrix);
                }
                else {
                    vector<int> orbitRow(possibleRow.begin() + f, possibleRow.end());
                    generate_orbit_matrix(matrix, row + 1, orbit_col_permutation(orbitRow, indexes));
                }
            }
        }
    }
}

void generate_fix_matrices(vector<vector<int>> matrix, int row, int fixDegree, int orbitDegree, vector<vector<int>> fixIndexes, vector<vector<int>> orbitIndexes) {
    if (row == 0) {
        vector<int> newRow;
        for (int i=0; i<(f-fixDegree); i++) newRow.push_back(0);
        for (int i=0; i<fixDegree; i++) newRow.push_back(1);
        for (int i=0; i<(o-orbitDegree); i++) newRow.push_back(0);
        for (int i=0; i<(orbitDegree); i++) newRow.push_back(1);

        matrix[row] = newRow;

        // transpose matrix
        for(int i=row+1; i<f; i++) matrix[i][row] = matrix[row][i];
        for(int i=f; i<f+o; i++) {
            if (matrix[row][i] == 1) matrix[i][row] = p;
            else matrix[i][row] = 0;
        }
        
        if (row + 1 == f) {
            vector<int> orbitRow(newRow.begin() + f, newRow.end());
            generate_orbit_matrix(matrix, f, fix_col_permutation(orbitRow, orbitIndexes));
        }
        else {
            vector<int> fixRow(newRow.begin(), newRow.begin() + f);
            vector<int> orbitRow(newRow.begin() + f, newRow.end());
            for (int i=0; i<linEqFixSolution.size(); i++) {
                generate_fix_matrices(matrix, row + 1, linEqFixSolution[i][1], linEqFixSolution[i][3], fix_col_permutation(fixRow, fixIndexes),
                 fix_col_permutation(orbitRow, orbitIndexes));
            }            
        }
    }
    else {
        // create lin. equations from matrix equation -> Behbahani
        // s_ir = mi n_i n_r + (l - m) c_ir n_r
        // s_ir = sum_k 1 to o {c_ik c_rk n_k}
        // number of lin. equations = row - 1 (s_12, s13)

        vector<vector<int>> linEquations;
        vector<int> vect;

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

            vect.push_back(value);
            linEquations.push_back(equation);
        }

        vector<vector<int>> sol = solve_linear_equations_with_constraints_fix(linEquations, vect, 0, 1, row, f, fixDegree-tmpSumFix, orbitDegree);

        for (int i=0; i<sol.size(); i++) {

            vector<int> possibleRow;
            for(int j=0; j<tmpSol.size(); j++) possibleRow.push_back(tmpSol[j]);
            for(int j=0; j<sol[i].size(); j++) possibleRow.push_back(sol[i][j]);

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
            if (testSeq == true) continue;

            for (int j=0; j<orbitIndexes.size(); j++) {
                int num = -1;
                if (testSeq == true) continue;
                for (int k=0; k<orbitIndexes[j].size(); k++) {
                    if (f + orbitIndexes[j][k] > row) {
                        if (possibleRow[f + orbitIndexes[j][k]] >= num) {
                            num = possibleRow[f + orbitIndexes[j][k]];
                            sumFix += possibleRow[f + orbitIndexes[j][k]];
                        }
                        else testSeq = true;
                    }
                }
            }
            if (testSeq == true) continue;

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
            if (testLM == true) continue;

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
                if (testIso == true) continue;

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
                if (testIso == true) continue;

                if (row + 1 == f) {
                    vector<int> orbitRow(possibleRow.begin() + f, possibleRow.end());
                    generate_orbit_matrix(matrix, f, fix_col_permutation(orbitRow, orbitIndexes));
                }
                else {
                    vector<int> fixRow(possibleRow.begin(), possibleRow.begin() + f);
                    vector<int> orbitRow(possibleRow.begin() + f, possibleRow.end());
                    for (int i=0; i<linEqFixSolution.size(); i++) {
                        generate_fix_matrices(matrix, row + 1, linEqFixSolution[i][1], linEqFixSolution[i][3], fix_col_permutation(fixRow, fixIndexes),
                         fix_col_permutation(orbitRow, orbitIndexes));
                    }            
                }           
            }
        }
    }
}

int main() {
    cout << "Program start " << n << " " << k << " " << l << " " << m << " " << p << " " << f << endl;
    cout << "Number of fixed points and orbits: " << f << " " << o << endl;
    cout << "Size of orbit matrix: " << (f + o) << "x" << (f + o) << endl;

    solve_linear_equations();

    cout << "Number of solutions on fixed rows: " << linEqFixSolution.size() << endl;
    cout << "Number of solutions on orbit rows: " << linEqOrbitSolution0.size() << " " << linEqOrbitSolution2.size() << endl;

    const_matrices();

    for (int i=0; i<linEqFixSolution.size(); i++) {
        generate_fix_matrices(create_matrix(f + o, 5), 0, linEqFixSolution[i][1], linEqFixSolution[i][3], create_indexes(f), create_indexes(o));
    }

    cout << "Number of solutions: " << solutions.size() << endl;

    return 0;
}
