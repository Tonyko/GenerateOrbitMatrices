#include<z3++.h>
#include<vector>

#include "libraryZ3.hpp"

using namespace z3;

expr mk_and(expr_vector args) {
    vector<Z3_ast> array;

    for (int i = 0; i < args.size(); i++)
      array.push_back(args[i]);

    return to_expr(args.ctx(), Z3_mk_and(args.ctx(), array.size(), &(array[0])));
}

expr mk_or(expr_vector args) {
    vector<Z3_ast> array;

    for (int i = 0; i < args.size(); i++)
      array.push_back(args[i]);

    return to_expr(args.ctx(), Z3_mk_or(args.ctx(), array.size(), &(array[0])));
}

expr mk_add(expr_vector args) {
    vector<Z3_ast> array;

    for (int i = 0; i < args.size(); i++)
      array.push_back(args[i]);

    return to_expr(args.ctx(), Z3_mk_add(args.ctx(), array.size(), &(array[0])));
}

vector<vector<int>> solveLinearEquations(vector<vector<int>> linEquations, vector<int> vect, int from, int to){
    context c;
    const unsigned N = linEquations[0].size();
    expr_vector x(c);
    for (unsigned i = 0; i < N; i++) { 
        std::stringstream x_name; 
        x_name << "x_" << i;
        x.push_back(c.int_const(x_name.str().c_str()));
    }

    solver s(c);

    for (unsigned i = 0; i < N; i++) {
        s.add(x[i] >= from);
        s.add(x[i] <= to);
    }

    for (int i = 0; i < linEquations.size(); i++) {
        expr_vector linVector(c);
        for (int j = 0; j< linEquations[i].size(); j++) {
            if (linEquations[i][j] != 0) {
                linVector.push_back(linEquations[i][j] * x[j]);
            }
        }
        if (!linVector.empty()) {
            s.add(mk_add(linVector) == vect[i]);
        }
    }

    vector<vector<int>> solutions;

    while(true) {
        if (s.check() == sat) {
            model m = s.get_model();
            expr_vector ve(c);
            vector<int> sol; sol.clear();
            for (unsigned i = 0; i < N; i++) {
                ve.push_back(x[i] != m.eval(x[i]));
                int val;
                Z3_get_numeral_int(c, m.eval(x[i]), &val);
                sol.push_back(val);
            }
            solutions.push_back(sol); 
            s.add(mk_or(ve));
        }
        else {
            break; 
        }
    }
    return solutions;
}