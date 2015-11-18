#include <ctime>
#include <iostream>

#include "main.hpp"

using namespace std;

void small_tests() {
    clock_t begin = clock();
    if (generate(15, 6, 1, 3, 3, 3, false) == 1) cout << "SRG(15,6,1,3,3,3) OK" << endl;
    else cout << "SRG(15,6,1,3,3,3) FAIL" << endl;

    if (generate(25, 12, 5, 6, 3, 1, false) == 9) cout << "SRG(25,12,5,6,3,1) OK" << endl;
    else cout << "SRG(25,12,5,6,3,1) FAIL" << endl;

    if (generate(25, 12, 5, 6, 3, 4, false) == 4) cout << "SRG(25,12,5,6,3,4) OK" << endl;
    else cout << "SRG(25,12,5,6,3,4) FAIL" << endl;

    if (generate(29, 14, 6, 7, 3, 5, false) == 1) cout << "SRG(29,14,6,7,3,5) OK" << endl;
    else cout << "SRG(29,14,6,7,3,5) FAIL" << endl;

    if (generate(26, 10, 3, 4, 3, 2, false) == 3) cout << "SRG(26,10,3,4,3,2) OK" << endl;
    else cout << "SRG(26,10,3,4,3,2) FAIL" << endl;

    if (generate(26, 10, 3, 4, 3, 8, false) == 0) cout << "SRG(26,10,3,4,3,8) OK" << endl;
    else cout << "SRG(26,10,3,4,3,8) FAIL" << endl;

    if (generate(28, 12, 6, 4, 3, 1, false) == 4) cout << "SRG(26,12,6,4,3,1) OK" << endl;
    else cout << "SRG(26,12,6,4,3,1) FAIL" << endl;

    if (generate(28, 12, 6, 4, 3, 4, false) == 0) cout << "SRG(26,12,6,4,3,4) OK" << endl;
    else cout << "SRG(26,12,6,4,3,4) FAIL" << endl;

    if (generate(28, 12, 6, 4, 3, 7, false) == 0) cout << "SRG(26,12,6,4,3,7) OK" << endl;
    else cout << "SRG(26,12,6,4,3,7) FAIL" << endl;

    if (generate(28, 12, 6, 4, 3, 10, false) == 2) cout << "SRG(26,12,6,4,3,10) OK" << endl;
    else cout << "SRG(26,12,6,4,3,10) FAIL" << endl;

    if (generate(35, 16, 6, 8, 3, 2, false) == 18) cout << "SRG(35,16,6,8,3,2) OK" << endl;
    else cout << "SRG(35,16,6,8,3,2) FAIL" << endl;

    if (generate(35, 16, 6, 8, 3, 5, false) == 3) cout << "SRG(35,16,6,8,3,5) OK" << endl;
    else cout << "SRG(35,16,6,8,3,5) FAIL" << endl;

    if (generate(35, 16, 6, 8, 3, 8, false) == 1) cout << "SRG(35,16,6,8,3,8) OK" << endl;
    else cout << "SRG(35,16,6,8,3,8) FAIL" << endl;

    if (generate(35, 18, 9, 9, 3, 2, false) == 18) cout << "SRG(35,18,9,9,3,2) OK" << endl;
    else cout << "SRG(35,18,9,9,3,2) FAIL" << endl;

    if (generate(35, 18, 9, 9, 3, 5, false) == 3) cout << "SRG(35,18,9,9,3,5) OK" << endl;
    else cout << "SRG(35,18,9,9,3,5) FAIL" << endl;

    if (generate(35, 18, 9, 9, 3, 8, false) == 1) cout << "SRG(35,18,9,9,3,8) OK" << endl;
    else cout << "SRG(35,18,9,9,3,8) FAIL" << endl;

    if (generate(36, 14, 4, 6, 3, 3, false) == 3) cout << "SRG(36,14,4,6,3,3) OK" << endl;
    else cout << "SRG(36,14,4,6,3,3) FAIL" << endl;

    if (generate(36, 14, 4, 6, 3, 6, false) == 0) cout << "SRG(36,14,4,6,3,6) OK" << endl;
    else cout << "SRG(36,14,4,6,3,3) FAIL" << endl;

    if (generate(36, 14, 4, 6, 3, 9, false) == 1) cout << "SRG(36,14,4,6,3,9) OK" << endl;
    else cout << "SRG(36,14,4,6,3,9) FAIL" << endl;

    if (generate(36, 15, 6, 6, 3, 3, false) == 30) cout << "SRG(36,15,6,6,3,3) OK" << endl;
    else cout << "SRG(36,15,6,6,3,3) FAIL" << endl;

    if (generate(36, 15, 6, 6, 3, 6, false) == 4) cout << "SRG(36,15,6,6,3,6) OK" << endl;
    else cout << "SRG(36,15,6,6,3,6) FAIL" << endl;

    if (generate(36, 15, 6, 6, 3, 9, false) == 1) cout << "SRG(36,15,6,6,3,9) OK" << endl;
    else cout << "SRG(36,15,6,6,3,9) FAIL" << endl;

    clock_t end = clock();
    double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    cout << "time: " << elapsed_secs << endl;
}

void all_tests() {
    clock_t begin = clock();
    small_tests();
    
    if (generate(40, 12, 2, 4, 3, 1, false) == 9) cout << "SRG(40,12,2,4,3,1) OK" << endl;
    else cout << "SRG(40,12,2,4,3,1) FAIL" << endl;

    if (generate(40, 12, 2, 4, 3, 4, false) == 5) cout << "SRG(40,12,2,4,3,4) OK" << endl;
    else cout << "SRG(40,12,2,4,3,4) FAIL" << endl;

    if (generate(40, 12, 2, 4, 3, 7, false) == 1) cout << "SRG(40,12,2,4,3,7) OK" << endl;
    else cout << "SRG(40,12,2,4,3,7) FAIL" << endl;

    if (generate(40, 12, 2, 4, 3, 13, false) == 1) cout << "SRG(40,12,2,4,3,13) OK" << endl;
    else cout << "SRG(40,12,2,4,3,13) FAIL" << endl;

    if (generate(45, 12, 3, 3, 3, 3, false) == 5) cout << "SRG(45,12,3,3,3,3) OK" << endl;
    else cout << "SRG(45,12,3,3,3,3) FAIL" << endl;

    if (generate(45, 12, 3, 3, 3, 6, false) == 6) cout << "SRG(45,12,3,3,3,6) OK" << endl;
    else cout << "SRG(45,12,3,3,3,6) FAIL" << endl;

    if (generate(45, 12, 3, 3, 3, 9, false) == 2) cout << "SRG(45,12,3,3,3,9) OK" << endl;
    else cout << "SRG(45,12,3,3,3,9) FAIL" << endl;

    if (generate(37, 18, 8, 9, 3, 1, false) == 18) cout << "SRG(37,18,8,9,3,1) OK" << endl;
    else cout << "SRG(36,18,8,9,3,1) FAIL" << endl;

    if (generate(41, 20, 9, 10, 3, 5, false) == 18) cout << "SRG(41,20,9,10,3,5) OK" << endl;
    else cout << "SRG(41,20,9,10,3,5) FAIL" << endl;

    if (generate(45, 22, 10, 11, 3, 9, false) == 7) cout << "SRG(45,22,10,11,3,9) OK" << endl;
    else cout << "SRG(45,22,10,11,3,9) FAIL" << endl;

    if (generate(49, 18, 7, 6, 3, 1, false) == 595) cout << "SRG(49,18,7,6,3,1) OK" << endl;
    else cout << "SRG(49,18,7,6,3,1) FAIL" << endl;

    if (generate(49, 18, 7, 6, 3, 4, false) == 107) cout << "SRG(49,18,7,6,3,4) OK" << endl;
    else cout << "SRG(49,18,7,6,3,4) FAIL" << endl;

    if (generate(49, 18, 7, 6, 3, 7, false) == 4) cout << "SRG(49,18,7,6,3,7) OK" << endl;
    else cout << "SRG(49,18,7,6,3,7) FAIL" << endl;

    if (generate(49, 24, 11, 12, 3, 1, false) == 5029) cout << "SRG(49,24,11,12,3,1) OK" << endl;
    else cout << "SRG(49,24,11,12,3,1) FAIL" << endl;

    if (generate(49, 24, 11, 12, 3, 4, false) == 124) cout << "SRG(49,24,11,12,3,4) OK" << endl;
    else cout << "SRG(49,24,11,12,3,4) FAIL" << endl;

    if (generate(49, 24, 11, 12, 3, 7, false) == 40) cout << "SRG(49,24,11,12,3,7) OK" << endl;
    else cout << "SRG(49,24,11,12,3,7) FAIL" << endl;

    if (generate(50, 21, 8, 9, 3, 2, false) == 5043) cout << "SRG(50,21,8,9,3,2) OK" << endl;
    else cout << "SRG(50,21,8,9,3,2) FAIL" << endl;

    if (generate(50, 21, 8, 9, 3, 5, false) == 106) cout << "SRG(50,21,8,9,3,5) OK" << endl;
    else cout << "SRG(50,21,8,9,3,5) FAIL" << endl;

    if (generate(50, 21, 8, 9, 3, 8, false) == 35) cout << "SRG(50,21,8,9,3,8) OK" << endl;
    else cout << "SRG(50,21,8,9,3,8) FAIL" << endl;

    if (generate(53, 26, 15, 16, 3, 5, false) == 1229) cout << "SRG(53,26,15,16,3,5) OK" << endl;
    else cout << "SRG(53,26,15,16,3,5) FAIL" << endl;

    clock_t end = clock();
    double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    cout << "time: " << elapsed_secs << endl;
}
