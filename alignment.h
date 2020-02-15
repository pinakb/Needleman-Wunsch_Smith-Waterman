#ifndef ALIGNMENT_H
#define ALIGNMENT_H
#include <iostream>
#include <fstream>
#include <string>
using namespace std;

class alignment
{

public:
    //cell Structure Declaration
    struct cell {
        short score;
        short diag;
        short left;
        short above;
    };

    //Function Prototype declaration
    alignment();
    void parameterReader(string fileName);
    void inputReader(string fileName);
    int globalAlignment(string outPut);
    int localAlignment(string outPut);



private:
    //Variable Declaration
    string seq1;
    string seq2;
    string alignA;
    string alignB;
    string localAlignA;
    string localAlignB;
    string delimiter;
    int seq1_len;
    int seq2_len;
    short m;
    short mm;
    short co_i;
    short co_j;
    short oGap;
    short gap_ext;
    short matchCount;
    short mismatchCount;
    short gapCount;
    short openingGapCount;
    int nucPerLine;
    unsigned long counter;
    ifstream parafile;
    ofstream outfile;




    //Function Prototype declaration

    void printMat(cell** scoreMat, string seq1, string seq2, int seq1_len, int seq2_len);
    void printSequences(string alignmentA, string alignmentB);

    //Global alignment prototypes
    void initGlobal(cell** scoreMat, int seq1_len, int seq2_len, short gap, short extraGap);
    int gAlign(cell** scoreMat, string seq1, string seq2, int seq1_len, int seq2_len, short m, short mm, short gap, short extraGap);
    int traceback(cell** scoreMat, string seq1, string seq2, string& alignA, string& alignB, int seq1_len, int seq2_len);

    //Local Alignment Prototypes
    void initlocal(cell** scoreMat, int seq1_len, int seq2_len, short oGap, short gap_ext);
    int lAlign(cell** scoreMat, string seq1, string seq2, int seq1_len, int seq2_len, short m, short mm, short oGap, short gap_ext);
    int maxMatrixValue(cell** scoreMat);
    int localTraceback(cell** scoreMat, string seq1, string seq2, string& alignA, string& alignB, int seq1_len, int seq2_len);

};

#endif // ALIGNMENT_H
#pragma once
