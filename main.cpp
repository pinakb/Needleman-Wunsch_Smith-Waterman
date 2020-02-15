#include "alignment.h"
#include <iostream>
#include <fstream>
#include <string>
using namespace std;

int main() {
    //enter the input file path
    string inputPath = "F:\\Abhyas\\WSU\\Sem 2\\Computational Genomics\\Homework\\dataCoding\\colorBlindness.fasta";
    //enter the parameter file path
    string paraPath = "F:\\Abhyas\\WSU\\Sem 2\\Computational Genomics\\Homework\\dataCoding\\parameter.txt";
    //enter the path for output file
    string outPut = "F:\\Abhyas\\WSU\\Sem 2\\Computational Genomics\\Homework\\dataCoding\\Report.txt";

    alignment A;
    A.parameterReader(paraPath);
    A.inputReader(inputPath);
    //A.globalAlignment(outPut);
    A.localAlignment(outPut);

}
