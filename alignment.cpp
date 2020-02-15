#include "alignment.h"
#include <iostream>
#include <fstream>
#include <string>
#include<algorithm>
#include <cstdlib>

using namespace std;

alignment::alignment()
{
    alignA = "";
    alignB = "";
    localAlignA = "";
    localAlignB = "";
    delimiter = "";
    matchCount = 0;
    mismatchCount = 0;
    gapCount = 0;
    openingGapCount = 0;
    nucPerLine = 60;
    counter = 1;
}

void alignment::parameterReader(string fileName) {
    parafile.open(fileName.c_str(), ios::in);
    if (!parafile) {
        cout << "error in opening parameter file" << endl;
    }
    else {
        parafile >> m;
        parafile >> mm;
        parafile >> oGap;
        parafile >> gap_ext;
    }
    parafile.close();
    cout << "match penalty= " << m << endl;
    cout << "mismatch penalty= " << mm << endl;
    cout << "gap penalty= " << oGap << endl;
    cout << "ext_gap penalty= " << gap_ext << endl;
}

void alignment::inputReader(string fileName) {

    ifstream file(fileName.c_str());
    if (!file.good()) {
        cerr << "Error opening '" << fileName << "'. Bailing out." << endl;
        //return -1;
    }

    string line, name, content;
    while (getline(file, line).good())
    {
        if (line.empty() || line[0] == '>')
        { // Identifier marker
            if (!name.empty()) { // Print out what we read from the last entry
                //cout<<content << endl;
                seq1 += content;
                name.clear();
            }
            if (!line.empty()) {
                name = line.substr(1);
            }
            content.clear();
        }
        else if (!name.empty())
        {

            if (line.find(' ') != string::npos) { // Invalid sequence--no spaces allowed
                name.clear();
                content.clear();
            }
            else {
                content += line;
            }
        }
    }
    if (!name.empty())
    { // Print out what we read from the last entry
        //cout <<content << endl;
        seq2 += content;
    }



    //file.close();
    seq1_len = seq1.length();
    seq2_len = seq2.length();

    //    cout << "seq1: " << seq1 << endl;
    cout << "seq1 length: " << seq1_len << endl;
    //    cout << "seq2: " << seq2 << endl;
    cout << "seq2 length: " << seq2_len << endl;
}

int alignment::globalAlignment(string outPut)
{
    cell** scoreMat = new cell * [seq1_len + 1];
    for (int i = 0; i <= seq1_len; i++)
    {
        scoreMat[i] = new cell[seq2_len];
    }

    initGlobal(scoreMat, seq1_len, seq2_len, oGap, gap_ext);
    gAlign(scoreMat, seq1, seq2, seq1_len, seq2_len, m, mm, oGap, gap_ext);
    //    printMat(scoreMat, seq1, seq2, seq1_len, seq2_len);
    traceback(scoreMat, seq1, seq2, alignA, alignB, seq1_len, seq2_len);
    int optimal_score = scoreMat[seq1_len][seq2_len].score;

    //report creation
    outfile.open(outPut.c_str(), ios::out);
    cout << "optimal score: " << optimal_score << endl;
    outfile << "optimal score: " << optimal_score << endl;
    cout << "length alignA: " << alignA.length() << endl;
    outfile << "length alignA: " << alignA.length() << endl;
    cout << "length alignB: " << alignB.length() << endl;
    outfile << "length alignB: " << alignB.length() << endl;
    cout << "Matches: " << matchCount << endl;
    outfile << "Matches: " << matchCount << endl;
    cout << "Mismatches: " << mismatchCount << endl;
    outfile << "Mismatches: " << mismatchCount << endl;
    cout << "opening gaps = " << openingGapCount << " gaps= " << gapCount << endl;
    outfile << "opening gaps = " << openingGapCount << " gaps= " << gapCount << endl;
    cout << "% identity " << matchCount << "/" << alignA.length() << " = " << ((100.0 * matchCount) / alignA.length()) << endl;
    outfile << "% identity " << matchCount << "/" << alignA.length() << " = " << ((100.0 * matchCount) / alignA.length()) << endl;
    cout << "% gap " << gapCount << "/" << alignA.length() << " = " << ((100.0 * gapCount) / alignA.length()) << endl;
    outfile << "% gap " << gapCount << "/" << alignA.length() << " = " << ((100.0 * gapCount) / alignA.length()) << endl;
    cout << "\n";
    cout << "Alignments: " << endl;
    outfile << "Alignments: " << endl;
    printSequences(alignA, alignB);
    outfile.close();
    cout << "Report is ready. Path: " << outPut << endl;
    delete[] scoreMat;

    return 0;
}

void alignment::initGlobal(cell** scoreMat, int seq1_len, int seq2_len, short gap, short extraGap)
{
    scoreMat[0][0].score = 0;

    for (int i = 1; i <= seq1_len; i++) {
        scoreMat[i][0].score = oGap + (i * gap_ext);
    }
    for (int j = 1; j <= seq2_len; j++) {
        scoreMat[0][j].score = oGap + (j * gap_ext);
    }
}

int alignment::gAlign(cell** scoreMat, string seq1, string seq2, int seq1_len, int seq2_len, short m, short mm, short gap, short extraGap)
{
    short value = 0;
    short diag, above, left;
    short score;

    for (int i = 1; i <= seq1_len; i++)
    {
        for (int j = 1; j <= seq2_len; j++)
        {

            int S = (seq1[i - 1] == seq2[j - 1]) ? m : mm;

            scoreMat[i][j].diag = scoreMat[i - 1][j - 1].score + S;
            diag = scoreMat[i][j].diag;

            if (value > 0)
            {
                scoreMat[i][j].above = scoreMat[i - 1][j].score + gap_ext;
                above = scoreMat[i][j].above;
                scoreMat[i][j].left = scoreMat[i][j - 1].score + gap_ext;
                left = scoreMat[i][j].left;
            }
            else {
                scoreMat[i][j].above = scoreMat[i - 1][j].score + oGap + gap_ext;
                above = scoreMat[i][j].above;
                scoreMat[i][j].left = scoreMat[i][j - 1].score + oGap + gap_ext;
                left = scoreMat[i][j].left;
            }

            score = max(diag, max(left, above));
            scoreMat[i][j].score = score;

            //affine gap value checker
            if (score == (diag)) {
                value = 0;
            }
            else {
                value += 1;
            }

        }
    }

    return 0;
}

int alignment::traceback(cell** scoreMat, string seq1, string seq2, string& alignA, string& alignB, int seq1_len, int seq2_len)
{
    int i = seq1_len;
    int j = seq2_len;
    int diag = 0;
    int left = 0;
    int above = 0;
    bool openFlag = false;

    while (i > 0 || j > 0)
    {

        diag = scoreMat[i][j].diag; //this score is post last penalty calculation
        left = scoreMat[i][j].left; // this score is post last penalty calculation
        above = scoreMat[i][j].above; //this score is post last penalty calculation

        if (i == 0)
        {
            alignB += seq2[j - 1];
            alignA += "-";
            j--;
        }
        else if (j == 0)
        {
            alignB += "-";
            alignA += seq1[i - 1];
            i--;
        }
        else
        {
            if (diag >= above && diag >= left)
            {
                alignA += seq1[i - 1];
                alignB += seq2[j - 1];
                i--;
                j--;
            }
            else if (above > left)
            {

                alignB += "-";
                alignA += seq1[i - 1];
                i--;
            }
            else
            {
                alignB += seq2[j - 1];
                alignA += "-";
                j--;
            }

        }
    }

    for (unsigned int x = 0; x < alignA.length(); x++)
    {
        if (alignA[x] != '-' && alignB[x] != '-')
        {
            if (alignA[x] == alignB[x])
            {
                matchCount++;
                delimiter += "|";
                if (openFlag == true) {
                    openingGapCount++;
                    openFlag = false;
                }
            }
            else
            {
                mismatchCount++;
                delimiter += ".";
                if (openFlag == true) {
                    openingGapCount++;
                    openFlag = false;
                }
            }
        }
        else
        {
            gapCount++;
            delimiter += " ";
            openFlag = true;
        }

    }

    reverse(alignA.begin(), alignA.end());
    reverse(alignB.begin(), alignB.end());
    reverse(delimiter.begin(), delimiter.end());

    cout << "trace complete" << endl;
    return 0;
}


//Local Alignment

void alignment::initlocal(cell** scoreMat, int seq1_len, int seq2_len, short oGap, short gap_ext)
{
    scoreMat[0][0].score = 0;

    for (int i = 1; i <= seq1_len; i++) {
        scoreMat[i][0].score = 0;
    }
    for (int j = 1; j <= seq2_len; j++) {
        scoreMat[0][j].score = 0;
    }
}

int alignment::localAlignment(string outPut)
{
    cell** scoreMat = new cell * [seq1_len + 1];
    for (int i = 0; i <= seq1_len; i++)
    {
        scoreMat[i] = new cell[seq2_len];
    }

    initlocal(scoreMat, seq1_len, seq2_len, oGap, gap_ext);
    lAlign(scoreMat, seq1, seq2, seq1_len, seq2_len, m, mm, oGap, gap_ext);
    localTraceback(scoreMat, seq1, seq2, alignA, alignB, seq1_len, seq2_len);
    //    printMat(scoreMat, seq1, seq2, seq1_len, seq2_len);

    int optimal_score = scoreMat[co_i][co_j].score;

    //Report Creation
    outfile.open(outPut.c_str(), ios::out);

    cout << "optimal score: " << optimal_score << endl;
    outfile << "optimal score: " << optimal_score << endl;
    cout << "length localAlignA " << localAlignA.length() << endl;
    outfile << "length localAlignA " << localAlignA.length() << endl;
    cout << "length localAlignB " << localAlignB.length() << endl;
    outfile << "length localAlignB " << localAlignB.length() << endl;
    cout << "Matches" << matchCount << endl;
    outfile << "Matches" << matchCount << endl;
    cout << "Mismatches" << mismatchCount << endl;
    outfile << "Mismatches" << mismatchCount << endl;
    cout << "opening gaps = " << openingGapCount << " gaps= " << gapCount << endl;
    outfile << "opening gaps = " << openingGapCount << " gaps= " << gapCount << endl;
    cout << "% identity " << matchCount << "/" << localAlignA.length() << " = " << ((100.0 * matchCount) / localAlignA.length()) << endl;
    outfile << "% identity " << matchCount << "/" << localAlignA.length() << " = " << ((100.0 * matchCount) / localAlignA.length()) << endl;
    cout << "% gap " << gapCount << "/" << localAlignA.length() << " = " << ((100.0 * gapCount) / localAlignA.length()) << endl;
    outfile << "% gap " << gapCount << "/" << localAlignA.length() << " = " << ((100.0 * gapCount) / localAlignA.length()) << endl;
    cout << "\n" << endl;
    cout << "Aligned sequences: " << "\n" << endl;
    outfile << "Aligned sequences: " << "\n" << endl;
    printSequences(localAlignA, localAlignB);

    outfile.close();
    cout << "Report is ready. Path: " << outPut << endl;

    delete[] scoreMat;

    return 0;
}

int alignment::lAlign(cell** scoreMat, string seq1, string seq2, int seq1_len, int seq2_len, short m, short mm, short oGap, short gap_ext)
{
    short value = 0;
    short diag, above, left;
    short zero = 0;
    short score;

    for (int i = 1; i <= seq1_len; i++)
    {
        for (int j = 1; j <= seq2_len; j++)
        {

            int S = (seq1[i - 1] == seq2[j - 1]) ? m : mm;

            scoreMat[i][j].diag = scoreMat[i - 1][j - 1].score + S;
            diag = scoreMat[i][j].diag;

            if (value > 0)
            {
                scoreMat[i][j].above = scoreMat[i - 1][j].score + gap_ext;
                above = scoreMat[i][j].above;
                scoreMat[i][j].left = scoreMat[i][j - 1].score + gap_ext;
                left = scoreMat[i][j].left;
            }
            else {
                scoreMat[i][j].above = scoreMat[i - 1][j].score + oGap + gap_ext;
                above = scoreMat[i][j].above;
                scoreMat[i][j].left = scoreMat[i][j - 1].score + oGap + gap_ext;
                left = scoreMat[i][j].left;
            }

            score = max(max(diag, zero), max(left, above));
            scoreMat[i][j].score = score;

            //affine gap value checker
            if (score == (diag)) {
                value = 0;
            }
            else {
                value += 1;
            }

        }
    }

    return 0;
}

int alignment::localTraceback(cell** scoreMat, string seq1, string seq2, string& alignA, string& alignB, int seq1_len, int seq2_len)
{
    int diag = 0;
    int left = 0;
    int above = 0;
    short x, y;
    bool openFlag = false;


    maxMatrixValue(scoreMat);
    x = co_i;
    y = co_j;
    while (scoreMat[x][y].score != 0)
    {

        diag = scoreMat[x][y].diag; //this score is post last penalty calculation
        left = scoreMat[x][y].left; // this score is post last penalty calculation
        above = scoreMat[x][y].above; //this score is post last penalty calculation

        if (x == 0)
        {
            localAlignB += seq2[y - 1];
            localAlignA += "-";
            y--;
        }
        else if (y == 0)
        {
            localAlignB += "-";
            localAlignA += seq1[x - 1];
            x--;
        }
        else
        {
            if (diag >= above && diag >= left)
            {
                localAlignA += seq1[x - 1];
                localAlignB += seq2[y - 1];
                x--;
                y--;
            }
            else if (above > left)
            {

                localAlignB += "-";
                localAlignA += seq1[x - 1];
                x--;
            }
            else
            {
                localAlignB += seq2[y - 1];
                localAlignA += "-";
                y--;
            }

        }
    }

    for (unsigned int z = 0; z < localAlignA.length(); z++)
    {
        if (localAlignA[z] != '-' && localAlignB[z] != '-')
        {
            if (localAlignA[z] == localAlignB[z])
            {
                matchCount++;
                if (openFlag == true)
                {
                    openingGapCount++;
                    openFlag = false;
                }
                delimiter += '|';
            }
            else
            {
                mismatchCount++;
                if (openFlag == true) {
                    openingGapCount++;
                    openFlag = false;
                }
                delimiter += '.';

            }
        }
        else
        {
            gapCount++;
            openFlag = true;
            delimiter += ' ';
        }

    }

    reverse(localAlignA.begin(), localAlignA.end());
    reverse(localAlignB.begin(), localAlignB.end());
    reverse(delimiter.begin(), delimiter.end());
    cout << "trace complete" << endl;
    return 0;
}

int alignment::maxMatrixValue(cell** scoreMat)
{
    int value = scoreMat[0][0].score;
    int x, y;
    for (x = 0; x <= seq1_len; x++)
    {
        for (y = 0; y <= seq2_len; y++)
        {
            if (scoreMat[x][y].score > value)
            {
                value = scoreMat[x][y].score;
                co_i = x;
                co_j = y;
            }
        }
    }
    cout << "max matrix value: " << value << endl;
    cout << "coordinates i j : " << co_i << "\t" << co_j << endl;
    return value;

}

void alignment::printMat(cell** scoreMat, string seq1, string seq2, int seq1_len, int seq2_len) {

    for (int j = 0; j < seq2_len; j++)
    {
        cout << seq2[j] << "   ";
    }
    cout << "\n  ";

    for (int i = 0; i <= seq1_len; i++)
    {
        if (i > 0)
        {
            cout << seq1[i - 1] << " ";
        }
        for (int j = 0; j <= seq2_len; j++)
        {
            cout.width(3);
            cout << scoreMat[i][j].score
                << " ";
        }
        cout << endl;
    }
    cout << endl;
}


void alignment::printSequences(string alignmentA, string alignmentB)
{
    string setA = "";
    string setB = "";
    string setD = "";

    if (alignmentA.length() == alignmentB.length())
    {
        do {

            for (int x = 1; x <= nucPerLine; x++) {
                setA += alignmentA[counter];
                setD += delimiter[counter];
                setB += alignmentB[counter];
                if (counter == alignmentA.length())
                    break;
                else
                    counter++;
            }
            cout << setA << endl;
            cout << setD << endl;
            cout << setB << endl;
            outfile << setA << endl;
            outfile << setD << endl;
            outfile << setB << endl;

            setA = "";
            setD = "";
            setB = "";

        } while (counter < alignmentA.length());
    }
    else
    {
        cout << "Alignment lengths don't match";
    }
}

