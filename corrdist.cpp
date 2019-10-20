/////////////////////////////////////////////////////////////////////////////////
// Program name     : CorrDist.cpp
//
// Version          : 1.00
//
// Author           : Lars S Jermiin
//
// Institution      : Australian National University
//                    Research School of Biology
//                    Acton, ACT 2601, Australia
//
//                    School of Biology and Environmental Sciences
//                    University College Dublin
//                    Belfield, Dublin 4, Ireland
//
// Date begun       : 14 December, 2018
//
// Date modified    : 20 October, 2019
//
// Copyright        : Copyright © 2019 Lars Sommer Jermiin. All rights reserved.
//
// Responsibility   : The copyright holder takes no legal responsibility for the
//                    correctness of results obtained using this program.
//
// Summary          : CorrDist computed the Corrected Distance between pairs
//                    of sequences of nucleotides, 10- and 14-state genotype
//                    sequences, and amino acids.
//
//                    Sequences must be stored in the FASTA format.
//
//                    Characters are translated to integers to speed up the program.
//
// Nucleotides      : Alphabet: [A,C.G,T/U,-] = [0,1,2,3,4].
//
//                    Ambiguous characters (i.e., ?, N, B, D, H, K, M, R, S, V, W and
//                    Y) are treated as if they were alignment gaps (-) (i.e., as
//                    missing data).
//
// Amino acids      : Alphabet: [A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,-] =
//                    [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]
//
//                    Ambiguous characters (i.e., ?, X and Z) are treated as if they
//                    were alignment gaps (-) (i.e., as missing data).
//
// SNPs (10 states) : Alphabet: [A,C,G,K,M,R,S,T/U,W,Y,-] = [0,1,2,3,4,5,6,7,8,9,10].
//
//                    Ambiguous characters (i.e., ?, N, B, D, G and V) are treated as
//                    if they were alignment gaps (-) (i.e., as missing data).
//
// SNPs (14 states) : Alphabet: [A,C,G,T/U,K,M,R,S,W,Y,B,D,H,V,-] =
//                    [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14].
//
//                    Ambiguous characters (i.e., ? and N) are treated as if they were
//                    alignment gaps (-) (i.e., as missing data).
//
// Distance metrics : F81 (Felsenstein 1981); LogDet (Lake 1994; Lockhart et al. 1994)
//
// Manuscript       : Jermiin LS, et al. (2020). Phylogenetic analysis of genotype
//                    sequence data. Syst. Biol. In prep.
//
// Reference        : Zhang...
/////////////////////////////////////////////////////////////////////////////////

#include <cctype>
#include <cmath>
#include <string>
#include <vector>
#include <random>
#include <iomanip>
#include <fstream>
#include <iostream>

#define SQR(a) ((a) * (a))
#define REP_MAX 1000

// The following are declared here because they are needed in different functions

std::vector<std::string> taxon;                // 2D container for sequence names
std::vector<std::vector<int> > alignment;      // 2D container for sequence data
std::vector<unsigned long> sites;
const unsigned FOUR(4);        // for 4-state alphabet (DNA)
const unsigned TEN(10);        // for 10-state alphabet (SNP data)
const unsigned FOURTEEN(14);   // for 14-state alphabet (SNP data)
const unsigned TWENTY(20);     // for 20-state alphabet (amino acids)
const unsigned max_array(21);  // for all alphabets plus gaps


// This function translates a string of characters into a vector of integers
std::vector<int> Translator(unsigned datatype, std::string seq) {
    std::vector<int> seq_data;
    
    switch (datatype) {
        case 1: // Nucleotides (A|C|G|T)
            for (std::string::size_type i = 0; i != seq.size(); ++i) {
                switch (toupper(seq[i])) {
                    case 'A': seq_data.push_back(0); break;
                    case 'C': seq_data.push_back(1); break;
                    case 'G': seq_data.push_back(2); break;
                    case 'T': seq_data.push_back(3); break;
                    case 'U': seq_data.push_back(3); break;
                    default : seq_data.push_back(4); break; // In case of other characters
                }
            }
            break;
        case 2: // 10-state genotype data
            for (std::string::size_type i = 0; i != seq.size(); ++i) {
                switch (toupper(seq[i])) {
                    case 'A': seq_data.push_back(0); break;
                    case 'C': seq_data.push_back(1); break;
                    case 'G': seq_data.push_back(2); break;
                    case 'T': seq_data.push_back(3); break;
                    case 'U': seq_data.push_back(3); break;
                    case 'K': seq_data.push_back(4); break;
                    case 'M': seq_data.push_back(5); break;
                    case 'R': seq_data.push_back(6); break;
                    case 'Y': seq_data.push_back(7); break;
                    case 'S': seq_data.push_back(8); break;
                    case 'W': seq_data.push_back(9); break;
                    default : seq_data.push_back(10);break; // In case of other characters
                }
            }
            break;
        case 3: // 14-state genotype data
            for (std::string::size_type i = 0; i != seq.size(); ++i) {
                switch (toupper(seq[i])) {
                    case 'A': seq_data.push_back(0); break;
                    case 'C': seq_data.push_back(1); break;
                    case 'G': seq_data.push_back(2); break;
                    case 'T': seq_data.push_back(3); break;
                    case 'U': seq_data.push_back(3); break;
                    case 'K': seq_data.push_back(4); break;
                    case 'M': seq_data.push_back(5); break;
                    case 'R': seq_data.push_back(6); break;
                    case 'Y': seq_data.push_back(7); break;
                    case 'S': seq_data.push_back(8); break;
                    case 'W': seq_data.push_back(9); break;
                    case 'B': seq_data.push_back(10);break;
                    case 'D': seq_data.push_back(11);break;
                    case 'H': seq_data.push_back(12);break;
                    case 'V': seq_data.push_back(13);break;
                    default : seq_data.push_back(14);break; // In case of other characters
                }
            }
            break;
        default: // amino acids (A|G|P|S|T|D|E|N|Q|H|K|R|M|I|V|L|W|F|Y|C)
            for (std::string::size_type i = 0; i != seq.size(); ++i) {
                switch (toupper(seq[i])) {
                    case 'A': seq_data.push_back(0); break;
                    case 'C': seq_data.push_back(1); break;
                    case 'D': seq_data.push_back(2); break;
                    case 'E': seq_data.push_back(3); break;
                    case 'F': seq_data.push_back(4); break;
                    case 'G': seq_data.push_back(5); break;
                    case 'H': seq_data.push_back(6); break;
                    case 'I': seq_data.push_back(7); break;
                    case 'K': seq_data.push_back(8); break;
                    case 'L': seq_data.push_back(9); break;
                    case 'M': seq_data.push_back(10);break;
                    case 'N': seq_data.push_back(11);break;
                    case 'P': seq_data.push_back(12);break;
                    case 'Q': seq_data.push_back(13);break;
                    case 'R': seq_data.push_back(14);break;
                    case 'S': seq_data.push_back(15);break;
                    case 'T': seq_data.push_back(16);break;
                    case 'V': seq_data.push_back(17);break;
                    case 'W': seq_data.push_back(18);break;
                    case 'Y': seq_data.push_back(19);break;
                    default : seq_data.push_back(20);break; // In case of other characters
                }
            }
            break;
    }
    return(seq_data);
}


// Function that reads input file and stores data in two 2D containers
void Read_Input(std::string inname, unsigned datatype){
    std::ifstream infile;
    std::string seq(""), str(""), tmp(""); // temporary string used to store input
    std::vector<int> sequence;     // temporary vector used to store input
    unsigned long counter(0), alignment_length(0);
    
    infile.open(inname.c_str());
    if (!infile) {
        std::cerr << "Input file not found" << std::endl;
        exit(1);
    }
    while (getline(infile, str)) {
        if (!str.empty()) {
            // remove blank space in string
            tmp.clear();
            for (std::string::size_type i = 0; i != str.size(); ++i) {
                if (!isblank(str[i])) {
                    tmp.push_back(str[i]);
                }
            }
            if (tmp[0] == '>') {
                if (seq.size() > 0) {
                    sequence = Translator(datatype, seq);
                    alignment.push_back(sequence); // stores sequence in vector
                    if (alignment_length == 0)
                        alignment_length = sequence.size();
                    sequence.clear();
                    seq.clear();
                }
                tmp.erase(tmp.begin()); // removes first character from name
                taxon.push_back(tmp); // stores sequence name in vector
            } else {
                seq += tmp;
            }
            str.clear();
        }
    }
    // Store last sequence in vector
    if (seq.size() > 0) {
        sequence = Translator(datatype, seq);
        alignment.push_back(sequence);
    } else {
        std::cerr << "Last sequence empty" << "\n" << std::endl;
        exit(1);
    }
    //Check whether the sequence names are unique
    for (std::vector<std::string>::const_iterator iter1 = taxon.begin(); iter1 != taxon.end(); ++iter1) {
        for (std::vector<std::string>::const_iterator iter2 = iter1 + 1; iter2 != taxon.end(); ++iter2) {
            if (*iter1 == *iter2) {
                std::cerr << "Program aborted due to error: Sequence name not unique -- look for " << *iter1 << "\n" << std::endl;
                exit(1);
            }
        }
    }
    // Check whether the sequences have the same length
    for (std::vector<std::vector<int> >::const_iterator iter = alignment.begin()+1; iter != alignment.end(); ++iter) {
        ++counter;
        sequence = *iter;
        if (sequence.size() != alignment_length) {
            std::cerr << "Program aborted due to error: sequences 1 and " << counter << " differ in length!\n" << std::endl;
            exit(1);
        }
    }
}


// Function that reads input file and stores data in two 2D containers

double Det(unsigned n, double mat[max_array][max_array]) {
    double l[max_array][max_array];
    double u[max_array][max_array];
    
    // LU decomposition of square matrix
    for (unsigned i = 0; i != n; i++) {
        for (unsigned j = 0; j != n; j++) {
            if (j < i)
                l[j][i] = 0.0;
            else {
                l[j][i] = mat[j][i];
                for (unsigned k = 0; k != i; k++) {
                    l[j][i] = l[j][i] - l[j][k] * u[k][i];
                }
            }
        }
        for (unsigned j = 0; j != n; j++) {
            if (j < i)
                u[i][j] = 0.0;
            else if (j == i)
                u[i][j] = 1;
            else {
                u[i][j] = mat[i][j] / l[i][i];
                for (unsigned k = 0; k != i; k++) {
                    u[i][j] = u[i][j] - ((l[i][k] * u[k][j]) / l[i][i]);
                }
            }
        }
    }
    // compute determinant from
    double detL = 1.0;
    for (unsigned i = 0; i != n; i++) {
        detL *= l[i][i];
    }
    double detU = 1.0;
    for (unsigned i = 0; i != n; i++) {
        detU *= u[i][i];
    }
    return detL * detU;
}


//Main program that calls the various functions above

int main(int argc, char** argv){
    std::ifstream infile;
    std::ofstream outfile;
    std::string inName, outName;
    std::string model, mode, nature_of_data;
    std::vector<int> sequence;     // temporary vector used to store input
    std::vector<double> row_of_double;
    std::vector<std::vector<double> > matDist;
    unsigned dataType(0), rows_columns (0), repeat_max(REP_MAX), seed;
    unsigned long sum_dm(0), site(0);
    unsigned long dm[max_array][max_array];   // 2D divergence matrix of integers
    double logdet(0.0), F81(0.0);
    double row_sumd[max_array];
    double col_sumd[max_array];
    double dmd[max_array][max_array];        // 2D divergence matrix of reals
    
// Start by checking format and content of input
    if(argc != 5){
        std::cerr << "\nCorrDist Copyright 2019, Lars Jermiin" << std::endl;
        std::cerr << " Contact: lars.jermiin [at] anu.edu.au / ucd.ie" << std::endl;
        std::cerr << "\nERROR -- use command: corrdist <infile> <F81|LD> <Y|N> <1|...|4>\n" << std::endl;
        std::cerr << "  infile   Fasta-formatted alignment" << std::endl;
        std::cerr << "     F81   Correction using Felsenstein (1981) model" << std::endl;
        std::cerr << "      LD   Correction using LogDet/Paralinear model" << std::endl;
        std::cerr << "     Y|N   Nonparametric boostrapping [Y|N] with " << REP_MAX << " replicates" << std::endl;
        std::cerr << "       1   Nucleotides;  4 states (A|C|G|T)" << std::endl;
        std::cerr << "       2   Genotypes;   10 states (A|C|G|T|K|M|R|Y|S|W)" << std::endl;
        std::cerr << "       3   Genotypes;   14 states (A|C|G|T|K|M|R|Y|S|W|B|D|H|V)" << std::endl;
        std::cerr << "       4   Amino acids; 20 states (A|G|P|S|T|D|E|N|Q|H|K|R|M|I|V|L|W|F|Y|C)" << std::endl;
        std::cerr << std::endl;
        exit(1);
    }
    inName = argv[1];
    model = argv[2];
    mode = argv[3];
    nature_of_data = argv[4];
    infile.open(inName.c_str());
    // Checking existence of infile
    if (!infile) {
        std::cerr << "\nInput file not found\n" << std::endl;
        exit(1);
    }
    infile.close();
    for (std::string::size_type i = 0; i != model.size(); ++i) {
        model[i] = toupper(model[i]);
    }
    // Checking model of correction
    if (model != "F81" && model != "LD") {
        std::cerr << "\nIncorrect choice of model of correction [F81|LD]\n" << std::endl;
        exit(1);
    }
    // Checking mode of analysis
    if (toupper(mode[0]) != 'Y' && toupper(mode[0]) != 'N') {
        std::cerr << "\nIncorrect choice of analysis: [Y|N]\n" << std::endl;
        exit(1);
    }
    dataType = stoi(nature_of_data);
    // Checking type of data
    if (dataType < 1 || dataType > 4) {
        std::cerr << "\nIncorrect choice of data: [1|...|4]\n" << std::endl;
        exit(1);
    }

// Prepare output file name
    outName.clear();
    for (std::string::size_type i = 0; i != inName.size() && inName[i] != '.'; ++i) {
        outName += inName[i];
    }
    if (toupper(mode[0]) == 'N') {
        outName = outName + "_real.dis";
        repeat_max = 1;
    } else {
        outName = outName + "_boot.dis";
        repeat_max = REP_MAX;
    }

// Set the sizes of matrices and vectors to hold interim data
    switch (dataType) {
        case 1: rows_columns = FOUR; break;
        case 2: rows_columns = TEN; break;
        case 3: rows_columns = FOURTEEN; break;
        default: rows_columns = TWENTY; break;
    }

// Read and check input file
    Read_Input(inName, dataType);

    std::cout << "\nData read ..." << std::endl;
// Retrieve first sequence so we know the length of the alignment
    sequence = alignment[0];

// Random number generator; required for bootstrapping
    // Get seed from system clock
    seed = (unsigned)std::chrono::system_clock::now().time_since_epoch().count();
    // Initialize randomisation in case it is n
    std::mt19937 generator(seed);
    // Generate a uniform distribution of random numbers between 0 and sequence.size() - 1
    std::uniform_int_distribution<unsigned long> distribution(0,sequence.size()-1);

// Create a 2D vector to hold distance values
    // Create and initialise 1D vector to hold distance values
    for (std::vector<double>::size_type i = 0; i != alignment.size(); i++) {
        row_of_double.push_back(0.0);
    }
    // Create and initialise 2D vector to hold distance values
    for (std::vector<std::vector<double> >::size_type i = 0; i != alignment.size(); i++) {
        matDist.push_back(row_of_double);
    }

// Open input file
    outfile.open(outName.c_str());

// Create vector sites so that its elements contain 1 in every cell
    for (std::vector<int>::size_type i = 0; i != sequence.size(); ++i) {
        sites.push_back(1);
    }

    std::cout << "\nDistances being computed ..." << std::endl;

// Start loop to compute distance matrix/matrices
    for (unsigned repeat = 0; repeat != repeat_max; ++repeat) {
        if (repeat_max > 1) {
            // Bootstrap mode
            for (std::vector<int>::size_type i = 0; i != sequence.size(); ++i) {
                sites[i] = 0;
            }
            // vector sites now contains 0 in every element
            for (std::vector<int>::size_type i = 0; i != sequence.size(); ++i) {
                site = distribution(generator);
                ++sites[site];
            }
        }
        
// Start doing calculations of corrected distances between pairs of sequences
        for (std::vector<std::vector<int> >::size_type iter1 = 0; iter1 != alignment.size(); ++iter1) {
            for (std::vector<std::vector<int> >::size_type iter2 = iter1 + 1; iter2 != alignment.size(); ++iter2) {
                // set all elements in the divergence matrix to zero
                for (size_t i = 0; i != max_array; ++i) {
                    for (size_t j = 0; j != max_array; ++j) {
                        dm[i][j] = 0;
                    }
                }
                // generate the divergence matrix
                for (std::vector<int>::size_type i = 0; i != sites.size(); ++i) {
                    dm[alignment[iter1][i]][alignment[iter2][i]] =
                    dm[alignment[iter1][i]][alignment[iter2][i]] + sites[i];
                }

                // generate the sum over all elements in the divergence matrix
                sum_dm = 0;
                for (size_t i = 0; i != rows_columns; ++i) {
                    for (size_t j = 0; j != rows_columns; ++j) {
                        sum_dm += dm[i][j];
                    }
                }
                // generate normalized divergence matrix
                for (size_t i = 0; i != rows_columns; ++i) {
                    for (size_t j = 0; j != rows_columns; ++j) {
                        dmd[i][j] = (double)dm[i][j]/sum_dm;
                    }
                }
                // set all elements in marginal vectors to zero
                for (size_t i = 0; i != max_array; ++i) {
                    row_sumd[i] = 0.0;
                    col_sumd[i] = 0.0;
                }
                // compute marginal frequencies
                for (size_t i = 0; i != rows_columns; ++i) {
                    for (size_t j = 0; j != rows_columns; ++j) {
                        row_sumd[i] += dmd[i][j];                }
                }
                for (size_t i = 0; i != rows_columns; ++i) {
                    for (size_t j = 0; j != rows_columns; ++j) {
                        col_sumd[i] += dmd[j][i];
                    }
                }
                if(model == "F81") {
                    // calculate F81 distance
                    double limit(1.0);
                    double hamming(1.0);
                    for (size_t i = 0; i != rows_columns; ++i) {
                        limit -= (SQR((row_sumd[i] + col_sumd[i])/2.0));
                        hamming -= dmd[i][i];
                    }
                    F81 = -1.0 * limit * log(1.0 - hamming/limit);
                    matDist[iter1][iter2] = F81;
                    matDist[iter2][iter1] = F81;
                } else {
                    // calculate LogDet distance
                    logdet = Det(rows_columns, dmd);
                    logdet = log(logdet);
                    for (size_t i = 0; i != rows_columns; ++i) {
                        logdet = logdet - 0.5 * (log(row_sumd[i]) + log(col_sumd[i]));
                    }
                    logdet = -1.0 * logdet/rows_columns;
                    matDist[iter1][iter2] = logdet;
                    matDist[iter2][iter1] = logdet;
                }
            }
        }

        // Printing output
        outfile << alignment.size() << std::endl;
        for (std::vector<std::vector<double> >::size_type i = 0 ; i != taxon.size(); i++) {
            outfile << std::left << std::setw(10) << taxon[i];
            for (std::vector<std::vector<double> >::size_type j = 0; j != taxon.size(); j++) {
                outfile << std::fixed << matDist[i][j];
                if (j != (taxon.size() - 1)) {
                    outfile << "\t";
                }
            }
            outfile << std::endl;
        }
    }
    outfile.close();
    std::cout << std::endl;
    std::cout << "--------------------------------------------------------------------" << std::endl;
    std::cout << "   RESULTS FROM ANALYSIS COMPLETE" << std::endl;
    std::cout << std::endl;
    if (repeat_max == 1) {
        std::cout << "   Matrix with " << model <<" distance estimates in .... " << outName << std::endl;
    } else {
        std::cout << "   Matrices with " << model <<" distance estimates in .. " << outName << std::endl;
    }
    std::cout << "--------------------------------------------------------------------" << std::endl;
    return 0;
}
