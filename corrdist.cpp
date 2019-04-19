/*--------------------------------------------------------------------------
 Program name     : CorrDist.cpp
 
 Version          : 0.9
 
 Author           : Lars S Jermiin
 
 Institution      : Australian National University
                    Research School of Biology
                    Acton, ACT 2601, Australia
 
                    School of Biology and Environmental Sciences
                    University College Dublin
                    Belfield, Dublin 4
                    Ireland
 
 Date begun       : 14 December, 2018
 
 Date modified    : 1 April, 2019
 
 Copyright        : Copyright © 2019 Lars Sommer Jermiin. All rights
                    reserved.
 
 Responsibility   : The copyright holder takes no legal responsibility for
                    the correctness of results obtained using this program.
 
 Summary          : CorDis computed the Corrected Distance between pairs of
                    sequences of nucleotides, amino acids and SNPs (10- and 
                    14-state alphabets).
 
                    Sequences must be stored in the FASTA format.
 
                    Characters are converted to integers to speed up the
                    program.
 
 Nucleotides      : Alphabet: [A,C.G,T/U,-] = [0,1,2,3,4].
 
                    Ambiguous characters (i.e., ?, N, B, D, H, K, M, R, S,
                    V, W and Y) are treated as if they were alignment gaps
                    (-) (i.e., as missing data).
 
 Amino acids      : Alphabet: [A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,-] =
                    [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]
 
                    Ambiguous characters (i.e., ?, X and Z) are treated as
                    if they were alignment gaps (-) (i.e., as missing data).
 
 SNPs (10 states) : Alphabet: [A,C,G,K,M,R,S,T/U,W,Y,-] =
                    [0,1,2,3,4,5,6,7,8,9,10].
 
                    Ambiguous characters (i.e., ?, N, B, D, G and V) are
                    treated as if they were alignment gaps (-) (i.e., as
                    missing data).
 
 SNPs (14 states) : Alphabet: [A,C,G,T/U,K,M,R,S,W,Y,B,D,H,V,-] =
                    [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14].
 
                    Ambiguous characters (i.e., ? and N) are treated as if
                    they were alignment gaps (-) (i.e., as missing data).
 
 Distance metrics : F81 (Felsenstein 1981)
                    LogDet (Lake 1994; Lockhart et al. 1994)
 
 Manuscript       : ...
 
 Reference        : Zhang...
 
 ----------------------------------------------------------------------------*/

#include <cctype>
#include <cmath>
#include <limits>
#include <string>
#include <vector>
#include <random>
#include <iomanip>
#include <fstream>
#include <iostream>

#define SQR(a) ((a) * (a))

// The following are declared here because they are needed in different functions

std::ifstream infile;
std::string inName, outName1, outName2;
std::vector<std::string> taxon;                // 2D container for sequence names
std::vector<std::vector<int> > alignment;      // 2D container for sequence data
std::vector<unsigned long> position;
std::vector<int> sites;
const unsigned FOUR(4);        // for 4-state alphabet (DNA)
const unsigned TEN(10);        // for 10-state alphabet (SNP data)
const unsigned FOURTEEN(14);   // for 14-state alphabet (SNP data)
const unsigned TWENTY(20);     // for 20-state alphabet (amino acids)
const unsigned max_array(21);  // for all alphabets plus gaps
unsigned dataType(0);
unsigned rows_columns (0), repeat_max(1);
unsigned long alignment_length(0), sets_of_sites(0), number_of_sites(0);
unsigned long sequence_number(0);
unsigned long total;


// Function used to communicate contact details

void Information(){
    std::cout << "\nProgram     LogDet - version 1.0\n\n";
    std::cout << "Author      Lars S. Jermiin\n\n";
    std::cout << "Address     Research School of Biology\n";
    std::cout << "            Australian National University\n";
    std::cout << "            Acton, ACT 2601, Australia\n\n";
    std::cout << "            Univerity College Dublin\n";
    std::cout << "            School of Biology & Environmental Science\n";
    std::cout << "            Belfield, Dublin 4, Ireland\n\n";
    std::cout << "E-mail      lars.jermiin [at] anu.edu.au\n";
    std::cout << "            lars.jermiin [at] ucd.ie\n\n";
    std::cout << "Manual      https://github.com/lsjermiin/...\n";
    std::cout << "____________________________________________________\n";
}


// Function used to introduce the program

int Welcome(void){
    int   counter(0);
    char  buffer[2];
    bool  success(false);
    
    do{
        counter++;
        std::cout << "\nInformation [Y|N] ? ";
        std::cin >> buffer;
        buffer[0] = toupper(buffer[0]);
        if((buffer[0] == 'Y' || buffer[0] == 'N') && buffer[1] == '\0'){
            success = true;
            if(buffer[0] == 'Y') Information();
        }
        else std::cout << "\nWRONG ENTRY -- try again\n\n";
    }	while(success == false && counter < 3);
    if (success == false && counter == 3) {
        std::cerr << "OPTIONS EXHAUSTED -- program aborted!\n";
        exit(1);
    }
    return success;
}


// Function used to prepare input and output files

int Prepare_Files() {
    int   counter(0);
    bool  success(false);
    
    do{
        counter++;
        std::cout << "\n\nEnter name of input file ................. [???.fst] : ";
        std::cin >> inName;
        infile.open(inName.c_str());
        if (infile) {
            success = true;
        }
        else std::cout << "\nWARNING -- file not found...\n\n\n";
    }while(success == false && counter < 3);
    if (success == false && counter == 3) {
        std::cerr << "OPTIONS EXHAUSTED -- program aborted!\n";
        exit(1);
    } else {
        outName1 = "";
        for (std::string::size_type i = 0; i != inName.size() && inName[i] != '.'; ++i) {
            outName1 += inName[i];
        }
        infile.close();
    }
    return success;
}


// Function used to determine data type

unsigned Data_Type(void){
    int   counter(0);
    char  buffer[2];
    bool  success(false);
    
    std::cout << "\n\nThe data should be treated as an alignment of\n\n";
    std::cout << "    nucleotides ............................... type 1\n";
    std::cout << "    amino acids ............................... type 2\n";
    std::cout << "    SNPs (10 states) .......................... type 3\n";
    std::cout << "    SNPs (14 states) .......................... type 4\n\n";
    do{
        counter++;
        std::cout << "Type your preference ............................... ";
        std::cin >> buffer;
        if(buffer[0] > '0' && buffer[0] < '5' && buffer[1] == '\0')
            success = true;
        else{
            if(success == false && counter < 3) std::cout << "\nWRONG ENTRY -- try again\n\n";
            else{
                std::cerr << "\nOPTIONS EXHAUSTED -- analysis aborted!\n";
                std::cerr << "\nEXPLANATION: You chose an unavailable option" << std::endl;
                exit(1);
            }
        }
    }while(success == false && counter < 3);
    return buffer[0] - 48;
}


// Function used to determine the mode of analysis.
// Returns three numbers: 1. mode of analysis
//                        2. repeat_max
//                        3. set_of_sites

unsigned Mode_of_analysis(void){
    int   counter(0);
    unsigned long pos;
    char  buffer[2];
    bool  success(false);
    
    std::cout << "\n\nSpecify mode of analysis\n\n";
    std::cout << "    Use the sites as they appear .............. type 1\n";
    std::cout << "    Use ordered sets of sites ................. type 2\n";
    std::cout << "    Use randomly chosen subsets of sites ...... type 3\n";
    std::cout << "    Use non-parametric bootstrapping .......... type 4\n\n";
    do{
        counter++;
        std::cout << "Type your preference ............................... ";
        std::cin >> buffer;
        if(buffer[0] > '0' && buffer[0] < '5' && buffer[1] == '\0')
            success = true;
        else{
            if(success == false && counter < 3) std::cout << "\nWRONG ENTRY -- try again\n\n";
            else{
                std::cerr << "\nOPTIONS EXHAUSTED -- analysis aborted!\n";
                std::cerr << "\nEXPLANATION: You chose an unavailable option" << std::endl;
                exit(1);
            }
        }
    }while(success == false && counter < 3);
    if (success == true) {
        switch (buffer[0]) {
            case '1':
                repeat_max = 1;
                break;
            case '2':
                repeat_max = 0;
                std::cout << "\nEnter first position of each set .. [e.g., 1 99 ...] ";
                do {
                    std::cin >> pos;
                    position.push_back(pos - 1);
                    ++repeat_max;
                } while (std::cin && std::cin.peek() != '\n');
                position.push_back(alignment_length);
                for (std::vector<unsigned>::size_type i = 0; i != position.size() - 1; i++) {
                    if (position[i+1] < position[i] || position[i+1] > alignment_length) {
                        std::cerr << "\n\nWRONG ENTRY -- analysis aborted!\n\n";
                        std::cerr << "EXPLANATION: A number was smaller than the previous one, or" << std::endl;
                        std::cerr << "             the last number exceeded the alignment length.\n" << std::endl;
                        exit(1);
                    }
                }
                break;
            case '3':
                std::cout << "\n\nNumber of randomly chosen sites ............. [+100] ";
                std::cin >> number_of_sites;
                if (number_of_sites < 10 || number_of_sites > alignment_length) {
                    std::cout << "\n\nWRONG ENTRY -- analysis aborted!\n\n";
                    std::cerr << "EXPLANATION: The number was too small or larger than" << std::endl;
                    std::cerr << "             the length of the alignment.\n" << std::endl;
                    exit(1);
                }
                std::cout << "\nNumber of replicates .................... [1, 10000] ";
                std::cin >> repeat_max;
                if (repeat_max < 1 || repeat_max > 10000) {
                    std::cout << "\n\nWRONG ENTRY -- analysis aborted!\n\n";
                    std::cerr << "EXPLANATION: The number was too small or larger than 10000." << std::endl;
                    exit(1);
                }
                break;
            default:
                std::cout << "\n\nNumber of replicates .................... [1, 10000] ";
                std::cin >> repeat_max;
                if (repeat_max < 1 || repeat_max > 10000) {
                    std::cout << "\nWRONG ENTRY -- analysis aborted!\n\n";
                    std::cerr << "EXPLANATION: The number was too small or larger than 10000.\n" << std::endl;
                    exit(1);
                }
                break;
        }
    }
    return buffer[0] - 48;
}


// Function that reads input file and stores data in two 2D containers

int Read_Input(void){
    unsigned long counter(0);
    std::string temp(""), str(""); // temporary string used to store input
    std::vector<int> sequence;     // temporary string used to store input
    bool success(true);
    
    infile.open(inName.c_str());
    while (getline(infile, str)) {
        if (!str.empty()) {
            if (str[0] == '>') {
                // Store sequence name in vector
                if (str.size() > 0) {
                    temp = str;
                    str.erase();
                    for (std::string::size_type i = 1; i != temp.size(); ++i) {
                        str.push_back(temp[i]);
                    }
                    taxon.push_back(str);
                    ++sequence_number; // Needed to control output
                    str = "";
                }
                // Store sequence in vector
                if (sequence.size() > 0){
                    alignment.push_back(sequence);
                    if (alignment_length == 0) {
                        alignment_length = sequence.size(); // Needed to control output
                    }
                    sequence.clear();
                }
            } else {
                if (dataType == 1) {
                    for (std::string::size_type i = 0; i != str.size(); ++i) {
                        switch (toupper(str[i])) {
                            case 'A': sequence.push_back(0); break;
                            case 'C': sequence.push_back(1); break;
                            case 'G': sequence.push_back(2); break;
                            case 'T': sequence.push_back(3); break;
                            case 'U': sequence.push_back(3); break;
                            default : sequence.push_back(4); break; // In case of other characters
                        }
                    }
                } else {
                    if (dataType == 2) {
                        for (std::string::size_type i = 0; i != str.size(); ++i) {
                            switch (toupper(str[i])) {
                                case 'A': sequence.push_back(0); break;
                                case 'C': sequence.push_back(1); break;
                                case 'D': sequence.push_back(2); break;
                                case 'E': sequence.push_back(3); break;
                                case 'F': sequence.push_back(4); break;
                                case 'G': sequence.push_back(5); break;
                                case 'H': sequence.push_back(6); break;
                                case 'I': sequence.push_back(7); break;
                                case 'K': sequence.push_back(8); break;
                                case 'L': sequence.push_back(9); break;
                                case 'M': sequence.push_back(10);break;
                                case 'N': sequence.push_back(11);break;
                                case 'P': sequence.push_back(12);break;
                                case 'Q': sequence.push_back(13);break;
                                case 'R': sequence.push_back(14);break;
                                case 'S': sequence.push_back(15);break;
                                case 'T': sequence.push_back(16);break;
                                case 'V': sequence.push_back(17);break;
                                case 'W': sequence.push_back(18);break;
                                case 'Y': sequence.push_back(19);break;
                                default : sequence.push_back(20);break; // In case of other characters
                            }
                        }
                    } else {
                        if (dataType == 3) {
                            for (std::string::size_type i = 0; i != str.size(); ++i) {
                                switch (toupper(str[i])) {
                                    case 'A': sequence.push_back(0); break;
                                    case 'C': sequence.push_back(1); break;
                                    case 'G': sequence.push_back(2); break;
                                    case 'T': sequence.push_back(3); break;
                                    case 'U': sequence.push_back(3); break;
                                    case 'K': sequence.push_back(4); break;
                                    case 'M': sequence.push_back(5); break;
                                    case 'R': sequence.push_back(6); break;
                                    case 'Y': sequence.push_back(7); break;
                                    case 'S': sequence.push_back(8); break;
                                    case 'W': sequence.push_back(9); break;
                                    default : sequence.push_back(10);break; // In case of other characters
                                }
                            }
                        } else {
                            for (std::string::size_type i = 0; i != str.size(); ++i) {
                                switch (toupper(str[i])) {
                                    case 'A': sequence.push_back(0); break;
                                    case 'C': sequence.push_back(1); break;
                                    case 'G': sequence.push_back(2); break;
                                    case 'T': sequence.push_back(3); break;
                                    case 'U': sequence.push_back(3); break;
                                    case 'K': sequence.push_back(4); break;
                                    case 'M': sequence.push_back(5); break;
                                    case 'R': sequence.push_back(6); break;
                                    case 'Y': sequence.push_back(7); break;
                                    case 'S': sequence.push_back(8); break;
                                    case 'W': sequence.push_back(9); break;
                                    case 'B': sequence.push_back(10);break;
                                    case 'D': sequence.push_back(11);break;
                                    case 'H': sequence.push_back(12);break;
                                    case 'V': sequence.push_back(13);break;
                                    default : sequence.push_back(14);break; // In case of other characters
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    // Store last sequence in vector
    if (sequence.size() > 0) {
        alignment.push_back(sequence);
    } else {
        std::cerr << "Last sequence empty -- check input file!\n" << std::endl;
        success = false;
    }
    //Check whether the sequence names are unique
    for (std::vector<std::string>::const_iterator iter1 = taxon.begin(); iter1 != taxon.end(); ++iter1) {
        for (std::vector<std::string>::const_iterator iter2 = iter1 + 1; iter2 != taxon.end(); ++iter2) {
            if (*iter1 == *iter2) {
                std::cerr << "Sequences have the same name -- look for " << *iter1 << "\n" << std::endl;
                success = false;
            }
        }
    }
    // Check whether the sequences have the same length
    
    for (std::vector<std::vector<int> >::const_iterator iter = alignment.begin(); iter != alignment.end(); ++iter) {
        ++counter;
        sequence = *iter;
        if (sequence.size() != alignment_length) {
            std::cerr << "Sequences 1 and " << counter << " differ in length!\n" << std::endl;
            success = false;
            exit(0);
            
        }
    }
    return success;
}

// Function used to introduce the program

bool Bayesian_Multiplicative_Treatment(void){
    int   counter(0);
    char  buffer[2];
    bool  success(false);
    
    do{
        counter++;
        std::cout << "\nBM treatment of count zeros .................. [Y|N] ";
        std::cin >> buffer;
        buffer[0] = toupper(buffer[0]);
        if((buffer[0] == 'Y' || buffer[0] == 'N') && buffer[1] == '\0'){
            success = true;
        }
        else std::cout << "\nWRONG ENTRY -- try again\n\n";
    }	while(success == false && counter < 3);
    if (success == false && counter == 3) {
        std::cerr << "OPTIONS EXHAUSTED -- program aborted!\n";
        exit(1);
    }
    if (buffer[0] == 'Y') {
        success = true;
    } else {
        success = false;
    }
    return success;
}


unsigned Select_Model(void){
    int   counter(0);
    char  buffer[2];
    bool  success(false);
    
    do{
        counter++;
        std::cout << "\n\nMarkov models of sequence evolution\n\n";
        std::cout << "    F81 model ................................. type 1\n";
        std::cout << "    General Markov model ...................... type 2\n\n";
        std::cout << "Type your preference ............................... ";
        std::cin >> buffer;
        if((buffer[0] == '1' || buffer[0] == '2') && buffer[1] == '\0'){
            success = true;
        }
        else std::cout << "\nWRONG ENTRY -- try again\n\n";
    }	while(success == false && counter < 3);
    if (success == false && counter == 3) {
        std::cerr << "OPTIONS EXHAUSTED -- program aborted!\n";
        exit(1);
    }
    return buffer[0] - 48;
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

int main(){
    std::ofstream outfile1;
    std::vector<double> row_of_double;
    std::vector<std::vector<double> > matDist;
    bool bmt(false);
    unsigned long sum_dm(0), site(0);;
    unsigned zeros(0);
    unsigned mode(0);
    unsigned metric(0);
    unsigned seed;
    double logdet(0.0), F81(0.0), zeroAdj(0.0);
    unsigned long dm[max_array][max_array];   // 2D divergence matrix of integers
    double dmd[max_array][max_array];        // 2D divergence matrix of reals
    double row_sumd[max_array];
    double col_sumd[max_array];
    
    Welcome();
    Prepare_Files();
    dataType = Data_Type();
    switch (dataType) {
        case 1: rows_columns = FOUR; break;
        case 2: rows_columns = TWENTY; break;
        case 3: rows_columns = TEN; break;
        case 4: rows_columns = FOURTEEN; break;
        default: break;
    }
    Read_Input();
    metric = Select_Model();
    switch (metric) {
        case 1:
            outName1 = outName1 + "_F81.dis";
            break;
        case 2:
            outName1 = outName1 + "_LD.dis";
            break;
        default:
            outName1 = outName1 + "_NT.dis";
            break;
    }
    mode = Mode_of_analysis();
    bmt = Bayesian_Multiplicative_Treatment();
    for (std::vector<int>::size_type i = 0; i != alignment_length; ++i) {
        sites.push_back(0);
    }
    std::cout << std::endl;
    /*
     for (std::vector<int>::size_type i = 0 ; i != sites.size(); ++i) {
     std::cout << sites[i];
     }
     */
    // create and initialise 1D vector to hold to hold LogDet values
    for (std::vector<double>::size_type i = 0; i != sequence_number; i++) {
        row_of_double.push_back(0.0);
    }
    // create and initialise 2D vector to hold to hold LogDet values
    for (std::vector<std::vector<double> >::size_type i = 0; i != sequence_number; i++) {
        matDist.push_back(row_of_double);
    }
    // Get seed from system clock
    seed = (unsigned)std::chrono::system_clock::now().time_since_epoch().count();
    // Initialize randomisation
    std::mt19937 generator(seed);
    // Generate a uniform distribution between 0 and sites.size() - 1
    std::uniform_int_distribution<unsigned long> distribution(0,sites.size() - 1);
    outfile1.open(outName1.c_str());
    for (unsigned repeat = 0; repeat != repeat_max; repeat++) {
        for (std::string::size_type i = 0; i != alignment_length; ++i) {
            sites[i] = 0;
        }
        switch (mode) {
            case 1:
                total = sequence_number * (sequence_number - 1)/2;
                for (std::vector<int>::size_type i = 0; i != sites.size(); ++i) {
                    sites[i] = 1;
                }
                /*
                 for (std::vector<int>::size_type i = 0 ; i != sites.size(); ++i) {
                 std::cout << sites[i];
                 }
                 */
                std::cout << std::endl;
                break;
            case 2:
                total = sequence_number * (sequence_number - 1)/2;
                for (std::vector<int>::size_type i = 0; i != sites.size(); ++i) {
                    if (i >= position[repeat] && i < position[repeat+1]) {
                        sites[i] = 1;
                    } else {
                        sites[i] = 0;
                    }
                }
                /*
                 for (std::vector<int>::size_type i = 0 ; i != sites.size(); ++i) {
                 std::cout << sites[i];
                 }
                 */
                std::cout << std::endl;
                break;
            case 3:
                for (std::vector<int>::size_type i = 0; i != number_of_sites; i++) {
                    // Get a random number between 0 and sites.size() - 1
                    site = distribution(generator);
                    if (sites[site] == 0) {
                        sites[site] = 1;
                    } else {
                        i--;
                    }
                }
                /*
                 for (std::vector<int>::size_type i = 0 ; i != sites.size(); ++i) {
                 std::cout << sites[i];
                 }
                 */
                std::cout << std::endl;
                break;
            default:
                for (std::vector<int>::size_type i = 0; i != sites.size(); i++) {
                    site = distribution(generator);
                    ++sites[site];
                }
                /*
                 for (std::vector<int>::size_type i = 0 ; i != sites.size(); ++i) {
                 std::cout << sites[i];
                 }
                 */
                std::cout << std::endl;
                break;
        }
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
                        //                        std::cout << fixed << dmd[i][j] << " " ;
                    }
                    //                    std::cout << std::endl;
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
                if (bmt == true) {
                    // count the number of cells with zero
                    zeros = 0;
                    for (size_t i = 0; i != rows_columns; ++i) {
                        for (size_t j = 0; j != rows_columns; ++j) {
                            if (dm[i][j] == 0) {
                                ++zeros;
                            }
                        }
                    }
                    if (zeros > 0) {
                        // use Bayesian-multiplicative treatment of count zeros to prevent
                        // multiplication by zero. Here I use a square-root prior: see
                        // Table 1 in Statistical Modelling 15: 134–158 [2015]
                        zeroAdj = 1.0/(((double)sqrt(sum_dm) + 1.0) * SQR(rows_columns));
                        for (size_t i = 0; i != rows_columns; ++i) {
                            for (size_t j = 0; j != rows_columns; ++j) {
                                if (dm[i][j] == 0) {
                                    dmd[i][j] = zeroAdj;
                                } else {
                                    dmd[i][j] = (1.0 - zeroAdj * zeros) * dmd[i][j];
                                }
                            }
                        }
                    }
                }
                // select distance metric
                switch (metric) {
                    case 1:
                    {
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
                        break;
                    }
                    case 2:
                    {
                        // calculate LogDet distance
                        logdet = Det(rows_columns, dmd);
                        logdet = log(logdet);
                        for (size_t i = 0; i != rows_columns; ++i) {
                            logdet = logdet - 0.5 * (log(row_sumd[i]) + log(col_sumd[i]));
                        }
                        logdet = -1.0 * logdet/rows_columns;
                        matDist[iter1][iter2] = logdet;
                        matDist[iter2][iter1] = logdet;
                        break;
                    }
                    default:
                        break;
                }
                switch (mode) {
                    case 1:
                        std::cout << "\rNumber of comparisons left = " << --total;
                        break;
                    case 2:
                        std::cout << "\rNumber of comparisons left = " << --total;
                        break;
                    case 3:
                        std::cout << "\rNumber of repeats left = " << repeat_max - repeat;
                        break;
                    default:
                        std::cout << "\rNumber of repeats left = " << repeat_max - repeat;
                        break;
                }
                fflush(NULL);
            }
        }
        // Printing output
        outfile1 << sequence_number << std::endl;
        for (std::vector<std::vector<double> >::size_type i = 0 ; i != taxon.size(); i++) {
            outfile1 << std::left << std::setw(10) << taxon[i];
            for (std::vector<std::vector<double> >::size_type j = 0; j != taxon.size(); j++) {
                outfile1 << "\t" << std::fixed << matDist[i][j];
            }
            outfile1 << std::endl;
        }
        outfile1 << std::endl;
    }
    outfile1.close();
    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << "--------------------------------------------------------------------" << std::endl;
    std::cout << "   RESULTS FROM ANALYSIS COMPLETE" << std::endl;
    std::cout << std::endl;
    std::cout << "   Matrix with estimates of LogDet distances .. " << outName1 << std::endl;
    std::cout << "--------------------------------------------------------------------" << std::endl;
    std::cout << std::endl;
    return 0;
}
