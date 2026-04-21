#include <iostream>
#include <iomanip>
#include <sstream>     // for stringstream
#include <fstream>
#include <vector>
#include <map>
#include <list>

using namespace std;

/******************************************************************************/
/******************   TOOL Functions from "tools.cpp"   ***********************/
/******************************************************************************/
string int_to_bstring(__int128_t bool_nb, unsigned int n);
unsigned int Bitset_count(__int128_t bool_nb);

/******************************************************************************/
/***********************    READ an MCM from a FILE  **************************/
/******************************************************************************/

// file where involved variables in each ICCs are given by an integer corresponding to their index:
map<unsigned int, __int128_t> read_MCM_fromfile_index(string Input_MCM_file, unsigned int r)
{
    map<unsigned int, __int128_t> Partition;

    string line, line2;
    __int128_t Op = 1;
    Op <<= r - 1;
    vector<int> comm;

    ifstream myfile(Input_MCM_file.c_str());
    if (myfile.is_open())
    {
        while (getline(myfile, line))
        {
            stringstream ss(line);
            while (getline(ss, line2, '\t'))
            {
                comm.push_back(stoi(line2));
            }
            Partition[comm[1]] += Op;
            Op >>= 1;

            comm.clear();
        }
        myfile.close();
    }
    return Partition;
}

// file in which ICCs are encoded by binary strings:
map<unsigned int, __int128_t> read_MCM_fromfile_bin(string MCM_file, unsigned int r)
{
    __int128_t one128 = 1;

    cout << "Read the MCM file: " << MCM_file << endl;
    cout << "Number of variables to read: n = " << r << endl;

    string line, line2;     
    char c = '1';

    __int128_t Op_bit = 0, ICC_bin = 0;

// ***** Store MCM in MCM_map:  **********************************************
    map<unsigned int, __int128_t> MCM_map;

    unsigned int i = 0;  // ICC counter --> quite useless here, but due to format of returned MCM "map<unsigned int, __int128_t>"

    ifstream myfile (MCM_file.c_str());
    if (myfile.is_open())
    {
        while (getline (myfile, line))
        {
            if (line.size() >= r && (line[0] == '0' || line[0] == '1'))  // ignored all the lines that don't start by '0' or '1'
            {
                line2 = line.substr(0,r);          //take the r first characters of line

                //convert string line2 into a binary integer:
                Op_bit = one128 << (r - 1);
                ICC_bin = 0;
                for (auto &elem: line2)
                {
                    if (elem == c) { ICC_bin += Op_bit; }
                    Op_bit = Op_bit >> 1;
                }
                //cout << "ICC = " << int_to_bstring(ICC_bin, r) << endl;

                MCM_map[i] = ICC_bin; 
                i++;
            }  
        }
        myfile.close();
    }
    else
    {
        cout << endl << "                     ########## Unable to open file ##########" << endl << endl;
    }

    return MCM_map;
}

/******************************************************************************/
/************************    PRINT MCM in TERMINAL  ***************************/
/******************************************************************************/
// r = number of variables in the dataset
void Print_MCM_Partition(map<unsigned int, __int128_t> partition, unsigned int r)
{
    map<unsigned int, __int128_t>::iterator it;
    int i = 1;

    for (it = partition.begin(); it != partition.end(); it++)
    {
        cout << "ICC" << setw(3) << setfill(' ') << right << i;
        cout << ":   " << int_to_bstring((*it).second, r) << "\t ICC size: " << Bitset_count((*it).second) << endl;
        i++;
    }
    cout << endl;
}


/******************************************************************************/
/**************************    PRINT MCM in FILE    ***************************/
/******************************************************************************/

void PrintFile_MCM_Info(list<__int128_t> Basis, map<unsigned int, __int128_t> MCM_Partition, unsigned int r, string filename) //, string filename = "Result")
{
    //***** PRINT BASIS: 
    fstream file_MCM_info;

    if (Basis.size() == 0) // MCM is in the original basis
    {
        cout << "--->> Print best Greedy MCM in the file: " << filename + "_MCM_inB0.dat" << endl;
        file_MCM_info.open(filename + "_MCM_inB0.dat", ios::out);
        file_MCM_info << "## The MCM Partition is defined in the original basis by the following partition:" << endl;
    }
    else
    {
        cout << "--->> Print best Greedy MCM in the file: " << filename + "_MCM_inBnew.dat" << endl;
        file_MCM_info.open(filename + "_MCM_inBnew.dat", ios::out);
        file_MCM_info << "## The MCM Partition is defined in the new basis by the following partition:" << endl;
    }

    //***** PRINT MCM: 
    int i = 1;
    file_MCM_info << " " << endl;

    for (map<unsigned int, __int128_t>::iterator it = MCM_Partition.begin(); it != MCM_Partition.end(); it++)
    {    
        __int128_t Part = (*it).second;

        file_MCM_info << int_to_bstring(Part, r) << "\t ICC " << i << endl; 
        i++;
    }

    if (Basis.size() != 0) // MCM is in the original basis
    {
        file_MCM_info << " " << endl;
        file_MCM_info << "## The ICC above are written in the new basis defined by the following operators," << endl;
        file_MCM_info << "## with the convention:" << endl;
        file_MCM_info << "## \t   -- sig_1 = RIGHTmost bit in the MCM above" << endl;
        file_MCM_info << "## \t   -- sig_n = LEFTmost bit  in the MCM above" << endl;
        file_MCM_info << " " << endl;

        i = 1;
        for (list<__int128_t>::iterator it = Basis.begin(); it != Basis.end(); it++)
        {
            file_MCM_info << "##\t sig_" << i << " = " << int_to_bstring((*it), r) << endl; i++;
        } 
        file_MCM_info << " " << endl;
    }

    file_MCM_info.close();
}

/******************************************************************************/
/****************    PRINT MCM in FILE DURING FULL MERGING   ******************/
/******************************************************************************/
void PrintFile_MCM_FullMerge(map<unsigned int, __int128_t> MCM_Partition, unsigned int r, fstream &file) 
{
    //cout << "--->> Print best Greedy MCM in the file: " << filename + "_MCM_inB0.dat" << endl;

    //***** PRINT MCM: 
    int i = 1;
    file << " " << endl;
    file << "## A = " << MCM_Partition.size() << endl;

    for (map<unsigned int, __int128_t>::iterator it = MCM_Partition.begin(); it != MCM_Partition.end(); it++)
    {    
        __int128_t Part = (*it).second;

        file << int_to_bstring(Part, r) << "\t ICC " << i << endl; 
        i++;
    }
}

