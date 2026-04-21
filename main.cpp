// To compile: make
// To run MCM Greedy Search: 
//         or: make run
//         or: ./GreedySearch.out data_filename n
// With options:
//      ./GreedySearch.out data_filename n [-b basis_filename] [--full] [--NoCheckPoint] [--proba]
//
// To sample from an MCM:
//         or: make sample
//         or: ./GreedySearch.out data_filename n --sample MCM_filename
// With options:
//      ./GreedySearch.out data_filename n --sample MCM_filename [-b basis_filename] [--proba]


#define _USE_MATH_DEFINES 
#include <iostream>
#include <fstream>
#include <sstream>
#include <list>
#include <map>
#include <vector>
#include <cmath>       /* tgamma */
#include <random>

#include <ctime> // for chrono
#include <ratio> // for chrono
#include <chrono> // for chrono

using namespace std;
using namespace std::chrono;

/*****************************************************************************************/
/*************************   CONSTANT VARIABLES  *****************************************/
/*****************************************************************************************/
#include "includes/default_data.h"

/*****************************************************************************************/
/****************   GREEDY SEARCH:   Useful functions and routines    ********************/
/*****************************************************************************************/
// **** Find the best MCM, Greedy Search:
map<unsigned int, __int128_t> MCM_GreedySearch(vector<pair<__int128_t, unsigned int>> Kset, unsigned int N, unsigned int r, bool print_info = true, bool Greedy_Full_merging = false);
map<unsigned int, __int128_t> MCM_GreedySearch_AND_printInfo(vector<pair<__int128_t, unsigned int>> Kset, unsigned int N, unsigned int r, bool print_info = true, bool Greedy_Full_merging = false);

// **** Find the best MCM, Greedy Search starting from the model MCM_0:
map<unsigned int, __int128_t> MCM_GreedySearch_MCM0(vector<pair<__int128_t, unsigned int>> Kset, unsigned int N, unsigned int r, map<unsigned int, __int128_t> MCM_0, bool print_info = true, bool Greedy_Full_merging = false);

// *** Greedy Search on Reduced dataset:
map<unsigned int, __int128_t> MCM_ReducedGreedySearch_AND_PrintInfo(vector<pair<__int128_t, unsigned int>> Kset, unsigned int K, unsigned int N, unsigned int r, bool print_it = false);

// *** Compare two MCMs:
void compare_two_MCMs_AND_printInfo(vector<pair<__int128_t, unsigned int>> Kset, unsigned int N, unsigned int r, map<unsigned int, __int128_t> fp1, map<unsigned int, __int128_t> fp2);


/*****************************************************************************************/
/*****************************   INFORMATION on a BASIS   ********************************/
/*****************************************************************************************/
bool Is_IndepModel(list<__int128_t> Basis_li, unsigned int n);
list<__int128_t> Invert_Basis(list<__int128_t> Basis_li, unsigned int n);


/*****************************************************************************************/
/*****************************   IMPORT an MCM from a FILE   *****************************/
/*****************************************************************************************/
// *** Read MCM from a file:
map<unsigned int, __int128_t> read_MCM_fromfile_bin(string MCM_file, unsigned int r);
map<unsigned int, __int128_t> read_MCM_fromfile_AND_printInfo(vector<pair<__int128_t, unsigned int>> Kset, unsigned int N, string Input_MCM_file, unsigned int r);


//map<__int128_t, unsigned int> build_Kset_ICC(vector<pair<__int128_t, unsigned int>> Kset, __int128_t Ai);
//vector<pair<__int128_t, unsigned int>> build_Kset_ICC_vect(vector<pair<__int128_t, unsigned int>> Kset, __int128_t Ai);


/*****************************************************************************************/
/**********************************   SAMPLING       *************************************/
/*****************************************************************************************/
void PrintTerm_samples(vector<pair<__int128_t, unsigned int>> Kset, map<unsigned int, __int128_t> MCM, list<__int128_t> Basis_li_invert, unsigned int n, unsigned int Nsample);
void PrintFile_samples(vector<pair<__int128_t, unsigned int>> Kset, map<unsigned int, __int128_t> MCM, list<__int128_t> Basis_li_invert, unsigned int n, unsigned int Nsample, string sample_filename);


/*****************************************************************************************/
/*************   PRINT FILE:    Model VS Data Probability Distribution   *****************/
/*****************************************************************************************/
void PrintFile_StateProbabilites_CurrentBasis(std::vector<std::pair<__int128_t, unsigned int>> Kset, std::map<unsigned int, __int128_t> MCM_Partition, unsigned int N, unsigned int r, std::string filename = "Result");
void PrintFile_StateProbabilites_OriginalBasis_NewBasis(std::vector<std::pair<__int128_t, unsigned int>> Nset, std::list<__int128_t> Basis, std::map<unsigned int, __int128_t> MCM_Partition, unsigned int N, unsigned int r, std::string filename = "Result");


/******************************************************************************/
/***************************   ADD OUTPUT FOLDER    ***************************/
/******************************************************************************/

//// ** location of the output folder:
// string OutputFile_Add_Location(string filename)
string Add_Output_Location(string filename)
{
    return (OUTPUT_directory + filename);
}

/****************************************************************************************************************************************************************************/
/****************************************************************************************************************************************************************************/
/**************************************************************************   """ TUTORIAL  """    **************************************************************************/
/****************************************************************************************************************************************************************************/
/****************************************************************************************************************************************************************************/

void tutorial(vector<pair<__int128_t, unsigned int>> Nset, unsigned int N,  unsigned int n, bool Greedy_Full_merging = false)
{
    cout << endl << "*******************************************************************************************";  
    cout << endl << "******************************  CHOICE OF THE BASIS:  *************************************";
    cout << endl << "*******************************************************************************************" << endl;

// original basis of the data: this is the most natural choice a priori:
    list<__int128_t> Basis_li = Original_Basis(n);

//// *** The basis can also be read from a file: Ex. the following files contain the best basis for the SCOTUS dataset:
//   list<__int128_t> Basis_li = Read_BasisOp_IntegerRepresentation(input_directory + basis_IntegerRepresentation_filename);
//   list<__int128_t> Basis_li = Read_BasisOp_BinaryRepresentation(n, input_directory + basis_BinaryRepresentation_filename);

    PrintTerm_Basis(Basis_li, n);

    cout << endl << "*******************************************************************************************";
    cout << endl << "**********************  TRANSFORM the DATA in the CHOSEN BASIS   **************************";
    cout << endl << "**********************************   Build Kset:   ****************************************";
    cout << endl << "*******************************************************************************************" << endl;

//// *** Transform the data in the specified in Basis_SCModel[];

    vector<pair<__int128_t, unsigned int>> Kset = build_Kset(Nset, Basis_li);
    cout << "\t Kset.size() = " << Kset.size() << endl;


    cout << endl << "*******************************************************************************************";
    cout << endl << "***********************  HIERARCHICAL GREEDY MERGING: BY STEPS:  **************************";
    cout << endl << "*******************************************************************************************" << endl;

    bool print_checkpoint = true;  

//// *** Finds the best MCM:
    map<unsigned int, __int128_t> mcm1 = MCM_GreedySearch(Kset, N, n, print_checkpoint, Greedy_Full_merging);

//// *** Print Log-Evidence:  
    double LogE_mcm1 = LogE_MCM(Kset, mcm1, N, n);
    cout << "Log-Evidence(MCM) = " << LogE_mcm1 << "\t = " << LogE_mcm1/((double) N)/log(2.) << " bits per datapoint \t" << endl;

//// *** Print max-Log-Likelihood:  
    double LogL_mcm1 = LogL_MCM(Kset, mcm1, N, n);
    cout << "Max Log-Likelihood(MCM) = " << LogL_mcm1 << "\t = " << LogL_mcm1/((double) N)/log(2.) << " bits per datapoint \t" << endl;

    // ***** Print info about the model -- Print Basis and MCM:  **************/
    PrintFile_MCM_Info(Basis_li, mcm1, n, "Best_MCM");

    cout << endl << "*******************************************************************************************";
    cout << endl << "******************************  READ an MCM from a FILE   *********************************";
    cout << endl << "*******************************************************************************************" << endl;

    cout << "#########  EX. READ a CHOSEN MCM:  #########" << endl;

    // the file MCM_ex = "INPUT/SCOTUS_Communities_inBestBasis.dat" contains the best MCM in the best basis:
    map<unsigned int, __int128_t> mcm2 = read_MCM_fromfile_bin(MCM_ex, n);
    Print_MCM_Partition(mcm2, n);

    cout << endl << "*******************************************************************************************";
    cout << endl << "*******************************  COMPARING TWO MCMs   *************************************";
    cout << endl << "*******************************************************************************************" << endl;

    compare_two_MCMs_AND_printInfo(Kset, N, n, mcm1, mcm2);

    cout << endl << "*******************************************************************************************";
    cout << endl << "***************************  Decomposition of Log-E   *************************************";
    cout << endl << "*******************************   over each ICC   *****************************************";
    cout << endl << "*******************************************************************************************" << endl;

    double LogE_final = LogE_MCM_infoICC(Kset, mcm1, N, n);
    //cout << "Log-Evidence(MCM) = " << LogE_final << "\t = " << LogE_final/((double) N)/log(2.) << " bits per datapoint \t" << endl;

    cout << endl << "*******************************************************************************************";
    cout << endl << "*************************  Decomposition of Max-Log-L   ***********************************";
    cout << endl << "*******************************   over each ICC   *****************************************";
    cout << endl << "*******************************************************************************************" << endl;

    double LogL_final = LogL_MCM_infoICC(Kset, mcm1, N, n);
    //cout << "Max-Log-Likelihood(MCM) = " << LogL_final << "\t = " << LogL_final/((double) N)/log(2.) << " bits per datapoint \t" << endl;

    cout << endl << "*******************************************************************************************";
    cout << endl << "***************************  Working with a Reduced Dataset   *****************************";
    cout << endl << "**********   Remove from Kset all the states that occur less than K times:   **************";
    cout << endl << "*******************************************************************************************" << endl;

    // All the states that occur less than K times will be removed from the dataset:
    unsigned int K=2;
    map<unsigned int, __int128_t> fp_reduced = MCM_ReducedGreedySearch_AND_PrintInfo(Kset, K, N, n);

    // ***** Print info about the model -- Print Basis and MCM:  **************/
    PrintFile_MCM_Info(Basis_li, fp_reduced, n, "Best_MCM_reduced");

    cout << endl << "*******************************************************************************************";
    cout << endl << "**********************  Print information about the found MCM:  ***************************";
    cout << endl << "*******************************************************************************************" << endl;

    // Prints 1) information about the MCM; 2) the state probabilities P(s) of observed states (in the Data VS MCM); 3) the probability P(k) of observing a state with k values "+1" (in the Data VS MCM) 
    PrintFile_StateProbabilites_OriginalBasis_NewBasis(Nset, Basis_li, mcm1, N, n, Add_Output_Location("Result"));

    // Print the state probabilities P(s) of observed states (in the Data VS MCM) using the data transformed in the bew basis:
    PrintFile_StateProbabilites_CurrentBasis(Kset, mcm1, N, n, Add_Output_Location("Result"));
}


/****************************************************************************************************************************************************************************/
/****************************************************************************************************************************************************************************/
/**************************************************************************     MAIN FUNCTION      **************************************************************************/
/****************************************************************************************************************************************************************************/
/****************************************************************************************************************************************************************************/
int main(int argc, char *argv[])
{

    cout << endl << "*******************************************************************************************";
    cout << endl << "************************************  CONFIGURATION:  *************************************";
    cout << endl << "*******************************************************************************************" << endl;

//// *** READ USER CONFIGURATION:
    cout << "--->> Input files are in the input directory: \"" << input_directory << "\"" << endl; 
    RunOptions options;

    if ( !(Read_argument(argc, argv, &datafilename, &n, &basis_filename, &options)) )  {   return 0;   }  // error flag --> quit


//// *** OUTPUT SETTINGS:
    cout << endl;
    cout << "--->> Create the \"OUTPUT\" Folder (if needed) ";
    system(("mkdir -p " + OUTPUT_directory).c_str());
    cout << endl;

    string prefix_datafilename = filename_remove_extension(datafilename); // For output specific to the Dataset

    cout << endl << "*******************************************************************************************";
    cout << endl << "***********************************  READ THE DATA:  **************************************";
    cout << endl << "*******************************************************************************************" << endl;

    //cout << "Read the dataset: " << datafilename << endl;
    //cout << "Number of variables to read: n = " << n << endl;

    unsigned int N = 0; // will contain the number of datapoints in the dataset
    vector<pair<__int128_t, unsigned int>> Nset = read_datafile(&N, input_directory + datafilename, n);

    if (N == 0)     // Terminate program if the file can't be found or read, or if it is empty:
        { 
        cout << "--->> Datafile cannot be read or is empty: terminate the program." << endl << endl;
        return 0; 
        }

    cout << endl << "###### Datafile has been read successfully:" << endl;
    cout << "\t Number of datapoints: N = " << N << endl;
    cout << "\t Number of different observed states = " << Nset.size() << endl;

    vector<pair<__int128_t, unsigned int>> Kset;

    list<__int128_t> Basis_li, Basis_li_invert;

    // ***** will comtain chosen MCM later:
    map<unsigned int, __int128_t> mcm3;

    if (options.change_basis)
    {
        cout << endl << "*******************************************************************************************";  
        cout << endl << "******************************  CHOICE OF THE BASIS:  *************************************";
        cout << endl << "*******************************************************************************************" << endl;

    //// *** The basis can also be read from a file: Ex. the following files contain the best basis for the SCOTUS dataset:

        Basis_li = Read_BasisOp_BinaryRepresentation(n, input_directory + basis_filename);

        if (Basis_li.size() == 0)     // Terminate program if the file can't be found or read, or if it is empty:
            { 
            cout << "--->> Basis file cannot be read or is empty: terminate the program." << endl << endl;
            return 0; 
            }

        cout << endl << "###### Basis file has been read successfully:" << endl;
        PrintTerm_Basis(Basis_li, n);

        if (!Is_IndepModel(Basis_li, n))
            { 
            cout << "--->> The Basis provided in " << input_directory + basis_filename << " is not an independent set: Terminate the program." << endl << endl;
            return 0; 
            }

        if(options.sampling)
        {
            cout << endl << "*******************************************************************************************";
            cout << endl << "*************************  PRINT INVERSE BASIS TRANSFORMATION:  ***************************";
            cout << endl << "*******************************************************************************************" << endl;

            //list<__int128_t> 
            Basis_li_invert = Invert_Basis(Basis_li, n);

            if (Basis_li_invert.size() == 0)  // Inverse not available
                {   cout << "Inverse transformation is not available." << endl;     }
            else
            {
                cout << "**** Inverse Basis Transformation: ****" << endl << endl;
                PrintTerm_Basis_invert(Basis_li_invert, n); 
            }
        }

        cout << endl << "*******************************************************************************************";
        cout << endl << "**********************  TRANSFORM the DATA in the CHOSEN BASIS   **************************";
        cout << endl << "**********************************   Build Kset:   ****************************************";
        cout << endl << "*******************************************************************************************" << endl;

        Kset = build_Kset(Nset, Basis_li);

        if (Kset.size() == 0)     // Terminate program if the file can't be found or read, or if it is empty:
        { 
            cout << "--->> Issue with conversion of the datafile in the new basis:" << endl;
            cout << " \t - Check that provided basis is actually a basis." << endl;
            cout << " \t Terminate the program." << endl << endl;
            return 0; 
        }

        cout << "###### Datafile has been converted successfully:" << endl;
        cout << "\t Number of datapoints: N = " << N << endl;
        cout << "\t Number of different observed states in the new basis = " << Kset.size() << endl;

    }
    else
    {
        Kset.swap(Nset);
    }

    if (!(options.sampling)) // DEFAULT MODE = GREEDY SEARCH:
    {
        if (options.change_basis)
        {
            cout << endl << "*******************************************************************************************";
            cout << endl << "*****************************  HIERARCHICAL GREEDY MERGING:  ******************************";
            cout << endl << "**********************************  in the NEW BASIS  *************************************";
            cout << endl << "*******************************************************************************************" << endl;
        }
        else 
        {
        //   Kset.swap(Nset);
    
            cout << endl << "*******************************************************************************************";
            cout << endl << "*****************************  HIERARCHICAL GREEDY MERGING:  ******************************";
            cout << endl << "********************************  in the ORIGINAL BASIS  **********************************";
            cout << endl << "*******************************************************************************************" << endl;
        }

    //// *** Finds the best MCM and print information about it in the terminal:
        map<unsigned int, __int128_t> mcm1 = MCM_GreedySearch_AND_printInfo(Kset, N, n, options.print_checkpoint, options.greedy_full_merging);

    //***** Print info about the model in file: Basis and MCM:  **************/
        PrintFile_MCM_Info(Basis_li, mcm1, n, Add_Output_Location(prefix_datafilename));

        cout << endl << "*******************************************************************************************";
        cout << endl << "***************************  Decomposition of Log-E   *************************************";
        cout << endl << "*******************************   over each ICC   *****************************************";
        cout << endl << "*******************************************************************************************" << endl;
        LogE_MCM_infoICC(Kset, mcm1, N, n);

        cout << endl << "*******************************************************************************************";
        cout << endl << "*************************  Decomposition of Max-Log-L   ***********************************";
        cout << endl << "*******************************   over each ICC   *****************************************";
        cout << endl << "*******************************************************************************************" << endl;

        LogL_MCM_infoICC(Kset, mcm1, N, n);

        mcm3 = mcm1;
    }
    else // MODE = SAMPLING:
    {
        cout << endl << "*******************************************************************************************";
        cout << endl << "******************************  READ MCM FROM FILE:  **************************************";
        cout << endl << "*******************************************************************************************" << endl;

        map<unsigned int, __int128_t> mcm2 = read_MCM_fromfile_bin(OUTPUT_directory + (options.MCM_file), n);

        if (mcm2.size() == 0)     // Terminate program if the file can't be found or read, or if it is empty:
            { 
            cout << "--->> MCM file cannot be read or is empty: terminate the program." << endl << endl;
            return 0; 
            }

        cout << endl << "###### MCM file has been read successfully:" << endl;
        Print_MCM_Partition(mcm2, n);


        cout << endl << "*******************************************************************************************";
        cout << endl << "*****************************  SAMPLE DATA FROM MCM:  *************************************";
        cout << endl << "*******************************************************************************************" << endl;

//      PrintTerm_samples(Kset, mcm2, Basis_li_invert, n, options.Nsample);
        PrintFile_samples(Kset, mcm2, Basis_li_invert, n, options.Nsample, OUTPUT_directory + prefix_datafilename);

    
// **** HISTOGRAM OF SAMPLED DATA:    
/*        int Nsample2 = 100000;

        vector<pair<__int128_t, unsigned int>> Kset_sample = histo_sample_bestMCM(Kset, mcm3, Basis_li_invert, n, Nsample2);
        Print_File_Kset(Kset_sample, Nsample2, n, OUTPUT_directory + "Kset_sample.dat");

        vector<pair<__int128_t, unsigned int>> Nset_sample = build_Kset(Kset_sample, Basis_li_invert);
        Print_File_Kset(Nset_sample, Nsample2, n, OUTPUT_directory + "Nset_sample.dat");
*/
        mcm3 = mcm2;
    }

    if (options.proba) // OPTION: print the empirical and model probability distributions
    {
        cout << endl << "*******************************************************************************************";
        cout << endl << "**********************  Print information about the found MCM:  ***************************";
        cout << endl << "*******************************************************************************************" << endl;


        if (Basis_li.size() != 0)  //  Kset != Nset:
        {
        // Prints:
        //    1) the state probabilities P(s) of observed states (in the Data VS MCM); 
        //    2) the probability P(k) of observing a state with k values "+1" (in the Data VS MCM) 
        PrintFile_StateProbabilites_OriginalBasis_NewBasis(Nset, Basis_li, mcm3, N, n, Add_Output_Location(prefix_datafilename)); // StateProba_Original_AND_NewBasis

        }
        else   // in this case Kset == Nset:
        {   
        //Print the state probabilities P(s) of observed states (in the Data VS MCM) using the data transformed in the bew basis:
        PrintFile_StateProbabilites_CurrentBasis(Kset, mcm3, N, n, Add_Output_Location(prefix_datafilename));       // StateProba_inBasis_of_Kset_only // StateProba_inKsetBasis() 
        // if Kset = Nset, then written in original basis variables
        // if Kset = transformed data, then only written in the new basis variables
        }
    }
    cout << endl;

/*
    cout << endl << "*******************************************************************************************************************";
    cout << endl << "*******************************************************************************************************************";
    cout << endl << "************************************************  TUTORIAL:  ******************************************************";
    cout << endl << "*******************************************************************************************************************";
    cout << endl << "*******************************************************************************************************************" << endl;

    bool greedy_full_merging = false;

    tutorial(Nset, N, n, greedy_full_merging);
*/

    return 0;
}
