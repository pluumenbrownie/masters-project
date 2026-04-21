#include <iostream>
#include <sstream>

using namespace std;


/******************************************************************************/
/**************************  ERROR MESSAGE  ***********************************/
/******************************************************************************/

void ERROR_message()
{
    cout << "       General call:   >> ./MCM_Greedy.out [datafilename] [n] (-b [basisfilename]) (--full) (--NoCheckPoint)" << endl;
    cout << "                       Commands in \'(...)\' are optional." << endl << endl; 
    cout << "       To get help:    >> ./MCM_Greedy.out -h" << endl << endl;
}

/******************************************************************************/
/**************************  HELP MESSAGE  ************************************/
/******************************************************************************/

void HELP_message()
{
    cout << endl << "*******************************************************************************************";
    cout << endl << "***************************************  HELP  ********************************************";
    cout << endl << "*******************************************************************************************" << endl;

    cout << endl << "*******************************************************************************************";
    cout << endl << "******************************  HOW TO RUN THE PROGRAM:  **********************************";
    cout << endl << "*******************************************************************************************" << endl;
    
    cout << "To perform the MCM Greedy search:" << endl;
    cout << "   1. Place the datafile(s) in the 'INPUT' folder (datafiles must be in the correct" << endl; 
    cout << "\t binary format -- see README file)" << endl;
    cout << "   2. Type in your terminal one of the commands below. Remember to:" << endl;
    cout << "\t replace [datafilename] by the name of your datafile" << endl;
    cout << "\t replace [n] by the number of variables in your data (must be an integer)" << endl;
    cout << "\t replace [basisfilename] by the name of the file containing your basis (in binary format)" << endl;
    cout << "   3. Options 1, 2, and 3 below can be combined in any order." << endl;
    cout << "\t Option 2 and 3 are incompatible:" << endl;
    cout << "\t Option 2 will always print the checkpoints & take over Option 3." << endl;

    cout << endl << "*******************************************************************************************";
    cout << endl << "***********************************  Default EXAMPLE  *************************************";
    cout << endl << "*******************************************************************************************" << endl << endl;

    cout << "This will run the example:" << endl << endl;
    cout << "\tRun: "; 
    cout << "   >> ./MCM_Greedy.out" << endl << endl;
    cout << "\tOR: "; 
    cout << "    >> make run" << endl;    
 
    cout << endl << "*******************************************************************************************";
    cout << endl << "********************************  Run the Greedy Search  **********************************";
    cout << endl << "**************************  in the ORIGINAL BASIS of the data  ****************************";
    cout << endl << "*******************************************************************************************" << endl << endl;

    cout << "This will run the greedy search for the selected binary datafile [datafilename]" << endl;
    cout << "\tover the [n] first variables:" << endl << endl;
    cout << "\tRun: ";
    cout << "   >> ./MCM_Greedy.out [datafilename] [n]" << endl << endl;

    cout << "Important: - The datafile must have at least [n] binary variables" << endl;

    cout << endl << "***************************************  OPTION 1  ****************************************";
    cout << endl << "********************************  Run the Greedy Search  **********************************";
    cout << endl << "********************************  in a chosen NEW BASIS  **********************************";
    cout << endl << "*******************************************************************************************" << endl << endl;

    cout << "OPTION:  -b [basisfilename]" << endl << endl;

    cout << "This will run the greedy search for the selected binary datafile [datafilename]" << endl;
    cout << "\tin the new basis [basisfilename]" << endl;
    cout << "\tover the [n] first basis variables:" << endl << endl;
    cout << "\tExample: ";
    cout << "   >> ./MCM_Greedy.out [datafilename] [n] -b [basisfilename]" << endl << endl;

    cout << "Important: - The datafile and the basis must have at least [n] binary variables" << endl;
    cout << "           - Make sure that your basis is a proper basis, as the program won't check this." << endl;

    cout << endl << "***************************************  OPTION 2  ****************************************";
    cout << endl << "********************************  Run the Greedy Search  **********************************";
    cout << endl << "******************************  until reaching FULL merging  ******************************";
    cout << endl << "*******************************************************************************************" << endl << endl;

    cout << "OPTION:  --full" << endl << endl;

    cout << "By default, the greedy algorithm will stop as soon as any additional merging starts" << endl;
    cout << "   decreasing the LogE." << endl;
    cout << "This option allows pursing merging until everything is merged, and will print the evolution" << endl;
    cout << "   of the LogE along the merging path." << endl;
    cout << "The best MCM along the greedy path is saved and returned at the end." << endl << endl;

    cout << "\tExample: ";
    cout << "   >> ./MCM_Greedy.out [datafilename] [n] --full" << endl << endl;    

    cout << endl << "***************************************  OPTION 3  ****************************************";
    cout << endl << "**********************************  Remove checkpoints  ***********************************";
    cout << endl << "*******************************************************************************************" << endl << endl;

    cout << "OPTION:  --NoCheckPoint" << endl << endl;

    cout << "By default, the program prints out the intermediate values of LogE along the Greedy path." << endl;
    cout << "   Using this option will stop printing this information." << endl;
    cout << "   This may save time, when printing is not needed, for datasets with a large number of" << endl;
    cout << "   variables and a long greedy path." << endl << endl;

    cout << "\tExample: ";
    cout << "   >> ./MCM_Greedy.out [datafilename] [n] --NoCheckPoint" << endl << endl; 

    cout << endl << "*******************************************************************************************" << endl;
    cout << endl;
}


/******************************************************************************/
/****************** STRUCTURE INFO Read Arguments *****************************/
/******************************************************************************/
// IMPORTANT: same structure must be defined in the file "library.hpp"

struct RunOptions
{
    bool change_basis = false;  // by default: Search in the original basis 
    bool print_checkpoint = true;   // by default: print all the checkpoints
    
    // by default: stop Greedy merging when LogE starts decreasing 
    bool greedy_full_merging = false; // if TRUE: keep on merging until everything is merged; save best MCM along the greedy path

    // MODE GREEDY:
    // by default in greedy merging mode:

    // MODE SAMPLING:
    // if sampling = true, then in sampling mode
    bool sampling = false;
    unsigned int Nsample = 1000; // default value
    std::string MCM_file = "";

    // MODE COMPUTE PROBABILITIES:
    bool proba = false;
};

/******************************************************************************/
/************************** Read Arguments ************************************/
/******************************************************************************/

int Read_argument(int argc, char *argv[], string *datafilename, unsigned int *n, string *basis_filename, RunOptions *options)
{
    int i = 3;

    if (argc == 2)
    {
        string help = argv[1];
        if (help == "-h")
            { HELP_message(); }
        else
            { cout << endl << "ERROR: The number of arguments is not correct." << endl << endl; ERROR_message(); } 
        return 0;
    }

    else if (argc >= 3)
    {
        (*datafilename) = argv[1];
        (*n) = stoul(argv[2]);

        i=3;
        while(i<argc)
        {   
            if(((string) argv[i]) == "-b")
            {
                if((i+1) < argc)
                {
                    (*basis_filename) = argv[i+1];
                    (*options).change_basis = true;
                    i+=2;
                }
                else
                { 
                    cout << endl << "ERROR: missing basis filename." << endl << endl; 
                    ERROR_message(); 
                    return 0;
                }
            } 
            else if(((string) argv[i]) == "--full")
            {
                (*options).greedy_full_merging = true; 
                i++;
            }
            else if(((string) argv[i]) == "--NoCheckPoint")
            {
                (*options).print_checkpoint = false;
                i++;
            }
            else if(((string) argv[i]) == "--sample")
            {
                if((i+1) < argc)
                {
                    (*options).MCM_file = argv[i+1];
                    (*options).sampling = true;
                    i+=2;
                }
                else
                { 
                    cout << endl << "ERROR: missing MCM filename after option \'--sample\'." << endl << endl; 
                    ERROR_message(); 
                    return 0;
                }
            }
            else if(((string) argv[i]) == "-N")
            {
                if((i+1) < argc)
                {
                    (*options).Nsample = stoul(argv[i+1]);
                    i+=2;
                }
                else
                { 
                    cout << endl << "ERROR: missing number of sample after option \'-N\'." << endl << endl; 
                    ERROR_message(); 
                    return 0;
                }
            }
            else if(((string) argv[i]) == "--proba")
            {
                (*options).proba = true;
                i++;
            }
            else
            {
                cout << endl << "ERROR: issue with the arguments of the program." << endl << endl; 
                ERROR_message(); 
                return 0; 
            }
        }
    }

    if((*options).greedy_full_merging) { (*options).print_checkpoint = true; }

    // **** Print information for the user:
    cout << endl;
    cout << "Dataset to read: " << (*datafilename) << " with n = " << (*n) << " variables" << endl << endl;
    cout << "OPTIONS: ";
    if((*options).sampling)
    {
        cout << " ********* MODE: SAMPLING ***************** " << endl << endl;

        cout << "\t - MCM to use:                    \"" << (*options).MCM_file << "\"" << endl;  
        cout << "\t - Number of samples:             " << ((*options).Nsample) << endl;      
    }
    else if((*options).proba)
    {
        cout << " ********* MODE: PROBABILITIES ************* " << endl;
        cout << "\t - Compute and compare the empirical probabilities and the model probabilities" << endl << endl;      

        cout << "\t - for the MCM in:                \"" << (*options).MCM_file << "\"" << endl;      
    }
    else 
    {
        cout << " ********* MODE: MCM GREEDY SEARCH ********* " << endl << endl;

        cout << "\t - Print checkpoints:             " << (((*options).print_checkpoint)?"Yes":"No") << endl;
        cout << "\t - Run Greedy until FULLY merged: " << (((*options).greedy_full_merging)?"Yes":"No") << endl;
    }

    if((*options).change_basis)
    { 
        cout << "\t - Change Basis:                  Yes" << endl;
        cout << "\t - New Basis in:                  \"" << (*basis_filename) << "\"" << endl; 
    }
    else
    {   cout << "\t - Change Basis:                  No" << endl; }

    return 1; 
}

