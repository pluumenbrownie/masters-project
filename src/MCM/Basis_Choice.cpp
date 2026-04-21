#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <list>

using namespace std;

/********************************************************************/
/**************************    CONSTANTS    *************************/
/********************************************************************/
// INPUT BASIS FILES (optional):
// const string basis_IntegerRepresentation_filename = "INPUT/SCOTUS_n9_BestBasis_Integer.dat";        // (optional) Input basis file 
// const string basis_BinaryRepresentation_filename = "INPUT/SCOTUS_n9_BestBasis_Binary.dat";      // (optional) Input basis file

/******************************************************************************/
/***************************   Constant variables   ***************************/
/******************************************************************************/
const __int128_t one128 = 1;

/******************************************************************************/
/******************   TOOL Functions from "tools.cpp"   ***********************/
/******************************************************************************/
string int_to_bstring(__int128_t bool_nb, unsigned int n);

/******************************************************************************/
/*****************   Read Basis Operators from file  **************************/
/******************************************************************************/

/******  1)  Operators should be written in one single column          ********/
/******  2)  Operators can be written in two different versions:       ********/
/***   a) as a binary representation of the spin involved in the Operator;    */
/***   b) as the integer value of that binary representation.                 */
/******************************************************************************/
/****  Ex. for a system with n=4 spin:  ***************************************/
/****      -- the operator s1 would be encoded as 0001,                      **/
/****      which has the integer representation 1  -->   0001 = 1      ********/
/****      -- the operator s1 s2 would be encoded as 0011,             ********/
/****      which has the integer representation 3  -->   0011 = 3      ********/
/******************************************************************************/


/******************************************************************************/
/*** VERSION a) Operators are written as the binary          ******************/
/****           representation of the interactions           ******************/
/****           Works for up to 128 bits           ****************************/
/******************************************************************************/
list<__int128_t> Read_BasisOp_BinaryRepresentation(unsigned int r, string Basis_binary_filename) //string Basis_binary_filename = basis_BinaryRepresentation_filename)
{
  __int128_t Op_bit = 0, Op = 0;

  cout << "Read the Basis file: " << Basis_binary_filename << endl;
  cout << "Number of variables to read: n = " << r << endl;

  string line, line2;   
  char c = '1';  

// ***** Store Basis in Basis_li:  **********************************************
  list<__int128_t> Basis_li;

  ifstream myfile (Basis_binary_filename.c_str());
  if (myfile.is_open())
  {
    while ( getline (myfile,line))
    {
      if (line.size() >= r && (line[0] == '0' || line[0] == '1'))
      {
        line2 = line.substr (0,r);          //take the r first characters of line

        //convert string line2 into a binary integer:
        Op_bit = one128 << (r - 1);
        Op = 0;
        for (auto &elem: line2)
        {
          if (elem == c) { Op += Op_bit; }
          Op_bit = Op_bit >> 1;
        }

        Basis_li.push_back(Op); 
      }  
    }
    myfile.close();
  }

  return Basis_li;
}

/******************************************************************************/
/*** VERSION b) Operators are written as the integer values of the binary *****/
/****           representation of the interactions           ******************/
/****           Only up to 64 bits           **********************************/
/******************************************************************************/
// Don't make sense to use for large number of bits 
// 2**50 ~ 10^15 --> printing and reading sucu large integers is not a very good ideas --> better to write it down in binary
// In any case, this function is not very usable for integer of more than 64 bits:
//         1) current function for printing basis doesn't integer with more than 64 bits (from the best basis algorithm)
//         2) "stoi(line)" will convert to a 64 bits integer, which will then be converted to 128 bits: so it doesn't work.

list<__int128_t> Read_BasisOp_IntegerRepresentation(string Basis_integer_filename) //string Basis_integer_filename = basis_IntegerRepresentation_filename)
{
  __int128_t Op = 0;

  cout << "Read the Basis file: " << Basis_integer_filename << endl;

  string line;    

// ***** Store Basis in Basis_li:  **********************************************
  list<__int128_t> Basis_li;

  ifstream myfile (Basis_integer_filename.c_str());
  if (myfile.is_open())
  {
    while ( getline (myfile,line))
    {
      Op = stoi(line);
      Basis_li.push_back(Op);
    }
    myfile.close();
  }

  return Basis_li;
}

/******************************************************************************/
/*************************    Original Basis     ******************************/
/******************************************************************************/
list<__int128_t> Original_Basis(unsigned int r)
{
  __int128_t Op = 1;
  list<__int128_t> Basis_li;

  for (int i=0; i<r; i++)
  {
    Basis_li.push_back(Op);
    Op = Op << 1;
  }

  return Basis_li;
}

/******************************************************************************/
/***************************    Print Basis     *******************************/
/******************************************************************************/
void PrintTerm_Basis(list<__int128_t> Basis_li, unsigned int r)
{
  cout << "## Number of basis operators = " << Basis_li.size() << endl;
  cout << "##" << endl;
  int i = 1;
  for (list<__int128_t>::iterator it = Basis_li.begin(); it != Basis_li.end(); it++)
  {
    cout << "##\t sig_" << setw(3) << setfill(' ') << left << i << " \t " << int_to_bstring((*it), r) << endl; i++;
  } 
  cout << endl;

  cout << "## Convention for converting datapoints in this basis in the following:" << endl;
  cout << "## \t   -- sig_1 = Rightmost bit" << endl;
  cout << "## \t   -- sig_n = Leftmost bit" << endl;
  cout << "##" << endl;
  cout << "## Besides, in the basis above:" << endl;
  cout << "## \t   -- Rightmost bit = Rightmost bit in the data = will be labeled s_1 in the following;" << endl;
  cout << "## \t   -- Leftmost bit  = Leftmost bit  in the data = will be labeled s_n in the following." << endl;
  cout << endl;
}

void PrintTerm_Basis_invert(list<__int128_t> Basis_li, unsigned int r)
{
  cout << "## Number of basis operators = " << Basis_li.size() << endl;
  cout << "##" << endl;
  int i = 1;
  for (list<__int128_t>::iterator it = Basis_li.begin(); it != Basis_li.end(); it++)
  {
    cout << "##\t   s_" << setw(3) << setfill(' ') << left << i << " \t " << int_to_bstring((*it), r) << endl; i++;
  } 
  cout << endl;

  cout << "Convention for reading this basis in comparison with the original basis of the data:" << endl;
  cout << "## \t   -- s_1 = Rightmost bit in the original data" << endl;
  cout << "## \t   -- s_n = Leftmost bit in the original data" << endl;
  cout << "##" << endl;
  cout << "## Besides, in the basis above:" << endl;
  cout << "## \t   -- Rightmost bit = first basis operator, labeled sig_1 above;" << endl;
  cout << "## \t   -- Leftmost bit  = last  basis operator, labeled sig_n above." << endl;
  cout << endl;
}



