#include <iostream>
#include <iomanip>
#include <map>
#include <list>
#include <vector>
#include <cmath>
#include <fstream>

using namespace std;


/******************************************************************************/
/**********************     CONSTANTS  and TOOLS  *****************************/
/******************************************************************************/
const __int128_t one128 = 1;

/*
string int_to_bstring(__int128_t bool_nb, unsigned int n);
void int_to_digits(__int128_t bool_nb, unsigned int n);
void int_to_digits_file(__int128_t bool_nb, unsigned int r, fstream &file);

unsigned int bitset_count(__int128_t bool_nb);
*/

/**************************************************************************************************************************************************/
/************************************************************   TOOLS MATRIX F2:    ***************************************************************/
/**************************************************************************************************************************************************/

/******************************************************************************/
/**********************   Using REF on Matrix over F2   ***********************/
/******************************************************************************/
void print_matrice(bool** M, int n, int m)
{
  for (int i=0; i<n; i++)
    { 
      for (int j=0; j<m; j++)  {  cout << M[i][j]; }
        cout << endl;
    }
}

/******************************************************************************/
/*******************   Convert  Basis  to  F2 Matrix   ************************/
/******************************************************************************/
// This function returns a binary matrix that has the operators of the basis 'Basis' as columns:
// Each Operator is a column of the matrix
// n = Number of spins = number of rows --> 1rst index          //     !! We placed the lowest bit (most to the right) in the bottom row !!
// m = Number of basis operators = number of columns --> 2nd index

bool** Basis_to_MatrixF2(list<__int128_t> Basis, unsigned int n)
{
  unsigned int m = Basis.size(); // number of operators
  __int128_t Op = 0;

  // Create a Boolean Matrix:
  bool** M = (bool**) malloc(n*sizeof(bool*));  // n rows --> 1rst index
  for (int i=0; i<n; i++)
    {   M[i] = (bool*) malloc(m * sizeof(bool));  }  // m columns --> 2nd index

  unsigned int i = 0; //iteration over the row;
  unsigned int j = 0; //iteration over the columns;

  // Copy the basis operators:
  for (auto& it_Op : Basis)
  {
    Op = it_Op;

    // *****  Filling in each column with a basis operator: ****************
      // -- Last bit (rightmost) at the bottom; first bit (leftmost) at the top.
      // -- First basis operator to the Left (i.e. placed first in the matrix) 
      //    Last basis operator to the Right (placed last)
      //       --> note that this is the opposite of the convention adopted in writing a state in the new basis
      //           (there sig_1 = Rightmost bit)
      // --> One must be careful when re-invert:Matrix to Basis !!!
    
    for (i=0; i<n; i++)
    { 
      M[(n-1-i)][j] = Op & one128;
      Op >>= 1;
    }

    j++; // next column
  } 

  return M;
}

/******************************************************************************/
/*******************   Convert  F2 Matrix  to  Basis   ************************/
/******************************************************************************/
list<__int128_t> MatrixF2_to_Basis(bool** M, unsigned int n)
{
  list<__int128_t> Basis;

  int i = 0; //iteration over the n row;
  int j = 0; //iteration over the m columns;

  __int128_t Op_bin = 1, state = 0;

  // Copy the basis operators from M to Basis:
  for (j=(n-1); j>=0; j--) // Copying each column into an operator: 
  // read operators from the right to the left!!!! 
  // as s1 = rightmost bit = first operator of the Basis_list_invert
  //    sn = left most bit
  { 
    Op_bin = 1;  // sig_1 (for i=0) = the lowest bit // sig_n (for i=(n-1)) = the highest bit
    //Op_bin = one128 << (n - 1);
    state = 0;
    for (i=0; i<n; i++) // Turn a column to an operator:
    {
      // first element of the column = sig_1 (see function "Basis_to_MatrixF2") 
            // should be placed as the Righmost bit for the user (cf. convention adopted in writing a state in the new basis)
      if(M[i][j] == 1)  { state += Op_bin; }    
      //cout << "Op_bin = " << int_to_bstring(Op_bin, n) << endl; 
      Op_bin = Op_bin << 1;    
    } 
    //j++; // next column

    Basis.push_back(state);
    //cout << int_to_bstring(state, n) << endl;
    //string int_to_bstring(__int128_t bool_nb, unsigned int n);
  }

  return Basis;
}

/**************************************************************************************************************************************************/
/************************************************************   GAUSSIAN ELIMINATION:    **********************************************************/
/*************************************************   Return LEAD positions for basis selection   **************************************************/
/**************************************************************************************************************************************************/
// Return Rank and invert matrix

/********************************************************************/
/*********    OPERATIONS for GAUSSIAN ELIMINATION on F2   ***********/
/********************************************************************/
void swap_row(bool** M, int i1, int i2, int n, int m)   //swap L_i1 and L_i2  in the matrix M
{
  if (i1>=n || i2>=n) { cout << "error swap" << endl; }
  else {
    bool* temp = M[i1]; //(bool*) malloc(m*sizeof(bool));
    M[i1]=M[i2];  M[i2]=temp;
    }
}

void add_row(bool** M, int i1, int i2, int n, int m)   //L_i2 <---- L_i1 XOR L_i2
{
  if (i1>=n || i2>=n) { cout << "error add" << endl; }
  else {
    for (int k = 0; k < m; k++)   {   M[i2][k] = ( M[i2][k] != M[i1][k] );   }
    }
}

/********************************************************************/
/**********************   Matrice REF:    ***************************/
/*********   return LEAD positions for basis selection   ************/
/********************************************************************/
// The matrix M is put in Reduced Row Echelon Form:
// And
// Returns a list of the column index where there is a lead position ==>> set of independent operators starting from the left

list<unsigned int> RREF_F2_rank(bool** M, int n, int m)    // (i_lead, j_lead) = positions of the lead
{
  int j_lead = 0;
  int rank=0;    //rank = final number of leads
  bool test_stop=false;

  list<unsigned int> lead_positions;  // list of the column indices in which there is a lead position

  int i=0;

  for (int i_lead = 0; i_lead < n  &&  j_lead < m; i_lead++) 
  {
    i = i_lead;

    // **** search for the 1rst non-zero element in the column
    while (! M[i][j_lead])  // while M[i][j_lead] = 0
    {
      i++;
      if (i == n)   //all elements are 0  --> no lead in this column, 1.go to the next column and 2.restart
      {
        j_lead++;   //1. go to the next column
        if (j_lead >= m)  
          { test_stop=true; break; }  // reduction finished 
        i = i_lead; //2. re-start the search for non-zero element from new (i_lead, j_lead)
      }
    }
    if(test_stop)   {  break;  }

    // **** has found a non-zero element --->  1. increase rank; and 2. swap row i and row i_lead
    rank++; 

    if (i != i_lead)
      { swap_row(M, i_lead, i, n, m); }

    //put to zero every element under the lead (starting from i = i_lead + 1)
    for (i=(i_lead + 1); i < n; i++)
    {
      if (M[i][j_lead])
        { add_row(M, i_lead, i, n, m); }  //L_i <---- L_i XOR L_i_lead
    }

    // Record lead position:
    lead_positions.push_back(j_lead);

    //go to the next column
    j_lead++;
  }

  // Reduced M:
//  cout << "Print reduced M:" << endl;
//  print_matrice(M, n, m);

//  cout << "\t Final rank = " << rank << endl; 
  return lead_positions;
}


/******************************************************************************/
/************   Check if a set of operators are all independent   *************/
/******************************************************************************/
// return 'True' if all the element in the matrix are independent

bool Is_IndepModel(list<__int128_t> Basis_li, unsigned int n)
{
  cout << "-->> Check if the set of operators are independent:" << endl;
  cout << "\t Number of operators analysed: " << Basis_li.size() << endl;
  // Convert the basis to a boolean matrix:
  bool** M_Basis = Basis_to_MatrixF2(Basis_li, n);

  // Row Reduction procedure:
  list<unsigned int> list_lead = RREF_F2_rank(M_Basis, n, Basis_li.size());
  cout << "\t Final rank = " << list_lead.size() << endl; 

  // Free memory:
  for (int i=0; i<n; i++)
    {   free(M_Basis[i]);   } 
  free(M_Basis);

  // Return if the operators are independent or not:
  if (list_lead.size() == Basis_li.size())
    { 
      cout << "\t -->> All the operators are independent." << endl << endl;
      return true;  
    }

  else
    { 
      cout << "\t -->> Not all the operators are independent:";
      cout << "   (number of independent operators = " << list_lead.size() << ") < (number of operators = " << Basis_li.size() << " ) " << endl << endl;
      return false; 
    }
}


/**************************************************************************************************************************************************/
/**************************************************************************************************************************************************/
/************************************************************    INVERT BASIS:    *****************************************************************/
/**************************************************************************************************************************************************/
/**************************************************************************************************************************************************/

// ******** INFO: ********
// The matrix M is put in Reduced Row Echelon Form
// And
// Returns the inverse of M
// ***********************

//bool** 
pair<int, bool**> RREF_F2_invert(bool** M, int n)    // (i_lead, j_lead) = positions of the lead
{
//  cout << "Initial matrix:  " << endl;
//  print_matrice(M, n, n);
//  cout << endl;

  // Create a Identity Matrix:
  bool** M_id = (bool**) malloc(n*sizeof(bool*));  // n rows --> 1rst index
  for (int i=0; i<n; i++)  // n columns --> 2nd index
  {   
    M_id[i] = (bool*) malloc(n * sizeof(bool));  
    for (int j=0; j<n; j++)  // n columns --> 2nd index
    {
      M_id[i][j] = 0;
    } 
    M_id[i][i] = 1;
  }

  // Reduction:
  int j_lead = 0;
  int rank=0;    //rank = final number of leads
  bool test_stop=false;

  //list<unsigned int> lead_positions;  // list of the column indices in which there is a lead position

  int i=0;

  for (int i_lead = 0; i_lead < n  &&  j_lead < n; i_lead++)   // 
  {
    i = i_lead;

    //search for the 1rst non-zero element in the column
    while (! M[i][j_lead])  // while M[i][j_lead] = 0
    {
      i++;
      if (i == n)   //all elements are 0  --> no lead in this column, 1.go to the next column and 2.restart
      {
        j_lead++;   //1. go to the next column
        if (j_lead >= n)  
        { 
          test_stop=true; break; 
        }  // reduction finished
        i = i_lead; //2. re-start the search for non-zero element from new (i_lead, j_lead)
      }
    }
    if(test_stop)   {  break;  }

    //has found a non-zero element --->  1. increase rank; and 2. swap row i and row i_lead
    rank++; 
    if (i != i_lead)
    { 
      swap_row(M, i_lead, i, n, n); 
      swap_row(M_id, i_lead, i, n, n); 
    }

    //put to zero all the element at i != i_lead:
    for (i=0; i < n; i++)
    {
      if (M[i][j_lead] && i!=i_lead)  //L_i <---- L_i XOR L_i_lead
      { 
        add_row(M, i_lead, i, n, n); 
        add_row(M_id, i_lead, i, n, n); 
      }
    }

    // Record lead position:
    //lead_positions.push_back(j_lead);  

    //go to the next column
    j_lead++;
  }
/*
  if(rank == n)
  {
    cout << "Rank = n = " << n << "\t: this is a basis." << endl;
  }
  else 
  {
    cout << "Rank = " << rank << " different from n = " << n << "\t: this is not a basis." << endl;
    cout << "Note that the inverse transformation provided is therefore incomplete." << endl;
  }    
*/
//  cout << "Inverted Matrice:" << endl;
//  print_matrice(M_id, n, n);

  return make_pair(rank, M_id); //M_id;
}

// ******** INFO: ********
// n = number of variables
// m = number of operators in the independent set
// if m > n : the set cannot be independent: stop the procedure;
// if n=m: check if rank=n, then everything is good and return the inverse basis
// if m < n:  even if it is an independent set, it is not a basis, and may not be invertable: stop the procedure.
// ***********************

list<__int128_t> Invert_Basis(list<__int128_t> Basis_li, unsigned int n)
{
  list<__int128_t> Basis_li_invert;

  if (Basis_li.size() == n)
  {
    bool** M_Basis = Basis_to_MatrixF2(Basis_li, n);

    // Row Reduction procedure to invert:
    pair<int, bool**> M_invert = RREF_F2_invert(M_Basis, n);

    if(M_invert.first == n) // M_invert.first = rank
    {
      cout << "Rank = n = " << n << "\t: this is a basis and can be inverted." << endl << endl;
      Basis_li_invert = MatrixF2_to_Basis(M_invert.second, n);
    }
    else 
    {
      cout << "The Rank = " << M_invert.first << " is smaller than the number of variables, n = " << n << "\t: this is not a basis." << endl << endl;
      //cout << "Note that the inverse transformation provided is therefore incomplete." << endl;
    } 

    // Free memory:
    for (int i=0; i<n; i++)
    {   
      free(M_Basis[i]);   free(M_invert.second[i]);   
    } 
    free(M_Basis);
    free(M_invert.second);
  }
  else
  {
    cout << "The number of operators in the basis provided, m=" << Basis_li.size() << ", is not equal to the number of spin variables, n=" << n << ":" << endl << endl;
  }

  return Basis_li_invert;
}

/*
vector<Operator128> Invert_Basis(vector<Operator128> Basis, unsigned int n)
{
  vector<Operator128> Basis_invert;

  if (Basis.size() == n)
  {
    bool** M_Basis = Basis_to_MatrixF2(Basis, n);

    // Row Reduction procedure to invert:
    pair<int, bool**> M_invert = RREF_F2_invert(M_Basis, n);

    if(M_invert.first == n) // M_invert.first = rank
    {
      cout << "Rank = n = " << n << "\t: this is a basis and can be inverted." << endl << endl;
      Basis_invert = MatrixF2_to_Basis(M_invert.second, n);
    }
    else 
    {
      cout << "The Rank = " << M_invert.first << " is smaller than the number of variables, n = " << n << "\t: this is not a basis." << endl << endl;
      //cout << "Note that the inverse transformation provided is therefore incomplete." << endl;
    } 

    // Free memory:
    for (int i=0; i<n; i++)
    {   
      free(M_Basis[i]);   free(M_invert.second[i]);   
    } 
    free(M_Basis);
    free(M_invert.second);  
  }
  else
  {
    cout << "The number of operators in the basis provided, m=" << Basis.size() << ", is not equal to the number of spin variables, n=" << n << ":" << endl << endl;
  }

  return Basis_invert;
}
*/





