#include <iostream>
#include <random>
#include <map>
#include <list>
#include <vector>
#include <fstream>

#include <array>

using namespace std;


/******************************************************************************/
/***************************   ADD OUTPUT FOLDER    ***************************/
/******************************************************************************/
//string Add_Output_Location(string filename);  // defined in main.cpp

/******************************************************************************/
/***************************   Constant variables   ***************************/
/******************************************************************************/
const __int128_t one128 = 1;

/********************************************************************/
/*********************    Proba Structure    ************************/
/********************************************************************/
struct Proba {
    __int128_t s;         // state in the original basis
    __int128_t sig;       // state in the new basis

    double P_D_s = 1.;    // empirical probability of s  
    double P_D_sig = 1.;  // empirical probability of sig --> this should be the same then s if r=n, i.e. if the MCM models all the n spins
    double P_MCM = 1.;  // model probability of s 
};

/******************************************************************************/
/******************   TOOL Functions from "tools.cpp"   ***********************/
/******************************************************************************/
unsigned int Bitset_count(__int128_t bool_nb);
string int_to_bstring(__int128_t bool_nb, unsigned int n);

//check if *Partition* is an actual partition of the basis elements, 
// i.e., that no basis element appears in more than 1 part of the partition.
// i.e., that each basis element only appears in a single part of the partition.
pair<bool, unsigned int> check_partition(map<unsigned int, __int128_t> Partition); 


/******************************************************************************/
/********************   Return Kset over a chosen ICC    **********************/
/******************************************************************************/
__int128_t transform_mu_basis(__int128_t mu, list<__int128_t> basis);

map<__int128_t, unsigned int> build_Kset_ICC(vector<pair<__int128_t, unsigned int>> Kset, __int128_t Ai); // Ai = integer indicated the spins included in the ICC

vector<pair<__int128_t, unsigned int>> build_Kset_ICC_vect(vector<pair<__int128_t, unsigned int>> Kset, __int128_t Ai);

vector<pair<__int128_t, unsigned int>> build_Kset(vector<pair<__int128_t, unsigned int>> Nset, list<__int128_t> Basis);

/******************************************************************************/
/*****************   Compute the contribution to P_MCM(s)   *******************/
/*****************  due to the sub-CM defined by Kset_ICC   *******************/
/******************************************************************************/

void update_proba_MCM(map<__int128_t, Proba> &all_P, map<__int128_t, unsigned int> Kset_ICC, __int128_t Ai, unsigned int N)
{
  map<__int128_t, Proba>::iterator it_P;

  __int128_t s, sig;        // states
//  unsigned int ks=0;      // number of time state s appear in the dataset
  double Nd = (double) N;

//Contribution of the basis "ba" of spins to the model probability P_MCM(s):
  for (it_P = all_P.begin(); it_P!=all_P.end(); ++it_P)
  {
    s = it_P->first;      // initial state s 
    sig = s & Ai;         // troncated state: take only the bits indicated by Ai
    all_P[s].P_MCM *= Kset_ICC[sig]/Nd;
  }
}

/******************************************************************************/
/******************************************************************************/
/*************      Compute and Print the model probability     ***************/
/*********************  for the MCM constructed from Kset  ********************/
/************************   with the given partition   ************************/
/******************************************************************************/
/******************************************************************************/

/*****************      Compute the model probabilities     *******************/
/*******    Computed only for the states that are observed in the data   ******/

// This function can be used directly on the original basis, by replacing Kset by Nset:
map<__int128_t, Proba> P_sig(vector<pair<__int128_t, unsigned int>> Kset, map<unsigned int, __int128_t> Partition, unsigned int N, unsigned int r) // Probabilities in the sigma basis
{
  // Fill in the data probability:
  map<__int128_t, Proba> all_P;
  double Nd = (double) N;

  __int128_t s;        // state
  //unsigned int ks=0; // number of time state s appear in the dataset

  // Check partition:
  pair<bool, unsigned int> Is_partition = check_partition(Partition);
  unsigned int rank = Is_partition.second;
//  cout << "Rank = " << rank << endl;

  if (!Is_partition.first) {cout << "Error, the argument is not a partition: the function returned an empty map for P[s]." << endl; }
  else
  { 
    double pre_factor = 1./((double) (one128 << (r-rank))); 

    for (auto const& it : Kset)
    {   
      s = (it).first;      // initial state s (in new basis)
      //ks = it->second;    // # of times s appears in the data set

      all_P[s].P_D_s = ((double) ((it).second))/Nd;  // (it).second = ks = # of times s appears in the data set
      all_P[s].P_MCM = pre_factor;
    }
  
    // Compute the Kset over each ICC: Kset_ICC:
    //map<__int128_t, unsigned int> Kset_ICC;
    //vector<pair<__int128_t, unsigned int>> Kset_ICC;

    map<unsigned int, __int128_t>::iterator Part;
    for (Part = Partition.begin(); Part != Partition.end(); Part++)
    {
      //Kset_ICC = build_Kset_ICC_vect(Kset, (*Part).second);         // (*Part).second) = Ai = integer indicated the spin elements included in b_a
      //Kset_ICC = build_Kset_ICC(Kset, (*Part).second);         // (*Part).second) = Ai = integer indicated the spin elements included in b_a
      //update_proba_MCM(all_P, Kset_ICC, (*Part).second, N);
      //Kset_ICC.clear();
      update_proba_MCM(all_P, build_Kset_ICC(Kset, (*Part).second), (*Part).second, N);
    }  
  }

  return all_P;
}

/*************      Print the model probabilities     ***************/

void PrintFile_StateProbabilites_CurrentBasis(vector<pair<__int128_t, unsigned int>> Kset, map<unsigned int, __int128_t> MCM_Partition, unsigned int N, unsigned int r, string filename) //, string filename = "Result"
{
  // Probabilities in the sigma basis:
  map<__int128_t, Proba> P_all = P_sig(Kset, MCM_Partition, N, r);
  map<__int128_t, Proba>::iterator it_P;

  //P_k:
  double *Pk_D = (double *)malloc((r+1)*sizeof(double)); 
  double *Pk_MCM = (double *)malloc((r+1)*sizeof(double)); 

  unsigned int k = 0;
  for(k=0; k<=r; k++)
  {
    Pk_D[k] = 0;
    Pk_MCM[k] = 0;
  }

  string Psig_filename = filename + "_DataVSMCM_Ps.dat";
  string Pk_filename = filename + "_DataVSMCM_Pk.dat";

  cout << "--->> Print the state probabilities P(s) in the file: \'" << Psig_filename << "\'" << endl << endl;
  cout << "--->> Print the probability of a state with k \'+1\' bits: \'" << Pk_filename << "\'" << endl << endl;

  // ***** Print P(s):  *****************************************************/
  fstream file_P_sig(Psig_filename, ios::out);

  file_P_sig << "## State probability: computed only for states that are observed in the datasets." << endl;
  file_P_sig << "## Other non-observed states may have a non-zero model probability, but are not reported here." << endl;
  file_P_sig << "## " << endl;
  file_P_sig << "## 1:s \t 2:P_D(s) \t 3:P_MCM(s)" << endl;

  __int128_t s;
  for (it_P = P_all.begin(); it_P!=P_all.end(); ++it_P)
  {   
    s = it_P->first;
    file_P_sig << int_to_bstring(s, r) << "\t " << (it_P->second).P_D_s << "\t " << (it_P->second).P_MCM << endl;

    k = Bitset_count(s);
    Pk_D[k] += (it_P->second).P_D_s;      // P[k] in the data
    Pk_MCM[k] += (it_P->second).P_MCM;    // P[k] from the MCM
  }

  file_P_sig.close();

  // ***** Print P(k):   ***************************************************/
  fstream file_Pk(Pk_filename, ios::out);

  file_Pk << "## 1:k \t 2:P_D(k) \t 3:P_MCM(k)" << endl;

  for(k=0; k<=r; k++)
  {
    file_Pk << k << "\t" << Pk_D[k] << "\t" << Pk_MCM[k] << endl;
  }
  file_Pk.close();
}

/******************************************************************************/
/******************************************************************************/
/*******************      Compute the model probability     *******************/
/**************************  for the MCM constructed   ************************/
/***************   in a given basis, with a given partition   *****************/
/******************************************************************************/
/******************************************************************************/

//// The partition must be guven in the new basis:
map<__int128_t, Proba> P_s(vector<pair<__int128_t, unsigned int>> Nset, list<__int128_t> Basis, map<unsigned int, __int128_t> Partition, unsigned int N, unsigned int r) // Probabilities in the sigma basis
{
  double Nd = (double) N;

  // Fill in the data probability:
  map<__int128_t, Proba> all_P;

  if (!check_partition(Partition).first) {cout << "Error, the argument is not a partition: the function returned an empty map for P[s]." << endl; }
  else
  { 
    // Build Kset:
    //map<__int128_t, unsigned int > Kset_map;
    __int128_t s;        // initial state
    __int128_t sig_m;    // transformed state and to the m first spins
    unsigned int ks=0; // number of time state s appear in the dataset

// ***** Compute P[s] from the original data: *************************************************************************************
    for (auto const& it : Nset)
    {
      s = (it).first;       // state s
      ks = (it).second;     // # of times s appears in the data set
      sig_m = transform_mu_basis(s, Basis);

      // Fill in P[s] empirical and value of transformed state:
      all_P[s].P_D_s = ((double) ks)/Nd;
      all_P[s].sig = sig_m;  

      // Fill in Kset:
      //Kset[sig_m] += ks;
    }

// ***** Build Kset: ***************************************************************************************************************
    vector<pair<__int128_t, unsigned int>> Kset = build_Kset(Nset, Basis);

// ***** Compute P_MCM[s] for the chosen MCM based on the new basis (using "Partition" on Kset): ***********************************
    //cout << "--->> Compute P[s] for the chosen MCM..." << endl << endl;
    map<__int128_t, Proba> all_P_sig = P_sig(Kset, Partition, N, r);   // Compute P[s] for the MCM in the new basis

  // Report the values of P_MCM in the original all_P:
    map<__int128_t, Proba>::iterator it_P;
    for (it_P = all_P.begin(); it_P!=all_P.end(); ++it_P)           // Report P[s] for the MCM in the original basis
    {
      sig_m = (it_P->second).sig;
      (it_P->second).P_MCM = all_P_sig[sig_m].P_MCM;
      (it_P->second).P_D_sig = all_P_sig[sig_m].P_D_s;  // not necessary a probability distribution anymore
    }
  
    all_P_sig.clear();
    Kset.clear();
  }

  return all_P;  
}

/******************************************************************************/
/*****************      PRINT FILE: INFO about an MCM     *********************/
/******************************************************************************/

void PrintFile_MCM_Info(list<__int128_t> Basis, map<unsigned int, __int128_t> MCM_Partition, unsigned int r, string filename); //, string filename = "Result")


/******************************************************************************/
/*************      Print the model probabilities in a file     ***************/
/******************************************************************************/
void PrintFile_StateProbabilites_OriginalBasis_NewBasis(vector<pair<__int128_t, unsigned int>> Nset, list<__int128_t> Basis, map<unsigned int, __int128_t> MCM_Partition, unsigned int N, unsigned int r, string filename) //, string filename = "Result")
{
  // Compute all the state probabilities:
  map<__int128_t, Proba> P_all = P_s(Nset, Basis, MCM_Partition, N, r);

  double *Pk_D = (double *)malloc((r+1)*sizeof(double)); 
  double *Pk_MCM = (double *)malloc((r+1)*sizeof(double)); 

  unsigned int k = 0;
  for(k=0; k<=r; k++)
  {
    Pk_D[k] = 0;
    Pk_MCM[k] = 0;
  }

  string Ps_filename = filename + "_DataVSMCM_Ps_sig.dat";
  string Pk_filename = filename + "_DataVSMCM_Pk.dat";

  //cout << "--->> Print information about the MCM in the file: \'" << filename << "_MCM_info.dat\'" << endl;
  cout << "--->> Print the state probabilities P(s) in the file: \'" << Ps_filename << "\'" << endl;
  cout << "--->> Print the probability of a state with k \'+1\' bits: \'" << Pk_filename << "\'" << endl << endl;

  // ***** Print info about the model -- Print Basis and MCM:  **************/
  //PrintFile_MCM_Info(Basis, MCM_Partition, r, filename);

  // ***** Print P(s):  *****************************************************/
  __int128_t s;

  fstream file_Ps(Ps_filename, ios::out);

  file_Ps << "## s = states in the original basis" << endl;
  file_Ps << "## sig = states in the chosen new basis (ideally the best basis)" << endl;
  file_Ps << "## The chosen \'sig\'-basis and the chosen MCM are printed in the file " << filename << "_MCM_info.dat" << endl;

  // Print P(s): 
  file_Ps << "## " << endl;
  file_Ps << "## 1:s \t 2:P_D(s) \t 3:P_MCM(s) \t 4:sig" << endl;

  for (map<__int128_t, Proba>::iterator it_P = P_all.begin(); it_P!=P_all.end(); ++it_P)
  {   
    s = it_P->first;
    
    file_Ps << int_to_bstring(s, r) << "\t" <<  (it_P->second).P_D_s << "\t" << (it_P->second).P_MCM << "\t" << int_to_bstring((it_P->second).sig, r) << endl;

    k = Bitset_count(s);
    Pk_D[k] += (it_P->second).P_D_s;      // P[k] in the data
    Pk_MCM[k] += (it_P->second).P_MCM;    // P[k] from the MCM
  }
  file_Ps.close();

  // ***** Print P(k):   ***************************************************/
  fstream file_Pk(Pk_filename, ios::out);

  file_Pk << "## 1:k \t 2:P_D(k) \t 3:P_MCM(k)" << endl;

  for(k=0; k<=r; k++)
  {
    file_Pk << k << "\t" << Pk_D[k] << "\t" << Pk_MCM[k] << endl;
  }
  file_Pk.close();
}

/******************************************************************************/
/***************************      SAMPLES         *****************************/
/******************************************************************************/
/*

std::string int_to_bstring(__int128_t bool_nb, unsigned int n);

void PrintTerm_samples(vector<pair<__int128_t, unsigned int>> Kset, map<unsigned int, __int128_t> MCM, list<__int128_t> Basis_li_invert, unsigned int n, unsigned int Nsample)
{
  // Create a random device and use it to generate a random seed
  std::random_device myRandomDevice;
  unsigned seed = myRandomDevice();

  // Initialize a default_random_engine with the seed
  default_random_engine generator(seed);


// ***** Vector of the weights from Kset  ********************************************
  vector<int> weights(Kset.size());

  int i=0;
  for (auto const& my_pair : Kset)
  {
      weights[i]=my_pair.second;
      i++;
  }

  //cout << weights.size() << " : " << Kset.size() << endl;

  cout << "## 1: sample i \t2: state_sig \t 3: state_s" << endl;

  discrete_distribution<int> distribution(weights.begin(), weights.end());

  int sample;
  __int128_t sig = 0, s = 0, Ai = 0;

  for (int i = 0; i < Nsample; i++)
  {
    sig=0;

    for (map<unsigned int, __int128_t>::iterator it = MCM.begin(); it != MCM.end(); it++) // for each ICC sample a new datapoint from the original dataset (in new basis)
    {    
      Ai = (*it).second;

      sample = distribution(generator);
      //cout << int_to_bstring(Kset[sample].first, n) << endl;

      //cout << "##\t ICC_sample_" << i << " = " << int_to_bstring(Ai, n) << endl;

      sig += (Kset[sample].first) & Ai;          // troncated state: take only the bits of s (=it.first) indicated by Ai

      //cout << "##\t ICC_sample_s = " << int_to_bstring(sig, n) << endl;
    }

    s = transform_mu_basis(sig, Basis_li_invert);
    cout << "\t sample " << i << " : " << int_to_bstring(sig, n) << " \t " << int_to_bstring(s, n) << endl;

  }

}


void PrintFile_samples(vector<pair<__int128_t, unsigned int>> Kset, map<unsigned int, __int128_t> MCM, list<__int128_t> Basis_li_invert, unsigned int n, unsigned int Nsample, string filename)
{
  // Create a random device and use it to generate a random seed
  std::random_device myRandomDevice;
  unsigned seed = myRandomDevice();

  // Initialize a default_random_engine with the seed
  default_random_engine generator(seed);


// ***** Vector of the weights from Kset  ********************************************
  vector<int> weights(Kset.size());

  int i=0;
  for (auto const& my_pair : Kset)
  {
      weights[i]=my_pair.second;
      i++;
  }

  string sample_filename = filename + "_samples_N" + to_string(Nsample) + ".dat"; 

  cout << "Print sampled data in:   " << sample_filename << endl;
  cout << "Number of samples: N = " << Nsample << endl;

// ***** create file:
  fstream file_sample(sample_filename, ios::out);

  bool newbasis_bool = (Basis_li_invert.size() != 0);

  if (newbasis_bool)
    { file_sample << "## 1:state_s \t 2: state_sig" << endl; }
  else 
    { file_sample << "## 1:state_s" << endl; }

  discrete_distribution<int> distribution(weights.begin(), weights.end());

  int sample;
  __int128_t sig = 0, s = 0, Ai = 0;

  for (int i = 0; i < Nsample; i++)
  {
    sig=0;

    for (map<unsigned int, __int128_t>::iterator it = MCM.begin(); it != MCM.end(); it++) // for each ICC sample a new datapoint from the original dataset (in new basis)
    {    
      Ai = (*it).second;

      sample = distribution(generator);
      //cout << int_to_bstring(Kset[sample].first, n) << endl;

      //cout << "##\t ICC_sample_" << i << " = " << int_to_bstring(Ai, n) << endl;

      sig += (Kset[sample].first) & Ai;          // troncated state: take only the bits of s (=it.first) indicated by Ai

      //cout << "##\t ICC_sample_s = " << int_to_bstring(sig, n) << endl;
    }

    s = transform_mu_basis(sig, Basis_li_invert);

    if (newbasis_bool)
    {
      file_sample << int_to_bstring(s, n) << "\t ";
    }
    file_sample << int_to_bstring(sig, n) << endl;

//    file_sample << "\t sample " << i << " : " << int_to_bstring(sig, n) << " \t " << int_to_bstring(s, n) << endl;
  }

  file_sample.close();
}


vector<pair<__int128_t, unsigned int>>  histo_sample_bestMCM(vector<pair<__int128_t, unsigned int>> Kset, map<unsigned int, __int128_t> MCM, list<__int128_t> Basis_li_invert, unsigned int n, unsigned int Nsample)
{
  default_random_engine generator;

  map<__int128_t, unsigned int> Kset_sample_map;

// ***** Vector of the weights from Kset  ********************************************
  vector<int> weights(Kset.size());

  int i=0;
  for (auto const& my_pair : Kset)
  {
      weights[i]=my_pair.second;
      i++;
  }

  //cout << weights.size() << " : " << Kset.size() << endl;

//  cout << "## 1: sample i \t2: state_sig \t 3: state_s" << endl;

  discrete_distribution<int> distribution(weights.begin(), weights.end());

  int sample;
  __int128_t sig = 0, s = 0, Ai = 0;

  for (i = 0; i < Nsample; i++)
  {
    sig=0;

    for (map<unsigned int, __int128_t>::iterator it = MCM.begin(); it != MCM.end(); it++) // for each ICC sample a new datapoint from the original dataset (in new basis)
    {    
      Ai = (*it).second;

      sample = distribution(generator);
      //cout << int_to_bstring(Kset[sample].first, n) << endl;

      //cout << "##\t ICC_sample_" << i << " = " << int_to_bstring(Ai, n) << endl;

      sig += (Kset[sample].first) & Ai;          // troncated state: take only the bits of s (=it.first) indicated by Ai

      //cout << "##\t ICC_sample_s = " << int_to_bstring(sig, n) << endl;
    }

    //s = transform_mu_basis(sig, Basis_li_invert);
    //cout << "\t sample " << i << " : " << int_to_bstring(sig, n) << " \t " << int_to_bstring(s, n) << endl;

    Kset_sample_map[sig] += 1;
  }

// ***** Convert map to a vector:  for faster reading later on ********************************
    vector<pair<__int128_t, unsigned int>> Kset_sample(Kset_sample_map.size());    //Nset.resize(Nset_map.size());

    i=0;
    for (auto& my_pair : Kset_sample_map)
    {
        Kset_sample[i]=my_pair;
        i++;
    }
    return Kset_sample;
}
*/

