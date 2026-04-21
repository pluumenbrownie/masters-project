#include <iostream>
#include <random>
#include <map>
#include <list>
#include <vector>
#include <fstream>

using namespace std;

/******************************************************************************/
/******************   TOOL Functions from "tools.cpp"   ***********************/
/******************************************************************************/
std::string int_to_bstring(__int128_t bool_nb, unsigned int n);

__int128_t transform_mu_basis(__int128_t mu, list<__int128_t> basis);


/******************************************************************************/
/***************************      SAMPLES         *****************************/
/******************************************************************************/
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

