/********************************************************************/
/**************************    CONSTANTS    *************************/
/********************************************************************/

// number of binary (spin) variables:
unsigned int n = 121;  

// INPUT FOLDER:
const std::string input_directory = "INPUT/";

// OUTPUT FOLDER:
const std::string OUTPUT_directory = "OUTPUT/";


// INPUT DATA FILES (optional):  // must be in binary representation: see Readme file
std::string datafilename = "MNIST11.sorted"; // "SCOTUS_n9_N895_Data.dat";

// INPUT BASIS FILES (optional):  // must be in binary representation: see Readme file
std::string basis_filename = "MNIST11.sorted_BestBasis_k4_Binary.dat";  //"SCOTUS_n9_BestBasis_Binary.dat";


// OTHER INPUT BASIS FILES (optional):
const std::string basis_IntegerRepresentation_filename = "SCOTUS_n9_BestBasis_Integer.dat";        // (optional) Input basis file 
const std::string basis_BinaryRepresentation_filename = "SCOTUS_n9_BestBasis_Binary.dat";      // (optional) Input basis file

// INPUT MCM FILES (optional): // For example: choice of a community to compare with the uncovered community:
const std::string MCM_ex = "INPUT/SCOTUS_Communities_inBestBasis.dat"; 

