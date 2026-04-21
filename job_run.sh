
### DATA:
## datafilename = datafile name, must be in INPUT folder
## n = number of binary variables in the datafile
## basisfilename = name of the basis-file, must be in INPUT folder

datafilename=MNIST11.sorted   #Big5PT.sorted  ## must be in INPUT folder
n=121  #50     # number of variables
basisfilename=MNIST11.sorted_BestBasis_k4_Binary.dat  #Big5PT.sorted_BestBasis_Binary.dat   ## must be in INPUT folder

### For Sampling:
## mcmfilename = name of the file containing the MCM you want to sample from, this file must be placed in the OUTPUT folder
## N = number of samples

mcmfilename=MNIST11_MCM_inBnew.dat  ## must be in OUTPUT folder
N=1000


##################################
### RUN in ORIGINAL BASIS:

#time ./MCM_Greedy.out $datafilename $n 

##################################
### RUN in chosen NEW BASIS:

time ./MCM_Greedy.out $datafilename $n -b $basisfilename

##################################
### RUN SAMPLING in ORIGINAL BASIS:

#time ./MCM_Greedy.out $datafilename $n --sample $mcmfilename -N $N

##################################
### RUN SAMPLING in chosen NEW BASIS:

#time ./MCM_Greedy.out $datafilename $n -b $basisfilename --sample $mcmfilename -N $N

