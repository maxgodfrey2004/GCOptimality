# GCOptimality (for CSB195)

Evaluate the efficiency of the Universal Genetic Code with respect to the similarity of adjacent codons.
We consider two codons to be adjacent if their sequences differ by exactly one nucleotide.

We use a genetic algorithm to evolve the existing universal genetic code with the intention of creating one with a greater similarity.

### Running `random_similarity.cpp`

In order to get `random_similarity.cpp` running, you will need to install R, C++, and GNU Make if you have not already done so.

* Modify the `Makefile` so that `REXE` and `CC` refer to your computer's way of running R and C++ programs respectively.
* Run `make data`.
  * This will produce a file `aadata.txt` containing eight different data points for each amino acid in the Universal Genetic Code.
* Uncomment the instance of the genetic algorithm you would like to run in `random_similarity.cpp` and then run `make run`.
