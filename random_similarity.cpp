#include <cassert>
#include <cmath>
#include <cstdlib>

#include <algorithm>
#include <fstream>
#include <functional>
#include <iostream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

using std::string;
using std::vector;
using fp_type = double;  // Type alias in case this needs changing later.

// Represents the universal genetic code.
const string AA = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";
const string B1 = "TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG";
const string B2 = "TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG";
const string B3 = "TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG";

// List of all one-letter amino acid codes for random selection (including stop codon).
const string acid_names = "*ACDEFGHIKLMNPQRSTVWY";

// Maps three-letter amino acid names to their corresponding one-letter names.
std::unordered_map<string, char> one_letter_codes = {
  {"Ala", 'A'}, {"Cys", 'C'}, {"Asp", 'D'}, {"Glu", 'E'}, {"Phe", 'F'},
  {"Gly", 'G'}, {"His", 'H'}, {"Ile", 'I'}, {"Lys", 'K'}, {"Leu", 'L'},
  {"Met", 'M'}, {"Asn", 'N'}, {"Pro", 'P'}, {"Gln", 'Q'}, {"Arg", 'R'},
  {"Ser", 'S'}, {"Thr", 'T'}, {"Val", 'V'}, {"Trp", 'W'}, {"Tyr", 'Y'}
};

// Maps codons to their one-letter amino acid name.
std::unordered_map<string, char> universal_gc;

// Stores numerical data for each amino acid.
std::unordered_map<char, vector<fp_type>> aa_data;

// Builds a mapping from the constants AA, B1, B2, and B3.
//
void load_universal_gc() {
  for (int i = 0; i < static_cast<int>(AA.size()); ++i) {
    universal_gc[string({B1[i], B2[i], B3[i]})] = AA[i];
  }
}

// Loads amino acid data passed to us by ./gen_aadata.r
//
void load_aadata(string infile_name) {
  std::ifstream fin(infile_name);
  string line;
  // Ignore the first line of data as it contains meaningless column headings.
  getline(fin, line);
  // Read in data for the 20 Amino Acids.
  for (int i = 0; i < 20; ++i) {
    getline(fin, line);
    std::stringstream data;
    data << line;
    string aa_name;
    data >> aa_name;
    vector<fp_type> values;
    for (fp_type val; data >> val; ) {
      values.push_back(val);
    }
    // Currently, `aa_name` has a leading and trailing quotation mark. We now
    // remove them.
    aa_name = string({aa_name[1], aa_name[2], aa_name[3]});
    aa_data[one_letter_codes[aa_name]] = values;
  }
  int num_categories = static_cast<int>(aa_data['A'].size());
  aa_data['*'] = vector<fp_type>(num_categories, fp_type(0));
}

// Computes the difference score of two amino acids.
//
fp_type difference(char a1, char a2) {
  fp_type sum_of_squared_diffs = 0;
  assert(aa_data[a1].size() == aa_data[a2].size());
  for (int i = 0; i < static_cast<int>(aa_data[a1].size()); ++i) {
    fp_type diff = aa_data[a1][i] - aa_data[a2][i];
    sum_of_squared_diffs += diff * diff;
  }
  return sqrt(sum_of_squared_diffs);
}

// Computes the total difference of a genetic code mapping.
//
fp_type compute_difference_sum(string aa, string b1 = B1, string b2 = B2, string b3 = B3) {
  fp_type diff = 0;
  for (int i = 0; i < static_cast<int>(aa.size()); ++i) {
    for (int j = 0; j < i; ++j) {
      // Calculate the number of different bases between the codons.
      int d = int(b1[i] != b1[j]) + int(b2[i] != b2[j]) + int(b3[i] != b3[j]);
      // If the codons differ by one base, add the difference of their
      // corresponding amino acids to the sum.
      if (d == 1) {
        diff += difference(aa[i], aa[j]);
      }
    }
  }
  return diff;
}

// ======================= Genetic Algorithm Utilities ========================

namespace mutation {
  // Returns a random integer in the range [low, high].
  //
  int randint(int low, int high) {
    return low + (rand() % static_cast<int>(high - low + 1));
  }

  // Returns a randomly chosen character from a string.
  //
  char rand_choice(string s) {
    return s[randint(0, static_cast<int>(s.size()) - 1)];
  }

  // Randomly replace a single character in an amino acid map with another
  // amino acid.
  //
  string random_replace(string aa) {
    // Strings are formatted in such a way that the i-th element of a string
    // corresponds to the codon B1[i]+B2[i]+B3[i].
    int replace_idx = randint(0, static_cast<int>(aa.size()));
    aa[replace_idx] = rand_choice(acid_names);
    return aa;
  }

  // Mutation where a genetic code mapping produces three offspring. One is
  // a replica, and the other two have one codon mapped to a different animo
  // acid.
  //
  vector<string> random_replace_offspring(string aa) {
    // Strings are formatted in such a way that the i-th element of a string
    // corresponds to the codon B1[i]+B2[i]+B3[i].
    vector<string> mutations({
      aa,
      mutation::random_replace(aa),
      mutation::random_replace(aa)
    });
    return mutations; 
  }
}  // namespace mutations

// Evolves a genetic code for a specified number of generations with a
// specified population size. A routine which generates offspring for a given
// member of the population is also passed to the routine.
//
void evolve_genetic_code(
  int num_generations,
  int population_size,
  std::function<vector<string>(string)> gen_offspring
) {
  // Create a starting pool.
  std::cout << "Universal Genetic Code:   " << AA << "\tScore: "
            << compute_difference_sum(AA) << std::endl;
  vector<string> population({AA});
  while (population.size() < population_size) {
    auto additions = gen_offspring(AA);
    population.insert(population.end(), additions.begin(), additions.end());
  }
  while (population.size() > population_size) {
    population.pop_back();
  }
  // Evolve the starting pool.
  for (int generation = 0; generation < num_generations; ++generation) {
    vector<std::pair<fp_type, string>> all_offspring;
    for (string parent : population) {
      vector<string> offspring = gen_offspring(parent);
      for (string child : offspring) {
        all_offspring.emplace_back(compute_difference_sum(child), child);
      }
    }
    std::sort(all_offspring.begin(), all_offspring.end());
    // Lowest difference scores will be at the front of the list.
    // Make the new generation the current population.
    for (int i = 0; i < population_size; ++i) {
      population[i] = all_offspring[i].second;
    }
    std::cout << "Generation " << generation + 1 << "\tBest map: "
              << all_offspring[0].second << "\tScore: "
              << all_offspring[0].first << std::endl;
  }
}

int main(int argc, char *argv[]) {
  if (argc != 2) {
    std::cerr << "Error: correct usage is " << argv[0] << " <data_filename>"
              << std::endl;
    return 1;
  }
  string infile_name = argv[1];
  load_universal_gc();
  load_aadata(infile_name);
  srand(100);  // Get consistent results every time by setting a random seed.
  
  // 1.1
  // evolve_genetic_code(10, 10, mutation::random_replace_offspring);

  // 1.2
  // evolve_genetic_code(1000, 100, mutation::random_replace_offspring);
}