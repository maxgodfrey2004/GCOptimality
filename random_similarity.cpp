#include <cassert>
#include <cmath>

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

using std::string;
using std::vector;
using fp_type = long double;

const string AA = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";
const string B1 = "TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG";
const string B2 = "TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG";
const string B3 = "TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG";

std::unordered_map<string, char> one_letter_codes = {
  {"Ala", 'A'}, {"Cys", 'C'}, {"Asp", 'D'}, {"Glu", 'E'}, {"Phe", 'F'},
  {"Gly", 'G'}, {"His", 'H'}, {"Ile", 'I'}, {"Lys", 'K'}, {"Leu", 'L'},
  {"Met", 'M'}, {"Asn", 'N'}, {"Pro", 'P'}, {"Gln", 'Q'}, {"Arg", 'R'},
  {"Ser", 'S'}, {"Thr", 'T'}, {"Val", 'V'}, {"Trp", 'W'}, {"Tyr", 'Y'}
};

std::unordered_map<string, char> universal_gc;
std::unordered_map<char, vector<fp_type>> aa_data;

void load_universal_gc() {
  for (int i = 0; i < static_cast<int>(AA.size()); ++i) {
    universal_gc[string({B1[i], B2[i], B3[i]})] = AA[i];
  }
}

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

fp_type difference(char a1, char a2) {
  fp_type sum_of_squared_diffs = 0;
  assert(aa_data[a1].size() == aa_data[a2].size());
  for (int i = 0; i < static_cast<int>(aa_data[a1].size()); ++i) {
    fp_type diff = aa_data[a1][i] - aa_data[a2][i];
    sum_of_squared_diffs += diff * diff;
  }
  return sqrt(sum_of_squared_diffs);
}

fp_type compute_difference_sum(string aa, string b1, string b2, string b3) {
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

int main(int argc, char *argv[]) {
  if (argc != 2) {
    std::cerr << "Error: correct usage is " << argv[0] << " <input_filename>"
              << std::endl;
    return 1;
  }
  string infile_name = argv[1];
  load_universal_gc();
  load_aadata(infile_name);
  fp_type diff = compute_difference_sum(AA, B1, B2, B3);
  std::cout << diff << std::endl;
}