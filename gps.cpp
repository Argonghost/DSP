#include <iostream>
#include <vector>
#include <string>
#include <iomanip>
#include <array>
#include <complex>
#include <cmath>
#include <fstream>
#include <cstdlib> // C - library to invoke rand()
#define PI 3.141592654
using namespace std;
using namespace std::complex_literals;

std::vector<int16_t> GEN_CA_CODE(int sv){

  int NUM_CODES =           37;             // Number of satellite id's
  int SR_LEN =              20;             
  int CA_PERIOD =           1023;
  int nchips =              CA_PERIOD;
  int LEN =                 nchips;
  std::vector<int16_t>      code(LEN,1);
  int FB_TAPS[8] =          {0, 1, 2, 4, 9, 12, 15, 18};
  int nShifts =             nchips - SR_LEN; // First SR_LEN bits are already in the shift reg

  // Init states of genL1_CA_INIT_TABLE
  int CA_INIT_STATES[37][20] = {
    { 1, 1, 0, 0, 1, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 1, 0, 1, 0, 0 },
    { 1, 1, 1, 0, 0, 1, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 1, 1 },
    { 1, 1, 1, 1, 0, 0, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0, 1, 0, 0, 0 },
    { 1, 1, 1, 1, 1, 0, 0, 1, 0, 0, 1, 1, 1, 0, 0, 0, 1, 1, 0, 1 },
    { 1, 0, 0, 1, 0, 1, 1, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 1, 1, 0 },
    { 1, 1, 0, 0, 1, 0, 1, 1, 0, 1, 0, 0, 0, 1, 1, 0, 1, 0, 1, 0 },
    { 1, 0, 0, 1, 0, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1 },
    { 1, 1, 0, 0, 1, 0, 1, 1, 0, 0, 0, 1, 1, 0, 1, 1, 0, 1, 1, 1 },
    { 1, 1, 1, 0, 0, 1, 0, 1, 1, 0, 1, 0, 1, 0, 0, 1, 0, 0, 1, 0 },
    { 1, 1, 0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 1 },
    { 1, 1, 1, 0, 1, 0, 0, 0, 1, 0, 1, 1, 0, 0, 0, 0, 1, 1, 0, 1 },
    { 1, 1, 1, 1, 1, 0, 1, 0, 0, 0, 0, 1, 1, 0, 1, 0, 1, 1, 1, 0 },
    { 1, 1, 1, 1, 1, 1, 0, 1, 0, 0, 1, 0, 1, 0, 0, 1, 1, 1, 1, 0 },
    { 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 1, 1, 0, 0, 0, 0, 0, 1, 1, 0 },
    { 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0, 0, 1, 0, 1, 0 },
    { 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 0, 1, 0, 1, 1, 0, 0 },
    { 1, 0, 0, 1, 1, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0 },
    { 1, 1, 0, 0, 1, 1, 0, 1, 1, 1, 1, 0, 0, 1, 0, 0, 1, 0, 1, 1 },
    { 1, 1, 1, 0, 0, 1, 1, 0, 1, 1, 0, 1, 0, 1, 1, 0, 1, 1, 0, 0 },
    { 1, 1, 1, 1, 0, 0, 1, 1, 0, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1 },
    { 1, 1, 1, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 1, 1, 0, 1, 1, 0 },
    { 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 0, 0, 1, 0, 1, 0, 0, 1, 0 },
    { 1, 0, 0, 0, 1, 1, 0, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 0, 0 },
    { 1, 1, 1, 1, 0, 0, 0, 1, 1, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0 },
    { 1, 1, 1, 1, 1, 0, 0, 0, 1, 1, 1, 1, 0, 1, 0, 0, 0, 0, 0, 1 },
    { 1, 1, 1, 1, 1, 1, 0, 0, 0, 1, 0, 1, 1, 1, 1, 0, 1, 0, 0, 1 },
    { 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 1, 0, 1, 1, 1, 1, 0, 1 },
    { 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 1, 0, 0, 0, 0, 1, 0, 1, 1, 1 },
    { 1, 0, 0, 1, 0, 1, 0, 1, 1, 1, 1, 0, 0, 1, 1, 0, 0, 1, 0, 1 },
    { 1, 1, 0, 0, 1, 0, 1, 0, 1, 1, 0, 1, 0, 1, 1, 1, 1, 0, 1, 1 },
    { 1, 1, 1, 0, 0, 1, 0, 1, 0, 1, 0, 0, 1, 1, 1, 1, 0, 1, 0, 0 },
    { 1, 1, 1, 1, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1 },
    { 1, 1, 1, 1, 1, 0, 0, 1, 0, 1, 1, 0, 0, 1, 0, 1, 0, 0, 0, 0 },
    { 1, 0, 0, 1, 0, 1, 1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0 },
    { 1, 1, 1, 1, 0, 0, 1, 0, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 0 },
    { 1, 1, 0, 0, 1, 0, 1, 1, 1, 0, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0 },
    { 1, 1, 1, 1, 0, 0, 1, 0, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 0 }
  };

  if (sv < 1 || sv > NUM_CODES){
    std::cout << " ERROR: sv must be in the range 1 - " << NUM_CODES << endl;
  } 
  else {
    // Initialise: Put the init state at the start of the code
    for ( int i = 0; i < 20; i++ ){
      code[i] = 1-2*CA_INIT_STATES[sv-1][i];
  }

    int *k = FB_TAPS;

    // Shifts through all the registers
    for (int i = 0; i < nShifts; i++){
      int row_prod = (code[k[0]]*code[k[1]]*code[k[2]]*code[k[3]]*code[k[4]]*code[k[5]]*code[k[6]]*code[k[7]]);
      code[i+SR_LEN] = row_prod;
      for (int j = 0; j < 8; j++){
        k[j] = k[j] + 1;
      }
    }
  }

  for (int i = 0; i < LEN; i++){
    code[i] = (1-code[i])/2;
  }

  return code;

};

class Base_Signal{
  private:
  double carrier_frequency;
  int num_samples;
  double code_duration;
  int chip_rate;
  int samples_per_chip;
  double sampling_frequency;

  public:
  Base_Signal(){
    carrier_frequency =             2.0e6;
    code_duration =                 10e-3;
    sampling_frequency =            20e6;
    num_samples =                   sampling_frequency * code_duration;
    chip_rate =                     1.023e6;
    samples_per_chip =              sampling_frequency / chip_rate;

  };

  vector<complex<double>> create_carrier(){
    vector<complex<double>> signal(num_samples);  
      const std::complex<double> i(0.0, 1.0);
      for(size_t j = 0; j < num_samples; j++){
        signal[j] = exp(i * (2.0 * PI * carrier_frequency * j / sampling_frequency));  // Create complex sinusoid
      };
      return signal;
  };


};

class BPSK: public Base_Signal{
  private:
  int sv;

  public:
  BPSK(){
    cout << "Creating Satellite SV" << endl;
  };
  BPSK(int s): sv{s}{
    cout << "Created Satellite with ID no : " << s << endl;
  };
  vector<int16_t> create_code(){
    vector<int16_t> vec = GEN_CA_CODE(this->sv);
    return vec;
  };

  void PRINT_RE(vector<double> CA);
  void PRINT_IM(vector<complex<double>> CA);

  vector<int16_t> REPEAT(vector<int16_t>& data, std::size_t count);
  vector<int16_t> create_nav();
  vector<double> gps_L1_modulate(Base_Signal& signal);
};
vector<int16_t> BPSK::create_nav() {
    int num_bits = 50; // Number of bits
    vector<int16_t> nav_data(num_bits);
    for (size_t i = 0; i < num_bits; i++) {
        nav_data[i] = rand() % 2;
    };

    return nav_data;

}
// vector<int16_t> BPSK::REPEAT(vector<int16_t>& data, std::size_t count) {
//     auto pattern_size = data.size();
//     data.resize(pattern_size * count);
//     const auto pbeg = data.begin();
//     const auto pend = std::next(pbeg, pattern_size);
//     auto it = std::next(data.begin(), pattern_size);
//     for(std::size_t k = 1; k < count; ++k) {
//         std::copy(pbeg, pend, it);
//         std::advance(it, pattern_size);
//     };

//     return data;
// };

vector<int16_t> BPSK::REPEAT(vector<int16_t>& data, std::size_t count) {
    vector<int16_t> repeated_data; // New vector to store repeated elements
    repeated_data.reserve(data.size() * count); // Reserve space for efficiency

    for (const auto& elem : data) {
        // Append `count` copies of `elem` to the repeated_data vector
        for (std::size_t i = 0; i < count; ++i) {
            repeated_data.push_back(elem);
        }
    }

    return repeated_data;
}

void BPSK::PRINT_RE(vector<double> CA){
    cout << "[";
    for(size_t i = 0; i < CA.size(); i++){
        cout << i << " " << CA[i]<< endl;
    };
    cout << "]";
};

void BPSK::PRINT_IM(vector<complex<double>> CA){
    cout << "[";
    for(size_t i = 0; i < CA.size(); i++){
        cout << i << " " << real(CA[i])<< endl;
    };
    cout << "]" << endl;
};

vector<double> BPSK::gps_L1_modulate(Base_Signal& signal){
  vector<complex<double>> complex_signal = signal.create_carrier();
  vector<int16_t> mod0 =                     GEN_CA_CODE(sv);
  vector<int16_t> mod1 = REPEAT(mod0, 8); // Creating CA code
  int code_size =                          mod1.size(); 

  vector<int16_t> mod2 =                     create_nav(); // Creating 50 bit Nav data
 
  /* Now lets XOR both modulation*/

  vector<double> modulated_signal(code_size);
  for(size_t i = 0; i < code_size; i++){
    int nav_index = (i * 50 / code_size);
    int mod_bits = mod1[i] ^ mod2[nav_index];
      // mod[i] = real(signal[i]) * static_cast<double>(code[i]);
    modulated_signal[i] = real(complex_signal[i]) * (2 * mod_bits - 1);
  };

  return modulated_signal;

};


int main(){

  /* Driver code ..... */
    int sv;
    cout << "Enter space vehicle number: " << endl;
    cin >> sv;
    
    Base_Signal base_signal;
    BPSK modulate(sv);

    modulate.gps_L1_modulate(base_signal);
    vector<int16_t> nav = modulate.create_nav();
    vector<int16_t> ca_code = modulate.create_code();
    vector<int16_t> new_code = modulate.REPEAT(ca_code, 8);

    vector<double> L1_signal = modulate.gps_L1_modulate(base_signal);
    ofstream outfile;
    outfile.open("gps_data_output.txt");
    for(size_t i = 0; i < L1_signal.size(); i++){
        int nav_index = (i * 50 / L1_signal.size()); 
        int CA_index = i * 1023 / L1_signal.size();
        outfile << i << std::setw(10) << CA_index <<std::setw(10) << (2*nav[nav_index] - 1)<< std::setw(10) << new_code[i] <<std::setw(10) <<  real(L1_signal[i])<< endl;
    };

    outfile.close();
    cout << "Tested successfully" << endl;    
    cout << "Done" << endl;
    cout << "Size of code is : " << new_code.size() << endl;
    for(auto i: new_code){
      cout << i << " ";
    };
      
    return 0;
};