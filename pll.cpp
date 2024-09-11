#include <iostream>
#include <fstream>
#include <cmath>
#include <complex>
#include <iomanip>
using namespace std;
using namespace std::complex_literals;

/*This is a simple second order Phase locked loop implementation*/
class SimulatePLL{
    private:
    float        phase_offset;  //     carrier phase offset
    float        frequency_offset;  // carrier frequency offset
    float        wn;  //               pll bandwidth
    float        zeta;              // pll damping factor
    float        K;   //               pll loop gain
    unsigned int n;    //              number of samples
    float t1; //                       tau_1
    float t2; //                       tau_2
    // b0,b1,b2 are feed forward coefficients

    float b0; 
    float b1;
    float b2;

    // a1, a2 are feedback coefficients
    float a1;
    float a2;

    // v0,v1,v2 are filter buffer
    float v0;
    float v1;
    float v2;

    float phi_hat; // PLL initial phase
    float phi;
    public:
    /* Default constructor that sets the simulation parameters upon calling*/
    SimulatePLL(){
        phase_offset = 0.00;
        frequency_offset = 0.30;
        wn = 0.01;
        zeta = 0.707;
        K = 1000;
        n = 400;
        t1 = K/(wn*wn); 
        t2 = 2*zeta/wn; 
        b0 = (4*K/t1)*(1.+t2/2.0f);
        b1 = (8*K/t1);
        b2 = (4*K/t1)*(1.-t2/2.0f);
        a1 = -2.0f;
        a2 =  1.0f;
        v0=0.0f, v1=0.0f, v2=0.0f;
        phi     = phase_offset;   
        phi_hat = 0.0f;      
        cout << "index" << std::setw(15) << "real(x)"<< std::setw(15) <<  "imag(x)"  << std::setw(15) << "real(y)"<< std::setw(15) <<  "imag(y)"<< std::setw(15) <<  "error" << endl;
    };

    void generate_loop_filter(){
        complex<float> x, y;
        for(size_t i = 0; i < n; i++){
            x = complex<float>(cosf(phi), sinf(phi)); 
            phi += frequency_offset;  // Update carrier phase

            // Estimated signal
            y = complex<float>(cosf(phi_hat), sinf(phi_hat));
            float delta_phi = arg(x * conj(y));
            cout << i << std::setw(15) << x.real() << std::setw(15) << x.imag() << std::setw(15) << y.real() << std::setw(15) <<y.imag() << std::setw(15) << delta_phi << endl;
            v2 = v1;  // shift center register to upper register
            v1 = v0;  // shift lower register to center register
            v0 = delta_phi - v1*a1 - v2*a2; // compute new lower register

            // compute new output
            phi_hat = v0*b0 + v1*b1 + v2*b2;
        };
    };


};


int main(){

    SimulatePLL pll;
    pll.generate_loop_filter();

    return 0;
};
