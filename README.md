## DSP
This is a micro repository that will contain common digital signal processing techniques often found in GNSS signal processing. Code will be in C and C++.  This includes Gold Code assignments for Pseudo-random number generation, cross correlation, and invoking use of Fast Fourier Transform.  Although there is a quick example ``FFT.c`` which uses the famous fftw3 library, I've preferred to use LiquidDSP, which is a compact end-to-end DSP & wireless communications toolbox in pure C, and my personal favorite.  ``liquiddsp`` also has its own wrapper for fftw3.

The ``gps.c`` example, is a complete example that simulates a transmitted GPS signal, that is received by an ADC that samples at 20 MSps, with attention paid to the Coarse Acquisition code (C/A) code only, Therefore, the digital input data (C/A code) has a frequency of 1.023 MHz, for simplicity, the carrier frequency is taken to be 2 MHz, and the CA code is made to be repeated (stretched) at 8 * 1.023e6 MSps  These parameters are obtained theoretically from the relation:

$f_o = f_i - n(f_s / 2) ~= f_s / 4$
This has to do with undersampling vs oversampling, which isn't crucial, since at the end carrier wipeoff and despreading will only leave the nav data.
Please see https://github.com/psas/gps/blob/master/docs/_01_Gold_Codes.ipynb for a better discussion on CA code.
The output is the plot PSD.pdf which is the power spectrum of a carrier wave of an GPS L1 signal BPSK modulated with the PRN code along with the 50 Hz nav data for a given a satellite, and PSD_AWGN.pdf is the recreated PSD with Complex Additive White Gaussian Noise typically introduced in wireless channel modeling.

For properly compiling ``gps.c``, simply invoke
```gcc gpc.c -o main -lm -lliquid```

To install ``liquiddsp``, check out the main repository of the author: https://github.com/jgaeddert/liquid-dsp.
The C++ version does not require any external library, and in my opinion is indeed easier to build off from.

### To do list:
- [x] Simulate power spectrum of GPS L1 signal
- [x] Add nav data to GPS Signal
- [x] Add complex Additive White Gaussian Noise
- [x] Simulate Phase Locked Loop for carrier tracking later
- [ ] Model BPSK constellation
- [ ] Carry out despreading, correlation, and tracking



