## DSP
This is a micro repository that will contain common digital signal processing techniques often found in GNSS signal processing. All code is in C.  This includes Gold Code assignments for Pseudo-random number generation, cross correlation, and invoking use of Fast Fourier Transform.  Although there is a quick example ``FFT.c`` which uses the famous fftw3 library, I've preferred to use LiquidDSP, which is a compact end-to-end DSP & wireless communications toolbox in pure C, and my personal favorite.  ``liquiddsp`` also has its own wrapper for fftw3.

The ``gps.c`` example, is a complete example that simulates a transmitted GPS signal, that is received by an ADC that samples at 8.184e3 MSps, with attention paid to the Coarse Acquisition code (C/A) code only, Therefore, the digital input data (C/A code) has a frequency of 1.023 MHz, for simplicity, the IF frequency is assumed to be 2.046 MHz and the sampling rate is assumed to be 8.184 MHz.  These parameters are obtained theoretically from the relation:

$f_o = f_i - n(f_s / 2) ~= f_s / 4$

The output is the plot PSD.pdf which is the power spectrum of a carrier wave of an GPS L1 signal BPSK modulated with the PRN code along with the 50 Hz nav data for a given a satellite, and PSD_AWGN.pdf is the recreated PSD with Complex Additive White Gaussian Noise typically introduced in wireless channel modeling.

For properly compiling ``gps.c``, simply invoke
```gcc gpc.c -o main -lm -lliquid```

To install ``liquiddsp``, check out the main repository of the author: https://github.com/jgaeddert/liquid-dsp.

### To do list:
- [x] Simulate power spectrum of GPS L1 signal
- [x] Add nav data to GPS Signal
- [x] Add complex Additive White Gaussian Noise
- [ ] Model BPSK constellation
- [ ] Carry out despreading, correlation, and tracking



