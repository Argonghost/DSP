### DSP
This is a micro repository that will contain common digital signal processing techniques often found in GNSS signal processing. All code is in C.  This includes Gold Code assignments for Pseudo-random number generation, cross correlation, and invoking use of Fast Fourier Transform.  Although there is a quick example ``FFT.c`` which uses the famous fftw3 library, I've preferred to use LiquidDSP, which is compact end-to-end DSP & wireless communications toolbox and my personal favorite.  ``liquiddsp`` also has its own wrapper for ffwt3.

The ``gps.c`` example, is a complete example that simulates a transmitted GPS signal, that is received by an ADC that samples at 8.184e3 MSps, with attention paid to the Coarse Acquisition code (C/A) code only, Therefore, the digital input data (C/A code) has a frequency of 1.023 MHz, for simplicity, the IF frequency is assumed to be 2.046 MHz and the sampling rate is assumed to be 8.184 MHz.  These parameters are obtained theoretically from the relation:

$\f_o = f_i - n(f_s / 2) ~= f_s / 4$

For properly compiling ``gps.c``, simply invoke
```gcc gpc.c -o main -lm -lliquid```

To install ``liquiddsp``, check out the main repository of the author: https://github.com/jgaeddert/liquid-dsp.




