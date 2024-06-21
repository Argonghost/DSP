### DSP
This is a micro repository that will contain common digital signal processing techniques often found in GNSS signal processing. All code is in C.  This includes Gold Code assignments for Pseudo-random number generation, cross correlation, and invoking use of Fast Fourier Transform.  Although there is a quick example ``FFT.c`` which uses the famous fftw3 library, I've preferred to use LiquidDSP, which is compact end-to-end DSP & wireless communications toolbox and my personal favorite.  ``liquiddsp`` also has its own wrapper for ffwt3.

For properly compiling ``gps.c``, simply invoke
```gcc gpc.c -o main -lm -lliquid```


