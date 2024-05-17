#include <stdio.h>
#include <complex.h>
#include <math.h>

#define N 100
#define Fs 20000.0
#define PI 3.141592653589793

int main() {
    double time[N];
    double complex IQ[N], IQ2[N];
    double complex correlated[2*N - 1];

    // Generate time axis
    for (int i = 0; i < N; i++) {
        time[i] = i / Fs;
    }

    // Generate I1 and Q1
    double f1 = 500.0;
    for (int i = 0; i < N; i++) {
        IQ[i] = cexp(I * 2 * PI * f1 * time[i]);
    }

    // Generate I2 and Q2
    for (int i = 0; i < N; i++) {
        IQ2[i] = cexp(I * 2 * PI * f1 * time[i] + PI / 2);
    }

    // Cross-correlation

    FILE* XCORR = fopen("xcorr.txt", "w");
    for (int i = 0; i < 2*N - 1; i++) {
        correlated[i] = 0.0 + 0.0 * I;
        for (int j = 0; j < N; j++) {
            int k = i - N + 1 + j;
            if (k >= 0 && k < 2*N - 1) {
                correlated[i] += conj(IQ[j]) * IQ2[k];
            }
        }
        correlated[i] /= N;
    }

    double lagTime;

    
    for (int i = 0; i < 2*N - 1; i++) {
        lagTime = (i - N + 1) / Fs;
        fprintf(XCORR,"%f %f\n", lagTime, creal(correlated[i]));
    }
    fclose(XCORR);



    FILE *pipe_gp = popen("gnuplot -persistent", "w");

    fprintf(pipe_gp, "set xlabel 'lag'\n");
    fprintf(pipe_gp, "set ylabel 'Re{FFT^-1}'\n");
    fprintf(pipe_gp, "set title 'Real-Time Plot'\n");
	fprintf(pipe_gp, "set grid\n");
    fprintf(pipe_gp, "plot 'xcorr.txt' using 1:2 with lines\n");
    pclose(pipe_gp);



    return 0;
}