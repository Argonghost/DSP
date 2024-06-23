#include <stdio.h>
#include <stdint.h>
#include <malloc.h>
#include <stdlib.h>
#include <math.h>
#include <liquid/liquid.h>
#include <complex.h>
#include "GoldCodelib.h"
#include <time.h>
#define PI 3.141592654

// The below function, get_code, takes in a given prn code, code length, and repeats the PRN code n times

int* get_code(int* code, int n, int code_len) {
    int* gc = malloc(n * code_len * sizeof(int));

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < code_len; j++) {
            gc[i*code_len + j] = code[j];
        }
    }
    return gc;
}

// The below function acts like np.repeat in python

void repeat_1d(int *src, int src_len, int repeats, int **dst, int *dst_len) {
    *dst_len = src_len * repeats;
    *dst = (int *)malloc((*dst_len) * sizeof(int));

    int j = 0;
    for (int i = 0; i < src_len; i++) {
        for (int k = 0; k < repeats; k++) {
            (*dst)[j++] = src[i];
        }
    }
}
// The below generates the complex sinusoid which serves as the carrier signal

double complex* generate_carrier(double carrier_freq, double code_duration, int num_samples) {
    double complex* carrier_signal = malloc(num_samples * sizeof(double complex));
    double time_step = code_duration / num_samples;
    double amplitude = sqrt(2 * pow(10, -6));
    for (int i = 0; i < num_samples; i++) {
        double t = i * time_step;
        carrier_signal[i] =cexp(I * 2 * PI * carrier_freq * t);
    }
    return carrier_signal;
}

int* generate_nav_data(int num_bits) {
    int* nav_data = (int*)malloc(num_bits * sizeof(int));
    
    // Generate random navigation data (0 or 1)
    for (int i = 0; i < num_bits; i++) {
        nav_data[i] = rand() % 2;
    }
    
    return nav_data;
}

// The below function computes the in-place FFT and Power Spectral Density of our signal (carrier + BPSK)

double* compute_power_spectrum(double complex* signal, int num_samples) {
    float complex* fft_signal = malloc(num_samples * sizeof(float complex));
    for (int i = 0; i < num_samples; i++) {
        fft_signal[i] = (float complex)signal[i];
    }
    
    fftplan p = fft_create_plan(num_samples, fft_signal, fft_signal, LIQUID_FFT_FORWARD, 0);
    fft_execute(p);
    fft_shift(fft_signal, num_samples);
    double* PSD = malloc(num_samples * sizeof(double));
    for (int i = 0; i < num_samples; i++) {
        PSD[i] = cabsf(fft_signal[i] / num_samples) / num_samples;
    }

    fft_destroy_plan(p);
    free(fft_signal);
    return PSD;
}

double randn() {
    static double U, V;
    static int phase = 0;
    double Z;

    if (phase == 0) {
        U = (double)rand() / RAND_MAX;
        V = (double)rand() / RAND_MAX;
        Z = sqrt(-2.0 * log(U)) * sin(2.0 * PI * V);
    } else {
        Z = sqrt(-2.0 * log(U)) * cos(2.0 * PI * V);
    }

    phase = 1 - phase;
    return Z;
}

float complex* generate_noise_matrix(int N, float SNR, float code_duration) {

    time_t t;
    // Seed the random number generator
    srand((unsigned) time(&t));

    float noise_scale = powf(10.0, SNR / 20.0) / sqrt(2.0);
    float complex* W = malloc(sizeof(float complex) * N);

    for (int i = 0; i < N; i++) {
        float real_part = randn();
        float imag_part = randn();
        W[i] = (real_part + I * imag_part) * noise_scale;   
    }
    return W;
}

int main() {
    int sat_24[1023];
    PRN(12, sat_24);
    printf("PRN code for satellite :\n");
    printf("------------------------------\n");
    printf("{");
    for (int i = 0; i < 1023; i++) {
        printf("%d", sat_24[i]);
    }

    printf("}\n");

    int *new = get_code(sat_24, 1023, 1023);
    int *code;
    int code_len;
    repeat_1d(new, 1023, 8, &code, &code_len);

    double carrier_freq = 1575.42e6; // GPS L1 carrier frequency (Hz)
    double code_duration = 1e-3; // 1 ms
    double fs = 8.184e6;
    int num_samples = fs * code_duration;
    float SNRdB = 40.0f; // signal-to-noise ratio [dB]
    float noise_floor = -40.0f; // Noise floor 
    float complex* W = generate_noise_matrix(num_samples, SNRdB, code_duration); // Generate noise matrix 

    // Generate the carrier signal
    double complex* carrier_signal = generate_carrier(carrier_freq, code_duration, num_samples);

    // Generate arbitrary 50 BPs nava data
    int* nav_data = generate_nav_data(50);

    // BPSK modulation
    double complex* bpsk_signal = malloc(num_samples * sizeof(double complex));
    for (int i = 0; i < num_samples; i++) {
        int chip_index = i / (num_samples / 1023);
        int nav_index = i / (num_samples / 50);
        bpsk_signal[i] = (2 * code[chip_index] - 1) *carrier_signal[i] * (2 * nav_data[nav_index] - 1);
    }
    double* PSD = compute_power_spectrum(bpsk_signal, num_samples);

    FILE *pipe_gp = popen("gnuplot -p", "w");
    FILE *output = fopen("output.txt", "w");

    // Write data to output file
    for (int i = 0; i < num_samples; i++) {
        double freq = (i < num_samples/2) ? i * fs / num_samples : (i - num_samples) * fs / num_samples;
        freq /= 1000;
        fprintf(output, "%lf %lf\n", freq, PSD[i]);
    }

    fclose(output);
    

    // Set up gnuplot commands
    fprintf(pipe_gp,"set xlabel 'Frequency [MHz]'\n");
    fprintf(pipe_gp,"set ylabel 'Power Spectrum [dB]'\n");
    fprintf(pipe_gp, "set grid linewidth 1\n");
    fprintf(pipe_gp, "set border linewidth 2\n");
    fprintf(pipe_gp, "set tics font 'Arial,10'\n");
    fprintf(pipe_gp, "set key font 'Arial,10'\n");
    fprintf(pipe_gp,"plot 'output.txt' using 1:2 with lines lc rgb 'blue' title 'Power Spectrum'\n");

    // Close the gnuplot pipe
    fflush(pipe_gp);
    pclose(pipe_gp);

    // Free memory
    free(new);
    free(code);
    free(carrier_signal);
    free(bpsk_signal);
    free(PSD);
    free(nav_data);
    free(W);

    return 0;
}