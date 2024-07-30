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

// Function to repeat PRN code
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

// Function to generate carrier signal
double complex* generate_carrier(double carrier_freq, double code_duration, int num_samples) {
    double complex* carrier_signal = malloc(num_samples * sizeof(double complex));
    for (int i = 0; i < num_samples; i++) {
        carrier_signal[i] = cexp(I * 2 * PI * carrier_freq * i / num_samples);
    }
    return carrier_signal;
}

// Function to generate navigation data
int* generate_nav_data(int num_bits) {
    int* nav_data = (int*)malloc(num_bits * sizeof(int));
    srand(time(NULL));
    for (int i = 0; i < num_bits; i++) {
        nav_data[i] = rand() % 2;
        if (i == num_bits - 1) {
            nav_data[i] = !nav_data[0];
        }
    }
    return nav_data;
}

// Function to compute Power Spectral Density
double* compute_power_spectrum(double complex* signal, int num_samples) {
    float complex* fft_signal = malloc(num_samples * sizeof(float complex));
    for (int i = 0; i < num_samples; i++) {
        fft_signal[i] = (float complex)signal[i];
    }
    
    fftplan p = fft_create_plan(num_samples, fft_signal, fft_signal, LIQUID_FFT_FORWARD, 0);
    fft_execute(p);
    double* PSD = malloc(num_samples * sizeof(double));
    for (int i = 0; i < num_samples; i++) {
        PSD[i] = cabsf(fft_signal[i]) / num_samples;
    }

    fft_destroy_plan(p);
    free(fft_signal);
    return PSD;
}

int main() {
    int sat_24[1023];
    PRN(1, sat_24);
    
    double carrier_freq = 2.0e6; // GPS L1 carrier frequency (Hz)
    double code_duration = 10e-3; // 10 ms
    double fs = 8 * 1.023e6; // Sampling frequency (Hz)
    int num_samples = fs * code_duration; // Total number of samples
    int chip_rate = 1.023e6;
    int samples_per_chip = fs / chip_rate;

    float SNRdB = -158.5f; // Signal-to-noise ratio [dB]
    float noise_floor = -20.0f; // Noise floor 
    float nstd = powf(10.0f, (SNRdB - noise_floor) / 20.0f);

    // Generate the carrier signal
    double complex* carrier_signal = generate_carrier(carrier_freq, code_duration, num_samples);

    // Upsample PRN code
    int *upsampled_code;
    int upsampled_code_len;
    repeat_1d(sat_24, 1023, samples_per_chip, &upsampled_code, &upsampled_code_len);

    // Ensure the upsampled code length matches num_samples
    if (upsampled_code_len < num_samples) {
        upsampled_code = realloc(upsampled_code, num_samples * sizeof(int));
        for (int i = upsampled_code_len; i < num_samples; i++) {
            upsampled_code[i] = upsampled_code[i % upsampled_code_len];
        }
        upsampled_code_len = num_samples;
    }

    // Generate fake nav data
    int* nav_data = generate_nav_data(50);

    // BPSK modulation
    double complex* bpsk_signal = malloc(num_samples * sizeof(double complex));
    for (int i = 0; i < num_samples; i++) {
        int nav_index = (i * 50 / num_samples);
        int data_bit = upsampled_code[i] ^ nav_data[nav_index]; // Modulo-2 addition (XOR)
        bpsk_signal[i] = (2 * data_bit - 1) * carrier_signal[i] + (randnf() + _Complex_I * randnf()) * 0.01;
    }

    // Compute PSD of the noisy signal
    double* PSD = compute_power_spectrum(bpsk_signal, num_samples);
    printf("Noise level is given by: %f\n", log10f(nstd));

    // Write PRN code and baseband signal to output.txt
    FILE *output = fopen("output.txt", "w");
    for (int i = 0; i < num_samples; i++) {
        fprintf(output, "%d %f\n", i, crealf(bpsk_signal[i]));
    }
    fclose(output);

    // Set up gnuplot commands
    FILE *pipe_gp = popen("gnuplot -p", "w");
        // Set up gnuplot commands
    fprintf(pipe_gp,"set xlabel 'Frequency [MHz]'\n");
    fprintf(pipe_gp,"set ylabel 'C/A Power Spectrum [dB]'\n");
    fprintf(pipe_gp, "set grid linewidth 1\n");
    fprintf(pipe_gp, "set border linewidth 2\n");
    fprintf(pipe_gp, "set tics font 'Arial,10'\n");
    fprintf(pipe_gp, "set key font 'Arial,10'\n");
    fprintf(pipe_gp,"plot 'output.txt' using 1:2 with lines lc rgb 'blue' title 'Baseband Signal'\n");
    fflush(pipe_gp);
    pclose(pipe_gp);

    // Free allocated memory
    free(carrier_signal);
    free(bpsk_signal);
    free(PSD);
    free(nav_data);
    free(upsampled_code);
    
    return 0;
}
