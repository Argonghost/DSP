#include <stdio.h>
#include <complex.h>
#include <math.h>
#include <time.h>
#include <liquid/liquid.h>
#include "GoldCodelib.h"
#include <complex.h>
#include "STAP_GNSS.h"
#include <lapacke.h>
#include <string.h>

#define GPS_L1_FREQ 1575.42e6
#define CODE_FREQ 1.023e6
#define NUM_CODE_CHIPS 1023
#define INTEGRATION_TIME 0.001 // 1 ms
#define PI 3.141592654

float complex *steering_vector(float *theta, int num_angles, int M, float spacing) {
    float complex *arr_sum = malloc(sizeof(float complex) * num_angles*M);
    for (int i = 0; i < num_angles; i++) {
        float complex sum = 0.0 + 0.0 * I;

        for (int j = 0; j < M; j++) {
            
            sum += cexp(I * PI * sin(theta[i] + 10* PI / 180) * j * spacing);
        }
        arr_sum[i] = sum;
    }
    return arr_sum;
}

float complex* generate_noise_matrix(int M /*Num antennas*/, int N /*Num samples*/) {
    time_t t;
    srand((unsigned) time(&t)); // Seed the random number generator
    float complex* W = malloc(sizeof(float complex) * M * N);

    for (int i = 0; i < M; i++) {
        for (int j = 0; j < N; j++) {
            float real_part = randnf();
            float imag_part = randnf();
            W[i*N + j] = (real_part + I * imag_part) * powf(10, 40 / 20);
        }
    }
    return W;
}

float complex *received_signal(float complex* W,  int M, int N) {
    int num_jam = 4;
    int num_SOI = 2;
    float complex* A = malloc(sizeof(float complex) * M * num_jam);
    float complex* final_signal = malloc(sizeof(float complex) * M * N );
    float angle[4] = {-15 *PI / 180, 0, 15*PI/180, 45*PI / 180};
    float phi[num_jam];

    for(int i = 0; i < num_jam; i++){
        phi[i] = PI * sin(angle[i]);
    };

    /* Define CW jam sources*/
    float complex* J1 = malloc(N * sizeof(float complex));
    float complex* J2 = malloc(N * sizeof(float complex));
    float complex* J3 = malloc(N * sizeof(float complex));
    float complex* J4 = malloc(N * sizeof(float complex));
    for (int i = 0; i < N; i++) {
        J1[i] =cexp(I * 2 * PI * 2.07e6 * i) *powf(10.0f, (80)/20.0f); 
        J2[i] =cexp(I * 2 * PI * 2.08e6 * i)*powf(10.0f, (80)/20.0f) ;
        J3[i] =cexp(I * 2 * PI * 2.09e6 * i)*powf(10.0f, (80)/20.0f) ;
        J4[i] =cexp(I * 2 * PI * 2.1e6 * i)*powf(10.0f, (80)/20.0f) ;
    };


    // Create Jam matrix

    float complex *J_tot = malloc(N * num_jam * sizeof(float complex));
    for (int j = 0; j < N; j++) {
        J_tot[0 * N + j] = J1[j];
        J_tot[1 * N + j] = J2[j];
        J_tot[2 * N + j] = J3[j];
        J_tot[3 * N + j] = J4[j];
    };

    for(int i = 0; i < M; i++){
        float complex a_1 = 0.0 + I*0.0;
        float complex a_2 = 0.0 + I*0.0;
        float complex a_3 = 0.0 + I*0.0;
        float complex a_4 = 0.0 + I*0.0;

        a_1+= cexp(I * phi[0] * i);
        a_2+= cexp(I * phi[1] * i);
        a_3+= cexp(I * phi[2] * i);
        a_4+= cexp(I * phi[3] * i);

        A[i * num_jam + 0] = a_1;
        A[i * num_jam + 1] = a_2;
        A[i * num_jam + 2] = a_3; 
        A[i * num_jam + 3] = a_4;       
    }

    matrixcf_trans(A, num_jam, M);
    float complex S_out[M*N];
    float complex * X = malloc(sizeof(float complex) * M * N);
    matrixcf_mul(A, M, num_jam, J_tot, num_jam, N, S_out, M, N); // Direction vector of Jam * Jam signals
    // matrixcf_mul(A_SOI, M, num_SOI, SOI, num_SOI, N, SOI_angle, M, N); // Direction vector of SOI * SOI 
    matrixcf_add(S_out, W, X, M, N);
    // matrixcf_add(SOI_angle, X, final_signal, M, N);
    free(A);
    free(J1);
    free(J2);
    free(J3);
    free(J4);
    free(J_tot);
    return X;
}

float complex* cov_matrix(float complex *X, int M, int N) {
    float complex* R = malloc(sizeof(float complex) * M * M);
    matrixcf_mul_hermitian(X, M, N, R);

    return R;
}

float complex* weight(float complex *R, int M, int N) {
    int num_SOI = 2;
    float complex g[] = {1.0 + I * 0.0, 1.0 + I * 0.0}; // Changed to have two constraints
    
    // Invert the covariance matrix R
    float complex *R_inv = malloc(sizeof(float complex) * M * M);
    memcpy(R_inv, R, sizeof(float complex) * M * M);
    matrixcf_inv(R_inv, M, M);


    // Create steering vector matrix (constraint matrix)
    float complex* C = malloc(sizeof(float complex) * M * num_SOI);
    float angle2[2] = {10*PI / 180, 65*PI / 180};
    float phi2[num_SOI];

    for (int i = 0; i < num_SOI; i++) {
        phi2[i] = PI * sin(angle2[i]);
    }

    for (int i = 0; i < M; i++) {
        float complex a_1 = 0.0 + I*0.0;
        float complex a_2 = 0.0 + I*0.0;

        a_1+= cexp(I * phi2[0] * i);
        a_2+= cexp(I * phi2[1] * i);

        C[i * num_SOI + 0] = a_1;
        C[i * num_SOI + 1] = a_2;
    }

    // Calculate C^H * R_inv * C
    float complex* C_H = malloc(sizeof(float complex) * num_SOI * M);
    memcpy(C_H, C, sizeof(float complex) * M * num_SOI);
    matrixcf_hermitian(C_H, M, num_SOI);

    float complex* R_inv_C = malloc(sizeof(float complex) * M * num_SOI);
    matrixcf_mul(R_inv, M, M, C, M, num_SOI, R_inv_C, M, num_SOI);

    float complex* C_H_R_inv_C = malloc(sizeof(float complex) * num_SOI * num_SOI);
    matrixcf_mul(C_H, num_SOI, M, R_inv_C, M, num_SOI, C_H_R_inv_C, num_SOI, num_SOI);

    // Invert C^H * R_inv * C
    matrixcf_inv(C_H_R_inv_C, num_SOI, num_SOI);

    // Calculate final weight vector
    float complex* temp = malloc(sizeof(float complex) * M * num_SOI);
    matrixcf_mul(R_inv_C, M, num_SOI, C_H_R_inv_C, num_SOI, num_SOI, temp, M, num_SOI);

    float complex* w = malloc(sizeof(float complex) * M * 1);
    matrixcf_mul(temp, M, num_SOI, g, num_SOI, 1, w, M, 1);

    // Free allocated memory
    free(R_inv);
    free(C);
    free(C_H);
    free(R_inv_C);
    free(C_H_R_inv_C);
    free(temp);
    return w;
}


float complex *arr_new(float *theta, int num_angles, int M, float spacing, float complex *w) {
    float complex *arr_sum = malloc(sizeof(float complex) * num_angles *M);
    for (int i = 0; i < num_angles; i++) {
        float complex sum = 0.0 + 0.0 * I;

        for (int j = 0; j < M; j++) {
            sum += conjf(w[j]) * cexp(I * PI * sin(theta[i]) * j * spacing);
        }
        arr_sum[i] = sum ;
    }
    return arr_sum;
}


int main(){

    int sat_1[1023];
    PRN(1, sat_1);
    printf("Anti Jam result for satellite 1:\n");
    double carrier_freq = 2.0e6; // GPS L1 carrier frequency (Hz)
    double code_duration = 1e-3; // 1 ms
    double fs = 8 * 1.023e6;
    int num_samples = fs * code_duration;
    int chip_rate = 1.023e6;
    int samples_per_chip = fs / chip_rate;


    float SNRdB = -158.5f; // signal-to-noise ratio [dB]
    float noise_floor = -20.0f; // Noise floor 
    float nstd = powf(10.0f, (SNRdB - noise_floor)/20.0f); 
    // Generate the carrier signal
    double complex* carrier_signal = generate_carrier(carrier_freq, code_duration, num_samples);

    int *upsampled_code;
    int upsampled_code_len;
    repeat_1d(sat_1, 1023, samples_per_chip, &upsampled_code, &upsampled_code_len);


    // Create sequence for satellite 12
    int sat_12[1023];
    PRN(12, sat_1);
    int *upsampled_code12;
    int upsampled_code_len12;
    repeat_1d(sat_12, 1023, samples_per_chip, &upsampled_code12, &upsampled_code_len12);

    // Ensure the upsampled code length matches num_samples
    if (upsampled_code_len < num_samples) {
        upsampled_code = realloc(upsampled_code, num_samples * sizeof(int));
        for (int i = upsampled_code_len; i < num_samples; i++) {
            upsampled_code[i] != upsampled_code[i % upsampled_code_len];
        }
    }

    // Do the same for satellite 12 .... Find a way to make this for efficient
    if (upsampled_code_len12 < num_samples) {
        upsampled_code12 = realloc(upsampled_code12, num_samples * sizeof(int));
        for (int i = upsampled_code_len12; i < num_samples; i++) {
            upsampled_code12[i] != upsampled_code12[i % upsampled_code_len12];
        }
    }

    // Generate fake nav data
    int* nav_data = generate_nav_data(50);

    // BPSK modulation
    float complex* bpsk_signal = malloc(num_samples * sizeof(double complex));
    for (int i = 0; i < num_samples; i++) {
        int nav_index = (i * 50 / num_samples);
        int data_bit = upsampled_code[i]^nav_data[nav_index]; // Modulo-2 addition (XOR)

        bpsk_signal[i] = nstd*(2 * data_bit - 1)* carrier_signal[i];
    }

    float complex* bpsk_signal12 = malloc(num_samples * sizeof(double complex));
    for (int i = 0; i < num_samples; i++) {
        int nav_index = (i * 50 / num_samples);
        int data_bit = upsampled_code12[i]^nav_data[nav_index]; // Modulo-2 addition (XOR)

        bpsk_signal12[i] = nstd*(2 * data_bit - 1)* carrier_signal[i];
    }

    // Now we construct array manifold with signal starting with steering vector
    int num_antennas = 8;
    int num_sources = 2;
    float spacing = 1.0f;
    int num_angles = 180;
    float theta[180];

    for (int i = 0; i < num_angles; i++) {
        theta[i] = (i - 90) * PI / 180.0;
    }

    float complex* S = malloc(sizeof(float complex) * 2 * num_samples);
    for (int i = 0; i < num_samples; i++) {
        float real_part1 = crealf(bpsk_signal[i]);
        float imag_part1 = cimagf(bpsk_signal[i]);
        float real_part2 = crealf(bpsk_signal12[i]);
        float imag_part2 = cimagf(bpsk_signal12[i]);
        
        S[0 * num_samples + i] = real_part1 + I * imag_part1;
        S[1 * num_samples + i] = real_part2 + I * imag_part2;
    }

    // Now we generate noise matrix
    float complex* W = generate_noise_matrix(num_antennas, num_samples);

    // Now we create received signal --> A_SOI(theta) * SOI(n) + A_J(n)*J(n) + N(n)
    float complex *X = received_signal(W, num_antennas,num_samples); 

    // Now we construct direction steering vector
    float complex *steer =  steering_vector(theta, num_angles, num_antennas, spacing);
    time_t start, stop;
    clock_t ticks;
    time(&start);
    float complex* R = cov_matrix(X,num_antennas, num_samples);

    // Now we do beamforming
    float complex* w = weight(R, num_antennas, num_samples); // Calling weight matrix
    matrixcf_trans(w, num_antennas, 1);
    float complex *arr_out = malloc(sizeof(float complex) *num_samples * 1);
    matrixcf_mul(w, 1, num_antennas, X, num_antennas, num_samples, arr_out, 1, num_samples);


    float complex* arr_w = arr_new(theta, num_angles, num_antennas, spacing, w);
    /*Plotting utils*/

    FILE *pipe_gp = popen("gnuplot -p", "w");
    FILE *output = fopen("output.txt", "w");
    double time_step = code_duration / num_samples;
    // for (int i = 0; i < num_samples; i++) {
    //     fprintf(output, "%d %lf\n",i, crealf(bpsk_signal[i]));
    // }
    for(int i = 0; i < num_angles; i++){
        fprintf(output,"%f %f %f\n", theta[i], 20.0 * log10f(cabsf(arr_w[i])), 20.0 * log10f(cabsf(steer[i])));
    }
 
    fclose(output);

    printf("Now we try Singular value decomposition of estimated Rx: \n");

    #define M num_antennas  // Number of rows
    #define N num_samples  // Number of columns
    #define MIN(a,b) (((a)<(b))?(a):(b))

    lapack_complex_float *U, *VT;
    float *Sigma;
    lapack_int info, lda, ldu, ldvt;
    float superb[MIN(M, N) - 1];
    Sigma = (float*)LAPACKE_malloc(MIN(M,N) * sizeof(float));
    U = (lapack_complex_float*)LAPACKE_malloc(M * M * sizeof(lapack_complex_float));
    VT = (lapack_complex_float*)LAPACKE_malloc(N * N * sizeof(lapack_complex_float));

    lda = M; // Changed from M to N
    ldu = M;
    ldvt = M;

    printf("Computing SVD...\n");

    // Enter here the SVD computation
    info = LAPACKE_cgesvd(LAPACK_ROW_MAJOR, 'A', 'A', M, M, R, lda, Sigma, U, ldu, VT, ldvt, superb);

    ticks = clock();
    printf("Performed SVD and used %0.2f seconds of CPU time.\n",(double)ticks/CLOCKS_PER_SEC);
    time(&stop);
    printf("Finished in about %.6f seconds.\n", difftime(stop, start));
    printf("SVD computation completed\n");

    printf("U matrix values:\n");
    matrixf_print(Sigma, MIN(M,N), 1);
    // Set up gnuplot commands
    fprintf(pipe_gp,"set xlabel 'Angle [degrees]'\n");
    fprintf(pipe_gp,"set ylabel 'BF Response [dB]'\n");
    fprintf(pipe_gp, "set grid linewidth 1\n");
    fprintf(pipe_gp, "set yrange [-50:40]\n");
    fprintf(pipe_gp, "set border linewidth 2\n");
    fprintf(pipe_gp, "set tics font 'Arial,10'\n");
    fprintf(pipe_gp, "set key font 'Arial,10'\n");
    fprintf(pipe_gp,"plot 'output.txt' using 1:2 with lines lc rgb 'blue', '' using 1:3 with lines lc rgb 'red'\n");

    // Close the gnuplot pipe
    fflush(pipe_gp);
    pclose(pipe_gp);

    free(carrier_signal);
    free(bpsk_signal12);
    free(bpsk_signal);
    free(nav_data);
    free(upsampled_code);
    free(upsampled_code12);
    free(W);
    free(X);
    free(S);
    free(steer);
    free(arr_w);
    free(w);
    free(R);
    LAPACKE_free(U);
    LAPACKE_free(VT);
    LAPACKE_free(Sigma);

    printf("Memory freed\n");
    return 0;
}
