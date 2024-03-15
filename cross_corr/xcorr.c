#include <stdio.h>
#include <fftw3.h>
#include <complex.h>
#include <time.h>
#include <stdlib.h>




int main(int argc, char* argv[]){
	int N = 5;
	fftw_complex *in1, *in2, *out1, *out2, *out_tot, *xcorr_out;
	fftw_plan my_plan1;

	fftw_plan my_plan2;
    fftw_plan my_plan3;


	in1 = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*N);
	out1 = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*N);

    in2 = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*N);
	out2 = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*N);

	out_tot =(fftw_complex*)fftw_malloc(sizeof(fftw_complex)*N);
    xcorr_out =(fftw_complex*)fftw_malloc(sizeof(fftw_complex)*N);
	srand(time(0));
    
    // Below we define subroutine for in1

	for(int i =0; i < N; i++){
		in1[i][0] = (double)(rand() %100) + 1;
		in1[i][1] = (double)(rand() %100) + 1;
	}

    // Below we define subroutine for in2

	for(int i =0; i < N; i++){
		in2[i][0] = (double)(rand() %100) + 1;
		in2[i][1] = (double)(rand() %100) + 1;
	}

	printf("\n");
	printf("Input FFT coefficients for input 1:\n");
	printf("\n");
	for(int i = 0; i < N; i++){
		printf("%3d %12f %12f\n", i , in1[i][0], in1[i][1]);
	}
    printf("Input FFT coefficients for input 2:\n");
	printf("\n");
	for(int i = 0; i < N; i++){
		printf("%3d %12f %12f\n", i , in2[i][0], in2[i][1]);
	}
	
	my_plan1 = fftw_plan_dft_1d(N, in1, out1, FFTW_FORWARD, FFTW_ESTIMATE);
	fftw_execute(my_plan1);
	printf("\n");
	printf("Output of FFT for input 1: \n");
	for(int i = 0; i < N; i++){
		printf("%3d %12f %12f\n", i, out1[i][0] / N, out1[i][1] / N);
	}

	printf("Output of  FFT for input 2: \n");
	my_plan2 = fftw_plan_dft_1d(N, in2, out2, FFTW_FORWARD, FFTW_ESTIMATE);

	fftw_execute(my_plan2);

	for(int i = 0; i < N; i ++){
		printf("%3d %12f %12f\n", i, out2[i][0], out2[i][1]);
	}

    printf("Output of  frequency domain cross correlation: \n");
    for (int i = 0; i < N ; i++){
        *out_tot[i] = *out1[i] * conj(*out2[i])/ N;
        printf("%3d %12f \n", i,*out_tot[i]);
    }

    my_plan3 = fftw_plan_dft_1d(N, out_tot, xcorr_out, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_execute(my_plan3);

    printf("Output of cross correlation is finally:\n");
    for(int i = 0; i < N; i ++){
		printf("%3d %12f %12f\n", i, xcorr_out[i][0]/N, xcorr_out[i][1]/N);
	}

	fftw_destroy_plan(my_plan1);
	fftw_destroy_plan(my_plan2);
    fftw_destroy_plan(my_plan3);
	fftw_free(in1);
	fftw_free(out1);
    fftw_free(in2);
	fftw_free(out2);
    fftw_free(out_tot);
    fftw_free(xcorr_out);
	return 0;
}
