#include <stdio.h>
#include <fftw3.h>
#include <time.h>
#include <stdlib.h>

int main(int argc, char* argv[]){
	int N = 100;
	fftw_complex *in, *out, *back_out;
	fftw_plan my_plan;

	fftw_plan back_plan;
	unsigned int seed = 123456789;

	in = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*N);
	out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*N);
	back_out =(fftw_complex*)fftw_malloc(sizeof(fftw_complex)*N);
	srand(seed);
	for(int i =0; i < N; i++){
		in[i][0] = (double)rand() / (double) RAND_MAX;
		in[i][1] = (double)rand() / (double) RAND_MAX;
	}
	printf("\n");
	printf("Input FFT coefficients:\n");
	printf("\n");
	for(int i = 0; i < N; i++){
		printf("%3d %12f %12f\n", i , in[i][0], in[i][1]);
	}
	
	my_plan = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
	fftw_execute(my_plan);
	printf("\n");
	printf("Output of FFT : \n");
	for(int i = 0; i < N; i++){
		printf("%3d %12f %12f\n", i, out[i][0], out[i][1]);
	}

	printf("Backward FFT: \n");
	back_plan = fftw_plan_dft_1d(N, out, back_out, FFTW_BACKWARD, FFTW_ESTIMATE);
	fftw_execute(back_plan);

	for(int i = 0; i < N; i ++){
		printf("%3d %12f %12f\n", i, back_out[i][0], back_out[i][1]);
	}

	fftw_destroy_plan(my_plan);
	fftw_destroy_plan(back_plan);
	fftw_free(in);
	fftw_free(out);
	fftw_free(back_out);
	return 0;
}
