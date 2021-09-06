/*
1. Monte Carlo estimator of X is [ sum(f(x)) / n ] and therefore we have:

	i) the estimator is consistent
	ii) variance of the estimator is (V^2) / (n-1) * { E[ f(x)^2 ] - ( E[ f(x) ]^2 }

2. Methods of Variance Reductions:

	i) Importance Sampling
	   See the other file for more details.

	ii) Stratified Sampling
		Stratified sampling is another approach to variance reduction by which the integration volume
		is divided into subdomains that are each evaluated separately. The estimates from each subdomain
		are then combined with a weight depending on their subdomain integration volume.

		The stratified sampling estimate of an integral over a function g(x) is given as:
			E[ g(x) ] = sum( V_j / n_j * sum( g(x_ij) ) ) where x_ij lies in V_j.

		The variance is therefore straght forward from the above expression.

*/

/*
Importance Sampling is suited for functions with a single peak globally and locally.
Stratified Sampling is suited for function with more than one (1) peak.

*/

/*
In this script we compare MC and MCSS on the function f(x) = exp(-(x-6)^4) + exp(-(x-14)^4), which has 2
peaks at 6 and 14.

*/

#include<iostream>
#include<stdlib.h>
#include<math.h>
#include<random>

double f(double x);
void MC(double lowerbound, double upperbound, int iterations, double mcstats[]);
void MCSS(double lowerbound, double upperbound, int iterations, double mcstats[],int subdomains);


int main() {

	double lb(0), ub(20);
	int iterations;
	
	double mc_stats[2];
	
	printf("Normal Monte Carlo Integration\n");
	for (int i = 1; i < 6; i++) {
		
		iterations = int(128 * pow(4, i));

		MC(lb, ub, iterations, mc_stats);

		printf("Estimate for %.1f -> %.1f is %.3f, std = %.4f, (%i iterations)\n",
			lb, ub, mc_stats[0], mc_stats[1], iterations);

	}

	printf("\nStratified Sampling Monte Carlo Integration\n");
	for (int i = 1; i < 6; i++) {

		iterations = int(128 * pow(4, i));

		MCSS(lb, ub, iterations, mc_stats, 4);

		printf("Estimate for %.1f -> %.1f is %.3f, STD = %.4f, (%i iterations)\n",
			lb, ub, mc_stats[0], mc_stats[1], iterations);

	}

	return 0;

	system("PAUSE");

}


double f(double x) {

	return exp(-1 * pow(x - 6, 4)) + +exp(-1 * pow(x - 14, 4));

}


void MC(double lowerbound, double upperbound, int iterations, double statsArray[]) {

	double sum = 0, sumsquared = 0;
	double randnum, funcval;

	int iter = 0;

	while (iter < iterations) {

		//Generate a random number within the boundary of the integration
		randnum = lowerbound + (float(rand()) / RAND_MAX) * (upperbound - lowerbound);

		//Value of the function at random points generated
		funcval = f(randnum);

		sum += funcval; //For future averaging
		sumsquared += pow(funcval, 2);

		iter++;

	}

	double estimate = (upperbound - lowerbound) * sum / iterations;

	double expected = sum / iterations;
	double expectedsquare = sumsquared / iterations;

	double std = (upperbound - lowerbound) * pow((expectedsquare - pow(expected, 2)) / (iterations - 1), 0.5);

	statsArray[0] = estimate;
	statsArray[1] = std;

}

void MCSS(double lowerbound, double upperbound, int iterations, double statsArray[], int subdomains) {
	
	//function to execute MC integral on predefined function.
	
	double* sum = new double[subdomains]; //WATCH OUT for the use of pointer here!!!
	double* sumsquared = new double[subdomains];

	int iter;
	iterations = int(float(iterations) / subdomains); //divide local iterations among subdomains.

	for (int i = 0; i < subdomains; ++i) {

		sum[i] = 0;
		sumsquared[i] = 0;

	}

	double increment = (upperbound - lowerbound) / float(subdomains);

	for (int seg = 0; seg < subdomains; ++seg) {

		iter = 0;

		double randnum, funcval;
		double rangestart = lowerbound + seg * increment;

		while (iter < iterations) {

			randnum = rangestart + (float(rand()) / RAND_MAX) * increment;
			funcval = f(randnum);
			sum[seg] += funcval;
			sumsquared[seg] += pow(funcval, 2);

			iter++;
		}

	}

	double* estimate = new double[subdomains];
	double* expected = new double[subdomains];
	double* expectedsquare = new double[subdomains];
	double* std = new double[subdomains];

	for (int i = 0; i < subdomains; ++i) {

		estimate[i] = increment * sum[i] / iterations;
		expected[i] = sum[i] / iterations;
		expectedsquare[i] = sumsquared[i] / iterations;

		std[i] = increment * pow((expectedsquare[i] - pow(expected[i], 2)) / (iterations - 1), 0.5);

	}
	
	double estimate_final = 0;
	double std_final = 0;

	for (int i = 0; i < subdomains; i++) {

		estimate_final += estimate[i];

		std_final += pow(increment, 2) * pow(std[i], 2) / iterations;

	}

	statsArray[0] = estimate_final;
	statsArray[1] = pow(std_final, 0.5);

	delete[] sum; sum = NULL;
	delete[] sumsquared; sumsquared = NULL;
	delete[] std; std = NULL;
	delete[] estimate; estimate = NULL;
	delete[] expected; expected = NULL;
	delete[] expectedsquare; expectedsquare = NULL;

}
