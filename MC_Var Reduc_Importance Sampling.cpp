/*
1. Monte Carlo estimator of X is [ sum(f(x)) / n ] and therefore we have:

	i) the estimator is consistent
	ii) variance of the estimator is (V^2) / (n-1) * { E[ f(x)^2 ] - ( E[ f(x) ]^2 }

2. Methods of Variance Reductions:

	i) Importance Sampling
	   This is a method that rejects samples in less important regions, e.g., regions of x where
	   f(x) is almost zeros (0).

	   Aim: to adjust the sampled distribution from a uniform distribution to a distribution with
			PDF that somewhat resembles the function of interest, then we can more quickly arrive
			at a statistically reliable integration estimation and reduce the sampling variance.
	   
	   Difficulty: choosing a companion function to sample from isn't always straightforward;
			and it sometimes requires the prior knowledge of the function being integrated.
	   
	   Algo: to estimate E[ f(X); p], the MC estimator is (1/n) * sum( p(X) * f(X) / q(X) ), where
			p(X) > 0 ==> q(X) > 0 for all x in Omega, i.e., we need to find such a PDF q.
	   
*/

/*
From the methods we could tell that Importance Sampling is suited for functions with a single peak globally
and locally, since cases where the function f has more than one (1) peak works better with Strtified method.

*/

/*
Therefore in this script, we compare standard MC with IS MC via the function f(x) = 10 * exp(-5*(x-3)^4)
over the interval [0,6]. It was then noticed that f(x) only significantly differs from zero (0) for x in [2,4].
We would like to have a sampling PDF that gives more likelihood to outcomes in [2,4] than a plain Unif[0,6].

N(3,1) is an ideal choice.

*/

#include<iostream>
#include<stdlib.h>
#include<math.h>
#include<random>

double f(double x);
double PDF(double x);
void MC(double lowerbound, double upperbound, int iterations, double mcstats[]);
void MCIS(double lowerbound, double upperbound, int iterations, double mcstats[]);


int main() {

	double lb(1), ub(5);
	int iterations;
	int counter(5);

	double mc_stats[2];
	double mcis_stats[2];

	//print table header.
	printf("MC Normal Sampling\t\tMC Importance Sampling\t\tIterations\n");
	printf("Estimate\tSTD\t\tEstimate\tSTD\n");

	for (int i = 0; i < counter; ++i) {
		iterations = 5 * pow(10, i);

		MC(lb, ub, iterations, mc_stats);
		MCIS(lb, ub, iterations, mcis_stats);

		printf("%.3f\t\t%.3f\t\t%.3f\t\t%.3f\t\t%i\n",
			mc_stats[0], mc_stats[1], mcis_stats[0], mcis_stats[1], iterations);

	}

	
	return 0;

	system("PAUSE");

}


double f(double x) {

	return 10 * exp(-5 * pow(x - 3, 4));

}

double PDF(double x) {

	return (1 / pow(2 * 3.14159, 0.5)) * exp(-(0.5) * pow(x - 3, 2));

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

void MCIS(double lowerbound, double upperbound, int iterations, double statsArray[]) {

	//Random number generator for samples from the companion distribution N(3,1)
	std::default_random_engine generator;
	std::normal_distribution<double> distribution(3, 1.0);

	double sum = 0, sumsquared = 0;

	int iter = 0;

	double randnum, funcval, weight;

	while (iter < iterations) {

		randnum = distribution(generator);

		weight = (1 / (upperbound - lowerbound)) / PDF(randnum);

		funcval = f(randnum) * weight;

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
