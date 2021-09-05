//This is 1D integral calculation

#include<iostream>
#include<stdlib.h>
#include<math.h>

double f(double x) {

	return pow(x, 4) * exp(-x);

}

double MonteCarloEstimate(double lowerbound, double upperbound, int iterations) {

	double sum = 0;
	double randnum, funcval;

	int iter = 0;

	while (iter < iterations) {

		//Generate a random number within the boundary of the integration
		randnum = lowerbound + (float(rand()) / RAND_MAX)*(upperbound - lowerbound);

		//Value of the function at random points generated
		funcval = f(randnum);

		sum += funcval; //For future averaging

		iter++;

	}

	double estimate = (upperbound - lowerbound) * sum / iterations;

	return estimate;
}

int main() {

	double lb(1), ub(5);
	int iteration(2000);

	double esti = MonteCarloEstimate(lb, ub, iteration);

	printf("Estimate for %.1f -> %.1f is %.2f, (%i iterations)\n",
		lb, ub, esti, iteration);

	return 0;

	system("PAUSE");

}
