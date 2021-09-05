#include <math.h>
#include <iostream>
#include <stdlib.h>
#include <ctime>

using namespace std;

double f(const double& x, const double& y) {
	
	double numer = exp(-2 * x - 3 * y);
	double denom = sqrt(x * x + y * y + 1);
	double retval = numer / denom;

	return retval;
}

int main() {

	const double area = atan(1); //atan(1) = arctan(1) = pi/4 = the area of the sector

	int n = 0;
	double sum = 0, x = 0, y = 0;

	srand((unsigned)time(0));

	cout << "Please enter the number of samples:" << endl;

	cin >> n;

	for (int i = 0; i < n; ++i) {
		while (true) {
			//Generate random points between 0 and 1
			x = static_cast<double>(rand()) / RAND_MAX;
			y = static_cast<double>(rand()) / RAND_MAX;

			//Check if the points are suitable for the denom
			if (x * x + y * y <= 1) {
				break;
			}
		}

		//update our sum with the given points
		sum += f(x, y);
	}

	//Integral = area * mean value of f simulated
	sum = area * sum / n;

	cout << "The integral is " << sum << endl;

	return 0;

	system("PAUSE");

}
