#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <cstdlib>
#include <algorithm>
#include <math.h>
#include <vector>
using namespace std;

double sigma_V = 0.3, sigma_U = 0.4, V0 = 50, U0 = 60, T = 5.0 / 12, r = 0.05, corr = 0.7;
int no_of_simulation=500, no_of_time = 48;

vector <vector <double>>  U; vector <vector <double>> V;
vector <double> price;

double max(double a, double b) {
	return (b < a) ? a : b;
}

double get_uniform()
{
	return (((double)rand()) / RAND_MAX);
}

double get_gaussian()
{
	return (sqrt(-2.0*log(get_uniform()))*cos(6.283185307999998*get_uniform()));
}

double N(const double& z) {
	if (z > 6.0) { return 1.0; };
	if (z < -6.0) { return 0.0; };
	double b1 = 0.31938153;
	double b2 = -0.356563782;
	double b3 = 1.781477937;
	double b4 = -1.821255978;
	double b5 = 1.330274429;
	double p = 0.2316419;
	double c2 = 0.3989423;
	double a = fabs(z);
	double t = 1.0 / (1.0 + a * p);
	double b = c2 * exp((-z)*(z / 2.0));
	double n = ((((b5*t + b4)*t + b3)*t + b2)*t + b1)*t;
	n = 1.0 - b * n;
	if (z < 0.0) n = 1.0 - n;
	return n;
};

void initialization() {
	for (int i = 0; i < no_of_simulation; i++) {
		U.push_back(vector <double> ());
		V.push_back(vector <double> ());
		price.push_back(0);
		for (int j = 0; j < no_of_time; j++) {
				U[i].push_back(0);
				V[i].push_back(0);
		}
	}
}

void simulation() {
	double delta_t = T / no_of_time;
	for (int i = 0; i < no_of_simulation; i++) {
		for (int j = 0; j < no_of_time; j++) {
			if (j == 0) {
				V[i][j] = V0; U[i][j] = U0;
			}else {
				double x = get_uniform();
				double y = get_uniform();
				double zv = sqrt(-2.0*log(x)) * cos(6.283185307999998*y);
				double zu = sqrt(-2.0*log(x)) * sin(6.283185307999998*y);
				while (isfinite(zv) != 1 && isfinite(zu) != 1) {
					x = get_uniform(); y = get_uniform();
					zv = sqrt(-2.0*log(x)) * cos(6.283185307999998*y);
					zu = sqrt(-2.0*log(x)) * sin(6.283185307999998*y);
				}
				V[i][j] = V0 * exp((r-0.5*sigma_V*sigma_V)*T + sigma_V*sqrt(T)*zv);
				U[i][j] = U0 * exp((r - 0.5*sigma_U*sigma_U)*T + sigma_U*sqrt(T)*(corr*zv + zu*sqrt(1 - corr*corr)));
				//	+ delta_t * r*(V[i][j - 1]) + sigma_V * V[i][j - 1] * sqrt(delta_t)*zv;
				//U[i][j] = U[i][j - 1] 
				//	+ delta_t * r*(U[i][j - 1]) + sigma_U * U[i][j - 1] * sqrt(delta_t)*(corr*zv + zu*(sqrt(1 - pow(corr, 2))));
			}
		}
	}
}

double get_price() {
	for (int i = 0; i < no_of_simulation; i++) {
		price[i] = exp(-r*T)*max(V[i][no_of_time - 1], U[i][no_of_time - 1]);
		
	}
	//cout << price[no_of_simulation -1] << endl;
	double sum = 0;
	for (int i = 0; i < no_of_simulation; i++) {
		sum = sum + price[i];
		
	}
	return  exp(-r * T)*(sum / no_of_simulation) ;
}

double BS() {
	double sigma_ = sqrt(pow(sigma_V, 2) + pow(sigma_U, 2) - 2 * corr*sigma_U*sigma_V);
	double d1 = (log(V0 / U0) + pow(sigma_, 2)*T*0.5) / (sigma_*sqrt(T));
	double d2 = d1 - sigma_ * sqrt(T);
	return V0*N(d1)  +U0 * N(-d2);
}



int main() {
	srand((unsigned)time(NULL));
	ofstream output_file("q5");
	vector <double> res;
	for (int i = 0; i < 500; i++) {
		no_of_simulation = 50 + 10 * i;
		initialization();
		simulation();
		output_file << get_price() << endl;
		res.push_back(get_price());
		//res[i] = 0;
	}

	
	/*for (int i = 0; i < 500; i++) {
		output_file << 50 + i * 10 << " " << res[i] << endl;
	}*/
	//initialization();
	//simulation();
	//cout << U[no_of_simulation-1][no_of_time-1] << endl;
	double aaa = get_price();
	cout << "monte carlo: " << aaa << endl;
	cout << "black_scholes: " << BS() << endl;

	

}