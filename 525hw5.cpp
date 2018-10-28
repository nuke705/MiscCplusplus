
#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <cstdlib>
#include <algorithm>
#include <math.h>
#include <vector>
#include <random>
#include <numeric>
using namespace std;
double get_uniform() {return ((double)rand() / RAND_MAX);}

float max(float a, float b) {
	return (b < a) ? a : b;
}

double N(const double& z) {
	if (z > 6.0) { return 1.0; }; // this guards against overflow 
	if (z < -6.0) { return 0.0; };
	double b1 = 0.31938153;
	double b2 = -0.356563782;
	double b3 = 1.781477937;
	double b4 = -1.821255978;
	double b5 = 1.330274429;
	double p = 0.2316419;
	double c2 = 0.3989423;
	double a = fabs(z);
	double t = 1.0 / (1.0 + a*p);
	double b = c2*exp((-z)*(z / 2.0));
	double n = ((((b5*t + b4)*t + b3)*t + b2)*t + b1)*t;
	n = 1.0 - b*n;
	if (z < 0.0) n = 1.0 - n;
	return n;
};

double betapdf(double x, double a, double b) {
	
	double B = tgamma(a)*tgamma(b) / tgamma(a + b);
	return 1 / B * pow(x, a - 1)*pow(1 - x, b - 1);
}

double gammarv(double a, double b) {
	double u1 = get_uniform();
	double u2 = get_uniform();
	double var = (a - 1) / (a + b - 2);
	double c;
	if (a < 1 && b < 1) {
		c = 0.99999;
	}else {
		c = betapdf(var, a, b);
	}
	/*if (a < 1 && b < 1) {
		if (u2 <= pow(4,a)*pow(u1*(1-u1),a)) {
			return u1;
		}else {
			return -10.0;
		}
	}else {*/
		
		if (c*u2 <= betapdf(u1, a, b)) {
			return u1;
		}
		else {
			return -10.0;
		}
	//}
}

vector <double> trial(double a, double b, int n) {
	vector <double> temp;
	for (int i = 0; i < n; i++) {
		//double res = gammarv(a, b);
		temp.push_back(gammarv(a, b));

	}
	
	return temp;
}

vector <double> euler_path() {
	int nt = 1000;
	//t fomr 0 to 10
	double a = 0.1, b = 0.4, r0 = 0.3, sigma = 2, h=10.0/1000.0;
	vector <double> path;
	for (int i = 0; i < nt; i++) {
		if (i == 0) {
			path.push_back(r0);
		}else {
			//default_random_engine generator;
			//normal_distribution<double> dist(0.0, sqrt(h));
			double x = get_uniform();
			double y = get_uniform();
			double z = sqrt(-2.0*log(x)) * cos(6.283185307999998*y);
			//double z = dist(generator);
			
			double prev = path[i - 1] + a*(b - path[i - 1])*h + sigma*sqrt(path[i - 1])*sqrt(h)*z;
			path.push_back(prev);
		}

	}

	return path;
}

vector <double> euler_path2() {
	int nt = 1000;
	//t fomr 0 to 10
	double a = 0.1, b = 0.4, r0 = 0.3, sigma = 2, h = 10.0 / 1000.0;
	vector <double> path;
	for (int i = 0; i < nt; i++) {
		if (i == 0) {
			path.push_back(r0);
		}
		else {
			//default_random_engine generator;
			//normal_distribution<double> dist(0.0, sqrt(h));
			double x = get_uniform();
			double y = get_uniform();
			double z = sqrt(-2.0*log(x)) * cos(6.283185307999998*y);
			//double z = dist(generator);
			double rt = max(path[i - 1], 0.0);
			double prev = path[i - 1] + a*(b - rt)*h + sigma*sqrt(rt)*sqrt(h)*z;
			path.push_back(prev);
		}

	}

	return path;
}

double variance(vector <double> price) {
	double mean = 0;
	double variance = 0;
	for (int i = 0; i < 10000; i++) {
		mean = mean + price[i];
	}
	mean = (double)mean / 10000;
	for (int i = 0; i < 10000; i++) {
		variance = variance + pow(price[i] - mean, 2);
	}
	variance = variance / (10000 * (10000 - 1));
	return variance;
}


void Q4() {
	int nt = 1000;
	//t fomr 0 to 10
	double mu = 0.1, s0 = 15, r = 0.0,K=16, sigma = 0.15, h = 1.0/24.0;
	double call_option_price = 0, put_option_price =0;
	double asian_call_price = 0;
	double expiration_time = 1;
	vector <double> calls; vector <double> puts; vector <double> asians;
	vector <double> path;
	for (int tr = 0; tr < 10000; tr++) {
		double S = 15; double sum = 15;
		for (int i = 0; i < 24; i++) {
			//if (i == 0) {
			//	path.push_back(s0);
			//}else {
				//default_random_engine generator;
				//normal_distribution<double> dist(0.0, sqrt(h));

				double x = get_uniform();
				double y = get_uniform();
				double z = sqrt(-2.0*log(x)) * cos(6.283185307999998*y);
				while (isfinite(z) != 1) {
					x = get_uniform();
					z = sqrt(-2.0*log(x)) * cos(6.283185307999998*y);
				}
				//double z = dist(generator);
				double R = (mu - 0.5*pow(sigma, 2))*h;
				double SD = sigma*sqrt(h);
				//path[i] = path[i - 1] * exp(R + SD*z);
				//sum += path[i];
				S = S* exp(R + SD*z); sum += S;
				//cout << sum << endl;
				//path.push_back(S1);
				//for (int )
			//}
			
		}
		double callthis = max(0.0, S - K), putthis= max(0.0, K - S);
		double asiancall = max(0.0, sum/25 - K);
		calls.push_back(callthis); puts.push_back(putthis);
		asians.push_back(asiancall);
		call_option_price += callthis ;
		put_option_price += putthis;
		asian_call_price += asiancall;
	}	
	double callprice = exp(-r)*call_option_price / 10000;
	double putprice = exp(-r)*put_option_price / 10000;
	double asianprice = exp(-r)*asian_call_price / 10000;
	double varcall = variance(calls);
	double varput = variance(puts); double varasian = variance(asians);

	cout << "call " <<callprice << endl; cout << "put  " << putprice << endl;
	cout << "asian  " << asianprice << endl;
	cout << "----------------" << endl;
	cout << "var call " << varcall << endl; cout << "var put " << varput << endl; 
	cout << "var as " << varasian << endl;
}

double option_price_put_black_scholes(const double& S,      // spot price
	const double& K,      // Strike (exercise) price,
	const double& r,      // interest rate
	const double& sigma,  // volatility
	const double& time) {
	double time_sqrt = sqrt(time);
	double d1 = (log(S / K) + r*time) / (sigma*time_sqrt) + 0.5*sigma*time_sqrt;
	double d2 = d1 - (sigma*time_sqrt);
	return K*exp(-r*time)*N(-d2) - S*N(-d1);
};

double option_price_call_black_scholes(const double& S,       // spot (underlying) price
	const double& K,       // strike (exercise) price,
	const double& r,       // interest rate
	const double& sigma,   // volatility 
	const double& time) {  // time to maturity 
	double time_sqrt = sqrt(time);
	double d1 = (log(S / K) + r*time) / (sigma*time_sqrt) + 0.5*sigma*time_sqrt;
	double d2 = d1 - (sigma*time_sqrt);
	return S*N(d1) - K*exp(-r*time)*N(d2);
};


int main(int argc, const char * argv[]) {
	//Q2
	/*ofstream out1("0505");
	ofstream out2("51");
	ofstream out3("25");
	ofstream out4("22");
	int nt = 10000;
	vector <double> outone = trial(0.5, 0.5, nt);
	for (int i = 0; i < nt; i++) {
		if (outone[i] > -10) {
			out1 << outone[i] << endl;
		}
	}

	vector <double> outtwo = trial(5, 1, nt);
	for (int i = 0; i < nt; i++) {
		if (outtwo[i] > -10) {
			out2 << outtwo[i] << endl;
		}
	}

	vector <double> outthree = trial(2, 5,nt); 
	for (int i = 0; i < nt; i++) {
		if (outthree[i] > -10) {
			out3 << outthree[i] << endl;
		}
	}

	vector <double> outfour = trial(2, 2, nt); 
	for (int i = 0; i < nt; i++) {
		if (outfour[i] > -10) {
			out4 << outfour[i] << endl;
		}
	}*/
	//Q3
	/*ofstream out3a1("3a1");
	vector <double> euler1 = euler_path();
	for (int i = 0; i < 1000; i++) {
			out3a1 << euler1[i] << endl;
	}
	ofstream out3a2("3a2");
	vector <double> euler2 = euler_path();
	for (int i = 0; i < 1000; i++) {
		out3a2 << euler2[i] << endl;
	}
	ofstream out3a3("3a3");
	vector <double> euler3 = euler_path();
	for (int i = 0; i < 1000; i++) {
		out3a3 << euler3[i] << endl;
	}
	ofstream out3a4("3a4");
	vector <double> euler4 = euler_path();
	for (int i = 0; i < 1000; i++) {
		out3a4 << euler4[i] << endl;
	}
	ofstream out3b1("3b1");
	euler1 = euler_path2();
	for (int i = 0; i < 1000; i++) {
		out3b1 << euler1[i] << endl;
	}
	ofstream out3b2("3b2");
	euler2 = euler_path2();
	for (int i = 0; i < 1000; i++) {
		out3b2 << euler2[i] << endl;
	}
	ofstream out3b3("3b3");
	euler3 = euler_path2();
	for (int i = 0; i < 1000; i++) {
		out3b3 << euler3[i] << endl;
	}
	ofstream out3b4("3b4");
	euler4 = euler_path2();
	for (int i = 0; i < 1000; i++) {
		out3b4 << euler4[i] << endl;
	}*/
	//Q4
	Q4();
	double theocall = option_price_call_black_scholes(15.0,16, 0,0.15, 1);
	double theoput = option_price_put_black_scholes(15.0,16, 0,0.15, 1);
	cout << endl;
	cout << "real call (or not) " <<theocall << endl;
	cout << "real put (or not) " <<theoput << endl;
}
