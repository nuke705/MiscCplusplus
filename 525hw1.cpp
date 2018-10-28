#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <fstream>
#include <cstdlib>
using namespace std;
// if you want to run the code, change the comment to use different methods 
float u, p, p_u, p_d, r, K;
float d, downtick_prob;
float S0, T, volatility,target = 9.472992 ;
int N;
float max(float a, float b) {
	return (b < a) ? a : b;
}
double Ncdf(const double& z) {
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

double option_price_put_black_scholes(const double& S,      // spot price
	const double& K,      // Strike (exercise) price,
	const double& r,      // interest rate
	const double& sigma,  // volatility
	const double& time) {
	double time_sqrt = sqrt(time);
	double d1 = (log(S / K) + r*time) / (sigma*time_sqrt) + 0.5*sigma*time_sqrt;
	double d2 = d1 - (sigma*time_sqrt);
	return K*exp(-r*time)*Ncdf(-d2) - S*Ncdf(-d1);
};

double americanprice (int N){
	
	double S0 = 100.0, K = 105.0, volatility = 0.2, r = 0.03, T = 1.0, deltaT; // N = 25;
	
		deltaT = T / N;
	
	u = exp(volatility*sqrt(deltaT)), d = 1 / u, p = (exp(r*deltaT) - d) / (u - d);
	double discount = exp(-r*deltaT);
	p_u = discount*p, p_d = discount*(1 - p);
	
	vector <vector <double>> option;// This is a N+1 x N+1 matrix
	for (int i = 0; i <= N; i++) {
		vector <double> temp;
		option.push_back(temp);
		for (int j = 0; j <= N; j++) {
			option[i].push_back(0);
		}
	}
	vector <vector <double>> stock;
	for (int i = 0; i <= N ; i++) {
		vector <double> temp;
		stock.push_back(temp);
		for (int j = 0; j <= N ; j++) {
			stock[i].push_back(0);
		}
	}
	
	for (int j = 0; j <= N; j++)
		for (int i = 0; i <= j; i++)
			stock[i][ j] = S0*pow(u, j - i)*pow(d, i);

	// Compute terminal payoffs
	for (int i = 0; i <= N; i++) {
		option[i][ N] = max(K - stock[i][N], 0.0);
	}

	// Backward through the tree
	for (int j = N - 1; j >= 0; j--){
		for (int i = 0; i <= j; i++) {
			//*********BBS for T-1*************
			if (j == N - 1) { double bsprice = option_price_put_black_scholes(stock[i][j], K, r, volatility, deltaT);
			option[i][j] = max(K - stock[i][j], bsprice);
			}else {
			//******************* below are normal binomial***********/
				option[i][j] = max(K - stock[i][j], exp(-r*deltaT)*(p*(option[i][j + 1]) + (1.0 - p)*(option[i + 1][j + 1])));
			}
		}
	}
	//ofstream outfile("BBSR");
	//double thedelta = (option[0][1] - option[1][1]) / (S0*u - S0 / u);
	//cout << thedelta << endl;
	//cout << option[0][0] << endl;
	return option[0][0];
}
int main(int argc, char* argv[]) {
	/*ofstream outfile("binomialBBS");
	for (int division = 25; division <= 2500; division+=50) {
		outfile <<americanprice(division) << endl;
	}*/
	ofstream outfilebbsr("BBSR");
	for (int division = 25; division <= 2500; division += 50) {
		//double errorBBSR =
		//abs(target - (2 * americanprice(division) - americanprice(division/2.0)) )  ;
		//cout << errorBBSR << endl;
		//outfile << errorBBSR << endl;
		double bbsrp = (2 * americanprice(division) - americanprice(division / 2.0));
		outfilebbsr << bbsrp << endl;
	}
}
/*
How does BBSR work?
the BBSR method computes the BBS price corresponding to a pair of options
with n and n/2 time steps (deltaT),Pn and Pn/2  respectively, and then sets the approximate price to P = 2*Pn-Pn/2 
*/
