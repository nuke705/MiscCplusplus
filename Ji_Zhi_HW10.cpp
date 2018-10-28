
#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <cstdlib>
#include "E:\study\hw10am\normdist.h"
#include <vector>
using namespace std;


float up_factor, uptick_prob, risk_free_rate, strike_price;
float downtick_prob, notick_prob;
float initial_stock_price, expiration_time, volatility, R;
int no_of_divisions;
double** memo_european_call; double** memo_european_put; double** memo_american_call; double** memo_american_put;
/*vector <vector <double> > memo_european_call; vector <vector <double> > memo_european_put; 
vector <vector <double> > memo_american_call; vector <vector <double> > memo_american_put;*/

float max(float a, float b) {
	return (b < a) ? a : b;
}




void initialize() {
	int size = 2 * no_of_divisions + 1;
	
	double ** temp = new double *[no_of_divisions + 1];
	for (int i = 0; i < no_of_divisions + 1; i++) {
		
		temp[i] = new double[size];
	}
	
	for (int i = 0; i < no_of_divisions + 1; i++) {
		for (int j = 0; j < size; j++) {
			
			temp[i][j] = -1;
		}
	}
	memo_american_call = temp; 
	memo_american_put = temp;
	memo_european_call = temp; 
	memo_european_put = temp;
	//cout << "reach2" << endl;
}



double european_call_option(int k, int i) {
	if (memo_european_call[k][no_of_divisions + i] >=0) {
		return memo_european_call[k][no_of_divisions + i];
	}
	else if (k == no_of_divisions) { 
		return memo_european_call[no_of_divisions][no_of_divisions + i] = max(0.0, initial_stock_price*pow(up_factor, (double)i) - strike_price);
	}else { 
		return memo_european_call[k][no_of_divisions + i] = (uptick_prob*european_call_option(k + 1, i + 1) +
			notick_prob*european_call_option(k + 1, i) +
			downtick_prob*european_call_option(k + 1, i - 1)) / R;
	}
}

double european_put_option(int k, int i) {
	if (memo_european_put[k][no_of_divisions + i] >= 0) {
		return memo_european_put[k][no_of_divisions + i];
	} else if (k == no_of_divisions) { 
		return memo_european_put[no_of_divisions][no_of_divisions + i] = max(0.0, strike_price - initial_stock_price*pow(up_factor, (double)i));
	} else { 
		return memo_european_put[k][no_of_divisions + i] = (uptick_prob*european_put_option(k + 1, i + 1) +
			notick_prob*european_put_option(k + 1, i) +
		downtick_prob*european_put_option(k + 1, i - 1)) / R;
	} 
}

double american_call_option(int k, int i, float current_stock_price) {
	if (memo_american_call[k][i + no_of_divisions] >= 0){
		return memo_american_call[k][i + no_of_divisions];
	}else if (k == no_of_divisions) {
		return  memo_american_call[no_of_divisions][no_of_divisions + i] = max(0.0, current_stock_price - strike_price);
	}else { 
		return  memo_american_call[k][no_of_divisions + i] = max(max(0.0, current_stock_price - strike_price),
			(uptick_prob*american_call_option(k + 1, i + 1, current_stock_price * up_factor) + 
				notick_prob*american_call_option(k + 1, i, current_stock_price) +
				downtick_prob*american_call_option(k + 1, i - 1, current_stock_price/ up_factor)) / R);
	}
}


double american_put_option(int k, int i, float current_stock_price) {
	if (memo_american_put[k][i + no_of_divisions] >= 0) {
		return memo_american_put[k][i + no_of_divisions];
	} else if (k == no_of_divisions) { 
		return  memo_american_put[no_of_divisions][i+no_of_divisions] = max(0.0, strike_price -  current_stock_price);
	} else { 
		return  memo_american_put[k][i+no_of_divisions] = max(max(0.0, strike_price -  current_stock_price),
			(uptick_prob*american_put_option(k + 1, i + 1, current_stock_price * up_factor) +
				notick_prob*american_put_option(k + 1, i, current_stock_price) +
				downtick_prob*american_put_option(k + 1, i - 1,  current_stock_price/ up_factor)) / R);
	} 
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

int main(int argc, char* argv[])
{
	
	sscanf_s(argv[1], "%f", &expiration_time);
	sscanf_s(argv[2], "%d", &no_of_divisions);
	sscanf_s(argv[3], "%f", &risk_free_rate);
	sscanf_s(argv[4], "%f", &volatility);
	sscanf_s(argv[5], "%f", &initial_stock_price);
	sscanf_s(argv[6], "%f", &strike_price);

	up_factor = exp(volatility*sqrt(2 * expiration_time / ((float)no_of_divisions)));

	R = exp(risk_free_rate*expiration_time / ((float)no_of_divisions));

	uptick_prob = pow((sqrt(R) - 1 / sqrt(up_factor)) / (sqrt(up_factor) - 1 / sqrt(up_factor)), 2);

	downtick_prob = pow(((sqrt(up_factor) - sqrt(R)) / (sqrt(up_factor) - 1 / sqrt(up_factor))), 2);
	//(R - (1 / up_factor)) / (up_factor - (1 / up_factor));
	notick_prob = 1 - uptick_prob - downtick_prob;
	

	cout << "(Memoized) Recursive Trinomial American-Asian Option Pricing" << endl;
	cout << "Expiration Time (Years) = " << expiration_time << endl;
	cout << "Number of Divisions = " << no_of_divisions << endl;
	cout << "Risk Free Interest Rate = " << risk_free_rate << endl;
	cout << "Volatility (%age of stock value) = " << volatility * 100 << endl;
	cout << "Initial Stock Price = " << initial_stock_price << endl;
	cout << "Strike Price = " << strike_price << endl;
	cout << "--------------------------------------" << endl;
	cout << "R = " << R << endl;
	cout << "Up Factor = " << up_factor << endl;
	cout << "Uptick Probability = " << uptick_prob << endl;
	cout << "Notick Probability = " << notick_prob << endl;
	cout << "Downtick Probability = " << downtick_prob << endl;
	cout << "--------------------------------------" << endl;
	initialize();
	double call_price = american_call_option(0, 0, initial_stock_price);
	cout << "Trinomial Price of an American Call Option = " << call_price << endl;
	initialize();
	double put_price = american_put_option(0, 0, initial_stock_price);
	cout << "Trinomial Price of an American Put Option = " << put_price << endl;
	cout << "--------------------------------------" << endl;
	initialize();
	double euro_call = european_call_option(0, 0);
	cout << "Trinomial Price of an European Call Option = " << euro_call << endl;
	cout << "Call Price according to Black-Scholes = " <<
		option_price_call_black_scholes(initial_stock_price, strike_price, risk_free_rate,
			volatility, expiration_time) << endl;
	cout << "--------------------------------------" << endl;
	initialize();
	double euro_put = european_put_option(0, 0);
	cout << "Trinomial Price of an European Put Option = " << euro_put << endl;
	cout << "Put Price according to Black-Scholes = " <<
		option_price_put_black_scholes(initial_stock_price, strike_price, risk_free_rate,
			volatility, expiration_time) << endl;
}