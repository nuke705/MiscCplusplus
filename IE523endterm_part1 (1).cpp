#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <cstdlib>
#include <vector>
#include <numeric>
#include <random>
//#include "E:\study\hw10am\normdist.h" 
#define UNIT_STEP(x) ((x)>0?(1):(0))
using namespace std;

double risk_free_rate, strike_price, initial_stock_price, expiration_time, volatility,bar_price;
int no_of_trials, no_of_divisions;

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




double max(double a, double b) {
	return (b < a) ? a : b;
}

/*
double get_uniform()
{
	uniform_real_distribution <double> distribution(0.0, 1.0);
	return  distribution(generator);
}*/

double get_uniform()
{
	return ((double)rand() / RAND_MAX);
}

double cdi_h_smaller_than_x (double s, double h, double x, double r, double sigma, double t) {
	double lambda = (r + (sigma*sigma) / 2) / (sigma*sigma);
	double y = log((h*h)/(s*x)) / (sigma*sqrt(t)) + lambda*sigma*sqrt(t);
	return s*pow(h / s, 2 * lambda)*N(y) - x*exp(-r*t)*pow(h / s, 2 * lambda - 2)*N(y - sigma*sqrt(t));

}

double cdo_h_samller_than_x() {
	return option_price_call_black_scholes(initial_stock_price, strike_price, risk_free_rate,
		volatility, expiration_time) - cdi_h_smaller_than_x(initial_stock_price, bar_price, strike_price, risk_free_rate,
			volatility, expiration_time);
}

double cdo_h_greater_than_x(double s, double h, double x, double r, double sigma, double t) {
	double lambda = (r + (sigma*sigma) / 2) / (sigma*sigma);
	double x1 = log(s / h) / (sigma*sqrt(t)) + lambda*sigma*sqrt(t);
	double y1 = log(h/s) / (sigma*sqrt(t)) + lambda*sigma*sqrt(t);
	return s*N(x1)- x*exp(-r*t)*N(x1- sigma*sqrt(t)) - s*pow(h / s, 2 * lambda)*N(y1)
		+ x*exp(-r*t)*pow(h / s, 2 * lambda - 2)*N(y1 - sigma*sqrt(t));

}

double cdi_h_greater_than_x() {
	return option_price_call_black_scholes(initial_stock_price, strike_price, risk_free_rate,
		volatility, expiration_time) - cdo_h_greater_than_x(initial_stock_price, bar_price, strike_price, risk_free_rate,
			volatility, expiration_time);
}

double pdo_h_greater_than_x() { return 0; }
double pdi_h_greater_than_x() {
	return option_price_put_black_scholes(initial_stock_price, strike_price, risk_free_rate,
		volatility, expiration_time);
}

double pdi_h_smaller_than_x(double s, double h, double x, double r, double sigma, double t) {
	double lambda = (r + (sigma*sigma) / 2) / (sigma*sigma);
	double x1 = log(s / h) / (sigma*sqrt(t)) + lambda*sigma*sqrt(t);
	double y1 = log(h / s) / (sigma*sqrt(t)) + lambda*sigma*sqrt(t);
	double y = log((h*h) / (s*x)) / (sigma*sqrt(t)) + lambda*sigma*sqrt(t);
	return -s*N(-x1) + x*exp(-r*t)*N(-x1 + sigma*sqrt(t)) + s*pow(h / s, 2 * lambda)*(N(y) - N(y1)) -
		x*exp(-r*t)*pow(h / s, 2 * lambda - 2)*(N(y - sigma*sqrt(t)) - N(y1 - sigma*sqrt(t)));

}

double pdo_h_smaller_than_x() {
	return option_price_put_black_scholes(initial_stock_price, strike_price, risk_free_rate,
		volatility, expiration_time) - pdi_h_smaller_than_x(initial_stock_price, bar_price, strike_price, risk_free_rate,
			volatility, expiration_time);

}

double get_pc(double ST) {
	if (initial_stock_price <= bar_price) {return 1;}
	if (ST <= bar_price) { return 1.0; }
	return exp(-(2 * log(initial_stock_price / bar_price)*log(ST / bar_price)) / (volatility*volatility*(double)expiration_time));
}

/*// perform one trial
vector <double> simulation(double delta_T, double delta_R, double delta_SD) {
	vector <vector<double>> price;
	vector <bool> breached;
	for (int i = 0; i < 4; i++) {
		vector <double> temp;
		price.push_back(temp);
		breached.push_back(false);
	}

	// by sharing random variables we create 4 paths 
	double current_stock_price1 = initial_stock_price;
	double current_stock_price2 = initial_stock_price;
	double current_stock_price3 = initial_stock_price;
	double current_stock_price4 = initial_stock_price;
	double k1 = strike_price, k2 = strike_price, k3 = strike_price, k4 = strike_price;

	for (int j = 0; j < no_of_divisions; j++) {
		// create the unit normal variates using the Box-Muller Transform
		double x = get_uniform();           double y = get_uniform();
		double a = sqrt(-2.0*log(x)) * cos(6.283185307999998*y);
		double b = sqrt(-2.0*log(x)) * sin(6.283185307999998*y);
		if (isfinite(a) == 1 && isfinite(b) == 1) {
			current_stock_price1 = current_stock_price1*exp(delta_R + delta_SD*a);
			current_stock_price2 = current_stock_price2*exp(delta_R - delta_SD*a);
			current_stock_price3 = current_stock_price3*exp(delta_R + delta_SD*b);
			current_stock_price4 = current_stock_price4*exp(delta_R - delta_SD*b);
			if (current_stock_price1 <= bar_price) {
				breached[0] = true; k1 = 0; current_stock_price1 = 0;
			}
			if (current_stock_price2 <= bar_price) {
				breached[1] = true; k2 = 0; current_stock_price2 = 0;
			}
			if (current_stock_price3 <= bar_price) {
				breached[2] = true; k3 = 0; current_stock_price3 = 0;
			}
			if (current_stock_price4 <= bar_price) {
				breached[3] = true; k4 = 0; current_stock_price4 = 0;
			}
		}
			if (j == no_of_divisions - 1) {
				price[0].push_back(current_stock_price1);
				price[1].push_back(current_stock_price2);
				price[2].push_back(current_stock_price3);
				price[3].push_back(current_stock_price4);
				
			}
			//cout << x<< " ------ " <<y << endl;
		
	}
	double put_option_price = 0.0, put_option_price_pc = 0.0, call_option_price_pc = 0.0, call_option_price = 0.0;
	
	// compute prices with pc regardless if knocked-out
	call_option_price_pc = (max(0.0, price[0][0] - strike_price)*(1 - get_pc(price[0][0])) +
		max(0.0, price[1][0] - strike_price)*(1 - get_pc(price[1][0])) +
		max(0.0, price[2][0] - strike_price)*(1 - get_pc(price[2][0])) +
		max(0.0, price[3][0] - strike_price)*(1 - get_pc(price[3][0]))) / 4.0;
	put_option_price_pc = (max(0.0, k1 - price[0][0])*(1 - get_pc(price[0][0])) +
		max(0.0, k2 - price[1][0])*(1 - get_pc(price[1][0])) +
		max(0.0, k3 - price[2][0])*(1 - get_pc(price[2][0])) +
		max(0.0, k4 - price[3][0])*(1 - get_pc(price[3][0]))) / 4.0;
	// Find out which path breaches barrier and set its ST = 0
	/*for (int i = 0; i < 4; i++) {
		if (breached[i] == true) {
			price[i][0] = 0.0;
		}
	}
	//cout << "reach1" << endl;
	// compute prices with simulation
	call_option_price = (max(0.0, price[0][0] - strike_price) +
		max(0.0, price[1][0] - strike_price) +
		max(0.0, price[2][0] - strike_price) +
		max(0.0, price[3][0] - strike_price)) / 4.0;
	put_option_price = (max(0.0, k1 - price[0][0]) +
		max(0.0, k2 - price[1][0]) +
		max(0.0, k3 - price[2][0]) +
		max(0.0, k4 - price[3][0])) / 4.0;
	
	// output result (it has order): call, put, call pc, put pc
	vector <double> call_and_put;
	call_and_put.push_back(call_option_price); 
	call_and_put.push_back(put_option_price); 
	call_and_put.push_back(call_option_price_pc); 
	call_and_put.push_back(put_option_price_pc);
	//cout << call_and_put[0] <<" "<< call_and_put[1]<< endl;
	return call_and_put;

}
*/
int main(int argc, char* argv[])
{
	sscanf_s(argv[1], "%lf", &expiration_time);
	sscanf_s(argv[2], "%lf", &risk_free_rate);
	sscanf_s(argv[3], "%lf", &volatility);
	sscanf_s(argv[4], "%lf", &initial_stock_price);
	sscanf_s(argv[5], "%lf", &strike_price);
	sscanf_s(argv[6], "%d", &no_of_trials);
	sscanf_s(argv[7], "%d", &no_of_divisions);
	sscanf_s(argv[8], "%lf", &bar_price);

	double R = (risk_free_rate - 0.5*pow(volatility, 2))*expiration_time;
	double SD = volatility*sqrt(expiration_time);

	cout << "--------------------------------" << endl;
	cout << "European Down-and-Out Continuous Barrier Options Pricing via Monte Carlo Simulation" << endl;
	cout << "Expiration Time (Years) = " << expiration_time << endl;
	cout << "Risk Free Interest Rate = " << risk_free_rate << endl;
	cout << "Volatility (%age of stock value) = " << volatility * 100 << endl;
	cout << "Initial Stock Price = " << initial_stock_price << endl;
	cout << "Strike Price = " << strike_price << endl;
	cout << "Barrier Price = " << bar_price << endl;
	cout << "Number of Trials = " << no_of_trials << endl;
	cout << "Number of Divisions = " << no_of_divisions << endl;
	cout << "--------------------------------" << endl;

	// chop expiration time into no_of_divisions many segments 
	// figure out the motion within each segment
	double delta_T = expiration_time / ((double)no_of_divisions);
	double delta_R = (risk_free_rate - 0.5*pow(volatility, 2))*delta_T;
	double delta_SD = volatility*sqrt(delta_T);


	double call = 0,put = 0, call_pc = 0, put_pc = 0;

	for (int i = 0; i < no_of_trials; i++) {
		double current_stock_price1 = initial_stock_price, current_stock_price2 = initial_stock_price;
		double current_stock_price3 = initial_stock_price, current_stock_price4 = initial_stock_price;
		double current_stock_price1_pc = initial_stock_price,current_stock_price2_pc = initial_stock_price; 
		double current_stock_price3_pc = initial_stock_price, current_stock_price4_pc = initial_stock_price;
		double k1 = strike_price, k2 = strike_price, k3 = strike_price, k4 = strike_price;
		for (int j = 0; j < no_of_divisions; j++) {
			double x = get_uniform(); double y = get_uniform();
			double a = sqrt(-2.0*log(x)) * cos(6.283185307999998*y);
			double b = sqrt(-2.0*log(x)) * sin(6.283185307999998*y);
			//check inf
			if (1 == isfinite(a) && 1 == isfinite(b)) {
				current_stock_price1 = current_stock_price1*exp(delta_R + delta_SD*a);
				current_stock_price2 = current_stock_price2*exp(delta_R - delta_SD*a);
				current_stock_price3 = current_stock_price3*exp(delta_R + delta_SD*b);
				current_stock_price4 = current_stock_price4*exp(delta_R - delta_SD*b);
				current_stock_price1_pc = current_stock_price1_pc*exp(delta_R + delta_SD*a);
				current_stock_price2_pc = current_stock_price2_pc*exp(delta_R - delta_SD*a);
				current_stock_price3_pc = current_stock_price3_pc*exp(delta_R + delta_SD*b);
				current_stock_price4_pc = current_stock_price4_pc*exp(delta_R - delta_SD*b);
				//check if breached
				if (current_stock_price1 <= bar_price) {
					current_stock_price1 = 0; k1 = 0;
				}
				if (current_stock_price2 <= bar_price) {
					current_stock_price2 = 0; k2 = 0;
				}
				if (current_stock_price3 <= bar_price) {
					current_stock_price3 = 0; k3 = 0;
				}
				if (current_stock_price4 <= bar_price) {
					current_stock_price4 = 0; k4 = 0;
				}
			}
		}//sum and avg
		call += (max(current_stock_price1 - strike_price, 0) +
			max(current_stock_price2 - strike_price, 0) +
			max(current_stock_price3 - strike_price, 0) +
			max(current_stock_price4 - strike_price, 0)) / 4.0;
		put += (max(k1 - current_stock_price1, 0) +
			max(k2 - current_stock_price2, 0) +
			max(k3 - current_stock_price3, 0) +
			max(k4 - current_stock_price4, 0)) / 4.0;
		// adjusted term
		call_pc += (max(current_stock_price1_pc - strike_price, 0)*(1 - get_pc(current_stock_price1_pc)) +
			max(current_stock_price2_pc - strike_price, 0)*(1 - get_pc(current_stock_price2_pc)) +
			max(current_stock_price3_pc - strike_price, 0)*(1 - get_pc(current_stock_price3_pc)) +
			max(current_stock_price4_pc - strike_price, 0)*(1 - get_pc(current_stock_price4_pc))) / 4.0;
		put_pc += (max(strike_price - current_stock_price1_pc, 0)*(1 - get_pc(current_stock_price1_pc)) +
			max(strike_price - current_stock_price2_pc, 0)*(1 - get_pc(current_stock_price2_pc)) +
			max(strike_price - current_stock_price3_pc, 0)*(1 - get_pc(current_stock_price3_pc)) +
			max(strike_price - current_stock_price4_pc, 0)*(1 - get_pc(current_stock_price4_pc))) / 4.0;
	}
	
	call = exp(-risk_free_rate*expiration_time)*(call / no_of_trials);
	put = exp(-risk_free_rate*expiration_time)*(put / no_of_trials);
	
	call_pc = exp(-risk_free_rate*expiration_time)*(call_pc / no_of_trials);
	put_pc = exp(-risk_free_rate*expiration_time)*(put_pc / no_of_trials);
	
	double theoretical_calldo = 0.0; double theoretical_putdo = 0.0;
	if (bar_price <= strike_price) {
		theoretical_calldo = cdo_h_samller_than_x(); theoretical_putdo = pdo_h_smaller_than_x();
	}
	else {
		theoretical_calldo =
			cdo_h_greater_than_x(initial_stock_price, bar_price, strike_price, risk_free_rate,volatility, expiration_time); 
		theoretical_putdo = pdo_h_greater_than_x();
	}
	
	cout << "--------------------------------" << endl;
	cout << "The average Call Price by explicit simulation  = " << call << endl;
	cout << "The call price using the (1-p)-adjustment term = " << call_pc << endl;
	cout << "Theoretical Call Price	                       = " << theoretical_calldo << endl;
	cout << "--------------------------------" << endl;
	cout << "The average Put Price by explicit simulation  = " << put << endl;
	cout << "The put price using the (1-p)-adjustment term = " << put_pc << endl;
	cout << "Theoretical Put Price                         = " << theoretical_putdo << endl;
	cout << "--------------------------------" << endl;
}