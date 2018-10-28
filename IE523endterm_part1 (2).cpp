// Down-and-out European Discrete Barrier Option Pricing Code written by Prof. Sreenivas
// The terminal payoff vector that corresponds to indices of the underlying 
// stock prices that are in the black vis-a-vis the barrier have been discounted 
// appropriately using the Brownian Bridge adjustment
#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <cstdlib>
#include <algorithm>
#include <vector>

using namespace std;

float up_factor, uptick_prob, risk_free_rate, strike_price;
float initial_stock_price, expiration_time, volatility, barrier_price;

int no_of_bar, no_of_trials;

float max(float a, float b) {
	return (b < a) ? a : b;
}


double get_uniform()
{
	return ((double)rand() / RAND_MAX);
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

double get_pd(double st) //P_sub_d, for Brownian-Bridge correction
{
	double pd = 1.0;
	if (st <= barrier_price) {
		return 0.0;
	}
	else {
		for (int i = 1; i < no_of_bar; i++) {
			double mu_i = (initial_stock_price + (i / no_of_bar)*(st - initial_stock_price));
			double var_i = ((i / no_of_bar)*expiration_time*(1 - (i / no_of_bar)));
			pd *= (1 - N((barrier_price - mu_i) / sqrt(var_i)));
		}
		return pd;
	}
}


/*
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
	double pd = 1.0, k1=strike_price, k2 = strike_price, k3 = strike_price, k4 = strike_price;
	for (int j = 1; j <= no_of_bar; j++) {
		// create the unit normal variates using the Box-Muller Transform
		double x = get_uniform();           double y = get_uniform();
		double a = sqrt(-2.0*log(x)) * cos(6.283185307999998*y);
		double b = sqrt(-2.0*log(x)) * sin(6.283185307999998*y);
		if (isfinite(a) == 1 && isfinite(b) == 1) {
			current_stock_price1 = current_stock_price1*exp(delta_R + delta_SD*a);
			current_stock_price2 = current_stock_price2*exp(delta_R - delta_SD*a);
			current_stock_price3 = current_stock_price3*exp(delta_R + delta_SD*b);
			current_stock_price4 = current_stock_price4*exp(delta_R - delta_SD*b);
			//if (j % (no_of_divisions / no_of_bar) == 0) {
				if (current_stock_price1 <= barrier_price) {
					breached[0] = true; k1 = 0; current_stock_price1 = 0;
				}
				if (current_stock_price2 <= barrier_price) {
					breached[1] = true; k2 = 0; current_stock_price2 = 0;
				}
				if (current_stock_price3 <= barrier_price) {
					breached[2] = true; k3 = 0; current_stock_price3 = 0;
				}
				if (current_stock_price4 <= barrier_price) {
					breached[3] = true; k4 = 0; current_stock_price4 = 0;
				}
			//}
		}
		if (j == no_of_bar) {
			price[0].push_back(current_stock_price1);
			price[1].push_back(current_stock_price2);
			price[2].push_back(current_stock_price3);
			price[3].push_back(current_stock_price4);
		}
		//cout << x<< " ------ " <<y << endl;
	}
	double put_option_price = 0.0, put_option_price_pd = 0.0, call_option_price_pd = 0.0, call_option_price = 0.0;

	// compute prices with pd regardless if knocked-out
	call_option_price_pd = (max(0.0, price[0][0] - strike_price)*get_pd(price[0][0]) +
		max(0.0, price[1][0] - strike_price)*get_pd(price[1][0]) +
		max(0.0, price[2][0] - strike_price)*get_pd(price[2][0]) +
		max(0.0, price[3][0] - strike_price)*get_pd(price[3][0])) / 4.0;
	put_option_price_pd = (max(0.0, k1 - price[0][0])*get_pd(price[0][0]) +
		max(0.0, k2 - price[1][0])*get_pd(price[1][0]) +
		max(0.0, k3 - price[2][0])*get_pd(price[2][0]) +
		max(0.0, k4 - price[3][0])*get_pd(price[3][0])) / 4.0;
	// Find out which path breaches barrier and set its ST = 0
	
	for (int i = 0; i < 4; i++) {if (breached[i] == true) {price[i][0] = 0.0;}}
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
	call_and_put.push_back(call_option_price_pd);
	call_and_put.push_back(put_option_price_pd);
	//cout << call_and_put[0] <<" "<< call_and_put[1]<< endl;
	return call_and_put;

}
*/
int main(int argc, char* argv[])
{

	sscanf_s(argv[1], "%f", &expiration_time);
	sscanf_s(argv[2], "%f", &risk_free_rate);
	sscanf_s(argv[3], "%f", &volatility);
	sscanf_s(argv[4], "%f", &initial_stock_price);
	sscanf_s(argv[5], "%f", &strike_price);
	sscanf_s(argv[6], "%d", &no_of_trials);
	sscanf_s(argv[7], "%d", &no_of_bar);
	sscanf_s(argv[8], "%f", &barrier_price);
	//no_of_divisions = 1000;
	double R = (risk_free_rate - 0.5*pow(volatility, 2))*expiration_time;
	double SD = volatility*sqrt(expiration_time);
	// chop expiration time into no_of_divisions many segments 
	// figure out the motion within each segment
	double delta_T = expiration_time / ((double)no_of_bar);
	double delta_R = (risk_free_rate - 0.5*pow(volatility, 2))*delta_T;
	double delta_SD = volatility*sqrt(delta_T);
	double call = 0, put = 0, call_pc = 0, put_pc = 0;

	for (int i = 0; i < no_of_trials; i++) {
		double current_stock_price1 = initial_stock_price, current_stock_price2 = initial_stock_price;
		double current_stock_price3 = initial_stock_price, current_stock_price4 = initial_stock_price;
		double current_stock_price1_pc = initial_stock_price, current_stock_price2_pc = initial_stock_price;
		double current_stock_price3_pc = initial_stock_price, current_stock_price4_pc = initial_stock_price;
		double k1 = strike_price, k2 = strike_price, k3 = strike_price, k4 = strike_price;
		for (int j = 0; j < no_of_bar; j++) {
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
				if (current_stock_price1 <= barrier_price) {
					current_stock_price1 = 0; k1 = 0;
				}
				if (current_stock_price2 <= barrier_price) {
					current_stock_price2 = 0; k2 = 0;
				}
				if (current_stock_price3 <= barrier_price) {
					current_stock_price3 = 0; k3 = 0;
				}
				if (current_stock_price4 <= barrier_price) {
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
		call_pc += (max(current_stock_price1_pc - strike_price, 0)*(get_pd(current_stock_price1_pc)) +
			max(current_stock_price2_pc - strike_price, 0)*(get_pd(current_stock_price2_pc)) +
			max(current_stock_price3_pc - strike_price, 0)*(get_pd(current_stock_price3_pc)) +
			max(current_stock_price4_pc - strike_price, 0)*(get_pd(current_stock_price4_pc))) / 4.0;
		put_pc += (max(strike_price - current_stock_price1_pc, 0)*( get_pd(current_stock_price1_pc)) +
			max(strike_price - current_stock_price2_pc, 0)*( get_pd(current_stock_price2_pc)) +
			max(strike_price - current_stock_price3_pc, 0)*( get_pd(current_stock_price3_pc)) +
			max(strike_price - current_stock_price4_pc, 0)*( get_pd(current_stock_price4_pc))) / 4.0;
	}

	call = exp(-risk_free_rate*expiration_time)*(call / no_of_trials);
	put = exp(-risk_free_rate*expiration_time)*(put / no_of_trials);

	call_pc = exp(-risk_free_rate*expiration_time)*(call_pc / no_of_trials);
	put_pc = exp(-risk_free_rate*expiration_time)*(put_pc / no_of_trials);
	
	cout << "European Down-and-Out Discrete Barrier Option Pricing via Monte Carlo Simulation" << endl;
	cout << "Expiration Time (Years) = " << expiration_time << endl;
	cout << "Risk Free Interest Rate = " << risk_free_rate << endl;
	cout << "Volatility (%age of stock value) = " << volatility * 100 << endl;
	cout << "Initial Stock Price = " << initial_stock_price << endl;
	cout << "Strike Price = " << strike_price << endl;
	cout << "Barrier Price = " << barrier_price << endl;
	cout << "Number of Trials = " << no_of_trials << endl;
	cout << "Number of Discrete Barriers = " << no_of_bar << endl;
	cout << "--------------------------------------" << endl;
	cout << "The average Call Price via explicit simulation of price paths		  = " << call << endl;
	cout << "The average Call Price with Brownian-Bridge correction on the final price = " << call_pc << endl;
	cout << "The average Put Price via explicit simulation of price paths		  = " << put << endl;
	cout << "The average Put Price with Brownian-Bridge correction on the final price  = " << put_pc<< endl;
	cout << "--------------------------------------" << endl;
}