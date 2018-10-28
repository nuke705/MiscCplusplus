//IE 523 HW9 AMERICAN OPTION BY TRINOMIAL with Recursion
//Zhi Ji
//Based on professor RS's code

#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <cstdlib>
using namespace std;

float up_factor, uptick_prob, risk_free_rate, strike_price;
float down_factor, downtick_prob, notick_prob;
float initial_stock_price, expiration_time, volatility, R;
int no_of_divisions;

float max(float a, float b) {
	return (b < a) ? a : b;
}

float american_call_option(int k, int i, float current_stock_price) {
	if (k == no_of_divisions)
		return max(0.0, (current_stock_price - strike_price));
	else
		return max((current_stock_price - strike_price),
		(uptick_prob*american_call_option(k + 1, i + 1, current_stock_price*up_factor) +
			(1- uptick_prob - downtick_prob)* american_call_option(k + 1, i, current_stock_price) + 
			downtick_prob *american_call_option(k + 1, i - 1, current_stock_price / up_factor)) / R);
}

float american_put_option(int k, int i, float current_stock_price) {
	if (k == no_of_divisions)
		return max(0.0, (strike_price- current_stock_price));
	else
		return max((strike_price - current_stock_price),
		(uptick_prob*american_put_option(k + 1, i + 1, current_stock_price*up_factor) +
			(1 - uptick_prob - downtick_prob)* american_put_option(k + 1, i, current_stock_price) +
			downtick_prob*american_put_option(k + 1, i - 1, current_stock_price / up_factor)) / R);
		//(uptick_prob*american_put_option(k + 1, i + 1, current_stock_price*up_factor) +
			//(1 - uptick_prob)*american_put_option(k + 1, i - 1, current_stock_price / up_factor)) / R);
}

int main(int argc, char* argv[])
{

	/*sscanf_s(argv[1], "%f", &expiration_time);
	sscanf_s(argv[2], "%d", &no_of_divisions);
	sscanf_s(argv[3], "%f", &risk_free_rate);
	sscanf_s(argv[4], "%f", &volatility);
	sscanf_s(argv[5], "%f", &initial_stock_price);
	sscanf_s(argv[6], "%f", &strike_price);*/

	double initial_stock_price = 100.0, strike_price = 105.0, volatility = 0.2, risk_free_rate = 0.03, expiration_time = 1.0;
	double no_of_divisions = 2500.0;

	up_factor = exp(volatility*sqrt(2*expiration_time / ((float)no_of_divisions)));

	R = exp(risk_free_rate*expiration_time / ((float)no_of_divisions));

	uptick_prob = pow((sqrt(R) - 1/sqrt( up_factor)) / (sqrt(up_factor) - 1 / sqrt(up_factor)), 2);

	downtick_prob = pow(((sqrt(up_factor) - sqrt(R) )/ (sqrt(up_factor) - 1 / sqrt(up_factor))),2);
		//(R - (1 / up_factor)) / (up_factor - (1 / up_factor));
	notick_prob = 1 - uptick_prob - downtick_prob;

	cout << "Recursive Trinomial American-Asian Option Pricing" << endl;
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
	cout << "Downtick Probability = " << downtick_prob << endl;
	cout << "Notick Probability = " << notick_prob << endl;
	cout << "--------------------------------------" << endl;
	double call_price = american_call_option(0, 0, initial_stock_price);
	cout << "Trinomial Price of an American Call Option = " << call_price << endl;
	double put_price = american_put_option(0, 0, initial_stock_price);
	cout << "Trinomial Price of an American Put Option = " << put_price << endl;
	/*
	cout << "--------------------------------------" << endl;
	cout << "Let us verify the Put-Call Parity: S+P-C = Kexp(-r*T) for American Options" << endl;
	cout << initial_stock_price << " + " << put_price << " - " << call_price;
	cout << " = " << strike_price << "exp(-" << risk_free_rate << " * " << expiration_time << ")" << endl;
	cout << initial_stock_price + put_price - call_price << " ?=? " << strike_price*exp(-risk_free_rate*expiration_time) << endl;
	if (abs(initial_stock_price + put_price - call_price - strike_price*exp(-risk_free_rate*expiration_time)) <= 1e-3)
		cout << "Looks like Put-Call Parity holds within three decimal places" << endl;
	else
		cout << "Looks like Put-Call Parity does NOT hold" << endl;
	cout << "--------------------------------------" << endl;*/

}