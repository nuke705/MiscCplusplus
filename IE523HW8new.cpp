// HW8new.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <fstream>
#include <algorithm>
#include <time.h>
#include "E:\study\newmat10\newmatap.h"           
#include "E:\study\newmat10\newmat.h"
#include "E:\study\newmat10\newmatio.h"

using namespace std;

Matrix repeated_squaring(Matrix A, int exponent, int no_rows) {
	Matrix B = IdentityMatrix(no_rows);
	if (exponent == 0) {
		return B;
	}
	if (exponent % 2 == 1) {
		return A*repeated_squaring(A*A, (exponent - 1) / 2, no_rows);
	}
	else {
		return repeated_squaring(A*A, exponent / 2, no_rows);
	}
}

Matrix bruteforce(Matrix A, int exponent, int no_rows) {
	Matrix C = A;
	for (int i = 1; i < exponent; i++) {
		C = A*C;
	}
	/*cout << "reach" << endl;
	cout << endl;
	cout << C << endl;*/
	return C;
}

double random_interval(int range)
{
	//float num = (((float) rand()) / (pow(2.0, 31.0) - 1.0));
	double num = (((float)rand()) / RAND_MAX);
	if (num < 0.5) { return -1 * range * num; }
	if (num >= 0.5) { return range * num; }
	//return 0;
}

Matrix initialize(int size) {
	Matrix A(size, size);
	for (int i = 1; i <= A.Nrows(); i++) {
		for (int j = 1; j <= A.Nrows(); j++) {
			A(i, j) = random_interval(5);
		}
	}
	return A;
}

int main(int argc, char* const argv[])
{

	int expo, size;
	sscanf_s(argv[1], "%d", &expo);
	sscanf_s(argv[2], "%d", &size);

	cout << "The number of rows/columns in the square matrix is: " << size << endl;
	cout << "The exponent is: " << expo << endl;

	Matrix randA = initialize(size);
	//cout << randA << endl;

	clock_t time_before = clock();
	Matrix result1 = repeated_squaring(randA, expo, size);
	clock_t time_after = clock();
	float diff = ((float)time_after - (float)time_before);
	cout << "Repeated squaring result: " << endl;
	cout << "It took " << setprecision(7) << diff / CLOCKS_PER_SEC << " seconds to complete " << endl;
	//cout << result1 << endl;


	clock_t time_before2 = clock();
	Matrix result2 = bruteforce(randA, expo, size);
	clock_t time_after2 = clock();
	float diff2 = ((float)time_after2 - (float)time_before2);
	cout << "Direct multiplication result: " << endl;
	cout << "It took " << setprecision(7) << diff2 / CLOCKS_PER_SEC << " seconds to complete " << endl;
	//cout << result2 << endl;

	int nn = 1000;
	vector <double> repeated_sqr_time, direct_time;
	for (int i = 1; i <= nn; i++) {


		clock_t t1 = clock();
		Matrix result_rptsqr = repeated_squaring(randA, i, size);
		clock_t t2 = clock();
		double diff_rptsqr = ((double)t2 - (double)t1);
		//double time_used_rptsqr = diff_rptsqr / CLOCKS_PER_SEC;
		double time_used_rptsqr = (diff_rptsqr) / CLOCKS_PER_SEC;
		//cout << before_rptsqr << endl;
		//cout << diff_rptsqr << endl;
		//printf("It took me (%f seconds).\n", before_rptsqr, ((float)before_rptsqr) / CLOCKS_PER_SEC);
		repeated_sqr_time.push_back(time_used_rptsqr);
	}

	for (int i = 1; i <= nn; i++) {
		double before_dir = clock();
		Matrix result_dir = bruteforce(randA, i, size);
		double after_dir = clock();
		double diff_dir = ((float)after_dir - (float)before_dir);
		double time_used_dir = diff_dir / CLOCKS_PER_SEC;
		//cout << time_used_dir << endl;
		direct_time.push_back(time_used_dir);
	}

	ofstream outf("repeated_sqr_time");
	for (int i = 0; i < nn; i++) {
		outf << repeated_sqr_time[i] << endl;
		//cout << repeated_sqr_time[i] << endl;
	}
	ofstream outf1("direct_time");
	for (int i = 0; i < nn; i++) {
		outf1 << direct_time[i] << endl;
		//cout << direct_time[i] << endl;
	}
	//cout << 1.000000011114643357 / 1.000011651321354131321<< endl;
   return 0;
}

