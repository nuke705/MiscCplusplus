// hw6.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <cstdlib>
#include <algorithm>
#include <complex>
#include "E:\study\newmat10\newmat.h"
#include "E:\study\newmat10\newmatap.h"  // Need this for the FFT algorithm

using namespace std;

// function that is used for sorting a 2-col array based on the value of the entries in the 2nd column.
// taken from http://stackoverflow.com/questions/3041897/sorting-a-2-dimensional-array-on-multiple-columns
bool compareTwoRows2(double* rowA, double* rowB) {
	return ((rowA[1]>rowB[1]) || ((rowA[1] == rowB[1]) && (rowA[0]>rowB[0])));
}

class Filtering_Instance
{
	// Private
	int no_of_terms, no_of_data_points;

	// Private member function that computes the mean of a data array
	double compute_mean(ColumnVector data)
	{
		// write the code to compute the mean of "data"
		double meantemp=0;
		for (int i = 1; i <= data.Nrows(); i++) {
			meantemp += data(i);
		}
		return meantemp / static_cast<double> (no_of_data_points);
	}

	// Private member function that computes the magnitude of (an array of) complex #s
	void compute_magnitude(ColumnVector &magnitude, ColumnVector real_part, ColumnVector imag_part)
	{
		// write the code to compute 
		for (int i = 1; i <= imag_part.Nrows(); i++) {
			magnitude(i) = sqrt(pow(real_part(i), 2) + pow(imag_part(i), 2));
		}
	}

	// Private member function that reads the data from the input file 
	// and stores it in "data" 
	void get_data(char* file_name, ColumnVector &data)
	{
		// write code that reads the ticker-data from the input file
		// and stores it in "data"
		double val;
		
		ifstream inputfile(file_name);
		if (inputfile.is_open()) {
			int i = 1;
			while (inputfile >> val) {
				data(i) = val;
				i++;
			}
		}
		//cout << "reachget" << endl;
	}

	// private member function that writes the data file into a file
	void write_data(char* file_name, ColumnVector &data)
	{
		// write code that writes "data" to file_name.
		ofstream outfile(file_name);
		for (int i = 1; i <= no_of_data_points; i++) {
			outfile << data(i) << endl;
		}
		//cout << "reachwrite" << endl;
	}

	// private member function that filters data using the FFT 
	// The filtered data is computed using the top "no_of_terms"-many 
	// magnitude-components of the orginal data
	void filter_the_data(ColumnVector &data, ColumnVector &filtered_data, int no_of_terms)
	{
		ColumnVector fft_real_part(data.Nrows()), fft_imag_part(data.Nrows());
		ColumnVector mean_adjusted_data(data.Nrows()), magnitude(data.Nrows());

		double mean = compute_mean(data);
		for (int i = 1; i <= data.Nrows(); i++)
			mean_adjusted_data(i) = data(i) - mean;

		RealFFT(mean_adjusted_data, fft_real_part, fft_imag_part);
		compute_magnitude(magnitude, fft_real_part, fft_imag_part);

		// creating a two dimensional array: first col is the index; second col is the 
		// magnitude.  The plan is to have this 2-D array sorted using the 2nd col (ie.
		// sorted based on magnitude). Then we pick the top "no_of_terms"-many of these
		// components to reconstitute/reconstruct the signal back. 

		double** two_dimensional_array = new double*[fft_imag_part.Nrows()];
		for (int i = 0; i < fft_imag_part.Nrows(); i++)
			two_dimensional_array[i] = new double[2];

		for (int i = 0; i < fft_imag_part.Nrows(); i++)
		{
			two_dimensional_array[i][0] = i;
			two_dimensional_array[i][1] = magnitude(i + 1);
		}
		std::sort(two_dimensional_array, two_dimensional_array + fft_imag_part.Nrows(), &compareTwoRows2);

		// if do_we_pick_this(i) == 1, then we keep that component for reconstruction
		// of the filtered signal.  The rest of the array-names should be self-explanatory
		ColumnVector do_we_pick_this(fft_imag_part.Nrows());
		ColumnVector filtered_fft_real_part(fft_imag_part.Nrows());
		ColumnVector filtered_fft_imag_part(fft_imag_part.Nrows());

		// write the code for picking the top "no_of_terms" many magnitudes
		// put the real-part in "filtered_fft_real_part" and imaginary-part in 
		// "filtered_fft_imag_part" -- and reconstruct the filtered signal as 
		// shown below. 
		//cout << "reachbefore picking" << endl;
		for (int i = 0; i < fft_imag_part.Nrows(); i++) {
				filtered_fft_real_part(i+1) = 0;
				filtered_fft_imag_part(i+1) = 0;
		}
		for (int i = 0; i < no_of_terms; i++) {
			int pick = two_dimensional_array[i][0]+1;
			filtered_fft_real_part(pick) = fft_real_part(pick);
			filtered_fft_imag_part(pick) = fft_imag_part(pick);
		}
		
		// reconstructed signal using just the "no_of_terms"-many top-magnitude 
		// components.
		RealFFTI(filtered_fft_real_part, filtered_fft_imag_part, filtered_data);

		// write code to add the mean-back to the reconstructed, filtered-signal
		for (int i = 1; i <= data.Nrows(); i++) { 
			filtered_data(i) = filtered_data(i) + mean; 
		}
		//cout << "reach4" << endl;
	}

public:
	// Public member function that reads the incomplete puzzle
	// we are not doing any checks on the input puzzle -- that is,
	// we are assuming they are indeed valid
	void run_the_filter(int argc, char * const argv[])
	{
		sscanf_s(argv[1], "%d", &no_of_terms);
		sscanf_s(argv[2], "%d", &no_of_data_points);

		std::cout << "Input File Name: " << argv[3] << std::endl;
		std::cout << "Number of data points in the input file = " << no_of_data_points << endl;
		std::cout << "Number of dominant terms in the FFT = " << no_of_terms << endl;
		std::cout << "Output File Name: " << argv[4] << std::endl;

		ColumnVector data(no_of_data_points), filtered_data(no_of_data_points);

		// get ticker data
		get_data(argv[3], data);

		// filter the ticker data
		filter_the_data(data, filtered_data, no_of_terms);

		// write the filtered data
		write_data(argv[4], filtered_data);
	}
};


int main(int argc, char* argv[])
{
	Filtering_Instance x;
	x.run_the_filter(argc, argv);
}

