#pragma once
#include <vector>
#include <iostream>
#include <random>
#include <stdio.h>

/* Structure to hold the current properties of the Ising Model */
struct Ising_Model
{
	Ising_Model(int n_rows, int n_cols, double T, long N_steps);  // signature for constructor, construct instance of class

	std::vector<unsigned long long> states; // 1d array to hold the name of all states
	std::vector<std::vector<int>> spin_matrix;	// 2d array for the current spin matrix
	long current_step; // integer denoting the number of steps taken
	int num_rows; // number of rows
	int num_cols; //number of columns
	long num_steps; // number of simulation steps
	double temp; // temperature of the model

	/* Define Methods */
	void evolve(std::ostream &time_series);
	void evolve(void);
};


/* Overload '<<' to easily write spin_matrix data to file */
std::ostream &operator<<(std::ostream &out, std::vector<std::vector<int>> spin_matrix);

/* Other Functions */
unsigned long long get_state(std::vector<std::vector<int>> spin_matrix, int num_rows, int num_cols);

//HEADER
// std::vector<std::vector<double>> get_TPM(std::vector<unsigned long long> states); // relevant only for EI.cpp 

// //MAIN:
// std::cout << "GET EI:" << std::endl;

// std::cout << "FIRST STATE: " << model.states[0] << std::endl;
// x = get_TPM(model.states);