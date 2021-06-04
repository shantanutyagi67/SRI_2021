
#include "header.h"
#include <fstream>


/* Main Function to Evolve Ising Model and Return list of States */
int main(int argc, char** argv)
{
	int n_rows = 3; // number of rows
	int n_cols = 3; // number of columns
	double T_c = 2.26918531421; // critical temp of ising model (coupling = +1, dimless temperature)
	double T = 10*T_c;  // Temperature of the ising model (Note, critical temp in 2d is 2.26918531421)
	long N_steps = 1000;  // Number of steps to take
	bool output_spins = false; // flag to output spin matrix at each time step

	std::vector<std::vector<double>> x;	// 2d array for the current spin matrix

	
	/* Make Sure Input Params are valid */
	try{
		if(n_rows != n_cols){
			throw std::runtime_error("SQUARE LATTICES ONLY - N_ROWS MUST MATCH N_COLS. CHECK VALUES IN MAIN.CPP");  // raise error if lattice isn't square
		}
		if(n_rows > 8){
			throw std::runtime_error("N_ROWS MUST BE LESS THAN 8!");  // raise error if we have a lattice bigger than 8x8
		}
		if(T <= 0.){
			throw std::runtime_error("TEMPERATURE MUST BE POSITIVE!");
		}
	}catch(std::runtime_error &e){
		std::cerr << "RUNTIME ERROR: " << e.what() << std::endl;
		return 1;
	}

	std::cout << "Starting...\n";
	
	/* Store parameters in file */
	std::ofstream params("params.txt");
	params << "n_rows\tn_cols\tTemperature\tSteps\t\n";
	params << n_rows << "\t" << n_cols << "\t" << T << "\t" << N_steps<< "\n";
	params.close();


	/* Initialize Simulation */
	Ising_Model model = Ising_Model(n_rows,n_cols,T,N_steps); // create Ising Model named 'model'
	
	/* Run Sim and Store Spin Matrices in File */
	if(output_spins == true){
		std::ofstream time_series("time_series.txt");
		model.evolve(time_series); // evolve it through N_steps
		time_series.close();		
	}else{
		model.evolve();
	}

	/* Write time series of states to output file */
	std::ofstream states("states.txt");
	for(long i=0;i<N_steps;i++){
		states << model.states[i] << std::endl;
	}
	states.close();

	std::cout << "Done.\n"; // End Main

}