
#include "header.h"
#include <fstream>
#include<cmath>
#include <time.h>
using namespace std;
/* Main Function to Evolve Ising Model and Return list of States */
int main(int argc, char** argv)
{
	int n_rows = 128; // number of rows
	int n_cols = 128; // number of columns
	double T_c = 2.26918531421; // critical temp of ising model (coupling = +1, dimless temperature)
	double T = (double)2.5*T_c;  // Temperature of the ising model (Note, critical temp in 2d is 2.26918531421)
	long N_steps = 100000000;  // Number of steps to take
	bool output_spins = false, output_states = false; // flag to output spin matrix at each time step
	double M = 0, avgM = 0;
	double monteCarlo = 1,include = 0;
	int skip = 1000;
	std::vector<std::vector<double>> x;	// 2d array for the current spin matrix

	
	/* Make Sure Input Params are valid */
	try{
		if(n_rows != n_cols){
			throw std::runtime_error("SQUARE LATTICES ONLY - N_ROWS MUST MATCH N_COLS. CHECK VALUES IN MAIN.CPP");  // raise error if lattice isn't square
		}
		if(n_rows > 2048){
			throw std::runtime_error("N_ROWS MUST BE LESS THAN 8!");  // raise error if we have a lattice bigger than 8x8
		}
		if(T < 0.){
			throw std::runtime_error("TEMPERATURE MUST BE POSITIVE!");
		}
	}catch(std::runtime_error &e){
		std::cerr << "RUNTIME ERROR: " << e.what() << std::endl;
		return 1;
	}

	std::cout << "Starting...\n";
	
	clock_t start_time = clock();

	for(int ii=0;ii<monteCarlo;ii++){
		M = 0;
		//cout<<"test1"<<endl;
		/* Initialize Simulation */
		Ising_Model model = Ising_Model(n_rows,n_cols,T,N_steps,skip); // create Ising Model named 'model'
		//cout<<"test2"<<endl;
		
		/* Run Sim and Store Spin Matrices in File */
		if(output_spins == true){
			std::ofstream time_series("time_series.txt");
			model.evolve(time_series); // evolve it through N_steps
			time_series.close();		
		}else{
		//cout<<"test3"<<endl;
			model.evolve();
		}

		/* Write time series of states to output file */
		if(output_states == true){
			std::ofstream states("states.txt");
			for(long i=0;i<N_steps;i++){
				states << model.states[i] << std::endl;
			}
			states.close();
		}
		
		//double diff = model.states[N_steps/skip -1] - model.states[N_steps/skip -2];
		//if(diff <0) diff*= -1;
		//if(diff<0.001){
			for(int i=0;i<n_rows;i++){
				for(int j=0;j<n_cols;j++){
					M+=model.spin_matrix[i][j];
					
				}
			}
			M = M / (double)(n_rows*n_cols);
			if (M < 0) M*=-1;
			avgM += M;
			include++;
		//}
	}
	
	avgM = (double)avgM/include;
	
	clock_t end_time = clock();
	double cpu_time = ((double)(end_time - start_time))/ CLOCKS_PER_SEC;
	
	/* Store parameters in file */
	std::ofstream params("params.txt");
	params << "n_rows\tn_cols\tT/Tc\tSteps(K)\tTime\tMagnetisation\t\n";
	params << n_rows << "," << n_cols << "," << T/T_c << "," << N_steps/1000<< "," << cpu_time << "," << avgM << "," <<"\n";
	params.close();
	
	std::cout << "Done in "<<cpu_time<<".\n"; // End Main

}
