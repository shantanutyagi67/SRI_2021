#include "header.h"
#include <cmath>
#include <stdexcept>

/* Define our constructor. Pass Args to this to create an instance of Ising Model Class */
/* Note: Its a constructor because there's no return type and method name 'Ising_Model' matches class name. */

Ising_Model::Ising_Model(int n_rows, int n_cols, double T, long N_steps){

	current_step = 0; // index for current step
	states = std::vector<unsigned long long>(N_steps,0); // initialize 1d array to hold sequence of integer states
	spin_matrix = std::vector<std::vector<int>>(n_rows,std::vector<int>(n_cols,0)); // initialize 2d matrix with null spins
	num_rows = n_rows;
	num_cols = n_cols;
	num_steps = N_steps;
	temp = T;

	/* Initialize Spin Matrix with +1 and -1 randomly drawn */
	std::random_device rd; // used to obtain seed for random number engine
	std::mt19937 gen(rd()); // standard mersenne_twister_engine seeded with rd()
	std::uniform_int_distribution<> uni_dis(0,1);
	for(int i=0;i<num_rows;i++){
		for(int j=0;j<num_cols;j++){
			int rand_int = uni_dis(gen);         // draw random int either 0 or 1
			if(rand_int == 1){
				 spin_matrix[i][j]= 1;
			}else{
				spin_matrix[i][j] = -1;  // if rand_int was 0 change it to -1
			}
		}
	}

	states[current_step] = get_state(spin_matrix,num_rows,num_cols);

}

// This method evolves and writes spins to file
void Ising_Model::evolve(std::ostream &time_series){
	int energy = 0;  // integer to hold the local energy value
	int row_index,col_index;
	int left_nn,right_nn,down_nn,up_nn; // integers to hold the values of the nearest neighbors

	std::random_device rd; // used to obtain seed for random number engine
	std::mt19937 gen(rd()); // standard mersenne_twister_engine seeded with rd()
	std::uniform_real_distribution<> real_dis(0,1); // real distribution to draw from (0 to 1)
	std::uniform_int_distribution<> uni_dis(0,num_rows-1); // uniform distribution to draw from (0 to n_rows-1 inclusive)

	//For each step get the spin matrix and convert to integer-named state
	for(int i=1;i<=num_steps;i++){
		//time_series << i << "\n";
		current_step = i;

		/* Choose a node at random */
		col_index = uni_dis(gen); // random int between 0 and n_cols-1
		row_index = uni_dis(gen); // random int between 0 and n_rows-1


		/* Get value of spin matrix at nearest neighbor sites using PBC*/
		/* Note, indexing is mixed up (i.e. rows and columns are switched when printing) */
		/* However, this does not adversely affect the results in any way */
		if((row_index-1) >= 0){
			left_nn = spin_matrix[row_index-1][col_index];
		}else{
			left_nn = spin_matrix[num_rows-1][col_index];
		}
		right_nn = spin_matrix[(row_index+1)%num_rows][col_index];
		up_nn = spin_matrix[row_index][(col_index+1)%num_cols];
		if((col_index -1) >= 0){
			down_nn = spin_matrix[row_index][col_index-1];
		}else{
			down_nn = spin_matrix[row_index][num_cols-1];
		}

		/* Calculate Energy and Decide whether or not to Flip */
		energy = -1*spin_matrix[row_index][col_index]*(left_nn+right_nn+down_nn+up_nn);
		if(energy > 0){
			spin_matrix[row_index][col_index] = -1*spin_matrix[row_index][col_index];
		}else{
			double r = real_dis(gen);
			if(r <= std::exp(2.*double(energy)/temp)){
				spin_matrix[row_index][col_index] = -1*spin_matrix[row_index][col_index];
			}
		}

		// save spin matrix in file
		for(int i=0;i<num_rows;i++){
			for(int j=0;j<num_rows;j++){
				if(j == num_rows-1 and i == num_rows-1){
					time_series << spin_matrix[i][j] << "\n"; // end of line
				}else{
					time_series << spin_matrix[i][j] << "\t"; // seperate entries by tab
				}
			}
		}	

		states[current_step] = get_state(spin_matrix,num_rows,num_cols);

	}
}

// This method evolves but does not write spins to file
void Ising_Model::evolve(void){
	int energy = 0;  // integer to hold the local energy value
	int row_index,col_index;
	int left_nn,right_nn,down_nn,up_nn; // integers to hold the values of the nearest neighbors

	std::random_device rd; // used to obtain seed for random number engine
	std::mt19937 gen(rd()); // standard mersenne_twister_engine seeded with rd()
	std::uniform_real_distribution<> real_dis(0,1); // real distribution to draw from (0 to 1)
	std::uniform_int_distribution<> uni_dis(0,num_rows-1); // uniform distribution to draw from (0 to n_rows-1 inclusive)

	//For each step get the spin matrix and convert to integer-named state
	for(int i=1;i<=num_steps;i++){
		//time_series << i << "\n";
		current_step = i;

		/* Choose a node at random */
		col_index = uni_dis(gen); // random int between 0 and n_cols-1
		row_index = uni_dis(gen); // random int between 0 and n_rows-1


		/* Get value of spin matrix at nearest neighbor sites using PBC*/
		/* Note, indexing is mixed up (i.e. rows and columns are switched when printing) */
		/* However, this does not adversely affect the results in any way */
		if((row_index-1) >= 0){
			left_nn = spin_matrix[row_index-1][col_index];
		}else{
			left_nn = spin_matrix[num_rows-1][col_index];
		}
		right_nn = spin_matrix[(row_index+1)%num_rows][col_index];
		up_nn = spin_matrix[row_index][(col_index+1)%num_cols];
		if((col_index -1) >= 0){
			down_nn = spin_matrix[row_index][col_index-1];
		}else{
			down_nn = spin_matrix[row_index][num_cols-1];
		}

		/* Calculate Energy and Decide whether or not to Flip */
		energy = -1*spin_matrix[row_index][col_index]*(left_nn+right_nn+down_nn+up_nn);
		if(energy > 0){
			spin_matrix[row_index][col_index] = -1*spin_matrix[row_index][col_index];
		}else{
			double r = real_dis(gen);
			if(r <= std::exp(2.*double(energy)/temp)){
				spin_matrix[row_index][col_index] = -1*spin_matrix[row_index][col_index];
			}
		}

		states[current_step] = get_state(spin_matrix,num_rows,num_cols);

	}
}


/* Define function to overload the operator "<<" such that if the input 
is the spin matrix it will print it as a table */
std::ostream &operator<<(std::ostream &out, std::vector<std::vector<int>> spin_matrix){
	
	int n_rows = spin_matrix.size();
	int n_cols = spin_matrix.size();
	
	out << "current_state\t" << get_state(spin_matrix,n_rows,n_cols) << "\n";
	for(int i=0;i<n_rows;i++){
		for(int j=0;j<n_rows;j++){
			if((j+1)%(n_rows) == 0){
				out << spin_matrix[i][j] << "\n";
			}else{
				out << spin_matrix[i][j] << "\t";	
			}
		}
	}

	out << std::endl;
	return out;
}

/* Function to convert 2d matrix into unique integer identifier via binary conversion */
unsigned long long get_state(std::vector<std::vector<int>> spin_matrix, int num_rows, int num_cols){

	double power = 0.;  // each node represents a power of 2
	unsigned long long state = 0;  // this will be the identifier of our state

	for(int i=0;i<num_rows;i++){
		for(int j=0;j<num_cols;j++){
			if(spin_matrix[i][j] == 1){
				state = state + std::pow(2.,power);
				//std::cout << "state: " << state << "\t";
			}
			power = power + 1.;
		}
	}

	return state;
}











