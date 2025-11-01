#include <NTL/ZZ.h>
#include <NTL/vector.h>
#include <NTL/matrix.h>
#include <NTL/mat_ZZ.h>
#include <NTL/vec_RR.h>
#include <NTL/RR.h>
#include <NTL/LLL.h>
#include <cmath>
#include <iostream>
#include <random>
#include <list>
#include <chrono>
using namespace std;
using namespace NTL;

random_device rand_dev;
mt19937       generator(rand_dev());

// I want to create a function that given a dimension gives me a unit vector with that many dimensions
// I need a function that will generate n random values from a range [-d, d]  
// I need a function that will generate an n dimensional vector with random values from a range
// I need a funciton that will generate k amount of n dimesnsional vectors with random values from a range
// I need a function that will compute the Hadamard ratio of a list of vectors.
// I need a function that will return a basis with a Hadamard ratio above a certain value
//
/*
 * Function that returns an n dimensional vector with coordinates initialized to 1 
*/

Vec<ZZ> GetAllOnesVec(unsigned int n){
	
	Vec<ZZ> v;
	v.SetLength(n);
	
	for(int i = 0; i < n; i++){
		v[i] = ZZ(1);
	}
	
	return v;
}


list<ZZ> GetRandValueByRange(unsigned int n, unsigned int d){

	list<ZZ> values;
	uniform_int_distribution<int> distr(-d, d);
	ZZ val;

	for (int i = 0; i < n; i++){
		val = ZZ(distr(generator)); 
		values.push_back(val);
	}

	return values;
}

Vec<ZZ> GetRandVec(unsigned int n, unsigned int d){	
	
	Vec<ZZ>  random_vec    = GetAllOnesVec(n);
	list<ZZ> random_values = GetRandValueByRange(n, d);
	
	int i = 0;
	for (ZZ val : random_values){
		random_vec[i] *= val;
		i++;
	}
	
	return random_vec;
}

Mat<ZZ> GetRandVectors(unsigned int amount, unsigned int dimension, unsigned int range){

	
	Mat<ZZ> rand_matrix;
	Vec<ZZ> rand_vec;

	rand_matrix.SetDims(amount,dimension);

	for(int i = 0; i < amount; i++){
		rand_vec = GetRandVec(dimension, range);
		rand_matrix[i] = rand_vec;
		generator.seed(rand_dev());        //reseed the generator or we will get the same values for all vectors
	}
	
	return rand_matrix;
}

RR GetVecNorm(Vec<ZZ>& vec){
 	
	RR sum = to_RR(0);
	RR norm;

	for(ZZ coord : vec){
		RR coordinate = conv<RR>(coord);
		sum += coordinate * coordinate;
	}	
	
	norm = sqrt(sum);
	return norm;
}

RR GetHadamardRatio(Mat<ZZ>& matrix){
	
	RR ratio;
	RR temp;
	RR n;
	RR det = conv<RR>(determinant(matrix));
	RR prod = to_RR(1);
	unsigned int dim = matrix.NumCols();

	for (int i=0; i < dim; i++){
		prod *= GetVecNorm(matrix[i]);
	}
	
	temp = abs(det/prod);
	n = to_RR(1)/to_RR(dim);
	
	ratio = pow(temp, n);
	return ratio;
}

Mat<ZZ> GetPrivKey(unsigned int dimension, unsigned int range,  float ratio){
	// Get random vectors 
	// Compute Hadamard ratio
	Mat<ZZ> random_M;
	RR computed_ratio;
	RR pre_lll_ratio;
	
	random_M.SetDims(dimension,dimension);
	
	
	
	while (computed_ratio < ratio){
		random_M =  GetRandVectors(dimension, dimension, range);
		BKZ_FP(random_M, 0.99, 30); // Sadly i must reduce if I want a higher ratio
		computed_ratio = GetHadamardRatio(random_M);
	}

	return random_M;	
}

int main(){
	
	auto start = chrono::high_resolution_clock::now();	
	
	unsigned int dimension = 20;
	unsigned int range     = 1000000;
	float hadamard_ratio   = 0.75;
	Mat<ZZ> Priv_key;
	
	Priv_key.SetDims(dimension, dimension);

	Priv_key = GetPrivKey(dimension, range, hadamard_ratio); 
	cout << Priv_key << "\n";

	auto end = chrono::high_resolution_clock::now();
	
	auto duration = chrono::duration_cast<chrono::seconds>(end - start);

    	cout << "Execution time: " << duration.count() << " seconds" << endl;

	return 0;
}
