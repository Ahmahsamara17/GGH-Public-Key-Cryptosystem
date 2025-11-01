#include <NTL/ZZ.h>
#include <NTL/vector.h>
#include <NTL/matrix.h>
#include <NTL/mat_ZZ.h>
#include <NTL/vec_RR.h>
#include <NTL/RR.h>
#include <NTL/LLL.h>
#include <cmath>
#include <iostream>
#include <cstdio>
#include <random>
#include <list>
#include <chrono>
using namespace std;
using namespace NTL;

random_device rand_dev;
mt19937       generator(rand_dev());
bool verbose = false;

// Private Key Gen:
// I want to create a function that given a dimension gives me a unit vector with that many dimensions
// I need a function that will generate n random values from a range [-d, d]  
// I need a function that will generate an n dimensional vector with random values from a range
// I need a funciton that will generate k amount of n dimesnsional vectors with random values from a range
// I need a function that will compute the Hadamard ratio of a list of vectors.
// I need a function that will return a basis with a Hadamard ratio above a certain value

// Public Key Gen:
//

/*
 * Function that returns an n dimensional vector with coordinates initialized to 1 
*/


Mat<ZZ> GetIdentityMatrix(unsigned int n){
	
	Mat<ZZ> I_M;

	I_M.SetDims(n, n);

	for(int i=0; i < n; i++){
		I_M[i][i] = ZZ(1);
	}
	
	return I_M;
}

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

	rand_matrix.SetDims(amount, dimension);
	
	if(verbose){
		printf("Creating a random basis with size %d and range [%d, %d]\n", dimension, -range, range);
		
		for(int i = 0; i < amount; i++){
			rand_vec = GetRandVec(dimension, range);
			printf("This is random vector #%d:\n", i);
			cout << rand_vec << "\n\n";
			rand_matrix[i] = rand_vec;
		}
	}

	else{
		for(int i = 0; i < amount; i++){
			rand_vec = GetRandVec(dimension, range);
			rand_matrix[i] = rand_vec;
		}
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


Mat<ZZ> CreateRowSwapMatrix(unsigned int n, int index1, int index2){
	
	Mat<ZZ> unit_M;
	
	unit_M = GetIdentityMatrix(n);

	if(index1 >= 0 && index1 < n && index2 >= 0 && index2 < n){
		swap(unit_M[index1], unit_M[index2]);
	}

	return unit_M;
}

Mat<ZZ> CreateRowScalingMatrix(unsigned int n, int index, int k){
	
	Mat<ZZ> unit_M;
	
	unit_M = GetIdentityMatrix(n);

	if(index >= 0 && index < n && k != 0){
		unit_M[index][index] = k;
	}

	return unit_M;
}

Mat<ZZ> CreateRowAdditionMatrix(unsigned int n, int index1, int index2, int k){
	
	Mat<ZZ> unit_M;
	
	unit_M = GetIdentityMatrix(n);

	if(index1 >= 0 && index1 < n && index2 >= 0 && index2 < n && index1 != index2){
		unit_M[index2][index1] = k;
	}

	return unit_M;
}

Mat<ZZ> GetElementaryMatrix(unsigned int n){
	
	Mat<ZZ> M;
 /*
  * I want to make a list of these elementary functions 
  * then pick a random amount and then continue to randomly pick funcs
  *
  *
  */
	uniform_int_distribution<int> dist(0, n-1);
	uniform_int_distribution<int> dist2(-10000, 10000);
	
	int rand1 = dist(generator);
	int rand2 = dist(generator);
	int rand3 = dist2(generator);

	printf("index1: %d, index2: %d, k: %d\n", rand1, rand2, rand3);
	M = CreateRowAdditionMatrix(n, rand1, rand2, rand3); 
	
	cout << M << endl;
	
	return M;
}

Mat<ZZ> GetBadMatrix(unsigned int dim){

	Mat<ZZ> bad_matrix;

	return bad_matrix;
}


Mat<ZZ> GetPrivKey(unsigned int dimension, unsigned int range,  float ratio){
	
	Mat<ZZ> random_M;
	RR computed_ratio;
	RR pre_lll_ratio;
	float improv;

	random_M.SetDims(dimension,dimension);
	
	if (verbose){
		
		printf("Creating a random square basis' of size %d and of range [%d, %d]\n\n", dimension, (-range), (range));	
		float best_improv = 0.0f;

		while (computed_ratio < ratio){
			random_M =  GetRandVectors(dimension, dimension, range);
			
			pre_lll_ratio = GetHadamardRatio(random_M);
			printf("This is the pre LLL ratio	: %f\n",conv<float>(pre_lll_ratio));
			
			LLL_FP(random_M, 0.99); // Sadly i must reduce if I want a higher ratio
			computed_ratio = GetHadamardRatio(random_M);
			
			printf("This is the post LLL ratio	: %f\n", conv<float>(computed_ratio));
			
			improv = conv<float>(computed_ratio - pre_lll_ratio);
			printf("improvment by %f\n\n", conv<float>(improv));

			if(improv > best_improv){
				best_improv = improv;
			}
		}

		printf("Best improvment was: %f\n\n", best_improv);
	}

	else{
	
		while (computed_ratio < ratio){
			random_M =  GetRandVectors(dimension, dimension, range);
			BKZ_FP(random_M, 0.99, 30); // Sadly i must reduce if I want a higher ratio
			computed_ratio = GetHadamardRatio(random_M);
		}
	}
	return random_M;	
}

int main(){
	
	GetElementaryMatrix(5);
	return 0;
	
}
