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
	return 0;
}
