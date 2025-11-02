#include <NTL/ZZ.h>
#include <NTL/vector.h>
#include <NTL/matrix.h>
#include <NTL/mat_ZZ.h>
#include <NTL/mat_RR.h>
#include <NTL/vec_RR.h>
#include <NTL/RR.h>
#include <NTL/LLL.h>
#include <cmath>
#include <iostream>
#include <cstdio>
#include <random>
#include <list>
#include <chrono>
#include <functional> 
using namespace std;
using namespace NTL;

random_device rand_dev;
mt19937       generator(rand_dev());

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

static inline RR DotRow(const Vec<RR>& a, const Vec<RR>& b){
	long n=a.length(); RR s=to_RR(0);
	for(long i=0;i<n;i++) s+=a[i]*b[i];
	return s;
}

static inline ZZ RoundRR(const RR& x){
	RR h=to_RR(0.5);
	return (x>=0? to_ZZ(floor(x+h)) : to_ZZ(ceil(x-h)));
}

static void GramSchmidtRows(const Mat<RR>& B, Mat<RR>& Bstar, Mat<RR>& mu, Vec<RR>& norm2){
	long n=B.NumRows(), m=B.NumCols();
	Bstar.SetDims(n,m); mu.SetDims(n,n); norm2.SetLength(n);
	for(long i=0;i<n;i++){
		Bstar[i]=B[i];
		for(long j=0;j<i;j++){
			mu[i][j]=DotRow(B[i],Bstar[j])/norm2[j];
			for(long k=0;k<m;k++) Bstar[i][k]-=mu[i][j]*Bstar[j][k];
		}
		norm2[i]=DotRow(Bstar[i],Bstar[i]);
	}
}

static Vec<ZZ> BabaiRows(const Mat<ZZ>& Bz, const Vec<ZZ>& e){
	Mat<RR> B; conv(B,Bz);
	Vec<RR> y; conv(y,e);
	long n=B.NumRows(), m=B.NumCols();
	Mat<RR> Bstar, mu; Vec<RR> n2;
	GramSchmidtRows(B,Bstar,mu,n2);
	Vec<ZZ> a; a.SetLength(n);
	for(long i=n-1;i>=0;i--){
		RR c=DotRow(y,Bstar[i])/n2[i];
		a[i]=RoundRR(c);
		for(long k=0;k<m;k++) y[k]-=to_RR(a[i])*B[i][k];
	}
	Vec<ZZ> v=a*Bz;
	return v;
}



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
	

	for(int i = 0; i < amount; i++){
		rand_vec = GetRandVec(dimension, range);
		rand_matrix[i] = rand_vec;
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

Mat<ZZ> GetBadMatrix(unsigned int n){
	Mat<ZZ> bad; 
	bad = GetIdentityMatrix(n);
	uniform_int_distribution<int> op(0, 2);
	uniform_int_distribution<int> idx(0, n - 1);
	uniform_int_distribution<int> factor(-17, 17);

	for(int t = 0; t < 20 * n; t++){
	
	int c = op(generator);
	int i = idx(generator);
	int j = idx(generator);
	
	if(i == j) continue;
	
	if(c == 0){
	    for(int k = 0; k < (int)n; k++) swap(bad[i][k], bad[j][k]);
	}

	else if(c == 1){
	    for(int k = 0; k < (int)n; k++) bad[i][k] = -bad[i][k];
		
	}	 

	else {
	    int kf = factor(generator);
	    if(kf == 0) continue;
	    for(int k = 0; k < (int)n; k++) bad[i][k] += bad[j][k] * kf;
	}

	}

	ZZ d = determinant(bad);
	if(!(d == 1 || d == -1)) return GetBadMatrix(n);

	return bad;
}


Mat<ZZ> GetPrivKey(unsigned int dimension, unsigned int range,  float ratio){
	
	Mat<ZZ> random_M;
	RR computed_ratio;
	RR pre_lll_ratio;
	float improv;

	random_M.SetDims(dimension, dimension);
	
	while (computed_ratio < ratio){
		random_M =  GetRandVectors(dimension, dimension, range);
		//BKZ_FP(random_M, 0.99, 30); // Sadly i must reduce if I want a higher ratio
		computed_ratio = GetHadamardRatio(random_M);
	}

	return random_M;	
}


Mat<ZZ> GetPublicKey(Mat<ZZ>& Priv_key){
	
	Mat<ZZ> Public_key;
	Mat<ZZ> M;
	unsigned int dimension = Priv_key.NumCols(); 
	
	M = GetBadMatrix(dimension);
	mul(Public_key, M, Priv_key); 

	return Public_key;
}


Vec<ZZ> EncryptGGH(Mat<ZZ> Public_key, Vec<ZZ> plain_text, unsigned int delta){
	Vec<ZZ> cipher_text;
	Vec<ZZ> ephemeral_key;
	ZZ 	dim;

	dim = conv<int>(Public_key.NumCols());

	ephemeral_key = GetRandVec(conv<int>(dim), delta);
	
	cipher_text = (plain_text * Public_key) + ephemeral_key;

	return cipher_text;
	
}




Vec<ZZ> DecryptGGH(Mat<ZZ> Priv_key, Mat<ZZ> Public_key, Vec<ZZ> cipher_text){
	Vec<ZZ> v=BabaiRows(Priv_key,cipher_text);
	Mat<RR> W; 
	conv(W,Public_key);
	Mat<RR> Winv; 
	inv(Winv,W);
	Vec<RR> vr; 
	conv(vr,v);
	Vec<RR> mr=vr*Winv;
	Vec<ZZ> m; 
	m.SetLength(mr.length());
	for(long i=0;i<m.length();i++) m[i]=RoundRR(mr[i]);
	return m;
}

int main(){
	auto start = chrono::high_resolution_clock::now();	

	Mat<ZZ> M;	
	unsigned int dimension = 1000;
	unsigned int range     = 10000;
	float hadamard_ratio   = 0.50;
	Mat<ZZ> Priv_key;
	Mat<ZZ> Public_key;
	Vec<ZZ> pt;
	Vec<ZZ> ct;
	Vec<ZZ> new_pt;

	RR good_ratio;
	RR bad_ratio;
	RR post_LLL_ratio;

	Priv_key.SetDims(dimension, dimension);

	Priv_key = GetPrivKey(dimension, range, hadamard_ratio); 
	cout << "Priv key acquired" << endl;
	
	Public_key = GetPublicKey(Priv_key);	
	cout << "Public key acquired" << endl;

	//bad_ratio  = GetHadamardRatio(Public_key);
	//good_ratio = GetHadamardRatio(Priv_key);
	
	//printf("Private Key Ratio: %f\n", conv<float>(good_ratio));
	//printf("Public Key Ratio : %f\n", conv<float>(bad_ratio));
	
	pt = GetRandVec(Public_key.NumCols(), 100);
	cout << "plaintext acquired" << endl;

	ct = EncryptGGH(Public_key, pt, 5);
	cout << "ciphertext acquired" << endl;

	//cout << ct << endl;
	new_pt = DecryptGGH(Priv_key, Public_key, ct);
	
	if(new_pt == pt){
		printf("plaintexts match\n");
	}
	auto end = chrono::high_resolution_clock::now();

	auto duration = chrono::duration_cast<chrono::seconds>(end - start);
    	cout << "Execution time: " << duration.count() << " seconds" << endl;
	return 0;
}
