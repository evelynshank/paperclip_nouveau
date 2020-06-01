#include <iostream>
using namespace std;

int main(){
	int numElements = 2664;
	
	double theta[numElements];
	double phi[numElements];
	double r[numElements];

	cout << "Enter r: "<< endl;

	for(int i=0; i< numElements; i++){
		cin >> r[i];
	}
	
// Populate the theta array
	for (int j = 0; j < 72; j++){
		for (int k=0; k<37; k++){
			theta[j*37+k]=5*k*3.14159/180;
		}
	}

// Populate the phi array
	for (int l = 0; l < 72; l++){
		for (int m=0; m < 37; m++){
			phi[l*37+m]=5*l*3.14159/180;
		}
	}

	for(int n=0; n< numElements; n++){
		cout <<"{"<<theta[n]<< ", "<<phi[n]<<", "<<r[n]<<"}, ";
	}
	
	return 0;
}
