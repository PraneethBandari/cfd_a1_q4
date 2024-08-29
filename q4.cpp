
#include <stdio.h>
#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>

using namespace std;


int main()
{
	double L = 0.13;          // Wall thickness (m)
	int N = 50;              // Number of nodes
	double dx = L / (N - 1); // Step size (m)
	double T_left = 600;
	double T_right = 400;
	//Plasterboard
	double K1 = 0.15;
	double L1 = 0.01;
	//Fiberglass
	double K2 = 0.038;
	double L2 = 0.11;//distance from starting
	//Plywood
	double K3 = 0.1;
	double L3 = 0.13;
	double accuracy = 1e-6;
	int n=0;//iteration count
	vector<double> x;
	for(int i=0; i<N; i++) {
		x.push_back(dx*i);
	}
	vector<double> K;
	//Analytical solution + assaigning values to k(heat conductivity)
	vector<double> T_analytical;
	double T1=T_left-(L1/K1)*((T_left-T_right)/((L1/K1)+(L2/K2)+(L3/K3)));
	double T2 = T_right+(L3/K3)*((T_left-T_right)/((L1/K1)+(L2/K2)+(L3/K3)));
	cout<<T2;
	for(int i=0; i<N; i++) {
		double T;
		if(x[i]<=L1) {
			T=T_left + (T1-T_left)*(x[i]/L1);
			T_analytical.push_back(T);
			K.push_back(K1);
		} else if(x[i]<=L2) {
			T=T1 + (T2-T1)*((x[i]-L1)/(L2-L1));
			T_analytical.push_back(T);
			K.push_back(K2);
		} else if(x[i]<=L3) {
			T=T2+(T_right-T2)*((x[i]-L2)/(L3-L2));
			T_analytical.push_back(T);
			K.push_back(K3);
		}
		cout << "Analytical temperature at node " << i + 1 << " is " << T_analytical[i] << "\n";

	}
	//fdm solution
	vector<double>Temp;
	for(int i=0; i<N; i++) {
		double T;
		T=(T_left+T_right)*0.5;
		Temp.push_back(T);
	}
	vector<double>T_new(N,0);
	T_new[0]=T_left;
	T_new[N-1]=T_right;
	double max_dt=1;
	while(max_dt>accuracy) {
		max_dt=0;
		for(int i=1; i<(N-1); i++) {
			T_new[i]=(K[i-1]*Temp[i-1]+K[i+1]*Temp[i+1])/(K[i+1]+K[i-1]);
		}
		// Check for convergence and update temperatures
		for (int i = 0; i < N; i++) {
			double dt = fabs(Temp[i] - T_new[i]);
			if (dt > max_dt) {
				max_dt = dt;
			}
			Temp[i] = T_new[i];
		}
		n++;
	}
	// Output the number of iterations and final temperatures
	cout << "Number of iterations: " << n << "\n";
	for (int i = 0; i < N; i++) {
		cout << "Node " << i + 1 << ": Analytical temp = " << T_analytical[i] << " FDM temp = " << Temp[i] << " \n";
	}
	ofstream myfile("fdmq4sol.txt");

	if (myfile.is_open()) {

		for (int i = 0; i < N; i++) {
			myfile<< i + 1 << ":" << Temp[i] << endl;
		}
		myfile.close();
		cout << "Temperature values have been written to fdmq4sol.txt" << endl;
	} else {
		cout << "Error opening file!" << endl;
	}
	ofstream outfile("fdmq4analyticalsol.txt");

	if (outfile.is_open()) {

		for (int i = 0; i < N; i++) {
			outfile << i + 1 << ":" << T_analytical[i] << endl;
		}
		outfile.close();
		cout << "Temperature values have been written to fdmq4analyticalsol.txt" << endl;
	} else {
		cout << "Error opening file!" << endl;
	}



	return 0;
}