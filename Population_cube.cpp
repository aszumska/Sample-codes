/*
 *  Ising type code to equilibrate MP13 lattice - Davide's project
 *  Created by Valerie on December 2014.
 *  All rights reserved.
 *
 */

#include <iostream>
#include <cmath>
#include <string>
#include "armadillo"
#include <fstream>
#include <stdio.h>
#include <math.h>
#include <vector>
using namespace arma;
using namespace std;

// Energy unit is eV, distance is Angstrom and the rest is in SI 
#define PI 3.14159265
#define h_bar 6.58211928e-16	//eV.s 
#define K_B 8.6173324e-5		//eV.K-1
#define T 300					//K
#define THRESHOLD 1e-5


#if (ARMA_VERSION_MAJOR > 0 || (ARMA_VERSION_MINOR >= 9 && ARMA_VERSION_PATCH >= 50))
#define RAND randu
#else
#define RAND rand
#endif

struct lattice_sites {
	int line;
	int col;
	int molecule;
	//int charge;
	int dihedral;
};


double coverage(double& concentration, double& binding_constant, int& n) {
	return n*n*concentration/(concentration+1/binding_constant);		//nb of molecules per square lattice param
}


int main ()
{
//Declaration of the parameters
	int nbNN=8;
    	double S=1;
	double Energy_binding_ti=0, Binding_sum=0;
	double A=1;
	double binding_a, binding_b;
	int n;
	double cov;
	char comment[5];
    	int count=0;
	double c=2e-4;											
	double kf=10e3;
	int nb_f=6;
	ifstream file1;
	int size_of_energy_file = 32;
	vector<double> Energies(size_of_energy_file);
	ofstream outfile,outfile2;
	srand(time(0));
	//Input parameters	
	cin >> comment;  //structure iterator - how you want to name output structure
	cin >> n;	//size of the facet
	cin >> S;	// S & A are probabilities that dye with diffuse to the facet (for 2 diff configurations)
	cin >> A;
	cin >> cov;	// coverage
	cin >> binding_a;	// binding energy is equal binding_a * <E_i> + binding_b*k_B*T;
	cin >> binding_b;
	string filetext="/home/aas215/Population_cube/Lattices/LatticePop_MP13_";
        string filetext2="/home/aas215/Population_cube/Lattices/MP13_Lattice_toWalk_";
        string filename, filename2;
        filename.append(filetext);
        filename.append(comment);
        filename.append(".txt");
        filename2.append(filetext2);
        filename2.append(comment);
        filename2.append(".txt");

//Read in the txt file where energies are stored 

	file1.open("/home/aas215/NBMP13energies.txt");
	if (file1.is_open()) {
                for (unsigned int i=0; file1; i++ ) {
                        file1 >> Energies[i];
			Binding_sum+=Energies[i];
                }
                file1.close();
		Energy_binding_ti=binding_a *Binding_sum/(4*nbNN)+binding_b*4*K_B*T/(27.2);
        } else {
                cerr << "Could not open file 1" << endl;
                return 1;
        }

	
	outfile.open(filename2.c_str());
	outfile2.open(filename.c_str());
//Creation of the lattice sites	
for (int f=0; f<nb_f; f++){
	vector <lattice_sites> sites(n*n, lattice_sites() );
		
//Add molecules on the sites	
		
	//double cov=coverage(c,kf,n)/double(n*n);					//Percentage of molecules to put on the n*n sites 
	//double cov=1;

		
    for (int i=0; i<n*n; i++){
        sites[i].col=i%n;
        sites[i].line=floor(i/n);
    }

		mat NN = zeros<mat>(n*n,nbNN);

                for (int m=0; m<=n*n-n-1; m++){
                   
                        NN(m,6)=m+n;
                        NN(m+n,2)=m;
                    

                    if ((m+1)%n!=0){
         
							NN(m,4)=m+1;
                            NN(m+1,0)=m;
							NN(m+n,3)=m+1;
                            NN(m+1,7)=m+n;
							NN(m,5)=m+n+1;
                            NN(m+n+1,1)=m;
							if (m<n){

								NN(m,3)=n*n-n+m%n+1;
								NN(n*n-n+m%n+1,7)=m;
							}
							if (m%n==0){
								NN(m,7)=m+2*n-1;
							}
						
                        
                    }
					

                    else {
                 
							NN(m,4)=m-n+1;
                            NN(m-n+1,0)=m;
							NN(m+n,3)=m-n+1;
                            NN(m+1,7)=m+n;
							NN(m,5)=m+1;
							NN(m+1,1)=m;
							
							if (m+1==n){
								NN(m,3)=n*n-n;
								NN(n*n-n,7)=m;
							}
                        
                    }

                } //end of for loop on m sites
                                
                for (int m=n*n-n; m<n*n; m++){ //Last row
                    if ((m+1)%n == 0){
							NN(m,4)=m-n+1;
                            NN(m-n+1,0)=m;
							NN(m,3)=m-2*n+1;
							NN(m,5)=m%n-n+1;
                            NN(m%n-n+1,1)=m;
                        
                    }	
					
					
                    else{
							NN(m,4)=m+1;
                            NN(m+1,0)=m;
							NN(m,3)=m-n+1;
                            NN(m-n+1,7)=m;
							NN(m,5)=m%n+1;
                            NN(m%n+1,1)=m;
							if (m%n==0){
								NN(m,7)=m%n+n-1;
							}
                        
                    }
                    
						NN(m,6)= m%n; // down of last row matrix is first row
                        NN(m%n,2)=m; // up first rwo is last row
                    
                }
		
    double dih=0.5; //ratio of dihedral -20 and 160
    int i=0;
    int nbMolecule=n*n*cov;

while (i !=nbMolecule){ // Redefine loop condition for site coverages
    
    int site=rand()%(n*n);
    while (sites[site].molecule!=0){
        site=rand()%(n*n);
    }
    
    double u_2=rand()/double(RAND_MAX);
    double dihedral=0;
    if (u_2<=dih){
        dihedral=-20;
    }else{
        dihedral=160;
    }
	int NNdihedral[2]={0,0};
//Calculating number of NN with different configurations NNdihedral[0] - the same dihedral as our, NNdihedral[1] - opposite
	for(int k=0; k<nbNN; k++){
		if(sites[NN(site,k)].dihedral==dihedral){
			NNdihedral[0]++;}
		else{ 
			if(sites[NN(site,k)].dihedral==0){;}
			else{NNdihedral[1]++;}
		}
	}
	double u_3=rand()/double(RAND_MAX);
	double Probability_dif=0;    // Probability dye will diffuse to the site
	Probability_dif=pow(S,NNdihedral[0])*pow(A,NNdihedral[1]);
//	cout << NNdihedral[0] <<"\t"<<NNdihedral[1]<<"\t"<<Probability_dif<<endl;
	if (u_3<=Probability_dif){
//Calculating energies for potentially new system with dye (specific dihedral) sticked to the pick side 

        double Energy_dye_site=0;
        vector<double> E_new(8);

        for (int m=0; m<8; m++){
        if (dihedral==-20){
                if (sites[NN(site,m)].dihedral==-20){
                        E_new[m]=Energies[m]-Energy_binding_ti;
                }else if (sites[NN(site,m)].dihedral==160){
                        E_new[m]=Energies[m+nbNN]-Energy_binding_ti;
                }
        }else if (dihedral==160){
            if (sites[NN(site,m)].dihedral==-20){
                E_new[m]=Energies[m+2*nbNN]-Energy_binding_ti;
            }else if (sites[NN(site,m)].dihedral==160){
                E_new[m]=Energies[m+3*nbNN]-Energy_binding_ti;
            }
        }else {E_new[m]=0;}
                Energy_dye_site+=0.5*E_new[m];
        }

// Delta_E
	double dE=Energy_dye_site;
	if(dE<0){
	sites[site].molecule=1;
            sites[site].dihedral=dihedral;
                        i++;
	}
	else
	{
	double R=rand()/double(RAND_MAX);
	double boltzmann=0;
	boltzmann=exp(-dE*27.2/K_B/T);
		if(R<boltzmann){
	        sites[site].molecule=1;
    	    sites[site].dihedral=dihedral;
//		cout<< i << '\t'<< (NNdihedral[0]+NNdihedral[1]) <<'\t'<< dE << '\t' << boltzmann <<endl;
			i++;
		}
	}
}
}
// Calutating avarage number of NN having the same dihedral
        double average_1, average_2, site, dihedral;
        double Sum[2]={0,0};
        double Dih[2]={0,0};
        for(int i=0; i<n*n; i++){
                int m=0;
                site=i;
                dihedral=sites[i].dihedral;
                if(dihedral==160){m=1;};
                Dih[m]++;
                for(int k=0; k<nbNN; k++){
                        if(sites[NN(site,k)].dihedral==dihedral){
                                Sum[m]++;
                        }
                }
        }
        average_1=Sum[0]/(nbNN*Dih[0]); //average for dihedral -20
        average_2=Sum[1]/(nbNN*Dih[1]);  //average for dihedral 160
        double ave=(Sum[0]+Sum[1])/(nbNN*(Dih[0]+Dih[1]));

	if(f==0){
	        outfile2<<comment<<"\t"<< n <<"\t"<<S <<"\t" << A << "\t" <<cov << "\t" << nb_f << "\t"<< binding_a << "\t" << binding_b << "\t" << average_1 << "\t" << average_2 << "\t" << ave << endl;}
	for (int y=0; y<n*n; y++){
		outfile2 << sites[y].line << "\t" << sites[y].col << "\t" << sites[y].dihedral <<  endl;
		outfile << sites[y].molecule << "\t" << sites[y].dihedral << endl;
	}
//	outfile.close();
}
outfile2.close();
outfile.close();
	

	return 0;
}
