//Here begins the file "main.cpp"
#include <ctime>
#include <cstdlib>
#include <math.h>
#include <fstream>
#include <chrono>
#include <iostream>

using namespace std;

//Size of the lattice
const int N=16;

//Taking Planck unit system, set Boltzmann constant
const float k_B=1;

//System Hamiltonian is in units of J, the exchange energy
const float J=1;

//Lattice properties+experimental conditions
class Lattice{
	public:
	int arr[N][N];	//Array to describe N by N lattice
	float tf;	//temperature we experiment on, as a factor of J/k_B
	
	float muH;	//product of magnetic moment of site with the applied magnetic field strength
	//...therefore the magnetic field is in units of J/mu
	
	float energy;	//Total energy of lattice
	int mag;	//Sum of magnetisation spins in lattice
};

//returns a 1D position along lattice
float randPos(){
	return  ((float)N*rand() / RAND_MAX);
}

//Take array to describe N by N lattice
void initUps(Lattice &lattice){
	
	//Choose starting spins as all +1 in lattice
	for(int i=0;i<N;i++){
		for(int j=0;j<N;j++){
			lattice.arr[i][j]=1;
		}

	}
	
	//Simple calculations for sum of all magnetisations and energies
	lattice.mag=N*N;
	lattice.energy=-2*J*N*N-lattice.muH*N*N;
}

void initRand(Lattice &lattice){
	//Choose starting spins as all random in lattice, with equal probability up and down
	for(int i=0;i<N;i++){
		for(int j=0;j<N;j++){
			if (((float)rand()/RAND_MAX)<0.5){
				lattice.arr[i][j]=1;
			}else{
				lattice.arr[i][j]=-1;
		}
	}
	int count=0;
	float interenergy=0;
	for(int i=0;i<N;i++){
		for(int j=0;j<N;j++){
			//Sum interaction energy terms
			interenergy+=-J*(lattice.arr[i][j])*(lattice.arr[(i+1)%N][j]+lattice.arr[i][(j+1)%N]);
			
			//Sum magnetisations
			count+=lattice.arr[i][j];
		}
	}
	lattice.energy=interenergy-lattice.muH*count;
	lattice.mag=count;
}
}

//Conditions to decide whether lattice site should flip
bool flipCond(float deltaE,float T){
	if(deltaE<0){
		return true;
	}else if(exp(-deltaE/k_B/T)>((float)rand()/RAND_MAX)){
		return true;
	}else{
		return false;
	}
}

//Function to flip one lattice site
void flip(Lattice &lattice){
	//pick random postion
	int x=(int)randPos()%N;
	int y=(int)randPos()%N;
	float self=lattice.arr[x][y];
	
	//Modulo for boundary conditions
	float neighbours=lattice.arr[(x+1)%N][y]+lattice.arr[x][(y+1)%N]+lattice.arr[(x-1+N)%N][y]+lattice.arr[x][(y-1+N)%N];
	
	//Only nearest neightbours relevant to energy change
	float deltaE=2*J*self*neighbours+2*lattice.muH*self;
	
	//Temperature given in terms of J/k_B (start near the order of room temp)
	float T=lattice.tf*J/k_B;
		
	if(flipCond(deltaE,T)){
		//update lattice energy and magnetisation
		lattice.energy=lattice.energy+deltaE;
		lattice.arr[x][y]=-lattice.arr[x][y];
		lattice.mag=lattice.mag+2*lattice.arr[x][y];
	}
}

//1 time step to determine flipping for N^2 random points
void timestep(Lattice (&lattice)){
	for(int i=1;i<=(N*N);i++){
		flip(lattice);		
	}
}

////////////////////////////////////////////////////////////////////
/*Everything before this point describes the core functions to run the simulations.
These could in fact be placed under a separate header file,
and exported for future projects.
The functions after this point test certain tasks under specific experimental conditions.
*/
////////////////////////////////////////////////////////////////////


//Print lattice as a matrix for all times
void Task1Draw(){
	Lattice lattice1;
	lattice1.muH=0;
	ofstream fout;
	
	//File name describes initial configuration 
	fout.open("drawRand.dat");
	
	for(int f=0; f<4;f++){
	initRand(lattice1);
	
	//Repeat over range of temperatures
	lattice1.tf=1+f*1.5;
		
	for(int a=1;a<=(1000);a++){
		for(int i=0;i<N;i++){
			for(int j=0;j<N;j++){
				fout<<lattice1.arr[i][j]<<" ";
			}
		fout<<"\n";
		}
		timestep(lattice1);
	}
	
	fout.close();
}
}

//Find how energy evolves with time to measure when it equilibrates
void Task1eq(){
	Lattice lattice1;
	lattice1.muH=0;
	ofstream fout;
	fout.open("Task1eq.dat");
	for(int i=0;i<5;i++){
		initRand(lattice1);
		
		//Repeat over range of temperatures
		lattice1.tf=0.5+i;
		for(int steps=0; steps<=1000;steps++){
			fout<<lattice1.energy<<" ";
			timestep(lattice1);
		}
		fout<<"\n";
	}
	fout.close();
}

//Prints total magnetisations after reaching equilibrium, to investigate magnetisation autocorrelation
void Task2(){
	Lattice lattice1;
	lattice1.muH=0;

	ofstream fout;
	fout.open("N15autocorr.dat");
	for(int f=0;f<50;f++){
		lattice1.tf=1.9+0.02*f ;
		initRand(lattice1);
		
		//evolve to equilibrium
		for(int a=0;a<1000;a++){
			timestep(lattice1);
		}
		for(int a=1;a<=5000;a++){
			timestep(lattice1);
			fout<<lattice1.mag<<" ";
		}
		fout<<"\n";
		
	}
	fout.close();
}

//Prints mean magnetisation at different temperatures
void Task3(){
	Lattice lattice1;
	lattice1.muH=0;
	
	//number of times over which we average 
	int longnumber=1000;
	ofstream fout;
	fout.open("Task3N128.dat");
	initRand(lattice1);
	for(int f=0;f<50;f++){
		//Repeat over range of temperatures
		lattice1.tf=5-0.1*f;
		
		//Variable to sum the magnetisations 
		int count=0;
		
		//evolve to equilibrium
		for(int a=0;a<500;a++){
			timestep(lattice1);
		}
		for(int a=1;a<=longnumber;a++){
			timestep(lattice1);
			count+=lattice1.mag;
		}
		//Calculate average magnetisation
		fout<<lattice1.tf<<"\t"<<(float)count/longnumber<<"\n";
		
	}
	fout.close();
}

//Heat capacity calculations from energy fluctuations
void Task4(){
	Lattice lattice1;
	lattice1.muH=0;

	ofstream fout;
	fout.open("Heatcapacity8.dat");
	for(int f=0;f<300;f++){
		//Repeat over range of temperatures
		lattice1.tf=0.8+0.05*f ;
		initRand(lattice1);
		
		//evolve to equilibrium
		for(int a=0;a<2000;a++){
			timestep(lattice1);
		}
		
		//Number of iterations over which we take 
		int iter=5000;
		
		//Variables to sum energies and energies squared
		float e=0 , e2=0;
		for(int a=1;a<=iter;a++){
			timestep(lattice1);
			e+=lattice1.energy;
			e2+=lattice1.energy*lattice1.energy;						
		}
		
		//Use values to calculate variance and subsequently the specific heat per lattice site
		fout<<lattice1.tf<<" "<<(e2/iter-(e/iter)*(e/iter))/((lattice1.tf*J)*(lattice1.tf*J/k_B)*N*N)<<"\n";
		
	}
	
}

//Run through a hysteresis loop, taking temperature and filename as argument
void hysteresis(float tf,string filename){
	Lattice lattice1;
	lattice1.tf=tf;
	
	//number of times over which we average 
	int longnumber=3000;
	ofstream fout;

	fout.open(filename.c_str());
	initRand(lattice1);
	lattice1.muH=-0.8;
	//evolve to equilibrium
		for(int a=0;a<1000;a++){
			timestep(lattice1);
		}
	
	for(int f=0;f<160;f++){
		//slowly increase H magnetic field from -0.8 to +0.8
		lattice1.muH=-0.8+0.01*f;
		
		//Variable to sum the magnetisations
		int count=0;
				
		for(int a=1;a<=longnumber;a++){
			timestep(lattice1);

			count+=lattice1.mag;
		}
		//Calculate average magnetisation
		fout<<lattice1.muH<<"\t"<<(float)count/longnumber<<"\n";
		
	}
	for(int f=0;f<160;f++){
		//turn H field back down to -0.8
		lattice1.muH=0.8-0.01*f;
		
		//Variable to sum the magnetisations
		int count=0;
				
		for(int a=1;a<=longnumber;a++){
			timestep(lattice1);
			count+=lattice1.mag;
		}
		//Calculate average magnetisation
		fout<<lattice1.muH<<"\t"<<(float)count/longnumber<<"\n";
		
	}
	fout.close();
}

//Produce necessary files for a 32 by 32 lattice
void Task6(){
	hysteresis(1,"hysteresismtf1000.dat");
	hysteresis(2,"hysteresismtf2000.dat");
	hysteresis(2.8,"hysteresismtf2800.dat");
	hysteresis(2.4,"hysteresismtf2400.dat");
	hysteresis(3.8,"hysteresismtf3800.dat");
}

//Basic function for timing the timesteps in lattice under certain conditions
void speedTest(){
	float totTime=0;
	int numLoops=50000;
	for(int i=0;i<numLoops;i++){
		// Use auto keyword to avoid typing long
		// type definitions to get the timepoint
		// at this instant use function now()
		auto start = chrono::high_resolution_clock::now();

		Lattice testLat;
		testLat.tf=3;
		testLat.muH=0;
		initUps(testLat);
		for(int i=0;i<80000;i++){
			timestep(testLat);
		}

		// After function call
		auto stop = chrono::high_resolution_clock::now();

		// Subtract stop and start timepoints and
		// cast it to required unit. Predefined units
		// are nanoseconds, microseconds, milliseconds,
		// seconds, minutes, hours. Use duration_cast()
		// function.
		auto duration = chrono::duration_cast<chrono::microseconds>(stop - start);
		
		// To get the value of duration use the count()
		// member function on the duration object
		float looptime=float(duration.count())/1000000;
		cout <<"Loop "<<i+1<<" executed in "<< looptime << " s"<<endl;
		totTime+=looptime;
	}
	cout<<numLoops<<" loops executed with an average time of "<<totTime/numLoops<<" s";
}

//susceptibility calculations from magnetisation fluctuations
void susceptibility(){
	Lattice lattice1;
	lattice1.muH=0;

	ofstream fout;
	fout.open("susceptibility12.dat");
	initRand(lattice1);
	for(int f=1;f<120;f++){
		//Repeat over a range of temperatures
		lattice1.tf=0.1+0.05*f ;
		
		//evolve to equilibrium
		for(int a=0;a<500;a++){
			timestep(lattice1);
		}
		
		//Number of iterations over which we take 
		int iter=2000;
		
		//Variables to sum magnetisations and magnetisation squared
		float m=0 , m2=0;
		for(int a=1;a<=iter;a++){
			timestep(lattice1);
			m+=lattice1.mag;
			m2+=lattice1.mag*lattice1.mag;						
		}
		
		//Use values to calculate variance and subsequently the susceptibility per lattice site
		fout<<lattice1.tf<<"\t"<<(m2/iter-(m/iter)*(m/iter))/(lattice1.tf*N*N)<<"\n";
		
	}
	
}

int main() {
	//Seed random number generator with time 0
	srand((unsigned) time(NULL));
	
	//Pick task to run
	//Task1Draw();
	//Task1eq();
	//Task3();
	// Task2();
	// Task4();
	//Task6();
	speedTest();
	//susceptibility();
}