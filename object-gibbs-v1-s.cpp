#include <utility>
#include <iostream>
#include <algorithm>
#include <cmath>  
#include <vector>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <ctime>
using namespace std;
struct spin{
    int val;
	int nb[4];
};

class IsingModel{
	public:
		std::vector<spin> matrix;
		double flipRatio[17][3];
		int L;
		double J,H,T;
		/*double posSum, negSum;
		int replica;
		int replicaCount[3];*/
		int magnet;
	public:
		void createSpin(int i, int j, int value);
		void initializeCongiguration(int seed, double thresh);
		void initializeFlipRatio();
		/*void initializeSum();
		void updateSum(int current);*/
		int nbSum(int current);
		double MHRatio(int current);
		//double switchRatio();
		int changeState(int loop);
		IsingModel(int N,double p1,double p2,double p3){
			matrix.reserve(N*N);
			L = N;
			J = p1;
			H = p2;
			T = p3;
			/*posSum=0;
			negSum=0;*/
			magnet=0;
		}
};
void IsingModel::createSpin(int i, int j, int value){
	matrix[i*L+j].val=value;
	/*int nb[0]= matrix[current].i ==0?L-1: i-1;
	int nb[1]= matrix[current].i ==L-1?:0:i+1;
	int nb[2]= matrix[current].j ==0?L-1: j-1;
	int nb[3]= matrix[current].j ==L-1?0: j+1;*/
	matrix[i*L+j].nb[0] = i == 0 ? L*(L-1)+j : (i-1)*L+j;
	matrix[i*L+j].nb[1] = i == L-1 ? j : (i+1)*L+j;
	matrix[i*L+j].nb[2] = j == 0 ? i*L+L-1 : i*L+j-1;
	matrix[i*L+j].nb[3] = j == L-1 ? 0 : i*L+j+1;
	//cout <<i*L+j<<endl;
	//cout <<matrix[i*L+j].nb[0] <<" "<<matrix[i*L+j].nb[1]<<" "<<matrix[i*L+j].nb[2]<< " "<< matrix[i*L+j].nb[3]<<endl;
}
void IsingModel::initializeCongiguration(int seed,double thresh){
	std::srand(seed);
	//cout << L <<J << H <<T << endl;
	for(int i=0; i<L; ++i){
		for(int j=0; j<L; ++j){
			double p= (double)rand()/(double)RAND_MAX;
			if(p<thresh){
				createSpin(i,j,1);
				magnet+=1;
				//cout<< p<< " " << magnet <<endl;
			}
			else{
				createSpin(i,j,-1);
				magnet+=-1;
				//cout<< p<< " " << magnet <<endl;
			}
			//cout << i << " " <<j <<endl;
		}
	}
	/*if(((double)rand()/(double)RAND_MAX)<1/2)
		replica = 1;
	else
		replica = -1;*/
}
void IsingModel::initializeFlipRatio(){
	for (int i = -8; i <= 8; i += 4) {
		flipRatio[i + 8][0] = 1<exp( (i * J + 2 * H) / T)?1:exp( (i * J + 2 * H) / T);
		flipRatio[i + 8][2] = 1<exp( (i * J - 2 * H) / T)?1:exp( (i * J - 2 * H) / T);
	}
}
int IsingModel::nbSum(int current){
	int total=0;
	for(int i =0; i<4; i++){
		int j= matrix[current].nb[i];
		total+=matrix[j].val;
	}
	return total;
}
double IsingModel::MHRatio(int current){
	int delta_ss =-2*nbSum(current)*matrix[current].val;
	return flipRatio[delta_ss+8][matrix[current].val+1];
}


int IsingModel::changeState(int loop){
	//std::srand(seed);	
	//srand(time(NULL));
	while(loop){
	int x = rand()%L;
	int y = rand()%L;
	int proposed=x*L+y;
	//cout << x << " "<< y <<" "<<matrix[proposed].val<<" " <<replica<<endl;
	//if(matrix[proposed].val==replica){
		double p=((double)rand()/(double)RAND_MAX);
		if(p<MHRatio(proposed)){
			cout << p << " " << MHRatio(proposed)<<endl;
			matrix[proposed].val=-matrix[proposed].val;
			/*replicaCount[replica+1]++;
			updateSum(proposed);*/
			//cout <<switchRatio()<<endl;
			magnet+=2*matrix[proposed].val;			
		}
	//}
	//double q=((double)rand()/(double)RAND_MAX);
	//cout<<q<< " "<<switchRatio()<<endl;
	/*if(q<switchRatio()){
		cout << q << " " << switchRatio() <<" "<<replica <<endl;
		replica=-replica;
		replicaCount[replica+1]++;
		
		return -1;
	}*/
	//else return 0;	
		loop--;	
	}
	return 0;
}

int main(){
	int Lat=32,iteration=100000,seed1=998,seed2=3394,step=1;
	//std::vector<int> magnetization;
	double p1=0.3,p2=0,p3=1;
	double thrh=0.5;
	/*cout << "Please input iteration,the value of lattice and J,H,T"<<endl;
	cin >>iteration>>L>>p1>>p2>>p3;
	cout <<L<<" "<<p1<<endl;*/
	IsingModel ismMH(Lat,p1,p2,p3);
	ismMH.initializeCongiguration(seed1,thrh);
	ismMH.initializeFlipRatio();
	//ism.initializeSum();
	ofstream file("ismMHsimulation-s.txt");
	srand(time(NULL));
	for(int i=0;i<iteration;i++){
		//ism.changeState(i);
		//magnetization.push_back(ism.magnet);
		//file << magnetization.back()<<endl;
		file << ismMH.changeState(step)<<" "<<ismMH.magnet<<endl;
		
	}
	//cout << magnetization[333]<< " "<<magnetization[iteration]<<endl;
	file.close();
	return  0;
}