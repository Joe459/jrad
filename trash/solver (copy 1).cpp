#include <iostream>
#include <boost/numeric/odeint.hpp>
#include <eigen3/Eigen/Dense>
#include <boost/math/distributions/normal.hpp>
#include <boost/numeric/odeint/external/eigen/eigen_algebra.hpp>
#include "mpi.h"
#include <random>
#include <ctime>

using namespace std;
using namespace boost::numeric::odeint;
using namespace Eigen;

typedef Vector<float, 4> four_vec;
typedef Vector<float, 3> three_vec;

#define N_TIME_STATES 200
#define N_PARTICLES 1e6
#define c 299.97

vector<three_vec> fields_;

class particle_states
{
public:

	std::vector<Vector<float,8>> data;

	particle_states(std::vector<Vector<float,8>>) : data(in) {}
	
	particle_states &operator+=(const particle_states &p)
	{
		data += p.data;
		return *this;
	}
	particle_states &operator*=(const particle_states &p)
	{
		data *= p.data;
		return *this;
	}
	
	
	
	
	four_vec position;
	four_vec momentum;
}
	
vector<vector<particle_state>> phase_space_;
vector<particle_state> temp_phase_space;

int rank_;
int size_;
int n_particles_;
int n_step;
MPI_Datatype MPI_PARTICLE_STATES;

void distribute_particles()
{
	default_random_engine generator(time(0) + rank_);
	vector<normal_distribution<double>> mom_dist = {normal_distribution<double>(0.0,1.0),normal_distribution<double>(0.0,1.0),normal_distribution<double>(0.0,1.0)};
	vector<normal_distribution<double>> pos_dist = {normal_distribution<double>(0.0,1.0),normal_distribution<double>(0.0,1.0),normal_distribution<double>(0.0,1.0)};
	
	n_particles_ = round(N_PARTICLES / size_)
	
	fields_.resize(n_particles_);
	phase_space_.resize(N_TIME_STATES);
	
	for (auto states : phase_space_){
		states.resize(n_particles_);
		for (particle_state state : states){
			state.momentum(seq(1,3)) << mom_dist[0](generator), mom_dist[1](generator), mom_dist[2](generator);
			state.momentum(0) = sqrt(state.momentum(seq(1,3)).squaredNorm() / pow(c,2) + 1);
			state.position(seq(1,3)) << pos_dist[0](generator), pos_dist[1](generator), pos_dist[2](generator);
		}	
	}
	
	cout << sizeof(particle_state)*n_particles_ << endl;
	
	MPI_Type_contiguous(sizeof(particle_state)*n_particles_,MPI_BYTE,&MPI_PARTICLE_STATES);
	MPI_Type_commit(&MPI_PARTICLE_STATES);
}


void push_particles( const four_vec &x , four_vec &dxdt , double t )
{
	n_step++;
	
	for (particle_state state : phase_space_[n_step % N_TIME_STATES])
		
	
}

void update_fields( const four_vec &x , const double t )
{
	for (int j = 0; j < size_; ++j){
	
		MPI_Sendrecv(phase_space_[n_step % N_TIME_STATES].data(),1,MPI_PARTICLE_STATES, (j + rank_) % size_ ,MPI_ANY_TAG,temp_phase_space.data(),1,MPI_PARTICLE_STATES, (j + rank_ - 1  + size_ ) % size_, MPI_ANY_TAG);
		
		for (particle_state state_1 : temp_phase_space){
			
			for (int time_idx = n_step; time_idx < n_step + N_TIME_STATES; time_idx++){
				
				for (particle_state state_2 : phase_space_[time_idx % N_TIME_STATES]){
					
					(state_1.position(seq(1,3)) - state_2.position(seq(1,3))).norm()
					
				}
				
			}
			
			
		}
		
		
	}
	
}

int main(int argc, char **argv)
{
	
	MPI_Init(&argc,&argv);

	MPI_Comm_rank(MPI_COMM_WORLD,&rank_);
    MPI_Comm_size(MPI_COMM_WORLD,&size_);
	
    four_vec x = { 10.0 , 1.0 , 1.0, 1.0}; // initial conditions
    integrate( push_particles , x , 0.0 , 25.0 , 0.1 , update_fields );
}
