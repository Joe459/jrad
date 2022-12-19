#include <iostream>
#include <eigen3/Eigen/Dense>
#include "mpi.h"
#include <stdlib.h>
#include <random>
#include <ctime>
#include "../include/input.h"

using namespace std;
using namespace Eigen;

typedef Vector<float, 4> four_vec;
typedef Vector<float, 3> three_vec;
const Matrix<float, 4, 4> I = Matrix<float, 4, 4>::Identity();
const Vector<float,4> metric_ = {1.0,-1.0,-1.0,-1.0};


vector<three_vec> fields_;

struct particle_state {
	
	three_vec position;
	three_vec momentum;
	
	three_vec E;
	three_vec B
	
	particle_state operator + (particle_state const &other_state) {
         particle_state res;
         res.position = position + other_state.position;
         res.momentum = momentum + other_state.momentum;
         return res;
    }
	

};
	
vector<vector<particle_state>> phase_space_;
vector<particle_state> temp_phase_space;

int rank_;
int size_;
int n_particles_;
int n_step = 0;
MPI_Datatype MPI_PARTICLE_STATES;
MPI_Status status;

void distribute_particles()
{
		
	default_random_engine generator(time(0) + rank_);
	vector<normal_distribution<double>> mom_dist = {normal_distribution<double>(0.0,1.0),normal_distribution<double>(0.0,1.0),normal_distribution<double>(0.0,1.0)};
	vector<normal_distribution<double>> pos_dist = {normal_distribution<double>(0.0,1.0),normal_distribution<double>(0.0,1.0),normal_distribution<double>(0.0,1.0)};
	
	n_particles_ = round(N_PARTICLES / size_);
	
	fields_.resize(n_particles_);
	phase_space_.resize(N_TIME_STATES);
	
	for (auto states : phase_space_){
		states.resize(n_particles_);
		for (particle_state state : states){
			state.momentum << mom_dist[0](generator), mom_dist[1](generator), mom_dist[2](generator);
			state.position << pos_dist[0](generator), pos_dist[1](generator), pos_dist[2](generator);
		}	
	}
	
	cout << sizeof(particle_state)*n_particles_ << endl;
	
	MPI_Type_contiguous(sizeof(particle_state)*n_particles_,MPI_BYTE,&MPI_PARTICLE_STATES);
	MPI_Type_commit(&MPI_PARTICLE_STATES);
}


void push_particles()
{
	n_step++;
	
	for (particle_state state : phase_space_[n_step % N_TIME_STATES])
	{
		const three_vec next_momentum = (I - 0.5 * dt * state.field_mat).inverse() * (I + 0.5 * dt * state.field_mat) * (metric_.cwiseProduct(state.momentum));
		state.position += 0.5 * (state.momentum + next_momentum) * dt;
		state.momentum = next_momentum;
	}
	
}

inline std::pair<three_vec,three_vec> lienard_wiechert_fields(float gamma, three_vec beta, three_vec accel,three_vec pos){
	
	const float R = pos.norm();
	const three_vec n = pos / R;
	
	const three_vec E = e * (n - beta) / (pow(gamma,2) * pow((1 - n.dot(beta)),3) * pow(R,2)) +
						e * n.cross((n - beta).cross(accel)) / (c * pow((1 - n.dot(beta)),3) * R);
						
	return std::make_pair(E, n.cross(E));
}


void update_fields()
{
	three_vec emitting_mom;
	three_vec emitting_accel;
	Vector<float, N_TIME_STATES> distances;


	// Iterate through neighboring ranks (including my own)
	for (int rank_idx = 0; rank_idx < size_; ++rank_idx){
	
		// "trade" particles with neighboring rank. The receiving particles' field matrices will be updated.
		MPI_Sendrecv(phase_space_[n_step % N_TIME_STATES].data(),1,MPI_PARTICLE_STATES, (rank_idx + rank_) % size_ ,MPI_ANY_TAG,temp_phase_space.data(),1,MPI_PARTICLE_STATES, (rank_idx + rank_ - 1  + size_ ) % size_, MPI_ANY_TAG, MPI_COMM_WORLD,&status);
		
		// Iterate through the particles acted UPON (particles that live in neighboring rank)
		for (particle_state &receiving_state : temp_phase_space){
		
			// Iterate through the particles EMITTING (particle_idx)
			for (int particle_idx = 0; particle_idx < n_particles_; ++particle_idx) {
							
				// Iterate through the history of the EMITTING particle
				for (int time_idx = 0; time_idx < N_TIME_STATES; ++time_idx) { 
					
					// The emitting state in between position & momentum update states, hence the averaging
					const three_vec emitting_pos = (phase_space_[(-time_idx + n_step + N_TIME_STATES) % N_TIME_STATES][particle_idx].position + phase_space_[(-time_idx - 1 + n_step + N_TIME_STATES) % N_TIME_STATES][particle_idx].position)/2;					
					distances(time_idx) = (emitting_pos - receiving_state.position.norm()  -  (time_idx * c * dt + c * dt / 2);
										
				}
				
				auto outer_ring = std::lower_bound(distances.begin(), distances.end(), 0);
				
				// The receiving particle is very close
				if (outer_ring == distances.begin()){
					const int idx = 0;
					const three_vec emitting_pos = (phase_space_[(-idx + n_step + N_TIME_STATES) % N_TIME_STATES][particle_idx].position + phase_space_[(-idx - 1 + n_step + N_TIME_STATES) % N_TIME_STATES][particle_idx].position)/2;	
					emitting_mom = (phase_space_[(-idx + n_step + N_TIME_STATES) % N_TIME_STATES][particle_idx].momentum + phase_space_[(-idx - 1 + n_step + N_TIME_STATES) % N_TIME_STATES][particle_idx].momentum)/2;				
					emitting_accel = (phase_space_[(-idx + n_step + N_TIME_STATES) % N_TIME_STATES][particle_idx].momentum - phase_space_[(-idx - 1 + n_step + N_TIME_STATES) % N_TIME_STATES][particle_idx].momentum/dt;
								
					std::pair<three_vec,three_vec> fields = lienard_wiechert_fields(
					  emitting_mom(0)/c, emitting_mom(seq(1,3)) / c, emitting_accel / c, receiving_state.position - emitting_pos);
					  
					receiving_state.update_field(fields.first,fields.second);
					  
				}
				
				// The receiving particle is very far
				else if (outer_ring == distances.end()){
					continue;
				}	
				
				// The receiving particle is just right
				else {
					
					// Compute the fields at the INNER RING
					int idx = (int)(outer_ring - distances.begin()) - 1;
						
					// the emitting state is in between neighboring particle-push states			
					const three_vec emitting_pos_inner = (phase_space_[(-idx + n_step + N_TIME_STATES) % N_TIME_STATES][particle_idx].position + phase_space_[(-idx - 1 + n_step + N_TIME_STATES) % N_TIME_STATES][particle_idx].position)/2;	
					emitting_mom = (phase_space_[(-idx + n_step + N_TIME_STATES) % N_TIME_STATES][particle_idx].momentum + phase_space_[(-idx - 1 + n_step + N_TIME_STATES) % N_TIME_STATES][particle_idx].momentum)/2;				
					emitting_accel = (phase_space_[(-idx + n_step + N_TIME_STATES) % N_TIME_STATES][particle_idx].momentum - phase_space_[(-idx - 1 + n_step + N_TIME_STATES) % N_TIME_STATES][particle_idx].momentum)/dt;
					
					const three_vec n = (receiving_state.position - emitting_pos_inner) / (receiving_state.position - emitting_pos_inner).norm();
					const float inner_arm_len = (idx * c * dt + c * dt / 2);
					
					std::pair<three_vec,three_vec> inner_fields = lienard_wiechert_fields(
					  sqrtf(emitting_mom.squaredNorm() + 1)/c, emitting_mom / c, emitting_accel / c, n * inner_arm_len);
					  
					// Compute the fields at the OUTER RING  
					idx++;
					const three_vec emitting_pos_outer = (phase_space_[(-idx + n_step + N_TIME_STATES) % N_TIME_STATES][particle_idx].position(seq(1,3)) + phase_space_[(-idx - 1 + n_step + N_TIME_STATES) % N_TIME_STATES][particle_idx].position(seq(1,3)))/2;	
					emitting_mom = (phase_space_[(-idx + n_step + N_TIME_STATES) % N_TIME_STATES][particle_idx].momentum + phase_space_[(-idx - 1 + n_step + N_TIME_STATES) % N_TIME_STATES][particle_idx].momentum)/2;				
					emitting_accel = (phase_space_[(-idx + n_step + N_TIME_STATES) % N_TIME_STATES][particle_idx].momentum(seq(1,3)) - phase_space_[(-idx - 1 + n_step + N_TIME_STATES) % N_TIME_STATES][particle_idx].momentum(seq(1,3)))/dt;
					
					const three_vec delta = emitting_pos_outer - emitting_pos_inner;	
					const float outer_arm_len = ( delta.dot(n) + sqrtf(powf(idx * c * dt + c * dt / 2,2) - delta.squaredNorm() + powf(delta.dot(n),2)) );
					 		 			  
					std::pair<three_vec,three_vec> outer_fields = lienard_wiechert_fields(
					  sqrtf(emitting_mom.squaredNorm() + 1), emitting_mom / c, emitting_accel / c, outer_arm_len * n - delta);
					  
					// linearly interpolate the fields and update the receiving state
					receiving_state.update_field( 
						(outer_fields.first - inner_fields.first) /  (outer_arm_len - inner_arm_len) * (receiving_state.position - emitting_pos_inner).norm() + inner_fields.first,
						(outer_fields.second - inner_fields.second) /  (outer_arm_len - inner_arm_len) * (receiving_state.position - emitting_pos_inner).norm() + inner_fields.second);
		  
				}	
			}	
		}
		
		// Send back particles with updated fields to neighboring rank
		MPI_Sendrecv(temp_phase_space.data(),1,MPI_PARTICLE_STATES, (rank_idx + rank_) % size_ ,MPI_ANY_TAG,phase_space_[n_step % N_TIME_STATES].data(),1,MPI_PARTICLE_STATES, (rank_idx + rank_ - 1  + size_ ) % size_, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
	
	}	
	
}

int main(int argc, char **argv)
{
	
	MPI_Init(&argc,&argv);

	MPI_Comm_rank(MPI_COMM_WORLD,&rank_);
    MPI_Comm_size(MPI_COMM_WORLD,&size_);
	
	while(true){
		push_particles();
		update_fields();
	}
}
