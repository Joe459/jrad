#include <iostream>
#include <eigen3/Eigen/Dense>
#include "mpi.h"
#include <stdlib.h>
#include <random>
#include <ctime>
#include "../include/input.h"
#include <iostream>
#include <fstream>
#include <string>
#include <eigen3/Eigen/StdVector>

using namespace std;
using namespace Eigen;

typedef Vector<float, 4> four_vec;
typedef Vector<float, 3> three_vec;
typedef Matrix<float,4,4> field_mat;
const Matrix<float, 4, 4> I = Matrix<float, 4, 4>::Identity();
const Vector<float,4> metric_ = {1.0,-1.0,-1.0,-1.0};

#define N_TIME_STATES 200
#define N_PARTICLES 2
#define c 299.97f //micrometers / picosecond
#define dt 0.01f //picoseconds
#define e 4.8032e-10f * 1e6f / 1e12f //esu converted to micrometers and picoseconds
#define m 9.1095e-28f //grams

//external electric field should be converted to these modified gaussian units

struct particle_state {
	four_vec x;
	four_vec U;
};

struct fields {
	three_vec E;
	three_vec B;
};

std::ofstream profile_fs;
 
float dtau;
 
vector<vector<particle_state>> phase_space_;
vector<particle_state> receiving_states;
vector<fields> acting_fields;
vector<fields> my_fields;

int rank_;
int size_;
int n_particles_;
int step = 0;
MPI_Datatype MPI_PARTICLE_STATES;
MPI_Datatype MPI_FIELDS;
MPI_Datatype MPI_THREE_VEC;
MPI_Status status;
MPI_Op FIELD_ADD_OP;


void update_field(field_mat &F, three_vec E, three_vec B){
	
	F(0,seq(1,3)) -= E;
	F(seq(1,3),0) += E;
	
	F(1,2) -= B(2);
	F(2,1) += B(2);
	F(3,1) -= B(1);
	F(1,3) += B(1);
	F(2,3) -= B(0);
	F(3,2) += B(0);
	
}
void set_field(field_mat &F,three_vec E, three_vec B){
	F <<  0 ,  -E(0), -E(1),-E(2),
	    E(0),    0  , -B(2), B(1),
	    E(1),   B(2),   0  ,-B(0),
	    E(2),  -B(1),  B(0),  0  ;
}

field_mat F(three_vec &E, three_vec &B){
	field_mat F;
	F <<  0 ,  -E(0), -E(1),-E(2),
	    E(0),    0  , -B(2), B(1),
	    E(1),   B(2),   0  ,-B(0),
	    E(2),  -B(1),  B(0),  0  ;
   return F;
}

field_mat F(fields &f){
	field_mat F;
	F <<  0 ,  -f.E(0),-f.E(1),-f.E(2),
	    f.E(0),    0  ,-f.B(2), f.B(1),
	    f.E(1), f.B(2),   0   ,-f.B(0),
	    f.E(2),-f.B(1), f.B(0),   0   ;

   return F;
}

void fields_add(fields *in, fields *inout, int *len, MPI_Datatype *dptr ){
	for (int i = 0; i < *len; i++){ 	
		inout[i].E += in[i].E;
		inout[i].B += in[i].B;
	}
} 

void print_states(int time_idx){
	
	if (profile_fs)
		for (int j = 0; j < n_particles_; j++){
			profile_fs << phase_space_[time_idx % N_TIME_STATES][j].x.transpose() << '\t';
			profile_fs << phase_space_[time_idx % N_TIME_STATES][j].U.transpose() << endl;
		}
}

void distribute_particles(){
		
	default_random_engine generator(time(0) + rank_);
	vector<normal_distribution<double>> mom_dist = {normal_distribution<double>(0.0,1e-2),normal_distribution<double>(0.0,1e-2),normal_distribution<double>(0.0,1e-2)};
	vector<normal_distribution<double>> pos_dist = {normal_distribution<double>(0.0,30.0),normal_distribution<double>(0.0,30.0),normal_distribution<double>(0.0,100.0)};
	
	n_particles_ = round(N_PARTICLES / size_);
	
	phase_space_.resize(N_TIME_STATES);

	for (int time_idx = 0; time_idx < N_TIME_STATES; time_idx++){
		
		if (time_idx == 0){	
			phase_space_[time_idx].resize(n_particles_);
			for (particle_state &state : phase_space_[0]){
				state.U(seq(1,3)) << mom_dist[0](generator), mom_dist[1](generator), mom_dist[2](generator);
				state.U(0) = sqrtf(state.U(seq(1,3)).squaredNorm() + 1);
				state.x(seq(1,3)) << pos_dist[0](generator), pos_dist[1](generator), pos_dist[2](generator);
				state.x(0) = 0;
				}		
		}
		else {	
			phase_space_[time_idx] = phase_space_[0];
			for (particle_state &state : phase_space_[time_idx]){
				state.x(0) = -(N_TIME_STATES-time_idx) * dt;
				state.x(seq(1,3)) += state.x(0) * state.U(seq(1,3)) * c / state.U(0);
			}	
			print_states(time_idx);
		}
	}

	receiving_states.resize(n_particles_);
	acting_fields.resize(n_particles_);
	my_fields.resize(n_particles_);
}


inline std::pair<three_vec,three_vec> lienard_wiechert_fields(float gamma, three_vec beta, three_vec accel,three_vec pos){
	
	const float R = pos.norm();
	const three_vec n = pos / R;
	
	const three_vec E = e * (n - beta) / (pow(gamma,2) * pow((1.f - n.dot(beta)),3) * pow(R,2)) +
						e * n.cross((n - beta).cross(accel)) / (c * pow((1 - n.dot(beta)),3) * R);
						
	return std::make_pair(E, n.cross(E));
}
inline void lienard_wiechert_fields(fields &f, float gamma, three_vec beta, three_vec accel,three_vec pos){
	
	const float R = pos.norm();
	const three_vec n = pos / R;
	
	f.E += e * (n - beta) / (pow(gamma,2) * pow((1.f - n.dot(beta)),3) * pow(R,2)) +
						e * n.cross((n - beta).cross(accel)) / (c * pow((1 - n.dot(beta)),3) * R);
	
	f.B += n.cross(f.E);

}

void update_fields()
{
	four_vec emitting_mom;
	three_vec emitting_accel;
	Vector<float, N_TIME_STATES> distances;

	// Iterate through ranks (including my own)
	for (int rank_idx = 0; rank_idx < size_; ++rank_idx){
	
		MPI_Scatter(phase_space_[step % N_TIME_STATES].data(),n_particles_,MPI_PARTICLE_STATES, receiving_states.data(),n_particles_,MPI_PARTICLE_STATES,rank_idx,MPI_COMM_WORLD);
		//MPI_Sendrecv(phase_space_[step % N_TIME_STATES].data(),1,MPI_PARTICLE_STATES, (rank_idx + rank_) % size_ , 0 ,receiving_states.data(),1,MPI_PARTICLE_STATES, (rank_idx + rank_ - 1  + size_ ) % size_, 0 , MPI_COMM_WORLD,&status);
		
		// Iterate through the particles acted UPON (Receiving particles)
		for (int j = 0; j < n_particles_; j++ ){
		
			// Iterate through the particles EMITTING (i)
			for (int i = 0; i < n_particles_; ++i) {
						
				// no self-fields allowed
				if (rank_ == rank_idx && i == j)
					continue;
							
				// Iterate through the history of the EMITTING particle
				for (int time_idx = 0; time_idx < N_TIME_STATES; ++time_idx) { 
					
					// The emitting state in between position & momentum update states, hence the averaging
					const four_vec emitting_x = (phase_space_[time_idx % N_TIME_STATES][i].x + phase_space_[(time_idx - 1 + N_TIME_STATES) % N_TIME_STATES][i].x)/2;					
					distances(time_idx) = (receiving_states[j].x(seq(1,3)) - (emitting_x(seq(1,3)) + c * (receiving_states[j].x(0) -  emitting_x(0)) * (receiving_states[j].x(seq(1,3)) - emitting_x(seq(1,3))).normalized()) ).norm();
										
				}
				
				auto outer_ring = std::lower_bound(distances.begin(), distances.end(), 0);
				
				// The receiving particle is very close
				if (outer_ring == distances.begin()){
						
					const int idx = 0;
					const four_vec emitting_x = (phase_space_[0][i].x + phase_space_[N_TIME_STATES - 1][i].x)/2;	
					emitting_mom = (phase_space_[0][i].U + phase_space_[N_TIME_STATES - 1][i].U)/2;				
					emitting_accel = (phase_space_[0][i].U(seq(1,3)) - phase_space_[N_TIME_STATES - 1][i].U(seq(1,3)))
					/(phase_space_[0][i].x(0) - phase_space_[N_TIME_STATES - 1][i].x(0));
								
					lienard_wiechert_fields(acting_fields[j],
					  emitting_mom(0), emitting_mom(seq(1,3)), emitting_accel, receiving_states[j].x(seq(1,3)) - emitting_x(seq(1,3)));

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
					const four_vec emitting_x_inner = (phase_space_[idx][i].x + phase_space_[(idx - 1 + N_TIME_STATES) % N_TIME_STATES][i].x)/2;	
					emitting_mom = (phase_space_[idx][i].U + phase_space_[(idx - 1 + N_TIME_STATES) % N_TIME_STATES][i].U)/2;				
					emitting_accel = (phase_space_[idx][i].U(seq(1,3)) - phase_space_[(idx - 1 + N_TIME_STATES) % N_TIME_STATES][i].U(seq(1,3)))
					 / (phase_space_[idx][i].x(0) - phase_space_[(idx - 1 + N_TIME_STATES) % N_TIME_STATES][i].x(0));
					
					const three_vec n = (receiving_states[j].x(seq(1,3)) - emitting_x_inner(seq(1,3))).normalized();
					const float inner_arm_len = c * (emitting_x_inner(0) - receiving_states[j].x(0));
					
					std::pair<three_vec,three_vec> inner_fields = lienard_wiechert_fields(
					  emitting_mom(0), emitting_mom(seq(1,3)), emitting_accel, n * inner_arm_len);
					
					// Compute the fields at the OUTER RING  
					idx++;
					const four_vec emitting_x_outer = (phase_space_[idx][i].x + phase_space_[(idx - 1 + N_TIME_STATES) % N_TIME_STATES][i].x)/2;	
					emitting_mom = (phase_space_[idx][i].U + phase_space_[(idx - 1 + N_TIME_STATES) % N_TIME_STATES][i].U)/2;				
					emitting_accel = (phase_space_[idx][i].U(seq(1,3)) - phase_space_[(idx - 1 + N_TIME_STATES) % N_TIME_STATES][i].U(seq(1,3)))
					/(phase_space_[idx][i].x(0) - phase_space_[(idx - 1 + N_TIME_STATES) % N_TIME_STATES][i].x(0));
					
					const three_vec delta = emitting_x_outer(seq(1,3)) - emitting_x_inner(seq(1,3));	
					const float outer_arm_len = ( delta.dot(n) + sqrtf(powf(c * (emitting_x_outer(0) - receiving_states[j].x(0)),2) - delta.squaredNorm() + powf(delta.dot(n),2)) );
					
					std::pair<three_vec,three_vec> outer_fields = lienard_wiechert_fields(
					  emitting_mom(0), emitting_mom(seq(1,3)) , emitting_accel, outer_arm_len * n - delta);	
					
					acting_fields[j].E += (outer_fields.first - inner_fields.first) /  (outer_arm_len - inner_arm_len) * (receiving_states[j].x(seq(1,3)) - emitting_x_inner(seq(1,3))).norm() + inner_fields.first;
					acting_fields[j].B += (outer_fields.second - inner_fields.second) /  (outer_arm_len - inner_arm_len) * (receiving_states[j].x(seq(1,3)) - emitting_x_inner(seq(1,3))).norm() + inner_fields.second;
				
				}	
			}	
		}
		
		MPI_Reduce(acting_fields.data(), my_fields.data(),  n_particles_, MPI_FIELDS, FIELD_ADD_OP, rank_idx, MPI_COMM_WORLD);
		
	}	
}

void push_particles()
{		
	particle_state prev_state;
	
	for (int j = 0; j < n_particles_; ++j){
		
		prev_state = phase_space_[(step - 1 + N_TIME_STATES) % N_TIME_STATES][j];
		
		dtau = dt / prev_state.U(0);
				
		phase_space_[step % N_TIME_STATES][j].U = (I - 0.5f * dtau * e / m / c * F(my_fields[j])).inverse() * (I + 0.5f * e / m / c * dtau * F(my_fields[j])) * (metric_.cwiseProduct(prev_state.U));
		phase_space_[step % N_TIME_STATES][j].x = prev_state.x + 0.5f * (prev_state.U + phase_space_[step % N_TIME_STATES][j].U) * dtau * c;
		phase_space_[step % N_TIME_STATES][j].x(0) /= c;
	
	}
}

int main(int argc, char **argv)
{
	MPI_Init(&argc,&argv);

	MPI_Comm_rank(MPI_COMM_WORLD,&rank_);
    MPI_Comm_size(MPI_COMM_WORLD,&size_);
	
	profile_fs.open("./bunch-prof/" + std::to_string(rank_) + ".txt",
	std::ios::trunc);  

	distribute_particles();
	
	MPI_Op_create((MPI_User_function*)fields_add, true, &FIELD_ADD_OP);
	
	MPI_Type_contiguous(sizeof(particle_state),MPI_BYTE, &MPI_PARTICLE_STATES);
	MPI_Type_commit(&MPI_PARTICLE_STATES);
	MPI_Type_contiguous(sizeof(fields),MPI_BYTE, &MPI_FIELDS);
	MPI_Type_commit(&MPI_FIELDS);
	MPI_Type_contiguous(sizeof(three_vec),MPI_BYTE, &MPI_THREE_VEC);
	MPI_Type_commit(&MPI_THREE_VEC);
	

	while(step < 1){
		update_fields();
		push_particles();
		print_states(step);
		step++;
	}
	
	profile_fs.close();

	MPI_Op_free(&FIELD_ADD_OP);
	MPI_Finalize();
	

}
