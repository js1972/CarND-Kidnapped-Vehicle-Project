/*
 * particle_filter.cpp
 *
 *  Created on: April 11, 2017
 *      Author: Jason Scott
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <time.h>

#include "particle_filter.h"

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).

	default_random_engine gen;
	normal_distribution<double> x_dist(x, std[0]);
	normal_distribution<double> y_dist(x, std[1]);
	normal_distribution<double> theta_dist(x, std[2]);

	num_particles = 100;
	weights.resize(num_particles);

	for (int i=0; i<num_particles; i++) {
		Particle p = {
			0,					//id
			x_dist(gen),		//x
			y_dist(gen), 		//y
			theta_dist(gen), 	//theta
			1 					//weight
		};

		particles.push_back(p);
	}

	is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/

	default_random_engine gen;
	normal_distribution<double> x_noise_dist(0, std_pos[0]);
	normal_distribution<double> y_noise_dist(0, std_pos[1]);
	normal_distribution<double> theta_noise_dist(0, std_pos[2]);

	// for (int i=0; i<num_particles; i++) {
	// 	Particle p = particles[i];

	// 	// predict next position of particle with bicycle motion equations:
	// 	// x_f = x_0 + v / yaw_rate * [sin(θ_0 + yaw_rate*dt) - sin(θ_0)]
	// 	// y_f = y_0 + v / yaw_rate * [cos(θ_0) - cos(θ_0 + yaw_rate*dt)]
	// 	// θ_f = θ_0 + yaw_rate * dt

	// 	p.x += velocity / yaw_rate * (sin(p.theta + yaw_rate*delta_t) - sin(p.theta)) + x_noise_dist(gen);
	// 	p.y += velocity / yaw_rate * (cos(p.theta) - cos(p.theta + yaw_rate*delta_t)) + y_noise_dist(gen);
	// 	p.theta += yaw_rate*delta_t + theta_noise_dist(gen);

	// 	particles[i] = p;
	// }

	// *** STILL NEED TO HANDLE WHEN YAW_RATE IS ZERO HERE ***

	for (auto &p: particles) {
		p.x += velocity / yaw_rate * (sin(p.theta + yaw_rate*delta_t) - sin(p.theta)) + x_noise_dist(gen);
		p.y += velocity / yaw_rate * (cos(p.theta) - cos(p.theta + yaw_rate*delta_t)) + y_noise_dist(gen);
		p.theta += yaw_rate*delta_t + theta_noise_dist(gen);
	}

	debug_iteration += 1;
	cout << "DEBUG ITER: " << debug_iteration << endl;
	if (debug_iteration == 100) {
		cout << "Output particles at step 100 for debug" << endl;
		debug_output_particles(particles);
	}
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.

	// Applying a rotation followed by a translation to change coordinate systems from the map that the
	// observations are given in to the particles
	//
	// Because our maps y-axis points down we have switch the signs below:
	// x: -xcosθ + ysinθ + x_t
	// y: -xsinθ - ycosθ + y_t
	// See http://planning.cs.uiuc.edu/node99.html for reference.

}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		std::vector<LandmarkObs> observations, Map map_landmarks) {
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation 
	//   3.33. Note that you'll need to switch the minus sign in that equation to a plus to account 
	//   for the fact that the map's y-axis actually points downwards.)
	//   http://planning.cs.uiuc.edu/node99.html

	// For each particle go through all the map landmarks and perform a rotational and
	// translational conversion to the particles viewpoint (coords).

	for (int i=0; i<num_particles; i++) {
		Particle &p = particles[i];

		// transform map_landmarks to particle space - place particle at the origin
		// pointing toward the positive y-axis.
		vector<LandmarkObs> landmarks_in_particle_space;
		for (int j=0; j<map_landmarks.landmark_list.size(); j++) {
			Map::single_landmark_s map_single_landmark = map_landmarks.landmark_list[j];

			LandmarkObs particle_space_landmark;
			particle_space_landmark.id = map_single_landmark.id_i;

			// re-point by 90 degrees not sure if this is necessary - need to plot this stuff???
			double new_theta = p.theta - M_PI / 2;

			float x_distance = map_single_landmark.x_f - p.x;
			float y_distance = map_single_landmark.y_f - p.y;

			particle_space_landmark.x = -(x_distance) * sin(new_theta) + (y_distance) * cos(new_theta);
			particle_space_landmark.y = -(x_distance) * cos(new_theta) - (y_distance) * sin(new_theta);

			landmarks_in_particle_space.push_back(particle_space_landmark);
		}

		// Sort the particle-space landmarks - this will place the landmarks closest
		// to the origin.
		// The observations are already sorted. Should we sort them anyway???
		// Allow the LandmarkObs to be sorted by the euclidean distance from the origin.

		// Perform the data association step
		vector<LandmarkObs> associated_landmarks;
		sort(landmarks_in_particle_space.begin(), landmarks_in_particle_space.end(), compare_distance);
		for (auto &lm: landmarks_in_particle_space) {
			double range = dist(lm.x, lm.y, 0, 0);
			if (range <= sensor_range) {
				associated_landmarks.push_back(lm);
			}
		}

		LandmarkObs associated_landmark;
		if (associated_landmarks.size() > 0) {
			sort(associated_landmarks.begin(), associated_landmarks.end(), compare_distance);
			associated_landmark = associated_landmarks[0];
		}

		//dataAssociation(landmarks_in_particle_space, observations);

		// Update the particles weights using a bivariate  gaussian distribution

		double posterior_observation = 1.0;

		for (int j=0; j<observations.size(); j++) {
			LandmarkObs observation = observations[j];

			// bi-variate Gaussian distribution with correlation (ρ) as zero
			double sigma_x = std_landmark[0];
			double sigma_y = std_landmark[1];
			double x_0 = observation.x;
			double y_0 = observation.y;
			double x = associated_landmark.x;
			double y = associated_landmark.y;

			double x_part = (x - x_0)*(x - x_0) / (sigma_x*sigma_x);
			double y_part = (y - y_0)*(y - y_0) / (sigma_y*sigma_y);
			double first_term = 1 / (2*M_PI*sigma_x*sigma_y*sqrt(1));

			posterior_observation *= first_term * exp(-1/2 * (x_part + y_part));
		}

		p.weight = posterior_observation;
		weights[i] = posterior_observation;
	}
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

	// C++ implementation of Sebastions python wheel-resampling logic - Particle Filters Lesson 20.

	// max element returns an iterator - deref with *
	// -> http://stackoverflow.com/questions/10158756/using-stdmax-element-on-a-vectordouble
	double max_weight = *max_element(begin(weights), end(weights));
	//srand(time(NULL)); //seed
	//int index = rand() % num_particles; // not a true uniform dist?!?

	// get random index
	default_random_engine gen;
	uniform_int_distribution<int> uni(0, num_particles);
	int index = uni(gen);

	uniform_real_distribution<double> uniform_dist(0, 2 * max_weight);

	double beta = 0.0;
	vector<Particle> resampled_particles;

	for (int i=0; i<num_particles; i++) {
		beta += uniform_dist(gen);
		while (beta > weights[index]) {
			beta -= weights[index];
			index = (index + 1) % num_particles;
		}

		resampled_particles.push_back(particles[index]);
	}

	particles = resampled_particles;
}

void ParticleFilter::write(std::string filename) {
	// You don't need to modify this file.
	std::ofstream dataFile;
	dataFile.open(filename, std::ios::app);
	for (int i = 0; i < num_particles; ++i) {
		dataFile << particles[i].x << " " << particles[i].y << " " << particles[i].theta << "\n";
	}
	dataFile.close();
}
