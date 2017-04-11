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

	for (auto & p: particles) {
		p.x += velocity / yaw_rate * (sin(p.theta + yaw_rate*delta_t) - sin(p.theta)) + x_noise_dist(gen);
		p.y += velocity / yaw_rate * (cos(p.theta) - cos(p.theta + yaw_rate*delta_t)) + y_noise_dist(gen);
		p.theta += yaw_rate*delta_t + theta_noise_dist(gen);
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

	// For each particle go through all the map landmarks

	for (int i=0; i<num_particles; i++) {
		Particle &p = particles[i];

		// transform map_landmarks to p's coordinate system
		vector<LandmarkObs> transformed_landmarks;

		for (int j=0; j<map_landmarks.landmark_list.size(); j++) {
			Map::single_landmark_s landmark = map_landmarks.landmark_list[j];

			LandmarkObs transformed_landmark;
			transformed_landmark.id = landmark.id_i;

			double cos_theta = cos(p.theta - M_PI / 2);
			double sin_theta = sin(p.theta - M_PI / 2);

			float x_distance = landmark.x_f - p.x;
			float y_distance = landmark.y_f - p.y;

			transformed_landmark.x = -(x_distance) * sin_theta + (y_distance) * cos_theta;
			transformed_landmark.y = -(x_distance) * cos_theta - (y_distance) * sin_theta;

			transformed_landmarks.push_back(transformed_landmark);
		}
	}
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

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
