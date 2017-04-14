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
#include <map>
#include <cctype>

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

	num_particles = 1000;
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

	// find nearest predicted landmark thats not been assigned
	for (auto prediction : predicted) {

		// calculate distance for the observation to all predicted
		double dist_prev;
		LandmarkObs *closest_observation;
		bool first_observation = true;

		for (auto &observation : observations) {

			// if this observation has already been matched - continue to the next
			if (observation.id > 0)
				continue;

			double d = dist(observation.x, observation.y, prediction.x, prediction.y);
			if (first_observation) {

				// first time through so just say this is the closest
				dist_prev = d;
				observation.id = prediction.id;
				closest_observation = &observation;
				first_observation = false;
				continue;
			}

			// found a closer observation to the predicted landmark
			if (d < dist_prev) {
				// reassign this predicted landmark id to the current
				closest_observation->id = -1;
				observation.id = prediction.id;

				// save current as closest
				dist_prev = d;
				closest_observation = &observation;
			}
		}

	}  // end for range over observations

	// remove unmatched observations
	observations.erase(
		std::remove_if(observations.begin(), observations.end(),
						[](const LandmarkObs lo) {return lo.id < 1;}),
		observations.end());
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

	// for each particle
	for (auto &p : particles) {
		// convert observations to map space http://planning.cs.uiuc.edu/node99.html
		vector<LandmarkObs> observations_map(observations);
		for (auto &o : observations_map) {
			o.x = o.x * cos(p.theta) - o.y * sin(p.theta) + p.x;
			o.y = o.y * sin(p.theta) + o.y * cos(p.theta) - p.y;
		}

		// predict landmarks within sensor range of this particle
		vector<LandmarkObs> predicted;

		LandmarkObs lo;
		for (auto landmark : map_landmarks.landmark_list) {
			if (dist(p.x, p.y, landmark.x_f, landmark.y_f) < sensor_range) {
				lo = {landmark.id_i, landmark.x_f, landmark.y_f};
				predicted.push_back(lo);
			}
		}

		// find nearest neighbor - updates observations map with landmark id
		dataAssociation(predicted, observations_map);

		// predicted landmark lookup
		map<int, LandmarkObs> predicted_lookup;
		for (auto prediction: predicted) {
			if (prediction.id > 0) {
				predicted_lookup[prediction.id]=prediction;
			}
		}

		// initialise measurement covariance matrix
		// MatrixXd measurementCovar;
		// measurementCovar << std_landmark[0], 0,
        // 	                0, std_landmark[1];

		// calculate how likely a measurement should be
		double weight_product;
		bool first_measurement = true;
		for (auto measurement : observations_map) {
			// nearest neighbor observation to the predicted landmark
			LandmarkObs predicted_measurement = predicted_lookup[measurement.id]; // id points to landmark
			//double weight = multiVariateGaussianWeight(predicted_measurement, measurement, measurementCovar);
			double weight = bivariate_gaussian(predicted_measurement, measurement, std_landmark[0], std_landmark[1]);
			if (first_measurement) {
				weight_product = weight;
				first_measurement = false;
			} else {
				weight_product *= weight;
			}
		}
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
