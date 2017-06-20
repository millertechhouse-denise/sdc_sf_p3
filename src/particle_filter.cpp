/*
* particle_filter.cpp
*
*  Created on: Dec 12, 2016
*      Author: Tiffany Huang
*/

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <math.h> 
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>

#include "particle_filter.h"

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {

	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
	num_particles = 100;

	std::default_random_engine gen;
	std::normal_distribution<double> N_x(x, std[0]);
	std::normal_distribution<double> N_y(y, std[1]);
	std::normal_distribution<double> N_theta(theta, std[2]);

	for (int i = 0; i < num_particles; i++) {

		Particle p;
		p.id = i;
		p.x = N_x(gen);
		p.x = N_y(gen);
		p.x = N_theta(gen);
		p.weight = 1.0;

		weights.push_back(1.0);
		particles.push_back(p);
	}

	is_initialized = true;
}


void ParticleFilter::prediction(double delta_t, double std[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/


	std::default_random_engine gen;
	std::normal_distribution<double> N_x(0.0, std[0]);
	std::normal_distribution<double> N_y(0.0, std[1]);
	std::normal_distribution<double> N_theta(0.0, std[2]);

	for (int i = 0;  i < num_particles; i++) {

		double phi = particles[i].theta + yaw_rate * delta_t;

		if (fabs(yaw_rate) < 1e-4) {

			particles[i].x += velocity * delta_t  * cos(particles[i].theta)+ N_x(gen);
			particles[i].y += velocity * delta_t  * sin(particles[i].theta) + N_y(gen);
			particles[i].theta += N_theta(gen);
			//std::cout << "straight" << std::endl;

		} else {

      
			particles[i].x += velocity / yaw_rate * (sin(phi) - sin(particles[i].theta)) + N_x(gen);
			particles[i].y += velocity / yaw_rate * (cos(particles[i].theta) - cos(phi)) + N_y(gen);
			particles[i].theta = phi + N_theta(gen);
			
			//	std::cout << "turn" << std::endl;
		}
	}
}


void ParticleFilter::dataAssociation(vector<LandmarkObs> predicted, vector<LandmarkObs>& observations){
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.
	

	for (int i = 0; i < observations.size(); i++) 
	{

		double min = 1e9;

		for (int j = 0; j < predicted.size(); j++) 
		{

			double distance = dist(predicted[j].x, predicted[j].y, observations[i].x, observations[i].y);

			if (distance< min) 
			{
				observations[i].id = j;
				min = distance;
			}
		}
	}
}



void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], std::vector<LandmarkObs> observations, Map map_landmarks) {
	
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation 
	//   3.33
	//   http://planning.cs.uiuc.edu/node99.html


	for (int  i = 0; i < num_particles; i++) {

		vector<LandmarkObs> landmarks_observations;
		vector<LandmarkObs> map_observations;

		//transform to map coordinates
		for (int j = 0; j < observations.size(); j++)
		{
			
			LandmarkObs obs;
			obs.x = particles[i].x + observations[j].x * cos(particles[i].theta) - observations[j].y * sin(particles[i].theta);
			obs.y = particles[i].y + observations[j].y * cos(particles[i].theta) + observations[j].x * sin(particles[i].theta);
			obs.id = j;

			map_observations.push_back(obs);
		}

		//find landmarks
		for (int j = 0;  j < map_landmarks.landmark_list.size(); j++) 
		{
			
			double distance = dist(map_landmarks.landmark_list[j].x_f, map_landmarks.landmark_list[j].y_f, particles[i].x, particles[i].y);	

			if (distance < sensor_range) 
			{

				LandmarkObs landmark;
				landmark.x = map_landmarks.landmark_list[j].x_f;
				landmark.y = map_landmarks.landmark_list[j].y_f;
				landmark.id = map_landmarks.landmark_list[j].id_i;

				landmarks_observations.push_back(landmark);
			}
		}

		//associate landmarks with observations
		dataAssociation(landmarks_observations, map_observations);

		//update weight
		double weight = 1.0;

		for (int j = 0; j < map_observations.size(); j++){

			double dx = map_observations[j].x - landmarks_observations[map_observations[j].id].x;
			double dy = map_observations[j].y - landmarks_observations[map_observations[j].id].y;
			
			double w = exp(-(0.5 / (std_landmark[0] * std_landmark[0]) * dx * dx + 0.5 / (std_landmark[1] * std_landmark[1]) * dy * dy)) / sqrt( 2.0 * M_PI * std_landmark[0] * std_landmark[1]);
			weight *= w;
		}

		particles[i].weight = weight;
		weights[i] = weight;
	}
}

void ParticleFilter::resample(){
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

	vector<Particle> resample_particles;
	std::default_random_engine gen;
	std::discrete_distribution<int> index(weights.begin(), weights.end());

	for (int c = 0; c < num_particles; c++) {

		int i = index(gen);

		Particle p;
		
		p.id = 1;
		p.x = particles[i].x;
		p.y = particles[i].y;
		p.theta = particles[i].theta;
		p.weight = 1.0;

		resample_particles.push_back(p);
	}

	particles = resample_particles;
}

Particle ParticleFilter::SetAssociations(Particle particle, std::vector<int> associations, std::vector<double> sense_x, std::vector<double> sense_y)
{
	//particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
	// associations: The landmark id that goes along with each listed association
	// sense_x: the associations x mapping already converted to world coordinates
	// sense_y: the associations y mapping already converted to world coordinates

	//Clear the previous associations
	particle.associations.clear();
	particle.sense_x.clear();
	particle.sense_y.clear();

	particle.associations= associations;
	particle.sense_x = sense_x;
	particle.sense_y = sense_y;

	return particle;
}

string ParticleFilter::getAssociations(Particle best)
{
	vector<int> v = best.associations;
	stringstream ss;
	copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
	string s = ss.str();
	s = s.substr(0, s.length()-1);  // get rid of the trailing space
	return s;
}
string ParticleFilter::getSenseX(Particle best)
{
	vector<double> v = best.sense_x;
	stringstream ss;
	copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
	string s = ss.str();
	s = s.substr(0, s.length()-1);  // get rid of the trailing space
	return s;
}
string ParticleFilter::getSenseY(Particle best)
{
	vector<double> v = best.sense_y;
	stringstream ss;
	copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
	string s = ss.str();
	s = s.substr(0, s.length()-1);  // get rid of the trailing space
	return s;
}
