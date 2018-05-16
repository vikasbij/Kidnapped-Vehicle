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
	// Number of particles 
	num_particles = 20;

	default_random_engine gen;
	// Create normal distributions for x 
	normal_distribution<double> N_x(x, std[0]);
	// Create normal distributions for y 
	normal_distribution<double> N_y(y, std[1]);
	// Create normal distributions for theta
	normal_distribution<double> N_theta(theta, std[2]);

	for(int i=0; i<num_particles; i++ )
	{
		Particle p;
		p.id = i;
		p.x = N_x(gen);
		p.y = N_y(gen);
		p.theta = N_theta(gen);
		p.weight = 1;

		particles.push_back(p);
		weights.push_back(1);
	}
	is_initialized = true;

}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/
	default_random_engine gen;
	//for i in random
	double x_new = std_pos[0];
	double y_new = std_pos[1];
	double theta_new = std_pos[2];
	fill(weights.begin(),weights.end(),1);
	for(int i=0; i<num_particles; i++)
	{
		double theta = particles[i].theta;
		if ( fabs(yaw_rate) < 0.001 ) 
		{ 
			particles[i].x += velocity * delta_t * cos( theta );
        	particles[i].y += velocity * delta_t * sin( theta );
        
        } 
		else 
		{
        	particles[i].x += velocity / yaw_rate * ( sin( theta + yaw_rate * delta_t ) - sin( theta ) );
        	particles[i].y += velocity / yaw_rate * ( cos( theta ) - cos( theta + yaw_rate * delta_t ) );
        	particles[i].theta += yaw_rate * delta_t;
        }
		normal_distribution<double> N_x(particles[i].x, x_new); 
		normal_distribution<double> N_y(particles[i].y, y_new);
		normal_distribution<double> N_theta(particles[i].theta, theta_new);

		// Adding noise.
        particles[i].x = N_x(gen);
        particles[i].y = N_y(gen);
        particles[i].theta = N_theta(gen);
        particles[i].weight = 1.0;
	}

}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.

}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		const std::vector<LandmarkObs> &observations, const Map &map_landmarks) {
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
	//vector<double> sense_x;
	//vector<double> sense_y;
	double normalizer = 2.0 * M_PI * std_landmark[0] * std_landmark[1];
	int index = 1;
	//vector<LandmarkObs> transformed_obs;
	//LandmarkObs obs;
	for (int i=0; i<num_particles; i++)
	{
		double part_x = particles[i].x;
		double part_y = particles[i].y;
		double part_theta = particles[i].theta;
		
		double total_prob = 1.0;
		//performing space transformation from vehicle to map
	for (int j=0; j < observations.size();j++)
	{	
		double transf_x = part_x + (observations[j].x * cos(part_theta)- observations[j].y * sin(part_theta));
		double transf_y = part_y + (observations[j].x * sin(part_theta)+ observations[j].y * cos(part_theta));	
		double min_dist = 100;
		
	for (int k =0; k<map_landmarks.landmark_list.size(); k++)
	{
		double distance = dist(transf_x,transf_y, map_landmarks.landmark_list[k].x_f, map_landmarks.landmark_list[k].y_f);
	    if (distance<min_dist)
		{
			min_dist = distance;
			index = k;
		}
	}
	
			double dx = transf_x - map_landmarks.landmark_list[index].x_f;
			double dy = transf_y - map_landmarks.landmark_list[index].y_f;
			double power = (dx * dx) / (2 * std_landmark[0] * std_landmark[0]) + (dy * dy) / (2 * std_landmark[1] * std_landmark[1]);
			total_prob *= exp(-power) / normalizer;
			}
			particles[i].weight = total_prob;
			weights[i] = total_prob;
			}
		}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

	vector<Particle> new_particles;

    random_device rd;
     mt19937 gen(rd());
     discrete_distribution<> distribution(weights.begin(), weights.end());

    for(int i = 0; i < num_particles; i++)
	{
        Particle p = particles[distribution(gen)];
        new_particles.push_back(p);
    }
    particles = new_particles;

}

Particle ParticleFilter::SetAssociations(Particle& particle, const std::vector<int>& associations, 
                                     const std::vector<double>& sense_x, const std::vector<double>& sense_y)
{
    //particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
    // associations: The landmark id that goes along with each listed association
    // sense_x: the associations x mapping already converted to world coordinates
    // sense_y: the associations y mapping already converted to world coordinates

    particle.associations= associations;
    particle.sense_x = sense_x;
    particle.sense_y = sense_y;
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
