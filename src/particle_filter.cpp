/**
 * particle_filter.cpp
 *
 * Created on: Dec 12, 2016
 * Author: Tiffany Huang
 */

#include "particle_filter.h"

#include <math.h>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <numeric>
#include <random>
#include <string>
#include <vector>

#include "helper_functions.h"

using std::string;
using std::vector;
using std::normal_distribution;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
  /**
   * TODO: Set the number of particles. Initialize all particles to 
   *   first position (based on estimates of x, y, theta and their uncertainties
   *   from GPS) and all weights to 1. 
   * TODO: Add random Gaussian noise to each particle.
   * NOTE: Consult particle_filter.h for more information about this method 
   *   (and others in this file).
   */
  num_particles = 1000;  // TODO: Set the number of particles
  
  weights.resize(num_particles);
  particles.resize(num_particles);
  
  std::default_random_engine gen;
  normal_distribution<double> dist_x(x, std[0]);
  normal_distribution<double> dist_y(y, std[1]);
  normal_distribution<double> dist_theta(theta, std[2]);
  
  for(int i = 0; i < num_particles; i++) {
    Particle particle;
    particle.id = i;
    particle.x = dist_x(gen);
    particle.y = dist_y(gen);
    particle.theta = dist_theta(gen);
    particle.weight = 1;
    
    particles[i] = particle;
    weights[i] = particle.weight;
    
  }
  
  is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], 
                                double velocity, double yaw_rate) {
  /**
   * TODO: Add measurements to each particle and add random Gaussian noise.
   * NOTE: When adding noise you may find std::normal_distribution 
   *   and std::default_random_engine useful.
   *  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
   *  http://www.cplusplus.com/reference/random/default_random_engine/
   */
  
//   std::cout << "PARAM: t "<< delta_t << " v "<< velocity << " yr "<< yaw_rate << std::endl;
  
  double epsilon = 0.0000000001;
  
  std::default_random_engine gen;
  for(int i = 0; i < particles.size(); i++) {
    Particle *p = &particles[i];
    

  	double new_x = particles[i].x + (velocity / (yaw_rate + epsilon)) * (sin(particles[i].theta + yaw_rate * delta_t) - sin(particles[i].theta));
    double new_y = particles[i].y + (velocity / (yaw_rate + epsilon)) * (cos(particles[i].theta) - cos(particles[i].theta + yaw_rate * delta_t));
    double new_theta = particles[i].theta + yaw_rate * delta_t;
    
    normal_distribution<double> dist_x(0, std_pos[0]);
  	normal_distribution<double> dist_y(0, std_pos[1]);
  	normal_distribution<double> dist_theta(0, std_pos[2]);
    p->x = new_x + dist_x(gen);
    p->y = new_y + dist_y(gen);
    p->theta = new_theta + dist_theta(gen);
    
//     std::cout << "PARAM: x "<< p->x << " y "<< p->y << " t "<< p->theta << std::endl;
  }

}

void ParticleFilter::dataAssociation(vector<LandmarkObs> predicted, 
                                     vector<LandmarkObs>& observations) {
  /**
   * TODO: Find the predicted measurement that is closest to each 
   *   observed measurement and assign the observed measurement to this 
   *   particular landmark.
   * NOTE: this method will NOT be called by the grading code. But you will 
   *   probably find it useful to implement this method and use it as a helper 
   *   during the updateWeights phase.
   */
	for(auto p : predicted) {
      	double dist_min = std::numeric_limits<double>::max();
    	for(auto o : observations) {
          double distance = dist(p.x, p.y, o.x, o.y);
          if(distance < dist_min) {
          	o.id = p.id; 
            dist_min = distance;
          }
          
        }
    }
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
                                   const vector<LandmarkObs> &observations, 
                                   const Map &map_landmarks) {
  /**
   * TODO: Update the weights of each particle using a mult-variate Gaussian 
   *   distribution. You can read more about this distribution here: 
   *   https://en.wikipedia.org/wiki/Multivariate_normal_distribution
   * NOTE: The observations are given in the VEHICLE'S coordinate system. 
   *   Your particles are located according to the MAP'S coordinate system. 
   *   You will need to transform between the two systems. Keep in mind that
   *   this transformation requires both rotation AND translation (but no scaling).
   *   The following is a good resource for the theory:
   *   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
   *   and the following is a good resource for the actual equation to implement
   *   (look at equation 3.33) http://planning.cs.uiuc.edu/node99.html
   */
  	double sum = 0.0;
  	for(int i = 0; i < num_particles; i++) {
      Particle *p = &particles[i];
      double weight = 1.0;
      for(auto o : observations) {
        double t_x = particles[i].x + cos(particles[i].theta) * o.x - sin(particles[i].theta) * o.y;
        double t_y = particles[i].y + sin(particles[i].theta) * o.x + cos(particles[i].theta) * o.y;
        
        Map::single_landmark_s nearest_l;
        double dist_min = sensor_range;
        for(auto l : map_landmarks.landmark_list) {
          double distance = dist(t_x, t_y, l.x_f,  l.y_f);
          if(distance < dist_min) {
            nearest_l = l;
            dist_min = distance;
          }
        }
        
        weight *= multiv_prob(std_landmark[0], std_landmark[1], t_x, t_y, nearest_l.x_f, nearest_l.y_f);

      }
      p->weight = weight;
      sum += weight;
      
    }
  	for(int i = 0; i < num_particles; i++) {
      	Particle *p = &particles[i];
       	p->weight /= sum;
      	weights[i] = p->weight;
      
    }
}

void ParticleFilter::resample() {
  /**
   * TODO: Resample particles with replacement with probability proportional 
   *   to their weight. 
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   */
  std::default_random_engine gen;
  std::discrete_distribution<int> discr(weights.begin(), weights.end());
  
  vector<Particle> resampled;
	for(int i = 0; i < num_particles; i++) {
      resampled.push_back(particles[discr(gen)]);
    }
  particles = resampled;
}

void ParticleFilter::SetAssociations(Particle& particle, 
                                     const vector<int>& associations, 
                                     const vector<double>& sense_x, 
                                     const vector<double>& sense_y) {
  // particle: the particle to which assign each listed association, 
  //   and association's (x,y) world coordinates mapping
  // associations: The landmark id that goes along with each listed association
  // sense_x: the associations x mapping already converted to world coordinates
  // sense_y: the associations y mapping already converted to world coordinates
  particle.associations= associations;
  particle.sense_x = sense_x;
  particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best) {
  vector<int> v = best.associations;
  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<int>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}

string ParticleFilter::getSenseCoord(Particle best, string coord) {
  vector<double> v;

  if (coord == "X") {
    v = best.sense_x;
  } else {
    v = best.sense_y;
  }

  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<float>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}