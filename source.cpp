#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include "kernel.hpp"


int main()
{
  float E = - 1.6 ;
  float ME = 9.1 * pow(10, -3);

  float mass = ME;
  float charge = E;


  std::vector<Particle> particles;
  //for(size_t x = -1; x <= 1 ; ++x)
    //particles.emplace_back(20, 1, 0, 0, x, x); //mass, charge, vx, vy, x, y


  particles.emplace_back(10,-1, 0.01, 0.05, 0.5, 0);
  particles.emplace_back(10, 1, 0, -0.05,-0.5, 0);
  /*particles.emplace_back(10,-1, -0.1, 0.1, 0, 0);
  particles.emplace_back(10,-1, 0.1, -0.7, 0.5, 2);
  particles.emplace_back(10,-1, 0.4, -0.7, -0.5, 2);
  particles.emplace_back(10,-1, 0, -1, 0, 3);
  particles.emplace_back(10,-1, -0.2, 0.2, -1, 1);
*/



  time(0.0005, 0, 15, particles);
  system("python3 path.py");
  return 0;

}
