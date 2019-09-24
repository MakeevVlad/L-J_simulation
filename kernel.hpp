#pragma ones

#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>

//Defines static charge
struct stch
{
  float q, x, y;
};

//Defines classical particle
class Particle
{
private:
  std::vector<float> state; //State vector (mass, charge, V_x, V_y, x, y)

public:
  Particle(float mass, float charge, float vx, float vy, float x, float y)
  {
      state.resize(6);
      state[0] = mass;
      state[1] = charge;
      state[2] = vx;
      state[3] = vy;
      state[4] = x;
      state[5] = y;
  }
//"Get" functions
float& m()      { return state[0]; }
float& q()      { return state[1]; }
float& vx()     { return state[2]; }
float& vy()     { return state[3]; }
float& x()      { return state[4]; }
float& y()      { return state[5]; }

};

//Saving data to the file
void collect(Particle p)
{
  std::ofstream file("data.txt", std::ios::app);
  file << p.x() << "  " << p.y( ) << "  " << p.vx() << "  " << p.vy()<<std::endl;
  file.close();
}

float axeleration(char axis, char type, size_t num, std::vector<Particle>& particles)
{
  //constants
  float k  = 1; //*9* pow(10, 9);
  float eps = 100;
  float min = 1;
  float r6 = pow(min, 6);

  float dx;
  float dy;
  float r;
  float a = 0;
  Particle p = particles[num];

  for(size_t i = 0; i < particles.size(); ++i)
  {
    if(i == num) continue;
    Particle cur = particles[i];

    dx = p.x() - cur.x();
    dy = p.y() - cur.y();
    r = sqrt(dx*dx + dy*dy);

    switch(axis)
    {
      case 'x':
        switch(type)
        {
          case 'k':
            a +=  cur.q() * dx / pow(r, 3);
            break;

          case 'l':
            a -= dx * ( 1 / pow(r, 8) - r6 / pow(r, 14));
            break;
        }
        break;
//*****************************************
      case 'y':
      switch(type)
      {
        case 'k':
          a +=  cur.q() * dy / pow(r, 3);
          break;

        case 'l':
          a -=  dy * ( 1 / pow(r, 8) - r6 / pow(r, 14));
          break;
      }
      break;
    }
  }
  switch(type)
  {
    case 'k':
      a *= k * p.q() / p.m();
      break;

    case 'l':
      a *= eps *12 * r6 / p.m();
      break;

  }
  return a;
}




//"Time" function (events recursion)
void time(float t0, float t, float maxt, std::vector<Particle>& particles)
{
  std::vector<float> ax(particles.size());
  std::vector<float> ay(particles.size());
  for(size_t num = 0; num < particles.size(); ++num)
  {
    ax[num] = axeleration('x', 'l', num, particles);
    ay[num] = axeleration('y', 'l', num, particles);
  }

  size_t num = 0;
  for(auto& p: particles)
  {

    p.x() = p.x() + p.vx() * t0 +  ax[num] * t0*t0 * 0.5;
    p.y() = p.y() + p.vy() * t0 +  ay[num] * t0*t0 * 0.5;

    p.vx() = p.vx() + t0 * ax[num];
    p.vy() = p.vy() + t0 * ay[num];

    collect(p);

    ++num;
  }

  if (t < maxt)
  {
    t += t0;
    time(t0, t, maxt, particles);
  }

}
