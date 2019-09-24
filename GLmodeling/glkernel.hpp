#include <GL/freeglut.h>
#include <GL/gl.h>
#include <GL/glu.h>

#include <omp.h>
#include <vector>
#include <cmath>
#include <fstream>
#include <iostream>


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
//constants
float k  = 1; //*9* pow(10, 9);
float eps = 10000;
float min = 20;
float r6 = pow(min, 6);
float t0 = 0.00001;
size_t curt = 0;


GLfloat rsize = 1.0f;


GLfloat windowWidth;
GLfloat windowHeight;

std::vector<Particle> particles;
void cr()
{
  int n = 3;

  for(int y = -n; y <= n; ++y)
    for(int x = (y % 2 ? -n : -n + 1); x<= (y % 2? n : n - 1); ++x)
    {
      if(x % 2 == 0)
        particles.emplace_back(1,-1, 0, 0, min * x, min * y * sqrt(3) * 0.5 );
      else
        particles.emplace_back(1,-1, 0, 0, min*0.5 + x*min, min * y * sqrt(3) * 0.5);
    }

      //particles.emplace_back(10,-1, 10, -500, 20, 0);

}



void RenderScene(void)
{
    glClear(GL_COLOR_BUFFER_BIT);

    #pragma omp parallel for ordered schedule(dynamic)
    for(auto& p: particles)
      {
        glColor3f(1.0f,1.0f,1.0f);
        //рисуем круг посредством соединения наружных точек с центральной
        glBegin( GL_TRIANGLE_FAN );
        glVertex2f(p.x(), p.y());
        #pragma omp parallel for ordered schedule(dynamic)
        for( size_t j = 0; j < 100; ++j )
        {
          float a = (float)j * 3.1415f * 2.0f/100;

          glVertex2f( cos( a ) * rsize + p.x(), sin( a ) * rsize + p.y());
        }
        glVertex2f(rsize + p.x(), p.y());
        glEnd();
      }

    glutSwapBuffers();
}



//Saving data to the file
void collect(Particle p)
{
  std::ofstream file("data.txt", std::ios::app);
  file << p.x() << "  " << p.y( ) << "  " << p.vx() << "  " << p.vy()<<std::endl;
  file.close();
}
void collectenergy()
{
  //std::ofstream Ken("KenVStime.txt", std::ios::app);
  //std::ofstream Wen("WenVStime.txt", std::ios::app);
  std::ofstream En ("EnVStime.txt", std::ios::app);

  float W = 0;
  float K = 0;
  float r;
  for(size_t i = 0; i < particles.size(); ++i)
  {
    K += (pow(particles[i].vx(), 2) + pow(particles[i].vy(), 2)) * particles[i].m();
    for(size_t j = 0; j < particles.size(); ++j)
    {
      if(i == j) continue;
      r = sqrt( pow( (particles[i].x() - particles[j].x() ), 2) +  pow( (particles[i].y() - particles[j].y() ), 2));
      W+= r6 / pow(r, 12) - 2 / pow(r, 6);
    }
  }
  //Ken << K / 2 << " " << curt << std::endl;
  //Wen << W * eps * r6<< " " << curt << std::endl;
  En << K/2 + W * eps * r6 << " " << curt << std::endl;
  //Ken.close();
  //Wen.close();
  En.close();
}


float axeleration(char axis, char type, size_t num)
{

  float dx;
  float dy;
  float r;
  float a = 0;
  Particle p = particles[num];
  #pragma omp parallel for ordered schedule(dynamic)
  for(size_t i = 0; i < particles.size(); ++i)
  {
    if(i == num) continue;
    Particle cur = particles[i];

    dx = p.x() - cur.x();
    dy = p.y() - cur.y();
    r = sqrt(dx*dx + dy*dy);

    if (r >= 2.23 * min)
      continue;

    switch(type)
    {
      case 'k':
        switch(axis)
        {
          case 'x':
            a +=  cur.q() * dx / pow(r, 3);
            break;

          case 'y':
          a +=  cur.q() * dy / pow(r, 3);
            break;
        }
        break;
      case 'l':
      switch(axis)
      {
        case 'x':
          a -= dx * ( 1 / pow(r, 8) - r6 / pow(r, 14));
          break;

        case 'y':
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


void TimerFunction(int value)
{
  std::vector<float> ax(particles.size());
  std::vector<float> ay(particles.size());

  #pragma omp parallel for ordered schedule(dynamic)
  for(size_t num = 0; num < particles.size(); ++num)
  {
    ax[num] = axeleration('x', 'l', num);
    ay[num] = axeleration('y', 'l', num);
  }

  size_t num = 0;
/*
  if(curt % 1000 == 0 ) {collectenergy();}
  curt += t0*1000;
*/
  #pragma omp parallel for ordered schedule(dynamic)
  for(auto& p : particles)
  {

    if (p.x() + rsize >= windowWidth || p.x() - rsize  <= -windowWidth)   /// проверка на достижение края экрана по OX
      p.vx() = -p.vx();
    if (p.y() + rsize >= windowHeight || p.y() - rsize <= -windowHeight )  /// ... по OY
      p.vy() = -p.vy();

    p.x() = p.x() + p.vx() * t0 +  ax[num] * t0*t0 * 0.5;
    p.y() = p.y() + p.vy() * t0 +  ay[num] * t0*t0 * 0.5;

    p.vx() = p.vx() + t0 * ax[num];
    p.vy() = p.vy() + t0 * ay[num];

    ++num;

/*
    if ( (p.x() + rsize) > (windowWidth + p.vx() * t0) )   /// проверяем не приведёт ли это к выходу за границу экрана
      p.x() = windowWidth - rsize - 1;
    else if ( (p.x() - rsize) < -windowWidth + p.vx() * t0)  // с другой стороны
      p.x() = - windowWidth + rsize + 1;

    if ( (p.y() + rsize) > (windowHeight + p.vy()) * t0)    ///по другим осям
      p.y() =  windowHeight - rsize - 1;

    else if( (p.y() - rsize) < - windowHeight + p.vy() * t0)
      p.y() = - windowHeight + rsize + 1;
*/
    glutPostRedisplay();  // перерисовываем экран
  }

  glutTimerFunc(t0 * 100,TimerFunction,1);  //запускаем таймер заново.
}



void SetupRC(void)
{
  glClearColor(0.0f, 0.0f, 0.0f, 0.0f);
  glClear(GL_COLOR_BUFFER_BIT);
}

void ChangeSize(GLsizei w, GLsizei h)
 {
    GLfloat aspectRatio;
    if (h == 0)
        h = 1;
    glViewport(0, 0, w, h);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    aspectRatio = (GLfloat)w / (GLfloat)h;
    if (w <= h){
        windowWidth = 100;
        windowHeight = 100 / aspectRatio;
        glOrtho(-100.0,100.0,-windowHeight,windowHeight,1.0,-1.0);
    } else {
        windowWidth = 100 * aspectRatio;
        windowHeight = 100;
        glOrtho(-windowWidth,windowWidth,-100.0,100.0,1.0,-1.0);
    }
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
}
