#include "glkernel.hpp"

#include <GL/freeglut.h>
#include <GL/gl.h>
#include <GL/glu.h>

#include <omp.h>
#include <cmath>
#include <iostream>



int main(int argc,char** argv) {

    glutInit(&argc,argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);
    glutCreateWindow("Bounce");
    cr();

    glutReshapeFunc(ChangeSize);
    glutDisplayFunc(RenderScene);

    glutTimerFunc(t0 * 100,TimerFunction,1);  //первый запуск таймера
    SetupRC();
    glutMainLoop();

    return 0;

}
