//
// Created by Wei Chen on 3/17/22
//

#include <Eigen/Eigen>

#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

int width = 720; // width and height of window
int height = 720;

void draw(){
    // glClearColor(0.1, 0.1, 0.1, 1.0);
    glClear(GL_COLOR_BUFFER_BIT); // clear buffer

    Eigen::MatrixXf tmp(width*3, height);
    for(int i=100; i<200; ++i){
        for(int j=600; j<700; ++j){
            tmp(Eigen::seq(i*3, (i+1)*3-1), j) = Eigen::Vector3f(0.6, 0.2, 0.2);
        }
    }
    float *pixel = tmp.data();
    // NOTE: pixel must be float, double will cause wrong!!!
    glDrawPixels(width, height, GL_RGB, GL_FLOAT, pixel); // draw pixels

    glutSwapBuffers(); // change two buffer
}

int main(int argc, char **argv){
    glutInit(&argc, argv);
    // setting up the display RGB color model
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);
    // configure window position and size
    glutInitWindowPosition(200, 200);
    glutInitWindowSize(width, height);
    // create window
    glutCreateWindow("draw picture");
    // call the drawing function
    glutDisplayFunc(draw);
    // keyboard control
    // glutKeyboardFunc(keyboard); // doesn't implement yet

    // loop require by opengl
    glutMainLoop();

    return 0;
}