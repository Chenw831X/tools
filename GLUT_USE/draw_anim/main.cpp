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
int x = 0, y = 360, step = 1; // move settings
Eigen::MatrixXf canvas(3*width, height);

void draw(){
    // glClearColor(0.1, 0.1, 0.1, 1.0);
    glClear(GL_COLOR_BUFFER_BIT); // clear buffer

    float *pixel = canvas.data();
    // NOTE: pixel must be float, double will cause wrong!!!
    glDrawPixels(width, height, GL_RGB, GL_FLOAT, pixel); // draw pixels

    glutSwapBuffers(); // change two buffer
}

void TimerFunction(int value){
    if(x==620 && step == 1){ // avoid collision with window boundary
        step = -step;
    }
    if(x==0 && step == -1){
        step = -step;
    }
    x += step;

    canvas.setZero();
    for(int i=0; i<100; ++i){
        for(int j=0; j<100; ++j){
            canvas(Eigen::seq((x+i)*3, (x+i+1)*3-1), y+j) =
                    Eigen::Vector3f(0.6, 0.2, 0.2);
        }
    }

    glutPostRedisplay(); // call draw()
    // call TimerFunction(...) after 0ms
    glutTimerFunc(0, TimerFunction, value);
}

int main(int argc, char **argv){
    glutInit(&argc, argv);
    // setting double buffer
    // setting up the display RGB color model
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);
    // configure window position and size 
    glutInitWindowPosition(200, 200); 
    glutInitWindowSize(width, height); 
    // create window 
    glutCreateWindow("draw animation");
    // call the drawing function 
    glutDisplayFunc(draw);
    glutTimerFunc(0, TimerFunction, 1);
    // keyboard control
    // glutKeyboardFunc(keyboard); // doesn't implement yet

    // loop require by opengl 
    glutMainLoop();

    return 0;
}