//Rachel Fernandes
//Assignment 1 - Triangles and Circles
//Computer Vision pd 4
//9/22/19

//run with following lines of code:
// gcc trianglesAndCircles.cpp -lm
// ./a/out > t.ppm
//display t.ppm

#include <stdio.h>      /* printf, scanf, puts, NULL */
#include <stdlib.h>     /* srand, rand, RAND_MAX */
#include <time.h>       /* time */
#include <math.h>


void fillRandom(double* arr) {
    for (int i=0; i<6; i++) arr[i] = (double)rand()/RAND_MAX;
    return;
}

double calculateDistance(double x1, double y1, double x2, double y2){
    return sqrt(pow((x1 - x2), 2) + pow((y1 - y2), 2));
}

double calculateSlope(double x1, double y1, double x2, double y2){
    return (y1 - y2) / (x1 - x2);
}

int scale(double xc)
{
    return (int)(600.0*xc);
}

void calculateMidpoint(double* arrr, double x1, double y1, double x2, double y2){
    arrr[0] = (double)(x1 + x2)/2;
    arrr[1] = (double)(y1 + y2)/2;
}

void plot(int myx, int myy, bool arr[1000][1000])
{
    if(myx<900 && myy<900 && myx > -200 && myy > -200)
        arr[myx+200][myy+200] = true;
}

void drawCircle(int cX, int cY, int r, bool myA[1000][1000]){
    int xmax = (int) (r * 0.64); // maximum x at radius/sqrt(2)
    int x;
    int y = r;
    int y2 = y * y;
    int ty = (2 * y) - 1;
    int y2_new = y2;
    for (x = 0; x <= xmax+(2+r/15); x++) {
        if ((y2 - y2_new) >= ty) {
            y2 -= ty;
            y -= 1;
            ty -= 2;
        }
        plot(x+cX, y+cY, myA);
        plot(x+cX, -y+cY, myA);
        plot(-x+cX, y+cY, myA);
        plot(-x+cX, -y+cY, myA);
        plot(y+cX, x+cY, myA);
        plot(y+cX, -x+cY, myA);
        plot(-y+cX, x+cY, myA);
        plot(-y+cX, -x+cY, myA);
        y2_new -= (2 * x) - 3;
    }
}


void bresenhamLine(int xi, int xii, int yi, int yii, int drivingAxis, bool myArr[1000][1000]){
    //I won't go into detail explaining this because its not my algorithm. Thanks Jack Elton Bresenham.
    int dx, dy, jj, ii, maxii;
    double m, ee;
        if(drivingAxis == 0)
        {
            dx = xii - xi;
            dy = yii - yi;
            m = (double) dy / (double) dx;
            jj = yi;
            ee = abs(m) - 1.0;
            maxii = xii - 1;
            for (ii = xi; ii < maxii+1; ii++) {
                plot(ii, jj, myArr);
                if (ee >= 0) {
                    if (dy < 0)
                        jj--;
                    else
                        jj++;
                    ee = ee - 1.0;
                }
                ee = ee + abs(m);
            }
        }
        if(drivingAxis == 1)
        {
            dx = xii - xi;
            dy = yii - yi;
            m = (double) dx / (double) dy;
            jj = xi;
            ee = abs(m) - 1.0;
            maxii = yii - 1;
            for (ii = yi; ii < maxii+1; ii++) {
                plot(jj, ii, myArr);
                if (ee >= 0) {
                    if (dx < 0)
                        jj--;
                    else
                        jj++;
                    ee = ee - 1.0;
                }
                ee = ee + abs(m);
            }
        }
}

int main ()
{
    srand (time(NULL));   /* set random seed once using clock */
    printf("P3 1000 1000 1\n");

//origin of coordinate plane is at (100, 700)
//allow 100 pixels on each side so circles don't get cut off
//subtract y coordinate from 700 to get coordinate on ppm file
//add 100 to x coordinate to get coordinate on ppm file

    bool coordPlane[1000][1000];
    for (int r=0; r < 1000; r++)
        for(int t=0; t < 1000; t++)
        {
            coordPlane[r][t] = false;
        }

    double randomDoubles[6];  /* declare the array,  */
    fillRandom(randomDoubles);  /* and pass it to be filled  */
    plot(scale(randomDoubles[0]), scale(randomDoubles[1]), coordPlane);
    plot(scale(randomDoubles[2]), scale(randomDoubles[3]), coordPlane);
    plot(scale(randomDoubles[4]), scale(randomDoubles[5]), coordPlane);

    //Bresenham algorithm

    if(abs(randomDoubles[0]-randomDoubles[2]) > abs(randomDoubles[1]-randomDoubles[3])){
        if (randomDoubles[2] > randomDoubles[0])
            bresenhamLine(scale(randomDoubles[0]), scale(randomDoubles[2]), scale(randomDoubles[1]), scale(randomDoubles[3]), 0, coordPlane);
        else
            bresenhamLine(scale(randomDoubles[2]), scale(randomDoubles[0]), scale(randomDoubles[3]), scale(randomDoubles[1]), 0, coordPlane);
    }
    else{
        if (randomDoubles[3] > randomDoubles[1])
            bresenhamLine(scale(randomDoubles[0]), scale(randomDoubles[2]), scale(randomDoubles[1]), scale(randomDoubles[3]), 1, coordPlane);
        else
            bresenhamLine(scale(randomDoubles[2]), scale(randomDoubles[0]), scale(randomDoubles[3]), scale(randomDoubles[1]), 1, coordPlane);
    }


    if(abs(randomDoubles[0]-randomDoubles[4]) > abs(randomDoubles[1]-randomDoubles[5])){
        if (randomDoubles[4] > randomDoubles[0])
            bresenhamLine(scale(randomDoubles[0]), scale(randomDoubles[4]), scale(randomDoubles[1]), scale(randomDoubles[5]), 0, coordPlane);
        else
            bresenhamLine(scale(randomDoubles[4]), scale(randomDoubles[0]), scale(randomDoubles[5]), scale(randomDoubles[1]), 0, coordPlane);
    }
    else{
        if (randomDoubles[5] > randomDoubles[1])
            bresenhamLine(scale(randomDoubles[0]), scale(randomDoubles[4]), scale(randomDoubles[1]), scale(randomDoubles[5]), 1, coordPlane);
        else
            bresenhamLine(scale(randomDoubles[4]), scale(randomDoubles[0]), scale(randomDoubles[5]), scale(randomDoubles[1]), 1, coordPlane);
    }

    if(abs(randomDoubles[2]-randomDoubles[4]) > abs(randomDoubles[3]-randomDoubles[5])){
        if (randomDoubles[4] > randomDoubles[2])
            bresenhamLine(scale(randomDoubles[2]), scale(randomDoubles[4]), scale(randomDoubles[3]), scale(randomDoubles[5]), 0, coordPlane);
        else
            bresenhamLine(scale(randomDoubles[4]), scale(randomDoubles[2]), scale(randomDoubles[5]), scale(randomDoubles[3]), 0, coordPlane);
    }
    else{
        if (randomDoubles[5] > randomDoubles[3])
            bresenhamLine(scale(randomDoubles[2]), scale(randomDoubles[4]), scale(randomDoubles[3]), scale(randomDoubles[5]), 1, coordPlane);
        else
            bresenhamLine(scale(randomDoubles[4]), scale(randomDoubles[2]), scale(randomDoubles[5]), scale(randomDoubles[3]), 1, coordPlane);
    }

    //side lengths
    double distAB = calculateDistance(randomDoubles[0], randomDoubles[1], randomDoubles[2], randomDoubles[3]);
    double distBC = calculateDistance(randomDoubles[0], randomDoubles[1], randomDoubles[4], randomDoubles[5]);
    double distAC = calculateDistance(randomDoubles[2], randomDoubles[3], randomDoubles[4], randomDoubles[5]);
    //x and y coordinates of midpoints of sides
    double midx1 = (randomDoubles[0] + randomDoubles[2])/2.0;
    double midy1 = (randomDoubles[1] + randomDoubles[3])/2.0;
    double midx2 = (randomDoubles[0] + randomDoubles[4])/2.0;
    double midy2 = (randomDoubles[1] + randomDoubles[5])/2.0;
    //slopes of lines perpendicular to sides (this includes altitudes and perpendicular bisectors
    double slope1 = -1.0/calculateSlope(randomDoubles[0], randomDoubles[1], randomDoubles[2], randomDoubles[3]);
    double slope2 = -1.0/calculateSlope(randomDoubles[0], randomDoubles[1], randomDoubles[4], randomDoubles[5]);
    //orthocenter will be used later to find nine point center
    double orthoCenterX = (randomDoubles[3]-(slope2*randomDoubles[2])+(slope1*randomDoubles[4])-randomDoubles[5])/(slope1-slope2);
    double orthoCenterY = (slope1*orthoCenterX)-(slope1*randomDoubles[4])+randomDoubles[5];

    double s = 0.5 * (distAB + distBC + distAC);
    double iR = sqrt(((s - distAB)*(s - distBC)*(s-distAC))/s);
    int inRad = scale(iR);
    int inCenterX = scale((distAB*randomDoubles[4] + distBC*randomDoubles[2] + distAC*randomDoubles[0])/(distAB+distBC+distAC));
    int inCenterY = scale((distAB*randomDoubles[5] + distBC*randomDoubles[3] + distAC*randomDoubles[1])/(distAB+distBC+distAC));
    int circumRad = scale((distAB*distBC*distAC)/(4.0*iR*s));
    double ccX = (midy2-(slope2*midx2)+(slope1*midx1)-midy1)/(slope1-slope2);
    double ccY = (slope1*ccX)-(slope1*midx1)+midy1;
    int circumCenterX = scale(ccX);
    int circumCenterY = scale(ccY);
    //nine point circle center lies on the midpoint of the circumCenter and orthocenter
    double ncX = (ccX + orthoCenterX)/2.0;
    double ncY = (ccY + orthoCenterY)/2.0;
    int nineCenterX = scale(ncX);
    int nineCenterY = scale(ncY);
    //nine point circle intersects midpoints of the edges so the radius should be the distance between the ninepoint center and one midpoint
    int nineRad = scale(calculateDistance(ncX, ncY, midx1, midy1));
    //draw all the circles
    drawCircle(inCenterX, inCenterY, inRad, coordPlane);
    drawCircle(circumCenterX, circumCenterY, circumRad, coordPlane);
    drawCircle(nineCenterX, nineCenterY, nineRad, coordPlane);

    //now that all the points are calculated, print them to the ppm file as pixels
    for (int j = 0; j < 1000; j++)
        for(int m = 0; m < 1000; m++)
        {
            if(coordPlane[j][m])
                printf("%d %d %d ", 0, 0, 0);
            else
                printf("%d %d %d ", 1, 1, 1);
        }
}