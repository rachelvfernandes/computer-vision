/* Rachel Fernandes pd 4
 * Convex Hull - Graham Scan Algorithm
 * Nov 25, 2019
 * I will uphold academic and personal integrity in the TJ community.
 */
#include <stdlib.h>     /* srand, rand, RAND_MAX */
#include <time.h>       /* time */
#include <vector>
#include <math.h>
#include <stdio.h>      /* printf, scanf, puts, NULL */
#include <algorithm>

using namespace std;

class Point {
public:
    double x, y, slopeFromP0;
    Point()
    {
        x = (double) rand() / RAND_MAX;
        y = (double) rand() / RAND_MAX;
        slopeFromP0 = 0;
    }
    Point(double mX, double mY)
    {
        x = mX;
        y = mY;
        slopeFromP0 = 0;
    }
    double lineDistance(double m, Point p, double mSquarePlusOne)
    {
        return abs(m*(x -p.x) - y + p.y)/mSquarePlusOne;
    }
    bool operator<(Point p)
    {
        return x < p.x;
    }
    double getArctanSlopeFromP0(Point myP0)
    {
        if (slopeFromP0 == 0)
            slopeFromP0 = atan((myP0.y-y)/(myP0.x-x));
        return slopeFromP0;
    }
};

vector<Point> grahamScanConvexHull;
int coordPlane[500][500] = {};
Point p0;

bool angleOperator(Point point1, Point point2)
{
    return atan((p0.y-point1.y)/(p0.x-point1.x)) < atan((p0.y-point2.y)/(p0.x-point2.x));
}

int plot(double x)
{
    return (int)(450 * x + 25.5);
}

void bresenhamLine(int xi, int xii, int yi, int yii)
{
    //I won't go into detail explaining this because its not my algorithm. Thanks Jack Elton Bresenham.
    int dx = xii - xi;
    int dy = yii - yi;
    int jj, ii;
    if(abs(dx) > abs(dy))
    {
        if(xi > xii) //make sure first point is the left most
        {
            swap(xi, xii);
            swap(yi, yii);
            dx = -dx;
            dy = -dy;
        }
        double e = dy - dx;
        jj = yi;
        for (ii = xi; ii < xii-1; ii++) {
            coordPlane[ii][jj] = 2;
            if (e >= 0)
            {
                if (dy < 0)
                    jj--;
                else
                    jj++;
                e = e - dx;
            }
            e = e + abs(dy);
        }
    }
    else
    {
        if(yi > yii) //make sure first point is the lower point
        {
            swap(xi, xii);
            swap(yi, yii);
            dx = -dx;
            dy = -dy;
        }
        double e = dy - dx;
        jj = xi;
        for (ii = yi; ii < yii-1; ii++) {
            coordPlane[jj][ii] = 2;
            if (e >= 0)
            {
                if (dx < 0)
                    jj--;
                else
                    jj++;
                e = e - dy;
            }
            e = e + abs(dx);
        }
    }
}

double ccw(Point point1, Point point2, Point point3)
{
    return (point2.x - point1.x)*(point3.y - point1.y) - (point2.y - point1.y)*(point3.x - point1.x);
}

void grahamScan(vector<Point> points)
{
    //find left most point with min y coordinate (always in the hull)
    sort(points.begin(), points.end());
    p0 = points.front();
    grahamScanConvexHull.emplace_back(p0);
    //sort the points in increasing polar angle from P0, remove the closer of any points with the same angle
    sort(points.begin()+1, points.end(), angleOperator);
    //grahamScanConvexHull is a vector which I am treating like a stack temporarily
    for(Point p: points)
    {
        while(grahamScanConvexHull.size() > 1 and ccw(grahamScanConvexHull[grahamScanConvexHull.size()-2], grahamScanConvexHull.back(), p) < 0)
            grahamScanConvexHull.pop_back();
        grahamScanConvexHull.emplace_back(p);
    }
    sort(grahamScanConvexHull.begin()+1, grahamScanConvexHull.end(), angleOperator);
}

int main()
{
    srand(time(nullptr)); /* set random seed once using clock */
    int N = 75;
    vector<Point> allPoints;

    for (int i = 0; i < N; i++)
    {
        allPoints.emplace_back(Point());
        coordPlane[plot(allPoints[i].x)][plot(allPoints[i].y)] = 1;
    }

    grahamScan(allPoints);
    sort(grahamScanConvexHull.begin()+1, grahamScanConvexHull.end(), angleOperator);
    for(int iter = 0; iter < grahamScanConvexHull.size()-1; iter++)
        bresenhamLine(plot(grahamScanConvexHull[iter].x), plot(grahamScanConvexHull[iter+1].x), plot(grahamScanConvexHull[iter].y), plot(grahamScanConvexHull[iter+1].y));
    bresenhamLine(plot(grahamScanConvexHull[0].x), plot(grahamScanConvexHull[grahamScanConvexHull.size()-1].x), plot(grahamScanConvexHull[0].y), plot(grahamScanConvexHull[grahamScanConvexHull.size()-1].y));


    for(Point inHull: grahamScanConvexHull)
    {
        //printf("( %f , %f )", inHull.x, inHull.y);
        coordPlane[plot(inHull.x)][plot(inHull.y)] = 2;
        coordPlane[plot(inHull.x)][1+plot(inHull.y)] = 2;
        coordPlane[plot(inHull.x)][-1+plot(inHull.y)] = 2;
        coordPlane[-1+plot(inHull.x)][plot(inHull.y)] = 2;
        coordPlane[1+plot(inHull.x)][plot(inHull.y)] = 2;
    }
    coordPlane[plot(grahamScanConvexHull[0].x)][plot(grahamScanConvexHull[0].y)] = 3;
    printf("P3 500 500 1\n");
    for (int j = 0; j < 500; j++)
        for (int m = 0; m < 500; m++)
            if (coordPlane[j][m] == 0)
                printf("%d %d %d ", 1, 1, 1);
            else if (coordPlane[j][m] == 2)
                printf("%d %d %d ", 1, 0, 0);
            else if (coordPlane[j][m] == 3)
                printf("%d %d %d ", 0, 1, 0);
            else
                printf("%d %d %d ", 0, 0, 0);
}