/* Rachel Fernandes pd 4
 * Convex Hull - Quick Hull Algorithm
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
    double x, y;
    Point()
    {
        x = (double) rand() / RAND_MAX;
        y = (double) rand() / RAND_MAX;
    }
    Point(double mX, double mY)
    {
        x = mX;
        y = mY;
    }
    double lineDistance(double m, Point p, double mSquarePlusOne)
    {
        return abs(m*(x -p.x) - y + p.y)/mSquarePlusOne;
    }
    bool operator<(Point p)
    {
        return x < p.x;
    }
};

vector<Point> quickConvexHull;
int coordPlane[500][500] = {};
Point p0;

int plot(double x) {
    return (int)(450* x + 25.5);
}

bool angleOperator(Point point1, Point point2)
{
    return atan((p0.y-point1.y)/(p0.x-point1.x)) < atan((p0.y-point2.y)/(p0.x-point2.x));
}

Point furthestPoint(vector<Point> points, double m, Point p)
{
    double mSquaredPlusOne = sqrt(m*m + 1);
    Point far = points.front();
    double farDist = far.lineDistance(m, p, mSquaredPlusOne);
    double tempDist;
    for(Point point: points)
    {
        tempDist = point.lineDistance(m, p, mSquaredPlusOne);
        if(tempDist > farDist)
        {
            farDist = tempDist;
            far = point;
        }
    }
    return far;
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

void findHull(vector<Point> points, Point p1, Point p2, bool above)
{
    if(points.empty())
        return;
    else
    {
        double m = (p1.y-p2.y)/(p1.x-p2.x); // delta y over delta x
        Point c = furthestPoint(points, m, p1);
        quickConvexHull.emplace_back(c);
        double mP1C = (p1.y-c.y)/(p1.x-c.x);
        double mP2C = (p2.y-c.y)/(p2.x-c.x);
        vector<Point> outsideLeft;
        vector<Point> outsideRight;
        if(above)
        {
            for(Point p: points)
                if(p.x < c.x and p.y > mP1C*(p.x - p1.x) + p1.y)
                    outsideLeft.emplace_back(p);
                else if(p.x > c.x and p.y > mP2C*(p.x - p2.x) + p2.y)
                    outsideRight.emplace_back(p);
        }
        else
        {
            for(Point p: points)
                if(p.x < c.x and p.y < mP1C*(p.x - p1.x) + p1.y)
                    outsideLeft.emplace_back(p);
                else if(p.x > c.x and p.y < mP2C*(p.x - p2.x) + p2.y)
                    outsideRight.emplace_back(p);
        }
        findHull(outsideLeft, p1, c, above);
        findHull(outsideRight, c, p2, above);
    }
}

void quickHull(vector<Point> points)
{
    //find points with min and max x (they are always in the hull)
    sort(points.begin(), points.end());
    p0 = points.front();
    Point p1 = points.front();
    Point p2 = points.back();
    quickConvexHull.emplace_back(p1);
    quickConvexHull.emplace_back(p2);
    //split the set along the line segment connecting those two points
    double m = (p1.y-p2.y)/(p1.x-p2.x);
    vector<Point> aboveSegment;
    vector<Point> belowSegment;
    for(Point p: points)
        if(p.y > m*(p.x - p1.x) + p1.y) //if the point is above the line
            aboveSegment.emplace_back(p);
        else
            belowSegment.emplace_back(p);
    //for each side of the line segments, cancel the points within the triangle formed by the line segment and the furthest point from the line segment
    //repeat the above step where the line segment is each of the two new sides of the triangle
    findHull(aboveSegment, p1, p2, true);
    findHull(belowSegment, p1, p2, false);
}

int main() {
    srand(time(nullptr)); /* set random seed once using clock */
    int N = 75;
    vector<Point> allPoints;


    for (int i = 0; i < N; i++)
    {
        allPoints.emplace_back(Point());
        coordPlane[plot(allPoints[i].x)][plot(allPoints[i].y)] = 1;
    }

    quickHull(allPoints);
    sort(quickConvexHull.begin()+1, quickConvexHull.end(), angleOperator);
    for(int iter = 0; iter < quickConvexHull.size()-1; iter++)
        bresenhamLine(plot(quickConvexHull[iter].x), plot(quickConvexHull[iter+1].x), plot(quickConvexHull[iter].y), plot(quickConvexHull[iter+1].y));
    bresenhamLine(plot(quickConvexHull[0].x), plot(quickConvexHull[quickConvexHull.size()-1].x), plot(quickConvexHull[0].y), plot(quickConvexHull[quickConvexHull.size()-1].y));


    for(Point inHull: quickConvexHull)
    {
        //printf("( %f , %f )", inHull.x, inHull.y);
        coordPlane[plot(inHull.x)][plot(inHull.y)] = 2;
        coordPlane[plot(inHull.x)][1+plot(inHull.y)] = 2;
        coordPlane[plot(inHull.x)][-1+plot(inHull.y)] = 2;
        coordPlane[-1+plot(inHull.x)][plot(inHull.y)] = 2;
        coordPlane[1+plot(inHull.x)][plot(inHull.y)] = 2;
    }


    printf("P3 500 500 1\n");

    for (int j = 0; j < 500; j++)
        for (int m = 0; m < 500; m++)
            if (coordPlane[j][m] == 0)
                printf("%d %d %d ", 1, 1, 1);
            else if (coordPlane[j][m] == 2)
                printf("%d %d %d ", 1, 0, 0);
            else
                printf("%d %d %d ", 0, 0, 0);

}