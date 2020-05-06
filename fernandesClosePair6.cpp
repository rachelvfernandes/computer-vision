/* Rachel Fernandes pd 4
 * Closest Pair Part 6 Final Improvement
 * I will uphold academic and personal integrity in the TJ community.
 */
#include <stdio.h>      /* printf, scanf, puts, NULL */
#include <stdlib.h>     /* srand, rand, RAND_MAX */
#include <time.h>       /* time */
#include <math.h>
#include <vector>
#include <map>
#include <fstream>

using namespace std;

class Point {
public:
    double x, y;
    Point() {
        x = (double) rand() / RAND_MAX;
        y = (double) rand() / RAND_MAX;
    }
    Point(double mX, double mY) {
        x = mX;
        y = mY;
    }
    bool operator<(Point p) {
        return x < p.x;
    }
};

double calcDistance(Point point1, Point point2)
{
    double distx = point1.x - point2.x;
    double disty = point1.y - point2.y;
    return (distx * distx) + (disty * disty);
}

void findClosestPair(int n, vector<Point> randomPoints)
{
    //assume point list is already randomized because it is made of randomly generated points
    // if not randomized, implement Knuth scrambling algorithm
    Point p1 = randomPoints[0];
    Point p2 = randomPoints[1];
    double currMinSquaredDist = calcDistance(p1, p2);
    double currMinDist = sqrt(currMinSquaredDist); //delta
    double halfMinDist = currMinDist/2.0;
    map <pair <int, int>, vector<Point>> subSquares;
    //keys are a pair of indexes of sub-squares and values are points from the randomized list in that sub-square.
    pair <int, int> tempSubSquare;
    double tempSquaredDist;
    int subSquareX;
    int subSquareY;
    bool foundCloserPair;
    Point currPoint;

    for(int it = 0; it < n; it++)
    {
        currPoint = randomPoints[it];
        subSquareX = (int)(currPoint.x/halfMinDist);
        subSquareY = (int)(currPoint.y/halfMinDist);
        foundCloserPair = false;
        for(int itx = -2; itx < 3; itx++)
            for(int ity = -2; ity < 3; ity++)
            {
                tempSubSquare = make_pair(subSquareX + itx, subSquareY + ity);
                if(subSquares.find(tempSubSquare) != subSquares.end()) //check that the subsquare is in the dictionary!!
                    for(Point subPoint : subSquares[tempSubSquare])
                    {
                        tempSquaredDist = calcDistance(currPoint, subPoint);
                        if (tempSquaredDist < currMinSquaredDist) //if a smaller distance is found
                        {
                            foundCloserPair = true;
                            currMinSquaredDist = tempSquaredDist;
                            p1 = currPoint;
                            p2 = subPoint;
                        }
                    }
            }
        if(foundCloserPair)
        {
            //remake the dictionary
            subSquares.clear();
            currMinDist = sqrt(currMinSquaredDist);
            halfMinDist = currMinDist/2.0;
            for(int e = 0; e < it+1; e++)
            {
                tempSubSquare = make_pair((int)(randomPoints[e].x/halfMinDist), (int)(randomPoints[e].y/halfMinDist));
                subSquares.insert(pair<pair<int, int>, vector<Point>>(tempSubSquare, {randomPoints[e]}));
            }
        }
        else
        {
            tempSubSquare = make_pair(subSquareX, subSquareY);
            //if(subSquares.find(tempSubSquare) != subSquares.end())
            //{
              //  subSquares[tempSubSquare].emplace_back(randomPoints[it]);
                //impossible?
            //}
            //else
                subSquares.insert(pair<pair<int, int>, vector<Point>>(tempSubSquare, {randomPoints[it]}));
        }
    }
    printf("(%f, %f) (%f, %f)", p1.x, p1.y, p2.x, p2.y);
}


int main() {
    srand(time(nullptr)); /* set random seed once using clock */
    vector <Point> randPoints;
    int N;
    ifstream inFile;
    inFile.open("points2.txt");
    if (inFile)
    {
        double value;
        double value2;
        if(inFile >> N)
            while ( inFile >> value)
                if(inFile >> value2)
                    randPoints.emplace_back(value, value2);

         findClosestPair(N, randPoints);
    }
}