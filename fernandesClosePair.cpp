/* Rachel Fernandes
 * Dr. Zacharias pd 4
 * Closest Pair Part 4
 * I will uphold academic and personal integrity in the TJ community.
 */
#include <stdio.h>      /* printf, scanf, puts, NULL */
#include <stdlib.h>     /* srand, rand, RAND_MAX */
#include <time.h>       /* time */
#include <math.h>
#include <bits/stdc++.h>
#include <vector>

using namespace std;
#define N 1000

class Point {
private:
    double x, y;
public:
    Point() {
        x = (double) rand() / RAND_MAX;
        y = (double) rand() / RAND_MAX;
    }

    Point(double mX, double mY) {
        x = mX;
        y = mY;
    }

    double getX() {
        return x;
    }

    double getY() {
        return y;
    }

    bool operator<(Point p) {
        return x < p.getX();
    }

};

static vector <Point> randomPoints;
static vector <Point> randomPointsY;

double calcDistances(Point point1, Point point2) {
    return sqrt(pow((point1.getX() - point2.getX()), 2) + pow((point1.getY() - point2.getY()), 2));
}

bool opY(Point point1, Point point2) {
    return point1.getY() < point2.getY();
}


vector <Point> recur(int startSearch, int endSearch) {
    double currMinValue = 5.0;
    int p1;
    int p2;
    vector <Point> closePair;
    if (endSearch - startSearch == 1) {
        currMinValue = calcDistances(randomPoints[startSearch], randomPoints[endSearch]);
        closePair.push_back(randomPoints[startSearch]);
        closePair.push_back(randomPoints[endSearch]);
    } else if (endSearch - startSearch == 2) {
        currMinValue = calcDistances(randomPoints[startSearch], randomPoints[startSearch + 1]);
        p1 = startSearch;
        p2 = startSearch + 1;
        double d1 = calcDistances(randomPoints[startSearch], randomPoints[startSearch + 2]);
        double d2 = calcDistances(randomPoints[startSearch + 1], randomPoints[startSearch + 2]);
        if (d1 < currMinValue) {
            currMinValue = d1;
            p1 = startSearch;
            p2 = startSearch + 2;
        }
        if (d2 < currMinValue) {
            currMinValue = d2;
            p1 = startSearch + 1;
            p2 = startSearch + 2;
        }
        closePair.push_back(randomPoints[p1]);
        closePair.push_back(randomPoints[p2]);
    } else {
        int split = startSearch + (endSearch - startSearch) / 2.0;
        vector <Point> dist1 = recur(startSearch, split);
        vector <Point> dist2 = recur(1 + split, endSearch);
        double temp1 = calcDistances(dist1[0], dist1[1]);
        double temp2 = calcDistances(dist2[0], dist2[1]);
        Point close1;
        Point close2;
        if (temp1 < temp2) {
            currMinValue = temp1;
            close1 = dist1[0];
            close2 = dist1[1];
        } else {
            currMinValue = temp2;
            close1 = dist2[0];
            close2 = dist2[1];
        }
        double xmid = (randomPoints[split].getX() + randomPoints[split + 1].getX()) / 2.0;
        double tempDist = 0;
        bool wasInStrip = false;
        for (int myp1 = 0; myp1 < N; myp1++)
            for (int myp2 = 0; myp2 < N; myp2++)
                if (randomPoints[myp1].getX() < xmid and randomPoints[myp2].getX() > xmid and
                    abs(randomPoints[myp1].getX() - xmid) < currMinValue and
                    abs(randomPoints[myp2].getX() - xmid) < currMinValue) {
                    tempDist = calcDistances(randomPoints[myp1], randomPoints[myp2]);
                    if (tempDist < currMinValue) {
                        wasInStrip = true;
                        currMinValue = tempDist;
                        p1 = myp1;
                        p2 = myp2;
                    }
                }
        if (wasInStrip) {
            closePair.push_back(randomPoints[p1]);
            closePair.push_back(randomPoints[p2]);
        } else {
            closePair.push_back(close1);
            closePair.push_back(close2);
        }
    }
    return closePair;
}


vector <Point> recurFast(int startSearch, int endSearch) {
    double currMinValue = 5.0;
    int p1;
    int p2;
    vector <Point> closePair;
    if (endSearch - startSearch == 1) {
        currMinValue = calcDistances(randomPoints[startSearch], randomPoints[endSearch]);
        closePair.push_back(randomPoints[startSearch]);
        closePair.push_back(randomPoints[endSearch]);
        return closePair;
    } else if (endSearch - startSearch == 2) {
        currMinValue = calcDistances(randomPoints[startSearch], randomPoints[startSearch + 1]);
        p1 = startSearch;
        p2 = startSearch + 1;
        double d1 = calcDistances(randomPoints[startSearch], randomPoints[startSearch + 2]);
        double d2 = calcDistances(randomPoints[startSearch + 1], randomPoints[startSearch + 2]);
        if (d1 < currMinValue) {
            currMinValue = d1;
            p1 = startSearch;
            p2 = startSearch + 2;
        }
        if (d2 < currMinValue) {
            currMinValue = d2;
            p1 = startSearch + 1;
            p2 = startSearch + 2;
        }
        closePair.push_back(randomPoints[p1]);
        closePair.push_back(randomPoints[p2]);
        return closePair;
    } else {
        int split = startSearch + (endSearch - startSearch) / 2.0;
        vector <Point> dist1 = recurFast(startSearch, split);
        vector <Point> dist2 = recurFast(1 + split, endSearch);
        double temp1 = calcDistances(dist1[0], dist1[1]);
        double temp2 = calcDistances(dist2[0], dist2[1]);
        Point close1;
        Point close2;
        if (temp1 < temp2) {
            currMinValue = temp1;
            close1 = dist1[0];
            close2 = dist1[1];
        } else {
            currMinValue = temp2;
            close1 = dist2[0];
            close2 = dist2[1];
        }
        double xmid = (randomPoints[split].getX() + randomPoints[split + 1].getX()) / 2.0;
        vector <Point> inStripY;
        bool wasInStrip = false;
        double tempDist = 0;
        for (int myPoint = 0; myPoint < N; myPoint++)
            if (abs(randomPointsY[myPoint].getX() - xmid) < currMinValue)
                inStripY.push_back(randomPointsY[myPoint]);
        for (int myp1 = 0; myp1 < inStripY.size(); myp1++)
            for (int myplus = 1; myplus <= 5; myplus++)
                if ((myp1 + myplus) < inStripY.size()) {
                    tempDist = calcDistances(inStripY[myp1], inStripY[myp1 + myplus]);
                    if (tempDist < currMinValue) {
                        wasInStrip = true;
                        currMinValue = tempDist;
                        p1 = myp1;
                        p2 = myp1 + myplus;
                    }
                }
        if (wasInStrip) {
            closePair.push_back(inStripY[p1]);
            closePair.push_back(inStripY[p2]);
        } else {
            closePair.push_back(close1);
            closePair.push_back(close2);
        }
        return closePair;
    }
}


int main() {
    srand(time(NULL)); /* set random seed once using clock */

    for (int i = 0; i < N; i++) {
        randomPoints.push_back(Point());
        randomPointsY.push_back(Point(randomPoints[i].getX(), randomPoints[i].getY()));
    }

    //brute force
    double bruteDistance;
    double currMin = 2.0;
    int currMinPoint1 = 0;
    int currMinPoint2 = 0;
    for (int r = 0; r < N; r++)
        for (int c = 0; c < N; c++) {
            if (r < c) {
                bruteDistance = calcDistances(randomPoints[r], randomPoints[c]);
                if (bruteDistance < currMin) {
                    currMin = bruteDistance;
                    currMinPoint1 = r;
                    currMinPoint2 = c;
                }
            }
        }

    printf("Brute Force O(n2): (%f, %f) (%f, %f) \n", randomPoints[currMinPoint1].getX(),
           randomPoints[currMinPoint1].getY(), randomPoints[currMinPoint2].getX(), randomPoints[currMinPoint2].getY());

    //recursive n2
    sort(randomPoints.begin(), randomPoints.end());
    vector <Point> f = recur(0, N - 1);
    printf("Recursive Solution O(n2): (%f, %f) (%f, %f) \n", f[0].getX(), f[0].getY(), f[1].getX(), f[1].getY());

    //recursive nlogn
    sort(randomPointsY.begin(), randomPointsY.end(), opY);
    vector <Point> g = recurFast(0, N - 1);
    printf("Recursive Solution O(nlogn): (%f, %f) (%f, %f) \n", g[0].getX(), g[0].getY(), g[1].getX(), g[1].getY());

    //printf("Minimum Distance: %f \n", currMin);
}
