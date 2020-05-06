/* Rachel Fernandes
 * Dr. Zacharias pd 4
 * Closest Pair Part 5
 * I will uphold academic and personal integrity in the TJ community.
 */
#include <stdio.h>      /* printf, scanf, puts, NULL */
#include <stdlib.h>     /* srand, rand, RAND_MAX */
#include <time.h>       /* time */
#include <math.h>
#include <bits/stdc++.h>
#include <vector>
#include <map>

using namespace std;
#define N 10000

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

static vector <Point> randomPoints;
static vector <Point> randomPointsY;

double calcDistances(Point point1, Point point2) {
    return sqrt(pow((point1.x - point2.x), 2) + pow((point1.y - point2.y), 2));
}

bool opY(Point point1, Point point2) {
    return point1.y < point2.y;
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
        int split = startSearch + (endSearch - startSearch) / 2.0; // >> 1 instead of divide by 2.0??
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
        double xmid = (randomPoints[split].x + randomPoints[split + 1].x) / 2.0;
        vector <Point> inStripY;
        bool wasInStrip = false;
        double tempDist = 0;
        for (int myPoint = 0; myPoint < N; myPoint++)
            if (abs(randomPointsY[myPoint].x - xmid) < currMinValue)
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


map<pair <int, int>, vector <Point>> makeDictionary(int lastIndex, double newMinValue){
    static map <pair <int, int>, vector<Point>> ss;
    pair <int, int> tempSS;
    for(int e = 0; e <= lastIndex; e++)
    {
        int ssX = (int)(randomPoints[e].x/(newMinValue/2));
        int ssY = (int)(randomPoints[e].y/(newMinValue/2));
        tempSS = make_pair(ssX, ssY);
        if(ss.find(tempSS) != ss.end())
            ss[tempSS].push_back(randomPoints[e]);
        else {
            vector<Point> randothing{randomPoints[e]};
            ss.insert(make_pair(tempSS, randothing));
        }

    }
    return ss;
}

pair <Point, Point> bigOn() {
    //assume point list is already randomized because it is made of randomly generated points
    // if not randomized, implement Knuth scrambling algorithm
    Point p1 = randomPoints[0];
    Point p2 = randomPoints[1];
    double currMinValue = calcDistances(p1, p2); //delta
    static map <pair <int, int>, vector<Point>> subSquares = makeDictionary(1, currMinValue);
    //keys are a pair of indexes of sub-squares and values are points from the randomized list in that sub-square.
    pair <int, int> tempSubSquare;
    double tempDistance;
    Point tempPoint;
    int subSquareX;
    int subSquareY;
    for(int it = 2; it < N; it++)
    {
        subSquareX = (int)(randomPoints[it].x/(currMinValue/2));
        subSquareY = (int)(randomPoints[it].y/(currMinValue/2));
        bool foundCloserPair = false;
        for(int itx = -2; itx <= 2; itx++)
            for(int ity = -2; ity <= 2; ity++)
                if(!(itx == 0 and ity == 0) and itx + subSquareX >= 0 and ity + subSquareY >= 0)
                { //check and make sure you're not comparing the point to itself and that you haven't found a closer pair yet
                    tempSubSquare = make_pair(subSquareX + itx, subSquareY + ity);
                    if(subSquares.find(tempSubSquare) != subSquares.end())
                    { //check that the subsquare is even in the dictionary!!
                        for(int thing = 0; thing < subSquares[tempSubSquare].size(); thing++)
                        {
                            tempPoint = subSquares[tempSubSquare][thing];
                            tempDistance = calcDistances(randomPoints[it], tempPoint);
                            if (tempDistance < currMinValue)
                            {
                                foundCloserPair = true;
                                currMinValue = tempDistance;
                                p1 = randomPoints[it];
                                p2 = tempPoint;
                            }
                        }
                    }
                }
        if(foundCloserPair)
            subSquares = makeDictionary(it, currMinValue);
        else
        {
            tempSubSquare = make_pair(subSquareX, subSquareY);
            if(subSquares.find(tempSubSquare) != subSquares.end())
                subSquares[tempSubSquare].push_back(randomPoints[it]);
            else {
                vector<Point> randothing{randomPoints[it]};
                subSquares.insert(make_pair(tempSubSquare, randothing));
            }
        }
    }
    return make_pair(p1, p2);
}


int main() {
    srand(time(NULL)); /* set random seed once using clock */
    time_t start_time;
    time_t current_time;
    double random_time = 0.0;
    double recur_time = 0.0;

    for (int i = 0; i < N; i++) {
        randomPoints.push_back(Point());
        randomPointsY.push_back(Point(randomPoints[i].x, randomPoints[i].y));
    }

    /*for (int j = 0; j < N; j++)
        printf("( %f, %f ) ", randomPoints[j].x, randomPoints[j].y);
    */
    start_time = time(NULL);
    pair <Point, Point> g = bigOn();
    current_time = time(NULL);
    random_time = random_time + current_time - start_time;
    //printf("\nRandom: (%f, %f) (%f, %f) : %f \n ", g.first.x, g.first.y, g.second.x, g.second.y, calcDistances(g.first, g.second));
    printf("Randomized time: %f \n", random_time);

    start_time = time(NULL);
    sort(randomPoints.begin(), randomPoints.end());
    sort(randomPointsY.begin(), randomPointsY.end(), opY);
    vector <Point> f = recurFast(0, N - 1);
    current_time = time(NULL);
    recur_time = recur_time + current_time - start_time;
    printf("Recursize time: %f \n", recur_time);
    printf("(%f, %f) (%f, %f) : %f \n ", f[0].x, f[0].y, f[1].x, f[1].y, calcDistances(f[0], f[1]));

}