#include "opencv2/imgproc/imgproc.hpp"
#include "opencv2/opencv.hpp"
#include "opencv2/core/core.hpp"
#include <vector>
#include <limits>

#define PI 3.14159265358979323846
#define PHI (1.0 + sqrt(5.0)) / 2.0
#define INVERSE_PHI 2.0/ (1.0 + sqrt(5.0))
#define CUBE 0
#define TETRAHEDRON 1
#define OCTAHEDRON 2
#define DODECAHEDRON 3
#define ICOSAHEDRON 4

using namespace std;
using namespace cv;

const vector<Point3f> centers{Point3f(90, 10, 120), Point3f(100, 5, 120), Point3f(110, 10, 120), Point3f(90, 0, 120), Point3f(110, 0, 120)};
const Mat camera = (Mat_<double>(4, 1) << 100.0, 4.0, -25.0, 1.0);
const int frameCount = 200;
const int fps = 60;

vector<pair<int, int>> getEdges(vector<Mat> points) {
    map<pair<int, int>, double> distances;
    int len = points.size();
    double minDist = numeric_limits<double>::max();
    for (int x = 0; x < len; x++)
        for (int y = x + 1; y < len; y++) {
            pair<int, int> xy(x, y);
            distances[xy] = ((points[x] - points[y]).dot((points[x] - points[y]))); //subtract to find difference between each coordinate, dot to square
            if (distances[xy] != 0)
                minDist = min(minDist, distances[xy]);
        }
    vector<pair<int, int>> edges;
    for (int x = 0; x < len; x++)
        for (int y = x + 1; y < len; y++) {
            pair<int, int> xy(x, y);
            if(abs(minDist - distances[xy]) < 0.0001)
                edges.emplace_back(x, y);
        }
    return edges;
}

vector<Mat> projectPoints(const vector<Mat>& points, Mat camera) {
    vector<Mat> projectedPoints;
    for (Mat point : points) {
        double projectionRatio = (-camera.at<double>(2, 0)) / (point.at<double>(2, 0) - camera.at<double>(2, 0));
        Mat projectedPoint = (point - camera) * projectionRatio;
        projectedPoint.at<double>(2, 0) = 0.0;
        projectedPoint.at<double>(3, 0) = 1.0;
        projectedPoints.emplace_back(projectedPoint);
    }
    return projectedPoints;
}

vector<Mat> drawCube(const Point3f& center, int sideLength) {
    vector<Mat> cubePoints;
    double rLength = sideLength / 2.0;
    for (double x : {-rLength, rLength})
        for (double y : {-rLength, rLength})
            for (double z : {-rLength, rLength})
                cubePoints.push_back((Mat_<double>(4, 1) << center.x + x, center.y + y, center.z + z, 1));
    return cubePoints;
}

vector<Mat> drawTetrahedron(const Point3f& center, int sideLength) {
    vector<Mat> tetrahedronPoints{ 
            (Mat_<double>(4, 1) << center.x + (sideLength / sqrt(3)), center.y, center.z, 1),
            (Mat_<double>(4, 1) << center.x, center.y + (2.0 * sideLength / sqrt(6)), center.z, 1),
            (Mat_<double>(4, 1) << center.x  - (sqrt(3) * sideLength/ 6.0), center.y, center.z + (sideLength / 2.0), 1),
            (Mat_<double>(4, 1) << center.x - (sqrt(3) *sideLength/ 6.0), center.y, center.z - (sideLength / 2.0), 1)
    };
    return tetrahedronPoints;
}

vector<Mat> drawOctahedron(const Point3f& center, int sideLength) {
    vector<Mat> octahedronPoints;
    double rLength = sideLength / 2.0;
    for (double x : {-rLength, 0.0, rLength})
        for (double y : {-rLength, 0.0, rLength})
            for (double z : {-rLength, 0.0, rLength})
                if ((!x + !y + !z) == 2) //two zeroes and one +- sidelength
                    octahedronPoints.push_back((Mat_<double>(4, 1) << center.x + x, center.y + y, center.z + z, 1));                 
    return octahedronPoints;
}

vector<Mat> drawIcosahedron(const Point3f& center, int sideLength) {
    vector<Mat> icosahedronPoints;
    double rLength = sideLength / 2.0;
    for (Point3f point : { Point3f(0.0, PHI, 1.0), Point3f(1.0, 0.0, PHI), Point3f(PHI, 1.0, 0.0), Point3f(0.0, -PHI, 1.0), Point3f(1.0, 0.0, -PHI), Point3f(-PHI, 1.0, 0.0) }){   
        icosahedronPoints.push_back((Mat_<double>(4, 1) << center.x + (rLength*point.x), center.y + (rLength * point.y), center.z + (rLength * point.z), 1));
        icosahedronPoints.push_back((Mat_<double>(4, 1) << center.x - (rLength * point.x), center.y - (rLength * point.y), center.z - (rLength * point.z), 1));
    }          
    return icosahedronPoints;
}

vector<Mat> drawDodecahedron(const Point3f& center, int sideLength) {
    vector<Mat> dodecahedronPoints;
    double rLength = sideLength / 2.0;
    for (double x : {-rLength, rLength})
        for (double y : {-rLength, rLength})
            for (double z : {-rLength, rLength}) 
                dodecahedronPoints.push_back((Mat_<double>(4, 1) << center.x + x, center.y + y, center.z + z, 1));
    for (Point3f point : { Point3f(0.0, INVERSE_PHI, PHI), Point3f(PHI, 0.0, INVERSE_PHI), Point3f(INVERSE_PHI, PHI, 0.0), Point3f(0.0, -INVERSE_PHI, PHI), Point3f(PHI, 0.0, -INVERSE_PHI), Point3f(-INVERSE_PHI, PHI, 0.0) }){
        dodecahedronPoints.push_back((Mat_<double>(4, 1) << center.x + (rLength * point.x), center.y + (rLength * point.y), center.z + (rLength * point.z), 1));
        dodecahedronPoints.push_back((Mat_<double>(4, 1) << center.x - (rLength * point.x), center.y - (rLength * point.y), center.z - (rLength * point.z), 1));
    }
    return dodecahedronPoints;
}

void drawShape(Mat& frame, int shapeType, int frameNum, int sideLength) {
    vector<Mat> shapePoints;
    Point3f center = centers[shapeType];
    switch (shapeType)
    {
    case CUBE:
        shapePoints = drawCube(center, sideLength);
        break;
    case TETRAHEDRON:
        shapePoints = drawTetrahedron(center, sideLength);
        break;
    case OCTAHEDRON:
        shapePoints = drawOctahedron(center, sideLength);
        break;
    case DODECAHEDRON:
        shapePoints = drawDodecahedron(center, sideLength);
        break;
    case ICOSAHEDRON:
        shapePoints = drawIcosahedron(center, sideLength);
        break;
    default:
        cerr << "Invalid shape, please try again" << endl;
        return;
    }
    Mat shift = (Mat_<double>(4, 4) <<
        1, 0, 0, -center.x,
        0, 1, 0, -center.y,
        0, 0, 1, -center.z,
        0, 0, 0, 1);
    Mat newShift = shift.inv();
    double rotationAngle = 2.0 * PI * frameNum / frameCount;
    double s = sin(rotationAngle), c = cos(rotationAngle);
    Mat rotate = (Mat_<double>(4, 4) <<
        c, 0, s, 0,
        0, 1, 0, 0,
        -s, 0, c, 0,
        0, 0, 0, 1);
    Mat scalar = newShift * rotate * shift;
    for (auto& shapePoint : shapePoints)
        shapePoint = scalar * shapePoint;  // Order matters
    vector<Mat> projectedPoints = projectPoints(shapePoints, camera);
    int cols = frame.cols, halfRows = frame.rows / 2, halfCols = cols / 2;
    for (pair<int, int> edgePoint : getEdges(shapePoints)) {
        Mat start = projectedPoints[edgePoint.first] * 100, end = projectedPoints[edgePoint.second] * 100;
        int x1 = int(start.at<double>(0, 0)) + halfRows, y1 = cols - int(start.at<double>(1, 0) + halfCols);
        int x2 = int(end.at<double>(0, 0)) + halfRows, y2 = cols - int(end.at<double>(1, 0) + halfCols);
        line(frame, Point(x1, y1), Point(x2, y2), Scalar(255, 255, 255), 3, LINE_AA);
    }
}

int main()
{
    Mat baseFrame(Size(1000, 1000), 16);
    VideoWriter outputVideo("FernandesProjectivePlatonicSolids.mov", VideoWriter::fourcc('m', 'p', '4', 'v'), fps, baseFrame.size(), true);    
    for (int frameNum = 0; frameNum < frameCount; frameNum++) {
        Mat frame = baseFrame.clone();
        drawShape(frame, CUBE, frameNum, 4);
        drawShape(frame, TETRAHEDRON, frameNum, 4);
        drawShape(frame, OCTAHEDRON, frameNum, 4);
        drawShape(frame, ICOSAHEDRON, frameNum, 4);
        drawShape(frame, DODECAHEDRON, frameNum, 4);
        outputVideo << frame;
        imshow("Rotating Solids", frame);
        waitKey(1);
    }
    outputVideo.release();
    destroyAllWindows();
    return 0;
}