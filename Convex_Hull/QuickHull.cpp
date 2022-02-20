//Lab for Computer Vision 
//Date: 12/2/2021
//Name: Tommy Williams 
//Task: Using the quick hull algorithm 
//create a convex hull for an image with n points

#include <iostream>
#include <fstream>
#include <math.h>
#include <iomanip>
#include <algorithm>
#include <vector>
#include <list>
#include <string>
#include <chrono>

using namespace std;

const int xlen = 400;
const int ylen = 400;
int pixels [ylen][xlen][3] = {};

int scaleUp(double n) {
    int newN = round(n*xlen);
    return newN;
}

void illuminate(int x, int y,string color) {
    if (x>=xlen || y>=xlen || y<=0 || x<=0) {
        return;
    }
    if (color == "r") {
        pixels[y][x][0]=100;
        pixels[y][x][1]=0;
        pixels[y][x][2]=0;
    }
    else if (color == "b") {
        pixels[y][x][0]=0;
        pixels[y][x][1]=0;
        pixels[y][x][2]=100;
    }
    else if (color == "g") {
        pixels[y][x][0]=0;
        pixels[y][x][1]=100;
        pixels[y][x][2]=0;
    }
    else {
        pixels[y][x][0]=120;
        pixels[y][x][1]=120;
        pixels[y][x][2]=120;
    }
}

void circle(double radius, double centerX, double centerY,string color) {
    int rad = scaleUp(radius);
    int cx = scaleUp(centerX);
    int cy = scaleUp(centerY);
    int x, y, xmax, y2, y2_new, ty;
    
    xmax = (int) (rad * 0.70710678); //45 degrees
    
    y=rad; //y is the radius 
    y2=y*y; //r squared 
    ty = (2*y)-1; 
    y2_new = y2;
    
    for (x=0; x<=xmax; x++) {
        if ((y2-y2_new) >= ty) { //if the radius is 
            y2 -= ty;
            y-=1;
            ty-=2;
        }
        illuminate(cx+x,cy+y,color);
        illuminate(cx+x,cy-y,color);
        illuminate(cx-x,cy+y,color);
        illuminate(cx-x,cy-y,color);
        illuminate(cx+y,cy+x,color);
        illuminate(cx+y,cy-x,color);
        illuminate(cx-y,cy+x,color);
        illuminate(cx-y,cy-x,color);
        
        y2_new -= (2*x);
    }
    
}

void Bresenham(double ex1,double ey1,double ex2,double ey2,string color) {
    //preconditions: x2>x1 and y2>y1 and driving axis is x axis 
    //one pixel on each value on driving axis from start to finish
    //driving axis is farther distance
    
    int x1 = scaleUp(ex1);
    int y1 = scaleUp(ey1);
    int x2 = scaleUp(ex2);
    int y2 = scaleUp(ey2);
    
    int dx = abs(x2-x1);
    int dy = abs(y2-y1);
    
    int xi = min(x1,x2);
    int xf = max(x1,x2);
    int yi = min(y1,y2);
    int yf = max(y1,y2);
    
    if (dy==0) {
        for (int i=xi; i<=xf; i++) {
            illuminate(i,yi,color);
        }
        return;
    }
    else if (dx==0) {
        for (int j=yi; j<=yf; j++) {
            illuminate(xi,j,color);
        }
        return;
    }
    
    int j = y1;
    int e;
    
    if (dy<=dx) { 
        e=dy-dx;
        for (int i=x1; (x2==xf) ? i<x2 : i>x2; (x2==xf) ? i++ : i--) {
            illuminate(i,j,color);
            if (e>=0) {
                y1==yi ? j++ : j--;
                e-=dx;
            }            
            e+=dy;
        }
    }
    else {
        e=dx-dy;  
        j=x1;
        for (int i=y1; (y2==yf) ? i<y2 : i>y2; (y2==yf) ? i++ : i--) {
            illuminate(j,i,color);
            //////////cout<< j << " " << i << "\n";
            if (e>=0) {
                x1==xi ? j++ : j--;
                e-=dy;
                
            }
            e+=dx;
        }
    }  
}

class point {
    private:
        double x;
        double y;
    public: 
        point(double xe, double ye) {
            x = xe;
            y = ye; 
        } 
        point() {
            x = 0;
            y = 0;
        }
        double getX() {
            return x;
        }
        double getY() {
            return y;
        }
        void set(double xe, double ye) {
            x = xe;
            y = ye;
        }
        double distanceTo(point o) {
            return sqrt(pow(x-o.getX(),2)+pow(y-o.getY(),2));
        }
        string toString() {
            ostringstream ste;
            ste << setprecision(17) << "(" << x << "," << y << ")";
            return(ste.str());
        }
        
        bool operator<(point & other) const {
            return x < other.getX();
        }
};



class Line {
    private:
        point refPoint;
        point refPoint2;
        bool undefined;
        bool isLine;
        double slope;
        string pointSlopeForm;
        string standardForm; //Ax+By=C
        double A=0;
        double B=0;
        double C=0;
    public:
        Line() {
            refPoint = point(0,0);
            refPoint2 = point(1,1);
        }
        Line(point start, point end) {
            refPoint = start;
            refPoint2 = end;
            slope = (start.getY()-end.getY())/(start.getX()-end.getX());
            if (isinf(slope)) { undefined = true; }
            else { undefined = false; }
            if (isnan(slope)) { isLine = false; }
            else { isLine = true; }
            
            //Setting up forms for equations//
            ostringstream stream;
            stream << "y-" << start.getY() << "=" << slope << "(x-" << start.getX() << ")";
            pointSlopeForm = stream.str();
            
            stream.str("");
            stream.clear();
            if (undefined) {
                A = 1;
                B = 0;
                C = start.getX();
            }
            else if (not isLine) {
                A = 0;
                B = 0; 
                C = 0;
            }
            else {
                A = slope;
                B = 1;
                C = (-1*start.getX()*slope)-(-1*start.getY());
            }
            stream << A << "x + " << B << "y = " << C;
            standardForm = stream.str();
        }
        double getSlope() {
            return slope;
        }
        point getRefPoint() {
            return refPoint;
        }
        point getRefPoint2() {
            return refPoint2;
        }
        bool isItLine() {
            return isLine;
        }
        bool undefinedSlope() {
            return undefined; 
        }
        vector<double> getStandardEquation() {
            vector<double> ABC;
            ABC.push_back(A);
            ABC.push_back(B);
            ABC.push_back(C);
            return ABC;
        }
        double solveAt(double x) {
            //Ax+By=C -> y = (C-Ax)/B
            return (C+A*x)/B;
        }
        double solveInverseAt(double y) {
            //y = (C-Ax)/B -> (BY-C)/-A=x
            return (B*y-C)/(A);
        }
        point intersection(Line other) {
            if (other.getSlope()==slope) {
                if (other.solveAt(refPoint.getX())==refPoint.getY()) {
                    //should never reach while completing this task 
                    return point(200,200); //to show that there are infinite intersections 
                }
                else {
                    //should also never reach while completing this task 
                    return point(400,400); //400 to show that there are no intersections 
                }
            }
            else {
                //(C1-A1x)/B1=(C2-A2x)/B2
                //(C1-A1x)B2=B1(C2-A2x)
                //(C1)B2-A1xB2=B1(C2)-B1A2x
                //B1(C2)-(C1)B2=B1(A2x)-A1x(B2)
                //B1(C2)-(C1)B2=(B1A2-B2A1)x
                //(B1(C2)-(C1)B2)/(B1A2-B2A1) = x 
                //y = solveAt(x)
                vector<double> sEOther = other.getStandardEquation();
                double A1 = sEOther.at(0);
                double B1 = sEOther.at(1);
                double C1 = sEOther.at(2);
                double xIntersect;
                double yUs;
                if (other.undefinedSlope()) {
                    xIntersect = C1/A1;
                    yUs = solveAt(xIntersect);
                }
               else if (!other.isItLine()) {
                    yUs = C1/B1;
                    xIntersect = solveInverseAt(yUs);
                }
                else if (undefined) {
                    xIntersect = C/A;
                    yUs = other.solveAt(xIntersect);
                }
                else if (!isLine) {
                    yUs = C/B;
                    xIntersect = solveInverseAt(yUs);
                }
                else {
                    xIntersect = (B*C1-C*B1)/(B1*A-B*A1);
                    yUs = other.solveAt(xIntersect);
                }
                
                return point(xIntersect,yUs);
            }
        }
        Line perp(point refPoint) {
            double newSlope;
            
            if (undefined) {
                point point2 = point(refPoint.getX()+0.01,refPoint.getY());
                return Line(refPoint,point2);
            }
            else if (!isLine) {
                point point2 = point(refPoint.getX(),refPoint.getY()+0.01);
                return Line(refPoint,point2);
            }
            else {
                newSlope = -1/slope;
                point point2 = point(refPoint.getX()+0.1,(refPoint.getY()+0.1*newSlope));
                return Line(refPoint,point2);
            }
        }
        vector<point> posPtAtDistFrom(point ref,double d) {
            vector<point> ePoints;
            double x = ref.getX();
            double y = ref.getY();
            point e1;
            point e2;
            if (undefined) {
                e1.set(x,y+d);
                e2.set(x,y-d);
            }
            else if (!isLine) {
                e1.set(x+d,y);
                e2.set(x-d,y);
            }
            else {
                e1.set(x+d*cos(atan(slope)),y+d*sin(atan(slope)));
         e2.set(x-d*cos(atan(slope)),y-d*sin(atan(slope)));
            }
            ePoints.push_back(e1);
            ePoints.push_back(e2);
            return ePoints;
        }
       string toString() {
            return pointSlopeForm;
       }
       pair<Line, Line> makeTriangle(point p) {
            return pair<Line,Line>(Line(p,refPoint),Line(p,refPoint2));
       }
};

bool xSort(point a, point b) {
    return (a.getX() < b.getX());
}
bool ySort(point a, point b) {
    return (a.getY() < b.getY());
}

double getArea(point p1, point p2, point p3) {
    double p1x = p1.getX(); double p1y = p1.getY();
    double p2x = p2.getX(); double p2y = p2.getY();
    double p3x = p3.getX(); double p3y = p3.getY();

    double a = sqrt(pow(p1x-p2x,2)+pow(p1y-p2y,2));
    double b = sqrt(pow(p2x-p3x,2)+pow(p2y-p3y,2));
    double c = sqrt(pow(p3x-p1x,2)+pow(p3y-p1y,2));
    double s = (a+b+c)/2;
    double area = sqrt(s*(s-a)*(s-b)*(s-c));
    return area;
}

double checkDistance(vector<point> two) {
    if (two.empty()) {
        return 1;
    }
    return two.at(0).distanceTo(two.at(1));
}

bool ttrol(point linePt1, point linePt2, point p) {
    //TO THE RIGHT OF LINE 
    if (linePt1.getY() > linePt2.getY()) { //if first point is higher than second point 
        return (p.getX() < Line(linePt1,linePt2).solveInverseAt(p.getY()));
        //return true if the x value of the point is lower than the x value of the line corresponding to the y value
        //this means that is 
    }
    else if (linePt1.getY() < linePt2.getY()) {
        return (p.getX() > Line(linePt1,linePt2).solveInverseAt(p.getY()));
    }
    else if (linePt1.getX() > linePt2.getX()) {
        return (p.getY() > linePt2.getY());
    }
    else {
        return (p.getY() < linePt1.getY());
    }
}


vector<point> ConvexHull;

void findHull(vector<point> &hull, Line side) {
    point p1 = side.getRefPoint(); 
    point p2 = side.getRefPoint2(); 
    point maxPt;
    double bigArea = 0;
    if (hull.size()==0) {
        Bresenham(p1.getX(), p1.getY(), p2.getX(), p2.getY(),"b");
        return;
    }
    
    double currArea = 0;
    int maxPtIndex;
    for (int i=0; i<hull.size(); i++) {
        currArea = getArea(p1,p2,hull.at(i));
        if (currArea>bigArea) {
            bigArea = currArea;
            maxPt = hull.at(i);
            maxPtIndex = i;
        }
    }
    ConvexHull.push_back(maxPt);
    Line p1Max = Line(p1,maxPt);
    Line p2Max = Line(p2,maxPt);
    
    double totalArea = 0;
    double area1; double area2; double area3;
    vector<point> s1;
    vector<point> s2;
    point p4;
    double eps = 0.000000000001;
    for (int i=0; i<hull.size(); i++) {
        if (i!=maxPtIndex) {
            p4 = hull.at(i);
            totalArea = getArea(p1,p2,p4) + getArea(p1,maxPt,p4) + getArea(p2,maxPt,p4);
            if (abs(totalArea-bigArea)<eps) {
                ;
            }
            else if (ttrol(p2,maxPt,p4)) {
                s2.push_back(p4);
            }
            else if (ttrol(maxPt,p1,p4)) {
                s1.push_back(p4);
            }
        }
    }
    findHull(s2,Line(maxPt,p2));
    findHull(s1,Line(p1,maxPt));
}

vector<point> part1() {
    vector<point> ptsVector;
    srand(time(NULL));
    for (int i=0; i<60; i++) {
        ptsVector.push_back(point((double) rand() / RAND_MAX,((double) rand() / RAND_MAX)));
    }
    
    fstream points;
    points.open("points.txt",ios::out);
    
    for (vector<point>::iterator pl = ptsVector.begin(); pl != ptsVector.end(); ++pl) {
        double lx = (*pl).getX();
        double ly = (*pl).getY();
        points << setprecision(17) << lx << "  " << ly << "\n";
        circle(2.0/400,lx,ly,"white");
    }
    
    points.close();
    
    //QUICK HULL//
    sort(ptsVector.begin(),ptsVector.end(),ySort);
    Line divide = Line(ptsVector.at(0),ptsVector.at(ptsVector.size()-1));
    ConvexHull.push_back(ptsVector.at(0));
    ConvexHull.push_back(ptsVector.at(ptsVector.size()-1));
    vector<point> s1;
    vector<point> s2;
    for (int i=0; i<ptsVector.size(); i++) {
        point pi = ptsVector.at(i);
        if (ttrol(ptsVector.at(0),ptsVector.at(ptsVector.size()-1),pi)) {
            s1.push_back(pi);
        }
        else {
            s2.push_back(pi);
        }
    }
    
    findHull(s1,Line(ptsVector.at(ptsVector.size()-1),ptsVector.at(0)));
    findHull(s2,divide);
    for (int i=0; i<ConvexHull.size(); i++) {
        circle(2.0/400,ConvexHull.at(i).getX(),ConvexHull.at(i).getY(),"g");
    }
    return ConvexHull;

    
}

int main() {
    
    part1();
    fstream imgFile;
    imgFile.open("quickhull.ppm",ios::out);
    imgFile << "P3 " << xlen << " " << ylen << " 120 \n";
    for (int ay=0; ay<ylen; ay++) {
        for (int ax=0; ax<xlen; ax++) {
            imgFile << pixels[ay][ax][0] << " " << pixels[ay][ax][1] << " " << pixels[ay][ax][2] << " ";
        }
        imgFile << "\n";
    }
    imgFile.close();
    return 0;
}
