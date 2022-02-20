//Name: Tommy Williams 
//Lab 3 Closest points 
//Part 1: Brute force

#include <iostream>
#include <fstream>
#include <math.h>
#include <iomanip>
#include <algorithm>
#include <vector>
#include <list>

using namespace std;

const int xlen = 800;
const int ylen = 800;
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
        pixels[y][x][0]=0;
        pixels[y][x][1]=0;
        pixels[y][x][2]=0;
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

void Bresenham(double ex1,double ey1,double ex2,double ey2) {
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
            illuminate(i,yi,"green");
        }
        return;
    }
    else if (dx==0) {
        for (int j=yi; j<=yf; j++) {
            illuminate(xi,j,"green");
        }
        return;
    }
    
    int j = y1;
    int e;
    
    if (dy<=dx) { 
        e=dy-dx;
        for (int i=x1; (x2==xf) ? i<x2 : i>x2; (x2==xf) ? i++ : i--) {
            illuminate(i,j,"green");
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
            illuminate(j,i,"green");
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
};

void part1() {
    
    for (int ye=0; ye<ylen; ye++) {
        for (int xe=0; xe<xlen; xe++) {
            pixels[ye][xe][0]=120;
            pixels[ye][xe][1]=120;
            pixels[ye][xe][2]=120;
        }
    }
    
    srand(time(NULL));
    list<point> pts;
    for (int i=0; i<60; i++) {
        pts.push_back(point((double) rand() / RAND_MAX,((double) rand() / RAND_MAX)));
    }
    double minDistance = 2.0;
    point minPoint1;
    point minPoint2;
    int count=0;
    
    fstream points;
    points.open("points.txt",ios::out);
    
    for (list<point>::iterator pl = pts.begin(); pl != pts.end(); ++pl) {
        double lx = (*pl).getX();
        double ly = (*pl).getY();
        points << setprecision(17) << lx << "  " << ly << "\n";
        circle(3.0/800,lx,ly,"black");
        list<point>::iterator ps = pl;
        for (ps = ++ps; ps != pts.end(); ++ps) {
            double d = (*pl).distanceTo(*ps);
            count++;
            ////cout << (*pl).getX() << " " << (*ps).getX() << "\n";
            if (d<minDistance) {
                minDistance = d;
                minPoint1 = *pl;
                minPoint2 = *ps;
            }
        }
    }
    circle(2.0/800,minPoint1.getX(),minPoint1.getY(),"r");
    circle(2.0/800,minPoint2.getX(),minPoint2.getY(),"r");
    circle(3.0/800,minPoint1.getX(),minPoint1.getY(),"r");
    circle(3.0/800,minPoint2.getX(),minPoint2.getY(),"r");
    //cout << minDistance << "\n";
    //cout << minPoint1.getX() << " " << minPoint1.getY() << "\n";
    //cout << minPoint2.getX() << " " << minPoint2.getY() << "\n";
    points.close();
    
    fstream pointsppm;
    pointsppm.open("points.ppm",ios::out);
    
    pointsppm << "P3 " << xlen << " " << ylen << " 120 \n";
    for (int ay=0; ay<ylen; ay++) {
        for (int ax=0; ax<xlen; ax++) {
            pointsppm << pixels[ay][ax][0] << " " << pixels[ay][ax][1] << " " << pixels[ay][ax][2] << " ";
        }
        pointsppm << "\n";
    }
    pointsppm.close();
    
    
}

int main() {
    part1();
    
    return 0;
}
