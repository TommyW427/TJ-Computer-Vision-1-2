
#include <iostream>
#include <fstream>
#include <math.h>
#include <iomanip>
#include <algorithm>
#include <vector>
using namespace std;

/*
TASK: Find the smallest square such that 
the extended edges of the square intersect each of the four 
given points. See sample output. 
*/

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
            ////////cout<< j << " " << i << "\n";
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
class Square {
    private:
        double area;
        Line l1;
        Line l2;
        Line l3;
        Line l4;
        point sect1;
        point sect2;
        point sect3; 
        point sect4;
    public: 
        Square() {
            area = 0;
        }
        Square(Line line1, Line line2, Line line3, Line line4) { 
            sect1 = line1.intersection(line2);
            sect2 = line1.intersection(line3);
            
            if (sect1.getX()==400) { //lines do not intersect (making sure lines are ordered correctly)
                sect1 = line1.intersection(line4);
                sect3 = line4.intersection(line2);
                sect4 = line3.intersection(line2);
            }
            else if (sect2.getX()==400) { //lines do not intersect (making sure lines are ordered correctly)
                sect2 = line1.intersection(line4);
                sect3 = line4.intersection(line3);
                sect4 = line2.intersection(line3);
            }
            else {
                sect3 = line2.intersection(line4);
                sect4 = line3.intersection(line4);
            }
            
            double sideLength = sect1.distanceTo(sect2);
            area = sideLength*sideLength;
            l1 = line1;
            l2 = line2;
            l3 = line3;
            l4 = line4;
            
        }
        double getArea() {
            return area;
        }
        void drawSquare() {
            vector<Line> sidesOfSquare;
            sidesOfSquare.push_back(l1);
            sidesOfSquare.push_back(l2);
            sidesOfSquare.push_back(l3);
            sidesOfSquare.push_back(l4);
            for (int i=0; i<4; i++) {
                Line lime = sidesOfSquare.at(i);
                double slope = lime.getSlope();
                if (isinf(slope)) {
                    double xRef = lime.getRefPoint().getX();
                    Bresenham(xRef,0,xRef,1);
                }
                else if (isnan(slope)) {
                    double yRef = lime.getRefPoint().getY();
                    Bresenham(0,yRef,1,yRef);
                }
                if (abs(slope)>=1) {
                    double xmin = lime.solveInverseAt(0);
                    double xmax = lime.solveInverseAt(1);
                    Bresenham(xmin,0,xmax,1);
                }
                else {
                    double ymin = lime.solveAt(0);
                    double ymax = lime.solveAt(1);
                    Bresenham(0,ymin,1,ymax);
                }
            }
            circle(2.0/800,sect1.getX(),sect1.getY(),"black");
            circle(2.0/800,sect2.getX(),sect2.getY(),"black");
            circle(2.0/800,sect3.getX(),sect3.getY(),"black");
            circle(2.0/800,sect4.getX(),sect4.getY(),"black");
            
        }
        string toString() {
            ostringstream sqs;
            sqs << setprecision(17) << sect1.toString() << " , ";
            sqs << setprecision(17) << sect2.toString() << " , ";
            sqs << setprecision(17) << sect3.toString() << ", ";
            sqs << setprecision(17) << sect4.toString();
            sqs << setprecision(17) << " Area = " << area;
            return(sqs.str());
        }
        
        
        
};

Square process(point Aa, point Cc, point Bb, point Dd, bool eDown) {
    //letters are purely for ordering
    Line Line1 = Line(Aa,Cc);
    double d1 = Aa.distanceTo(Cc);
    Line Line2 = Line1.perp(Bb);
    vector<point> ePts = Line2.posPtAtDistFrom(Bb,d1);
    point Ee;
    if (eDown) {
       Ee = ePts[0]; 
    }
    else {
       Ee = ePts[1];
    }
    Line Line3 = Line(Dd,Ee);
    Line Line4 = Line3.perp(Aa);
    Line Line5 = Line3.perp(Cc);
    Line Line6 = Line4.perp(Bb);
    return Square(Line3,Line4,Line5,Line6);
}

void part2() {
    fstream txtFile;
    txtFile.open("output.txt",ios::out);
    
    ifstream four ("points.txt");
    string pts;
    getline(four,pts);
    four.close();
    
    txtFile << pts << "\n";
    
    pts.erase(remove(pts.begin(),pts.end(),'('),pts.end());
    pts.erase(remove(pts.begin(),pts.end(),')'),pts.end());
    pts.erase(remove(pts.begin(),pts.end(),' '),pts.end());
    stringstream pou (pts);
    
    vector<double> ptsList;
    string temp;
    while (getline(pou,temp,',')) {
        ptsList.push_back(stod(temp));
    }
    
    point a = point(ptsList[0],ptsList[1]);
    point b = point(ptsList[2],ptsList[3]);
    point c = point(ptsList[4],ptsList[5]);
    point d = point(ptsList[6],ptsList[7]);
    
    vector<Square> results;
    results.push_back(process(a,b,c,d,true));
    results.push_back(process(a,b,c,d,false));
    results.push_back(process(a,c,b,d,true));
    results.push_back(process(a,c,b,d,false));
    results.push_back(process(a,d,b,c,true));
    results.push_back(process(a,d,b,c,false));
    double min = 10;
    Square minSquare;
    
    for (int i=0; i<results.size(); i++) {
        Square curSquare = results.at(i);
        if (curSquare.getArea()<min) {
            minSquare = curSquare;
            min = curSquare.getArea();
        }
        txtFile << curSquare.toString() << "\n";
    }
    
    
    
    circle(2.0/800,a.getX(),a.getY(),"blue");
    circle(2.0/800,b.getX(),b.getY(),"blue");
    circle(2.0/800,c.getX(),c.getY(),"blue");
    circle(2.0/800,d.getX(),d.getY(),"blue");
    
    
    
    minSquare.drawSquare();
    
    fstream imgFile;
    imgFile.open("output.ppm",ios::out);
    imgFile << "P3 " << xlen << " " << ylen << " 120 \n";
    for (int ay=0; ay<ylen; ay++) {
        for (int ax=0; ax<xlen; ax++) {
            imgFile << pixels[ay][ax][0] << " " << pixels[ay][ax][1] << " " << pixels[ay][ax][2] << " ";
        }
        imgFile << "\n";
    }
    imgFile.close();
}

int main() {
    
    for (int ye=0; ye<ylen; ye++) {
        for (int xe=0; xe<xlen; xe++) {
            pixels[ye][xe][0]=120;
            pixels[ye][xe][1]=120;
            pixels[ye][xe][2]=120;
        }
    }
    part2();
    
    
    return 0;
}