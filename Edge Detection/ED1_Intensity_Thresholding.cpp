//Lab for Computer Vision 
//Name: Tommy Williams 

#include <sstream>
#include <iostream>
#include <fstream>
#include <math.h>
#include <iomanip>
#include <algorithm>
#include <vector>
#include <list>
#include <string>
#include <chrono>
#include <stack>

using namespace std;

int xlen = 800;
int ylen = 600;

int scaleUp(double n) {
    int newN = round(n*xlen);
    return newN;
}

class Image {
    private:
        string link;
    public:
        int xLength = 800;
        int yHeight = 800;
        vector<vector<array<int,3>>> pixels;
        Image(const int xSize,const int ySize, string lonk) {
            pixels.resize(ySize, vector<array<int,3>>(xSize));
            xLength = xSize;
            yHeight = ySize;
            link = lonk;
        }
        void set(int y, int x, int c, int value) {
            pixels[y][x][c] = value;
        }
        void setAll(int y, int x, int value) {
            pixels[y][x][0] = value;
            pixels[y][x][1] = value;
            pixels[y][x][2] = value;
        }
        int get(int y, int x, int c) {
            return pixels[y][x][c];
        }
        void drawImageToFile() {
            fstream imgFile;
            imgFile.open(link,ios::out);
            imgFile << "P3" << "\n" << xLength << " " << yHeight << "\n" << 120 << "\n";
            for (int ay=0; ay<yHeight; ay++) {
                for (int ax=0; ax<xLength; ax++) {
                imgFile << pixels[ay][ax][0] << " " << pixels[ay][ax][1] << " " << pixels[ay][ax][2] << " ";
                }
                imgFile << "\n";
            }
            imgFile.close();
        }
};

void illuminate(int x, int y,string color, Image img) {
    if (x>=xlen || y>=xlen || y<=0 || x<=0) {
        return;
    }
    if (color == "r") {
        img.set(y,x,0,100);
        img.set(y,x,1,0);
        img.set(y,x,2,0);
    }
    else if (color == "b") {
        img.set(y,x,0,0);
        img.set(y,x,1,0);
        img.set(y,x,2,100);
    }
    else if (color == "g") {
        img.set(y,x,0,0);
        img.set(y,x,1,100);
        img.set(y,x,2,0);
    }
    else {
        img.setAll(y,x,120);
    }
}

void circle(double radius, double centerX, double centerY,string color, Image img) {
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
        illuminate(cx+x,cy+y,color,img);
        illuminate(cx+x,cy-y,color,img);
        illuminate(cx-x,cy+y,color,img);
        illuminate(cx-x,cy-y,color,img);
        illuminate(cx+y,cy+x,color,img);
        illuminate(cx+y,cy-x,color,img);
        illuminate(cx-y,cy+x,color,img);
        illuminate(cx-y,cy-x,color,img);
        
        y2_new -= (2*x);
    }
    
}

void Bresenham(double ex1,double ey1,double ex2,double ey2,string color, Image img) {
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
            illuminate(i,yi,color,img);
        }
        return;
    }
    else if (dx==0) {
        for (int j=yi; j<=yf; j++) {
            illuminate(xi,j,color,img);
        }
        return;
    }
    
    int j = y1;
    int e;
    
    if (dy<=dx) { 
        e=dy-dx;
        for (int i=x1; (x2==xf) ? i<x2 : i>x2; (x2==xf) ? i++ : i--) {
            illuminate(i,j,color,img);
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
            illuminate(j,i,color,img);
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
        double angleTo(point p1) {
            double lenP1 = sqrt(x*x + y*y);
            double lenP0 = sqrt(p1.getX()*p1.getX() + p1.getY()*p1.getY());
            double dot = x*p1.getX() + y*p1.getY();
            double cosineOf = dot / (lenP1 * lenP0);
            return cosineOf;
        }
        double getSlopeTo(point p1) {
            double r = ((p1.getY()-y)/((p1.getX()-x)));
            return (-1/r);
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

void part1() {
    //make sure x and ys are different in all functions so that they are compatible with rectangles not just squares
    ifstream inputPPM ("input.ppm");
    string version;
    string xSize; string ySize;
    string colorRange;
    getline(inputPPM,version);
    getline(inputPPM,xSize,' ');
    getline(inputPPM,ySize);
    //cout << xSize << " " << ySize;
    
    xlen = stoi(xSize); ylen = stoi(ySize);
    getline(inputPPM,colorRange);
    
    Image gray = Image(xlen,ylen,"greyimagem.ppm");
    //int gray [ylen][xlen][3] = {};
    string r; string g; string b;
    for (int y=0; y<ylen; y++) {
        for (int x=0; x<xlen; x++) {
            getline(inputPPM,r,' ');
            getline(inputPPM,g,' ');
            if (x==xlen-1) {
                getline(inputPPM,b);
            }
            else {
                getline(inputPPM,b,' ');
            }
            int greyIntensity = (stoi(r)+stoi(g)+stoi(b))/3;
            gray.setAll(y,x,greyIntensity);
            //gray[ylen][xlen][0] = greyIntensity;
            //gray[ylen][xlen][1] = greyIntensity;
            //gray[ylen][xlen][2] = greyIntensity;
        }
    }
    
    gray.drawImageToFile();
    
    //Image gradient = Image(xlen,ylen,"gradient.ppm");
    int gradient[ylen][xlen] = {};
    
    for (int y=1; y<ylen-1; y++) {
        for (int x=1; x<xlen-1; x++) {
            //cout << "\nprogress";
            //xSobel[y].push_back(sobelX(y,x,gray));
            //ySobel[y].push_back(sobelY(y,x,gray));
            int sobelYe = -1 * (gray.get(y+1,x+1,0) + 2*(gray.get(y+1,x,0)) + gray.get(y+1,x-1,0)) + gray.get(y-1,x-1,0) + 2*gray.get(y-1,x,0) + gray.get(y-1,x+1,0);
            int sobelXe = gray.get(y+1,x+1,0) - gray.get(y+1,x-1,0) + 2*gray.get(y,x+1,0) - 2*gray.get(y,x-1,0) + gray.get(y-1,x+1,0) - gray.get(y-1,x-1,0);
            gradient[y][x] = pow(sobelXe,2) + pow(sobelYe,2);
            //gradient.setAll(y,x,pow(sobelX(y,x,gray),2) + pow(sobelY(y,x,gray),2));
        }
    } 
    
    Image edgeThingy = Image(xlen,ylen,"imagem.ppm");
    int threshold = 20000;
    
    for (int y=0; y<ylen; y++) {
        for (int x=0; x<xlen; x++) {
            if (gradient[y][x] >= threshold) {
                edgeThingy.setAll(y,x,120);
            }
        }
    }
    edgeThingy.drawImageToFile();
}


int main() {

    part1();
    
    return 0;
}
