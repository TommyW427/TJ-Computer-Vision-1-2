//Lab for Computer Vision 
//Date: 12/2/2021
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

int xlen = 4000;
int ylen = 4000;

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
        int m;
        vector<vector<array<int,3>>> pixels;
        Image(const int xSize,const int ySize, string lonk, int maxIntensity) {
            pixels.resize(ySize, vector<array<int,3>>(xSize));
            /*for (int y=0; y<ySize; y++) {
                pixels[y] = vector<array<int,3>>(xLength);
            }*/
            m = maxIntensity;
            xLength = xSize;
            yHeight = ySize;
            link = lonk;
            
        }
        pair<int,int> getDimensions() {
            return make_pair(xLength,yHeight);
        }
        double getAverageMax() {
            return (xLength+yHeight)/2;
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
            return pixels.at(y).at(x).at(c);
        }
        void drawImageToFile() {
            fstream imgFile;
            imgFile.open(link,ios::out);
            imgFile << "P3" << "\n" << xLength << " " << yHeight << "\n" << m << "\n";
            for (int ay=0; ay<yHeight; ay++) {
                for (int ax=0; ax<xLength; ax++) {
                    imgFile << pixels[ay][ax][0] << " " << pixels[ay][ax][1] << " " << pixels[ay][ax][2] << " ";
                }
                imgFile << "\n";
            }
            imgFile.close();
        }
        vector<pair<int,int>> getEdges(int y, int x) {
            vector<pair<int,int>> edges;
            for (int hor=x-1; hor<=x+1; hor++) {
                if (hor < xLength && hor >= 0) {
                    for (int ver=y-1; ver<=y+1; ver++) {
                        if (hor!=x || ver!=y) {
                            if (ver >=0 && ver < yHeight) {
                                edges.push_back(make_pair(ver,hor));
                            }
                        }
                    }
                }
            }
            return edges;
        }
        void setMaxValue(int mex) {
            m = mex;
        }
};

void illuminate(int x, int y,string color, Image* img) {
    if (x>=xlen || y>=ylen || y<=0 || x<=0) {
        return;
    }
    if (color == "r") {
        img->set(y,x,0,100);
        img->set(y,x,1,0);
        img->set(y,x,2,0);
        //pixels[y][x][0]=100;
        //pixels[y][x][1]=0;
        //pixels[y][x][2]=0;
    }
    else if (color == "b") {
        img->set(y,x,0,0);
        img->set(y,x,1,0);
        img->set(y,x,2,1);
        //pixels[y][x][0]=0;
        //pixels[y][x][1]=0;
        //pixels[y][x][2]=100;
    }
    else if (color == "g") {
        img->set(y,x,0,0);
        img->set(y,x,1,1);
        img->set(y,x,2,0);
        //pixels[y][x][0]=0;
        //pixels[y][x][1]=100;
        //pixels[y][x][2]=0;
    }
    else {
        img->setAll(y,x,120);
        //pixels[y][x][0]=120;
        //pixels[y][x][1]=120;
        //pixels[y][x][2]=120;
    }
    return;
}

void circle(int rad, int cx, int cy,string color, Image* img) {
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
        //cout << "hello there";
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
    return;
}

void Bresenham(double ex1,double ey1,double ex2,double ey2,string color, Image& img) {
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
            illuminate(i,yi,color,&img);
        }
        return;
    }
    else if (dx==0) {
        for (int j=yi; j<=yf; j++) {
            illuminate(xi,j,color,&img);
        }
        return;
    }
    
    int j = y1;
    int e;
    
    if (dy<=dx) { 
        e=dy-dx;
        for (int i=x1; (x2==xf) ? i<x2 : i>x2; (x2==xf) ? i++ : i--) {
            illuminate(i,j,color,&img);
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
            illuminate(j,i,color,&img);
            if (e>=0) {
                x1==xi ? j++ : j--;
                e-=dy;
                
            }
            e+=dx;
        }
    }  
}

void addVote(int y, int x, Image img) {
    int prevValue = img.get(y,x,0);
    img.setAll(y,x,prevValue+1);
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
        double distanceToXY(int ye, int xe) {
            //yes sqrt
            return sqrt(pow(y-ye,2)+pow(x-xe,2));
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

int addFunction(int d, double avgLen) {
    //cout << d << " " << avgLen << "\n";
    if (d>175) return 0;
    if (d>avgLen/2) {
        return 1;
    }
    else if (d>avgLen/4) {
        return 2;
    }
    else if (d>avgLen/10 && d>150) {
        return 3;
    }
    else if (d>avgLen/22 && d>150) {
        return 4;
    }
    else if (d>150) {
        return 5;
    }
    else if (d>100) {
        return 8;
    }
    else if (d>75) {
        return 6;
    }
    else {
        return 3;
    }
}

void voteBresenham(int y1, int x1, int y2, int x2, point edge, Image* votingMatrix, Image* edgeMatrix) {
    
    double averageMax = votingMatrix->getAverageMax();
    
    int maxVote = 0;
    
    int dx = abs(x2-x1);
    int dy = abs(y2-y1);
    
    int xi = min(x1,x2);
    int xf = max(x1,x2);
    int yi = min(y1,y2);
    int yf = max(y1,y2);
    
    //without edge-breaking//
    
    if (dy==0) {
        for (int i=xi; i<=xf; i++) {
            double toAdd = addFunction(edge.distanceToXY(yi,i),averageMax);
            int prevValue = votingMatrix->get(yi,i,0);
            votingMatrix->setAll(yi,i,prevValue+1+toAdd);
            //if (i!=xi && edgeMatrix->get(yi,i,0)==1) {
            //    break;
            //}
        }
        return;
    }
    else if (dx==0) {
        for (int j=yi; j<=yf; j++) {
            int toAdd = addFunction(edge.distanceToXY(j,xi),averageMax);
            int prevValue = votingMatrix->get(j,xi,0);
            votingMatrix->setAll(j,xi,prevValue+toAdd);
            
            //if (j!=yi && edgeMatrix->get(j,xi,0)==1) {
            //    break;
            //}
        }
        return;
    }
    
    int j = y1;
    int e;
    
    if (dy<=dx) { 
        e=dy-dx;
        for (int i=x1; (x2==xf) ? i<x2 : i>x2; (x2==xf) ? i++ : i--) {
            int toAdd = addFunction(edge.distanceToXY(j,i),averageMax);
            int prevValue = votingMatrix->get(j,i,0);
            votingMatrix->setAll(j,i,prevValue+1+toAdd);
            //if (i!=x1 && j!=y1 && edgeMatrix->get(j,i,0)==1) {
            //        break;
            //}
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
            double toAdd = addFunction(edge.distanceToXY(i,j),averageMax);
            int prevValue = votingMatrix->get(i,j,0);
            votingMatrix->setAll(i,j,prevValue+1+toAdd);
            //if (i!=y1 && j!=x1 && edgeMatrix->get(i,j,0)==1) {
            //       break;
            //}
            if (e>=0) {
                x1==xi ? j++ : j--;
                e-=dy;
            }
            e+=dx;
        }
    }
    return;
}

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

void recurFlow(int currentY, int currentX, Image& visitedImage, Image& thresholdImage) {
    int value = thresholdImage.get(currentY,currentX,1);
    if (value==2 || value==1) {
        visitedImage.setAll(currentY,currentX,1);
        vector<pair<int,int>> edges = thresholdImage.getEdges(currentY,currentX);
        for (int i=0; i<edges.size(); i++) {
            if (visitedImage.get(edges.at(i).first,edges.at(i).second,1)==0) {
                recurFlow(edges.at(i).first,edges.at(i).second,visitedImage,thresholdImage);
            };
        }
    }
}



void part1(string input, int thres1, int thres2, int thres3, int subd) {
    ifstream inputPPM (input);
    string version;
    string xSize; string ySize;
    string colorRange;
    getline(inputPPM,version);
    getline(inputPPM,xSize,' ');
    while (xSize[0]=='#') {
        getline(inputPPM,xSize,' ');
    }
    getline(inputPPM,ySize);
    cout << "Image sized " << xSize << " by " << ySize << "\n";
    
    xlen = stoi(xSize); ylen = stoi(ySize);
    getline(inputPPM,colorRange);
    Image gray = Image(xlen,ylen,"imageg.ppm",120);
    Image finalWCircles = Image(xlen,ylen,"imageCC.ppm",stoi(colorRange));
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
            int greyIntensity = stoi(r)*0.2126+stoi(g)*0.7152+stoi(b)*0.0722;
            gray.setAll(y,x,greyIntensity);
            finalWCircles.set(y,x,0,stoi(r));
            finalWCircles.set(y,x,1,stoi(g));
            finalWCircles.set(y,x,2,stoi(b));
        }
    }
    
    gray.drawImageToFile();
    vector<vector<int>> gradient(ylen);
    vector<vector<int>> angular(ylen);
    cout << ylen << " " << xlen;
    
    for (int x=0; x<xlen; x++) {
        gradient.at(0).push_back(0);
        gradient.at(ylen-1).push_back(0);
        angular.at(0).push_back(0);
        angular.at(ylen-1).push_back(0);
    }
    
    for (int y=1; y<ylen-1; y++) {
        gradient.at(y).push_back(0);
        angular.at(y).push_back(0);
        for (int x=1; x<xlen-1; x++) {
            //cout << x << " " << y << "\n"; 
            int sobelYe = -1 * (gray.get(y+1,x+1,0) + 2*(gray.get(y+1,x,0)) + gray.get(y+1,x-1,0)) + gray.get(y-1,x-1,0) + 2*gray.get(y-1,x,0) + gray.get(y-1,x+1,0);
            int sobelXe = gray.get(y+1,x+1,0) - gray.get(y+1,x-1,0) + 2*gray.get(y,x+1,0) - 2*gray.get(y,x-1,0) + gray.get(y-1,x+1,0) - gray.get(y-1,x-1,0);
            gradient.at(y).push_back(pow(sobelXe,2) + pow(sobelYe,2));
            int degAngle = atan2(sobelYe,sobelXe) * 180/3.14159265;
            angular.at(y).push_back(degAngle);
            
            //gradient.setAll(y,x,pow(sobelX(y,x,gray),2) + pow(sobelY(y,x,gray),2));
        }
        gradient.at(y).push_back(0);
        angular.at(y).push_back(0);
    } 
    Image angularBoogwa = Image(xlen,ylen,"image2.ppm",1);
    pair<int,int> otherCoord1;
    pair<int,int> otherCoord2;
    for (int y=1; y<ylen-1; y++) {
        for (int x=1; x<xlen-1; x++) {
            int angleProx = 0;
            double mod45 = angular.at(y).at(x) % 45;
            if (mod45 > 22.5) {
                angleProx = angular.at(y).at(x)+(45-mod45);
            }
            else {
                angleProx = angular.at(y).at(x)-mod45;
            }
            if (angleProx == 45 || angleProx == -135) {
                otherCoord1 = make_pair(y-1,x+1); otherCoord2 = make_pair(y+1,x-1);
            }
            else if (angleProx == 90 || angleProx == -90) {
                otherCoord1 = make_pair(y-1,x); otherCoord2 = make_pair(y+1,x);
            }
            else if (angleProx == 0 || angleProx >= 180 || angleProx <= -180) {
                otherCoord1 = make_pair(y,x+1); otherCoord2 = make_pair(y,x-1);
            }
            else if (angleProx == -45 || angleProx == 135) {
                otherCoord1 = make_pair(y+1,x+1); otherCoord2 = make_pair(y-1,x-1);
            }
            else {
                cout << angleProx;
            }
            
            if (gradient.at(y).at(x) > max(gradient.at(otherCoord1.first).at(otherCoord1.second),gradient.at(otherCoord2.first).at(otherCoord2.second))) {
                angularBoogwa.setAll(y,x,1);
            }
        }
    }
    angularBoogwa.drawImageToFile();
    
    
    Image edgeThingy = Image(xlen,ylen,"imagem.ppm",2);
    int threshold1 = pow(thres1,2);
    int threshold2 = pow(thres2,2);
    vector<pair<int,int>> ooga;
    
    for (int y=1; y<ylen-1; y++) {
        for (int x=1; x<xlen-1; x++) {
            if (gradient.at(y).at(x) >= threshold2) {
                edgeThingy.setAll(y,x,2);
                ooga.push_back(make_pair(y,x));
            }
            else if (gradient.at(y).at(x) >= threshold1) {
                edgeThingy.setAll(y,x,1);
            }
        }
    }
    edgeThingy.drawImageToFile();
    
    Image visitedImage = Image(xlen,ylen,"image1.ppm",1);
    for (int i=0; i<ooga.size(); i++) {
        int startIndexY = ooga.at(i).first;
        int startIndexX = ooga.at(i).second;
        recurFlow(startIndexY,startIndexX,visitedImage,edgeThingy);
    }
    
    Image finalImage = Image(xlen,ylen,"imagef.ppm",1);
    for (int y=1; y<ylen-1; y++) {
        for (int x=1; x<xlen-1; x++) {
            if (visitedImage.get(y,x,0)==1 and angularBoogwa.get(y,x,0)==1) {
                finalImage.setAll(y,x,1);
            }
        }
    }
    
    visitedImage.drawImageToFile();
    finalImage.drawImageToFile();
    
    Image voteImage = Image(xlen,ylen,"imagev.ppm",200);
    for (int y=1; y<ylen-1; y++) {
        for (int x=1; x<xlen-1; x++) {
            point edge = point(x,y);
            if (finalImage.get(y,x,0) == 1) {
                double radianAngle = angular.at(y).at(x) * 3.14159265/180;
                //cout << x << " " << y << "\n";
                //cout << radianAngle << "\n";
                int topX = x + (y)/tan(radianAngle);
                int botX = x - (ylen-1-y)/tan(radianAngle);
                int topY = y-(xlen-1-x)*tan(radianAngle);
                int botY = y+x*tan(radianAngle);
                //cout << topX << " " << botX << "\n";
                //cout << topY << " " << botY << "\n";
                //cout << "\n\n";
                if (topX < xlen && topX >= 0) {
                    if (topY < ylen && topY>=0) {
                        voteBresenham(0,topX,topY,xlen-1,edge,&voteImage,&finalImage);
                    }
                    else if (botY < ylen && botY>=0) {
                        voteBresenham(0,topX,botY,0,edge,&voteImage,&finalImage);
                    }
                    else {
                        voteBresenham(0,topX,ylen-1,botX,edge,&voteImage,&finalImage);
                    }
                }
                else if (botX < xlen && botX >=0) {
                    if (topY < ylen && topY >= 0) {
                        voteBresenham(ylen-1,botX,topY,xlen-1,edge,&voteImage,&finalImage);
                    }
                    else if (botY < ylen && botY>=0) {
                        voteBresenham(ylen-1,botX,botY,0,edge,&voteImage,&finalImage);
                    }
                    else {
                        voteBresenham(0,topX,ylen-1,botX,edge,&voteImage,&finalImage);
                    }
                }
                else {
                    voteBresenham(topY,xlen-1,botY,0,edge,&voteImage,&finalImage);
                }    
            }
        }
    }
    
    //100 seperate regions 
    //in each region get density of center matches
    //if density is high increase threshold accordingly and reapply threshold in that region 
    
    vector<pair<int,int>> glob;
    vector<pair<int,int>> local;
    int threshold3 = thres3;
    if (xlen>800 && ylen>800) {
    int yincr = ylen/subd;
    int xincr = xlen/subd;
    cout << "\n";
    cout << xincr << " " << subd << " " << xincr*subd << "\n";
    cout << yincr << " " << subd << " " << yincr *subd << "\n";
    for (int j=0; j<subd; j++) {
        for (int i=0; i<subd; i++) {
            int count=0;
            local.clear();
            for (int p=j*yincr; p<(j+1)*yincr; p++) {
                for (int q=i*xincr; q<(i+1)*xincr; q++) {
                    if (voteImage.get(p,q,0)>=threshold3) {
                        count++;
                        //if (j*yincr>3024) {cout << i*xincr << " " << j*yincr << "\n";}
                        local.push_back(make_pair(p,q));
                    }
                }
            }
            double density = count/((float) xincr*yincr);
            //cout << "Density in section " << j << " " << i << " is " << density << " with " << count << "matches" << "\n";
            if (density>0.15) {
                double newThreshold3 = 0;
                if (density>0.8) {newThreshold3 = threshold3 * (1+0.8*density);}
                else if (density<0.4) {newThreshold3 = threshold3 * 0.6;}
                else {newThreshold3 = threshold3 * (1+0.5*density);}
                for (int i=0; i<local.size(); i++) {
                    int p = local.at(i).first; int q = local.at(i).second;
                    if (voteImage.get(p,q,0)>=newThreshold3) {
                        //if (j*yincr>3024) {cout << i*xincr << " " << j*yincr << " " << p << " " << q << "\n";}
                        glob.push_back(make_pair(p,q));
                    }
                }
            }
            
            else {
                glob.insert(glob.end(),local.begin(), local.end());
            }
        }
    }
    }
    
    for (int u=0; u<glob.size(); u++) {
        int y=glob.at(u).first; int x=glob.at(u).second;
        //if (y>xlen) cout << y << "\n";
        circle(4,x,y,"r",&finalWCircles);
    }

    cout << "Number of Centers: " << glob.size() << "\n";
    /*for (int y=1; y<ylen-1; y++) {
        for (int x=1; x<xlen-1; x++) {
            if (voteImage.get(y,x,0)>=threshold3) {
                //cout << finalWCircles.get(y,x+2,0) << " ";
                circle(4,x,y,"r",&finalWCircles);
                //cout << finalWCircles.get(y,x+2,0) << "\n";
            }
        }
    }*/
    voteImage.drawImageToFile();
    finalWCircles.drawImageToFile();
}

int main(int argc, char** argv) {
    string image = "image.ppm";
    int threshold1 = 110; 
    int threshold2 = 200;
    int threshold3 = 300;
    int subd = 10;
    for (int i=0; i<argc; i++) {
        
        if (string(argv[i])=="-L" && argc>i+1) {
            threshold1 = atoi(argv[i+1]);
        }
        else if (string(argv[i])=="-H" && argc>i+1) {
            threshold2 = atoi(argv[i+1]);
        }
        else if (string(argv[i])=="-F" && argc>i+1) {
            cout << "hi";
            image = string(argv[i+1]);
        }
        else if (string(argv[i])=="-TC" && argc>i+1) {
            threshold3 = atoi(argv[i+1]);
        }
        else if (string(argv[i])=="-SD" && argc>i+1) {
            subd = atoi(argv[i+1]);
        }
    }
    cout << "Checking for centers on " << image << "\n";
    cout << "Lower Edge Threshold: " << threshold1 << " Upper Edge Threshold: " << threshold2 << "\n";
    cout << "Angular Threshold: " << threshold3 << "\n";
    part1(image,threshold1,threshold2,threshold3,subd);
    
    return 0;
}
