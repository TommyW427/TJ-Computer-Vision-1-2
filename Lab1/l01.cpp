//First Lab for Computer Vision
//Name: Tommy Williams 
//Date: 9/7/2021
#include <iostream> 
#include <math.h>
#include <fstream>
#include <stdio.h>
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
    if (color == "red") {
        pixels[y][x][0]=100;
        pixels[y][x][1]=0;
        pixels[y][x][2]=0;
    }
    else if (color == "blue") {
        pixels[y][x][0]=0;
        pixels[y][x][1]=0;
        pixels[y][x][2]=100;
    }
    else if (color == "green") {
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
        
        y2_new -= (2*x); //make sure this is right 
    }
    
}

double getInRadius(double x1, double y1, double x2, double y2, double x3, double y3) {
    double l1,l2,l3;
    l1 = sqrt(pow(x2-x1,2)+pow(y2-y1,2));
    l2 = sqrt(pow(x3-x2,2)+pow(y3-y2,2));
    l3 = sqrt(pow(x1-x3,2)+pow(y1-y3,2));
    double semi = 0.5*(l1+l2+l3);
    double area = sqrt(semi * (semi-l1) * (semi-l2) * (semi-l3));
    return (double) area/semi;
    //return (double) sqrt(((semi-l1)*(semi-l2)*(semi-l3))/semi);
}

double getOutRadius(double x1, double y1, double x2, double y2, double x3, double y3) {
    double l1,l2,l3;
    l1 = sqrt(pow(x2-x1,2)+pow(y2-y1,2));
    l2 = sqrt(pow(x3-x2,2)+pow(y3-y2,2));
    l3 = sqrt(pow(x1-x3,2)+pow(y1-y3,2));
    double semi = 0.5*(l1+l2+l3);
    double smarl = getInRadius(x1,y1,x2,y2,x3,y3);
    return (l1*l2*l3)/(4*semi*smarl);
}

double* getOutCenter(double x1, double y1, double x2, double y2, double x3, double y3,double circumcenter[]) {
    double slope1 = -1 * 1/((y2-y1)/(x2-x1));
    double slope2 = -1 * 1/((y3-y2)/(x3-x2));
    
    double midY1 = (y2+y1)/2;
    double midY2 = (y3+y2)/2;
    double midX1 = (x2+x1)/2;
    double midX2 = (x3+x2)/2;
   
    double xValue1 = slope1;
    double const1 = slope1*(-1)*midX1+midY1; 
    double xValue2 = slope2;
    double const2 = slope2*(-1)*midX2+midY2;
    
    double centerX,centerY;

    if (isnan(slope1)) {
        centerX = midX1;
        centerY = midX1*xValue2+const2;
    }
    else if (isinf(slope1)) {
        centerX = (midX1-const2)/xValue2; 
        centerY = centerX*xValue1+const1;
    }
    else if (isnan(slope2)) {
        centerX = midX2;
        centerY = midX2*xValue1+const1;
    }
    else if (isinf(slope2)) {
        centerX = (midX2-const1)/xValue1;
        centerY = centerX*xValue1+const1;
    }
    else {
       centerX = (const2-const1)/(xValue1-xValue2); 
       centerY = centerX*xValue1+const1;
    }
    
    cout << "CIRCUMCENTER:\n";
    cout << centerX << " " << centerY << "\n";
    circumcenter[0] = centerX;
    circumcenter[1] = centerY;
    return circumcenter;
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
            if (e>=0) {
                x1==xi ? j++ : j--;
                e-=dy;
                
            }
            e+=dx;
        }
    }  
}
double* getInCenter(double x1, double y1, double x2, double y2, double x3, double y3,double incenter[]) {
    
    double c = sqrt(pow(x2-x1,2)+pow(y2-y1,2));
    double a = sqrt(pow(x3-x2,2)+pow(y3-y2,2));
    double b = sqrt(pow(x1-x3,2)+pow(y1-y3,2));
    
    double incenterX = (c*x3+a*x1+b*x2)/(a+b+c);
    double incenterY = (c*y3+a*y1+b*y2)/(a+b+c);
    
    incenter[0] = incenterX;
    incenter[1] = incenterY;
    
    return incenter;
}

void euler(double x1, double y1, double x2, double y2, double x3, double y3,double circCent[]) {
    double centroidX = (x1+x2+x3)/3;
    double centroidY = (y1+y2+y3)/3;
    double slope = (circCent[1]-centroidY)/(circCent[0]-centroidX);
   
    if (isinf(slope)) {
        Bresenham(centroidX,0,centroidX,1,"red");
    }
    else if (isnan(slope)) {
        Bresenham(0,centroidY,1,centroidY,"red");
    }
    if (abs(slope)>=1) {
        double xValue = slope;
        double constValue = slope*(-1)*centroidX + centroidY;
        double xmin = -1*constValue/xValue;
        double xmax = (1-constValue)/xValue;
        Bresenham(xmin,0,xmax,1,"red");
        cout << "Euler Line: y = " << xValue << "x + " << constValue << "\n";
    }
    else {
        cout << "Euler Line: y = " << slope << "x + " << slope*(-1)*centroidX + centroidY << "\n";
        double ymin = slope*(0-centroidX)+centroidY;
        double ymax = slope*(1-centroidX)+centroidY;
        Bresenham(0,ymin,1,ymax,"red");
    }
    return;
}


void nineCircle(double x1, double y1, double x2, double y2, double x3, double y3, double circumcenter[]) {
    double centroidX = (x1+x2+x3)/3;
    double centroidY = (y1+y2+y3)/3;
    double orthX = 3*centroidX-2*circumcenter[0];
    double orthY = 3*centroidY-2*circumcenter[1];
    
    double x9 = (orthX+circumcenter[0])/2;
    double y9 = (orthY+circumcenter[1])/2;
    double rad = getOutRadius(x1,y1,x2,y2,x3,y3)/2;
    
    circle(rad,x9,y9,"green");
    
    cout << "Orthocenter: [" << orthX << ", " << orthY << "]\n";
    cout << "Centroid: [" << centroidX << ", " << centroidY << "]\n";
    cout << "9 point circle: Center = [" << x9 << ", " << y9 << "]; Radius = " << rad << "\n";
    
    illuminate(orthX,orthY,"blue");
    illuminate(centroidX,centroidY,"blue");
    illuminate(circumcenter[0],circumcenter[1],"blue");
    
    //NOTES:
    //Radius = 0.5circumradius 
    //Center = 1/4 of the way from centroid to orthocenter on the euler line 
}

int main() {
    srand(time(NULL));
    double p1 [2] = {((double) rand() / RAND_MAX),((double) rand() / RAND_MAX)};
    double p2 [2] = {((double) rand() / RAND_MAX),((double) rand() / RAND_MAX)};
    double p3 [2] = {((double) rand() / RAND_MAX),((double) rand() / RAND_MAX)};
    //double p1[2] = {, 10/800};
    //double p2[2] = {1/8, 1/8};
    //double p3[2] = {10/800, 1/8};
    cout << "Point 1: [" << scaleUp(p1[0]) << ", " << scaleUp(p1[1]) << "] \n";
    cout << "Point 2: [" << scaleUp(p2[0]) << ", " << scaleUp(p2[1]) << "] \n";
    cout << "Point 3: [" << scaleUp(p3[0]) << ", " << scaleUp(p3[1]) << "] \n";
    //Creating file and setup//
    fstream imgFile;
    imgFile.open("triangle.ppm",ios::out);
    
    for (int ye=0; ye<ylen; ye++) {
        for (int xe=0; xe<xlen; xe++) {
            pixels[ye][xe][0]=0;
            pixels[ye][xe][1]=0;
            pixels[ye][xe][2]=0;
        }
    }
    
    double inRad = getInRadius(p1[0],p1[1],p2[0],p2[1],p3[0],p3[1]);
    double outRad = getOutRadius(p1[0],p1[1],p2[0],p2[1],p3[0],p3[1]);
    
    Bresenham(p1[0],p1[1],p2[0],p2[1],blue);
    Bresenham(p2[0],p2[1],p3[0],p3[1],blue);
    Bresenham(p3[0],p3[1],p1[0],p1[1],blue);
    
    double *circumcenter = new double[2];
    double *incenter = new double[2];
    circumcenter = getOutCenter(p1[0],p1[1],p2[0],p2[1],p3[0],p3[1],circumcenter);
    incenter = getInCenter(p1[0],p1[1],p2[0],p2[1],p3[0],p3[1],incenter);
    
    cout << "CIRCUMCIRCLE\n";
    cout << outRad << " " << circumcenter[0] << " " << circumcenter[1] << "\n";
    cout << "INCIRCLE: " << inRad << " " << incenter[0] << " " << incenter[1] << "\n";
    
    circle(outRad,circumcenter[0],circumcenter[1],"green");
    circle(inRad,incenter[0],incenter[1],"green");
    euler(p1[0],p1[1],p2[0],p2[1],p3[0],p3[1],circumcenter);
    nineCircle(p1[0],p1[1],p2[0],p2[1],p3[0],p3[1],circumcenter);
    
    illuminate(p1[0]*xlen,p1[1]*ylen,"b");
    illuminate(p2[0]*xlen,p2[1]*ylen,"b");
    illuminate(p3[0]*xlen,p3[1]*ylen,"b");
    
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
