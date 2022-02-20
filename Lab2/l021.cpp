


#include <iostream>
#include <fstream>
#include <math.h>
#include <iomanip>
using namespace std;



double getArea(double p1x,double p1y,double p2x,double p2y,double p3x,double p3y) {
    double a = sqrt(pow(p1x-p2x,2)+pow(p1y-p2y,2));
    double b = sqrt(pow(p2x-p3x,2)+pow(p2y-p3y,2));
    double c = sqrt(pow(p3x-p1x,2)+pow(p3y-p1y,2));
    double s = (a+b+c)/2;
    double area = sqrt(s*(s-a)*(s-b)*(s-c));
    return area;
}

void generatePoints() {
    srand(time(NULL));
    double p1 [2] = {((double) rand() / RAND_MAX),((double) rand() / RAND_MAX)};
    double p2 [2] = {((double) rand() / RAND_MAX),((double) rand() / RAND_MAX)};
    double p3 [2] = {((double) rand() / RAND_MAX),((double) rand() / RAND_MAX)};
    
    
    
    cout << "(" << p1[0] << "," << p1[1] << ") , \n";
    cout << "(" << p2[0] << "," << p2[1] << ") , \n";
    cout << "(" << p3[0] << "," << p3[1] << ") , \n";
    
    bool validFourth = false;
    double p4 [2] = {0.97370024902263463, 0.13582933867278912};
    double totalArea = getArea(p1[0],p1[1],p2[0],p2[1],p3[0],p3[1]);
    double area1;
    double area2;
    double area3;
    double bigArea;
    double eps = 0.000000000001;;
    
    while (validFourth==false) {
        p4[0] = ((double) rand() / RAND_MAX);
        p4[1] = ((double) rand() / RAND_MAX);
        
        area1 = getArea(p1[0],p1[1],p2[0],p2[1],p4[0],p4[1]);
        area2 = getArea(p1[0],p1[1],p3[0],p3[1],p4[0],p4[1]);
        area3 = getArea(p2[0],p2[1],p3[0],p3[1],p4[0],p4[1]);
        
        bigArea = area1+area2+area3;
        eps = 0.000000000001;
        if (abs(totalArea-bigArea)>eps) {
            validFourth=true;
        }
    }
    
    cout << "(" << p4[0] << "," << p4[1] << ") , \n";
    
    fstream points;
    points.open("points.txt",ios::out);
    points << setprecision(17) << "(" << p1[0] << "," << p1[1] << ") , ";
    points << setprecision(17) << "("  << p2[0] << "," << p2[1] << ") , ";
    points << setprecision(17) << "("  << p3[0] << "," << p3[1] << ") , ";
    points << setprecision(17) << "("  << p4[0] << "," << p4[1] << ")";
    points.close();
    
}


int main() {
    generatePoints();
    return 0;
} 