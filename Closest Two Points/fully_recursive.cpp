//Name: Tommy Williams 
//Lab 3 Closest Two Points
//Approach 3: Fully recursive 
//Time complexity: O(nlg^2n)


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

const int xlen = 800;
const int ylen = 800;
int pixels [ylen][xlen][3] = {};
chrono::duration<double> t1;
chrono::duration<double> t2;
chrono::duration<double> t3;

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
        double nosqrtDistanceTo(point o) {
            return pow(x-o.getX(),2)+pow(y-o.getY(),2);
        }
        string toString() {
            ostringstream ste;
            ste << setprecision(17) << "(" << x << "," << y << ")";
            return(ste.str());
        }
        
};

class twoPoint {
    private:
        double dist;
        point p1;
        point p2;
    public:
        twoPoint(point point1, point point2) {
            p1 = point1;
            p2 = point2;
            dist = point1.distanceTo(point2);
        }
        twoPoint(point point1, point point2, double d) {
            p1 = point1; 
            p2 = point2; 
            dist = d;
        }
        point gp1() {
            return p1;
        }
        point gp2() {
            return p2;
        }
        double getDist() {
            return dist;
        }
        double getRealDist() {
            return sqrt(dist);
        }
        bool operator<(twoPoint other) const {
            return dist < other.getDist();
        }
        bool operator>(twoPoint other) const {
            return dist > other.getDist();
        }
};

bool xSort(point a, point b) {
    return (a.getX() < b.getX());
}
bool ySort(point a, point b) {
    return (a.getY() < b.getY());
}

vector<point> part1() {
    
    for (int ye=0; ye<ylen; ye++) {
        for (int xe=0; xe<xlen; xe++) {
            pixels[ye][xe][0]=120;
            pixels[ye][xe][1]=120;
            pixels[ye][xe][2]=120;
        }
    }
    
    srand(time(NULL));
    list<point> pts;
    //vector<points> psv;
    for (int i=0; i<6000; i++) {
        pts.push_back(point((double) rand() / RAND_MAX,((double) rand() / RAND_MAX)));
    }
    //recur(ptsv,)
    double minDistance = 2.0;
    point minPoint1;
    point minPoint2;
    int count=0;
    
    fstream points;
    points.open("points.txt",ios::out);
    
    chrono::time_point<chrono::system_clock> start, end;
    start = chrono::system_clock::now();
    for (list<point>::iterator pl = pts.begin(); pl != pts.end(); ++pl) {
        double lx = (*pl).getX();
        double ly = (*pl).getY();
        points << setprecision(17) << lx << "  " << ly << "\n";
        circle(3.0/800,lx,ly,"black");
        list<point>::iterator ps = pl;
        for (ps = ++ps; ps != pts.end(); ++ps) {
            double d = (*pl).distanceTo(*ps);
            count++;
            if (d<minDistance) {
                minDistance = d;
                minPoint1 = *pl;
                minPoint2 = *ps;
            }
        }
    }
    end = chrono::system_clock::now();
    t1 = end-start;
    
    circle(2.0/800,minPoint1.getX(),minPoint1.getY(),"r");
    circle(2.0/800,minPoint2.getX(),minPoint2.getY(),"r");
    circle(3.0/800,minPoint1.getX(),minPoint1.getY(),"r");
    circle(3.0/800,minPoint2.getX(),minPoint2.getY(),"r");
    points.close();
    
    vector<point> toReturn;
    toReturn.push_back(minPoint1);
    toReturn.push_back(minPoint2);
    return toReturn;
}



twoPoint recur(vector<point> &section, int i, int j) {   
    if (j-i==1) {
        return twoPoint(section.at(i),section.at(j));
    }
    
    else if (j-i==2) {
        double d1 = section.at(i).distanceTo(section.at(i+1));
        double d2 = section.at(i).distanceTo(section.at(j));
        double d3 = section.at(i+1).distanceTo(section.at(j));
        
        if (d2<d1) {
            if (d3<d2) {
                return twoPoint(section.at(i+1),section.at(j));
            }
            else {
                return twoPoint(section.at(i),section.at(j));
            }
        }
        else {
            if (d3<d1) {
                return twoPoint(section.at(i+1),section.at(j));
            }
            else {
                return twoPoint(section.at(i),section.at(i+1));
            }
        }
    }
    
    else {
        int mid = (j+i)/2;
        double medianX = (section.at(mid).getX()+section.at(mid+1).getX())/2;
        
        twoPoint minClass = min(recur(section,i,mid),recur(section,mid+1,j));
        double minDist = minClass.getDist();
        
        int n=0;
        vector<point> lefts;
        vector<point> rights;
        point left = section.at(mid);
        point right = section.at(mid+1);
        double rlDist;
        
        //double realDist = sqrt(minDist);
        while (medianX-left.getX()<minDist and mid-n>=i) {
            left = section.at(mid-n);
            lefts.push_back(left);
            n++;
        }
        n=1;
        while (right.getX()-medianX<minDist and mid+n<=j) {
            right = section.at(mid+n);
            rights.push_back(right);
            n++;
        }
        
        for (vector<point>::iterator l = lefts.begin(); l!=lefts.end(); ++l) {
            for (vector<point>::iterator r = rights.begin(); r!=rights.end(); ++r) {
                double rl = (*l).distanceTo(*r);
                if (rl<minDist) {
                    minDist = rl;
                    minClass = twoPoint(*l,*r);
                }
            }
        }
        
        return minClass;
    }
    
}

twoPoint fullRecur(vector<point> &section, int i, int j) {
    if (j-i==1) {
        return twoPoint(section.at(i),section.at(j));
    }
    
    else if (j-i==2) {
        double d1 = section.at(i).distanceTo(section.at(i+1));
        double d2 = section.at(i).distanceTo(section.at(j));
        double d3 = section.at(i+1).distanceTo(section.at(j));
        
        if (d2<d1) {
            if (d3<d2) {
                return twoPoint(section.at(i+1),section.at(j));
            }
            else {
                return twoPoint(section.at(i),section.at(j));
            }
        }
        else {
            if (d3<d1) {
                return twoPoint(section.at(i+1),section.at(j));
            }
            else {
                return twoPoint(section.at(i),section.at(i+1));
            }
        }
    }
    
    else {
        int mid = (j+i)/2; 
        double medianX = (section.at(mid).getX()+section.at(mid+1).getX())/2;
        twoPoint minClass = min(fullRecur(section,i,mid),fullRecur(section,mid+1,j));
        double minDist = minClass.getDist();
        point mp1 = minClass.gp1();
        point mp2 = minClass.gp2();
        
        
        int n=0;
        vector<point> stripPts;
        while (medianX-section.at(mid-n).getX()<minDist and mid-n>i) {
            stripPts.push_back(section.at(mid-n));
            n++;
        }
        int leftIndex = mid-n;
        n=1;
        while (section.at(mid+n).getX()-medianX<minDist and mid+n<j) {
            stripPts.push_back(section.at(mid+n));
            n++;
        }
        int rightIndex = mid+n;
        //vector<point> stripPts(section.begin()+leftIndex,section.begin()+rightIndex);
        
        sort(stripPts.begin(),stripPts.end(),ySort);
        
        //((*b).getY()-(*l).getY())<minDist
        double tempDist;
        for (vector<point>::iterator l = stripPts.begin(); l<stripPts.end(); ++l) {
            for (vector<point>::iterator b = l+1; b<l+7 and b<stripPts.end(); ++b) {
                //cout << (*l).toString() << " " << (*b).toString() << "\n";
                tempDist = (*l).distanceTo(*b);
                if (tempDist<minDist) {
                    minDist = tempDist;
                    mp1 = *l;
                    mp2 = *b;
                }
            }
        }
        
        return twoPoint(mp1,mp2);
    }
    
}

twoPoint part2() {
    ifstream pointReader ("points.txt");
    string ptsString;
    vector<point> ptsVector;
    while (getline(pointReader,ptsString)) {
        stringstream pro (ptsString);
        string X;
        string Y; 
        getline(pro,X,' ');
        pro.ignore();
        getline(pro,Y,' ');
        double x = stod(X);
        double y = stod(Y);
        ptsVector.push_back(point(x,y));
    }
    chrono::time_point<chrono::system_clock> start, end;
    start = chrono::system_clock::now();
    sort(ptsVector.begin(),ptsVector.end(),xSort);
    twoPoint finale = recur(ptsVector,0,ptsVector.size()-1);
    end = chrono::system_clock::now();
    t2 = end-start;
    return finale;
}

twoPoint part3() {
    ifstream pointReader ("points.txt");
    vector<point> ptsVector;
    string ptsString;
    while (getline(pointReader,ptsString)) {
        stringstream pro (ptsString);
        string X;
        string Y; 
        getline(pro,X,' ');
        pro.ignore();
        getline(pro,Y,' ');
        double x = stod(X);
        double y = stod(Y);
        ptsVector.push_back(point(x,y));
    }
    chrono::time_point<chrono::system_clock> start, end;
    start = chrono::system_clock::now();
    sort(ptsVector.begin(),ptsVector.end(),xSort);
    twoPoint finale = fullRecur(ptsVector,0,ptsVector.size()-1);
    end = chrono::system_clock::now();
    t3 = end-start;
    return finale;
}

int main() {
    //vector<points> p1 = part1();
    /*fstream points("points.txt",ios::out);
    for (int i=0; i<1000000; i++) {
        points << setprecision(17) << (double) rand() / RAND_MAX << "  " << (double) rand() / RAND_MAX << "\n";
    }
    points.close();*/
    
    
    twoPoint p2 = part2(), p3 = part3();
    double dp2 = p2.getDist(), dp3 = p3.getDist();
    string p12 = (p2.gp1()).toString(); string p22 = (p2.gp2()).toString();
    string p13 = (p3.gp1()).toString(); string p23 = (p3.gp2()).toString();
    
    fstream results; 
    results.open("results.txt",ios::out);
     
    results << "PART 2 (RECURSIVE) \n";
    results << "Closest Points: " << p12 << " " << p22 << "\n";
    results << setprecision(17) << "Smallest Distance: " << dp2 << "\n";
    results << setprecision(17) << "Time taken for part 2: " << t2.count() << "\n";
    results << "PART 3 (FULL RECURSIVE) \n";
    results << "Closest Points: " << p13 << " " << p23 << "\n";
    results << setprecision(17) << "Smallest Distance: " << dp3 << "\n";
    results << setprecision(17) << "Time taken for part 3: " << t3.count() << "\n";
    
    cout << "PART 2 (RECURSIVE) \n";
    cout << "Closest Points: " << p12 << " " << p22 << "\n";
    cout << setprecision(17) << "Smallest Distance: " << dp2 << "\n";
    cout << setprecision(17) << "Time taken for part 2: " << t2.count() << "\n";
    cout << "PART 3 (FULL RECURSIVE) \n";
    cout << "Closest Points: " << p13 << " " << p23 << "\n";
    cout << setprecision(17) << "Smallest Distance: " << dp3 << "\n";
    cout << setprecision(17) << "Time taken for part 3: " << t3.count() << "\n";
    
    results.close();
    
    return 0;
}
