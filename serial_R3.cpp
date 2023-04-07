#include <iostream>
#include <bits/stdc++.h>
#include <fstream>

using namespace std;
const double INF_ = 1e9;
const double epsilon = 0.0000000001;

typedef vector<double> vd;
typedef vector<vd> vvd;
typedef vector<bool> vb;
typedef vector<vb> vvb;
typedef vector<int> vi;
typedef vector<vi> vvi;
typedef pair<int,int> ii;
typedef vector<ii> vii;

#define pb push_back
#define F first
#define S second
#define MP make_pair

void printgrid(vvi &grid){
    for(int i=0;i<grid.size();i++){
        for(int j=0;j<grid[i].size();j++){
            cout<<grid[i][j]<<" ";
        }
        cout<<endl;
    }
    cout<<endl;
}

void readfile_and_make_grid(vvi &grid){
    //std::ifstream infile("/Users/avirajbevli/Desktop/Sem10/HPPP/cuda_term project/data/cdnf43w/heights.txt");
    cout<<"Reading from input file..."<<endl;
    std::ifstream infile("heights.txt");
    if(!infile.is_open()) {
        std::cerr << "Error: Could not open read file" << std::endl;
        return  ;
    }

    string line;
    while(std::getline(infile, line)) {
        vi temp;
        string str; 
        int i=0;
        while(i<line.size()){
            if(line[i]==','){
                temp.pb(stoi(str));
                str="";
                i++;
            }
            else{
                str+=line[i];
            }
            i++;
        }
        grid.pb(temp);
    }
    infile.close();
    cout<<"Done reading"<<endl;
}

// function for line generation
vii plotPixel(int x1, int y1, int x2, int y2, int dx, int dy, int decide){
    // pk is initial decision making parameter
    // Note:x1&y1,x2&y2, dx&dy values are interchanged
    // and passed in plotPixel function so
    // it can handle both cases when m>1 & m<1
    int pk = 2 * dy - dx;
    vii line;
    for (int i = 0; i <= dx; i++) {
        // cout << x1 << "," << y1 << endl;
        // checking either to decrement or increment the
        // value if we have to plot from (0,100) to (100,0)
        x1 < x2 ? x1++ : x1--;
        if (pk < 0) {
            // decision value will decide to plot
            // either  x1 or y1 in x's position
            if (decide == 0) {
                // putpixel(x1, y1, RED);
                line.pb(MP(x1,y1));
                pk = pk + 2 * dy;
            }
            else {
                //(y1,x1) is passed in xt
                // putpixel(y1, x1, YELLOW);
                line.pb(MP(y1,x1));
                pk = pk + 2 * dy;
            }
        }
        else {
            y1 < y2 ? y1++ : y1--;
            if (decide == 0) {
                line.pb(MP(x1,y1));
                // putpixel(x1, y1, RED);
            }
            else {
                line.pb(MP(y1,x1));
                //  putpixel(y1, x1, YELLOW);
            }
            pk = pk + 2 * dy - 2 * dx;
        }
    }
    return line;
}

vii bresenham(int x1, int y1, int x2, int y2){
    int dx = abs(x2 - x1);
    int dy = abs(y2 - y1);
    // If slope is less than one
    if (dx > dy) {
        // passing argument as 0 to plot(x,y)
        return plotPixel(x1, y1, x2, y2, dx, dy, 0);
    }
 
    // if slope is greater than or equal to 1
    else {
        // passing argument as 1 to plot (y,x)
        return plotPixel(y1, x1, y2, x2, dy, dx, 1);
    }
}

double calc_dist(ii p1, ii p2){
    double x1=p1.F; double y1=p1.S;
    double x2=p2.F; double y2=p2.S;
    double dx = x2-x1; double dy = y2-y1;
    return sqrt((dx*dx) + (dy*dy));
}

bool isvalid(ii p, int n, int m, double dist_from_origin, double range){
    if(p.F < 0 || p.S < 0) return 0;
    if(p.F >= n || p.S >= m) return 0;
    if(dist_from_origin > range) return 0;
    return 1;
}

vii find_points_in_range(ii observer, vvi &grid, double range, vvi& viewshed){
    int n = grid.size();
    int m = grid[0].size();

    // bfs
    int dx[4] = {0,1,0,-1};
    int dy[4] = {1,0,-1,0};
    queue<ii> q;
    q.push(observer);
    
    vii pts_inrange;
    //cout<<"Points in range: "<<endl;
    while(!q.empty()){
        ii curr = q.front();
        pts_inrange.pb(curr);
        viewshed[curr.F][curr.S]=1;
        q.pop();
        //cout<<curr.F<<","<<curr.S<<" ";
        for(int i=0;i<4;i++){
            ii next = MP(curr.F + dx[i], curr.S + dy[i]);
            double dist_from_origin = calc_dist(observer, next);
            if(isvalid(next, n, m, dist_from_origin, range)){
                if(viewshed[next.F][next.S]==-1){
                    q.push(next);
                }
            }
        }
    }
    return pts_inrange;
}

double find_slope(vvi& grid, ii observer, ii curr, int observer_ht){
    double dist = calc_dist(observer, curr);
    double dh = double(grid[curr.F][curr.S]) - (double(grid[observer.F][observer.S]) + (double)observer_ht);
    double m_curr = dh/dist;
    return m_curr;
}

// -1: out of range, 0: not visible, 1: visible
void calculate_viewshed_and_writefile_r3(vvi &grid){
    // parameters to the algo
    ii observer = MP(1000,3000);
    int observer_ht = 10; // height of observer from the ground at <observer_x, observer_y>

    int n = grid.size();
    int m = grid[0].size();

    vvi viewshed(n,vi(m,-1));

    // double range = 500;
    // vii pts_inrange = find_points_in_range(observer,grid,range,viewshed);
    // cout<<"Number of points in range: "<<pts_inrange.size()<<endl;

    // for(int i=0;i<pts_inrange.size();i++){
    //     ii curr = pts_inrange[i];
    //     vii line = bresenham(observer.F, observer.S, curr.F, curr.S);
    //     // do something with the obtained line
    //     double m_max = -INF_;
    //     ii observer = line[0];
    //     for(int i=1;i<line.size();i++){
    //         double dist = calc_dist(observer, curr);
    //         double dh = grid[curr.F][curr.S] - grid[observer.F][observer.S];
    //         double m_curr = dist/dh;
    //         if(m_curr >= m_max){
    //             m_max = m_curr;
    //         }
    //         else{
    //             viewshed[curr.F][curr.S]=0;
    //         }
    //     }
    // }
    int range = 400;
    int li = max(0, (observer.F) - range);
    int ri = min(n,(observer.F) + range);
    int lj = max(0, (observer.S) - range);
    int rj = min(m,(observer.S) + range);
    cout<<"li:"<<li<<",ri:"<<ri<<", lj:"<<lj<<",rj:"<<rj<<endl;
    for(int i = li; i < ri; i++){
        for(int j = lj; j < rj; j++){
            ii curr = MP(i,j);
            viewshed[curr.F][curr.S]=1;
        } 
    }
    for(int i = li; i < ri; i++){
        for(int j = lj; j < rj; j++){
            ii target = MP(i,j);
            double m_target = find_slope(grid, observer, target, observer_ht);
            vii line = bresenham(observer.F, observer.S, target.F, target.S);
            // do something with the obtained line

            viewshed[target.F][target.S]=1; // target is visible by default
            double m_max = -INF_;
            int linesz = line.size();
            for(int lid=1;lid<linesz-1;lid++){
                // ii curr = line[lid];
                double m_curr = find_slope(grid, observer, line[lid], observer_ht);
                m_max = max(m_max, m_curr);
                // if(m_curr >= m_max){
                //     m_max = m_curr;
                // }
                // else{
                //     viewshed[curr.F][curr.S]=0;
                // }
            }
            if(m_max > m_target + epsilon){
                viewshed[target.F][target.S]=0;
            }
        } 
    }

    cout<<"Writing to output viewshed file ...."<<endl;
    std::ofstream outFile("viewshed_r3.txt");   
    if (!outFile.is_open()) {
        std::cout << "Unable to open write file" << std::endl;
        return;
    }
    for(int i=0;i<viewshed.size();i++){
        for(int j=0;j<viewshed[i].size();j++){
            outFile<<viewshed[i][j]<<" ";
        }
        outFile<<endl;
    }
    // Close the file
    outFile.close();
    cout<<"Done writing"<<endl;
}

int main(){
    vvi grid;   
    readfile_and_make_grid(grid);
    calculate_viewshed_and_writefile_r3(grid);

    // vii testline = bresenham(10,5,3,20);
    // for(int i=0;i<testline.size();i++){
    //     cout<<testline[i].F<<","<<testline[i].S<<" ";
    // }cout<<endl;

    // vii testline2 = bresenham(5,10,20,3);
    // for(int i=0;i<testline2.size();i++){
    //     cout<<testline2[i].F<<","<<testline2[i].S<<" ";
    // }cout<<endl;

    // ii observer = MP(10,5);
    // double m_max = -INF_;
    // for(int i=1;i<testline.size();i++){
    //     double dist = calc_dist(observer, curr);
    //     double dh = double(grid[curr.F][curr.S]) - double(grid[observer.F][observer.S]);
    //     double m_curr = dh/dist;
    //     if(m_curr >= m_max){
    //         m_max = m_curr;
    //         viewshed[curr.F][curr.S]=1;
    //     }
    //     else{
    //         viewshed[curr.F][curr.S]=0;
    //     }
    // }

    return 0;
}
