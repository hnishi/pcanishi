#include"nlib.h"
#include<cmath>
#include<vector>
using namespace std;

/*#ifndef _INCLUDE_AVE_
#define _INCLUDE_AVE_
        inline double ave(float *all, int i){         // to determine averagei
                float total=0;
                int u;
                for(u=0;u<i;u++){
                        total+=all[u];
                }
                return total/i;
        }
#endif*/

#ifndef _INCLUDE_RMSD_
#define _INCLUDE_RMSD_
//inline double rmsd(vector<double> &ax,vector<double> &ay,vector<double> &az,
double rmsd(vector<double> &ax,vector<double> &ay,vector<double> &az,
		vector<double> &bx,vector<double> &by,vector<double> &bz); // to calculate R.M.S.D.
#endif

#ifndef _INCLUDE_RMSD_2_
#define _INCLUDE_RMSD_2_
double rmsd2(vector<double> &a,vector<double> &b); // to calculate R.M.S.D.
#endif
