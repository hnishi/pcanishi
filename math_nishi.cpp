#include"math_nishi.h"

double rmsd(vector<double> &ax,vector<double> &ay,vector<double> &az,
		vector<double> &bx,vector<double> &by,vector<double> &bz) // to calculate R.M.S.D.
{
	double rmsd,tot_x=0,tot_y=0,tot_z=0;
	//double sum=0;
	if(ax.size() != bx.size()){
		cerr<<"ERROR: both structures should have the same num. of atoms\n";
		return -1;
	}
	for(unsigned int u=0;u<ax.size();u++){
		tot_x += ( ax[u]- bx[u])*( ax[u]- bx[u]);
		tot_y += ( ay[u]- by[u])*( ay[u]- by[u]);
		tot_z += ( az[u]- bz[u])*( az[u]- bz[u]);
	}
	//sum = tot_x + tot_y + tot_z;
	rmsd=sqrt((tot_x + tot_y + tot_z)/ ax.size());
        return rmsd;
}
double rmsd2(vector<double> &a,vector<double> &b) // to calculate R.M.S.D.
{
	double rmsd,tot_x=0,tot_y=0,tot_z=0;
	if(a.size() != b.size()){
		cerr<<"ERROR: both structures should have the same num. of atoms\n";
		return -1;
	}
	double sum = tot_x + tot_y + tot_z;
	for(unsigned int u=0;u<a.size();u=u+3){
		//cout<<"ax ay az = "<<a[u]<<" "<<a[u+1]<<" "<<a[u+2]<<endl;
		//cout<<"bx by bz = "<<b[u]<<" "<<b[u+1]<<" "<<b[u+2]<<endl;
		tot_x += ( a[u] - b[u])*( a[u] - b[u]);
		tot_y += ( a[u+1] - b[u+1])*( a[u+1] - b[u+1]);
		tot_z += ( a[u+2] - b[u+2])*( a[u+2] - b[u+2]);
		sum = tot_x + tot_y + tot_z;
		if( sum < 0) cerr<<"!!! sum = "<<sum<<endl;
	}
	//sum = tot_x + tot_y + tot_z;
	//cout<<"!!!! sum = "<<sum<<endl;
	rmsd=sqrt((tot_x + tot_y + tot_z)/ (a.size()/3));
        return rmsd;
}
