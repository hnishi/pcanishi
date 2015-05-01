#include<iostream>
#include<vector> 
#include<string> // for vector<string>
#include<stdlib.h> // for exit(1); error end
#include<stdio.h>
#include<string.h>
#include<fstream>
#include<cmath> // for sqrt(); square root
//#include<iomanip> //setw()
//#include"math_nishi.h"

#define DEBUG 0	// for debugging, instead of comment out
//	#if DEBUG == 1
//	........................
//	#endif

//#define SCAN_FORMAT "%6s%5d%4s%4s%1s%4d%7f%7f%7f%5f%5f%s" // no missing
#define WRITE_FORMAT_1 "%-6s%5d%5s %-4s%c %-5d%10.3f%8.3f%8.3f%6.2f%6.2f%12s\n"
#define WRITE_FORMAT_2 "%-6s%5d  %-3s %-4s%c %-5d%10.3f%8.3f%8.3f%6.2f%6.2f%12s\n"
//#define SCAN_FORMAT "%6s%5d%4s%4s%5d%7f%7f%7f%5f%5f" // missing chain element
//#define WRITE_FORMAT_1 "%-6s%5d%5s %-4s%5d%12.3f%8.3f%8.3f%6.2f%6.2f\n"
//#define WRITE_FORMAT_2 "%-6s%5d  %-3s %-4s%5d%12.3f%8.3f%8.3f%6.2f%6.2f\n"

using namespace std;

inline double ave(float *all, int i); // for average

// #########################################################################
//                             pdb_nishi
// #########################################################################
// pdbnishi.cpp
#ifndef _INCLUDE_PDB_NISHI_
#define _INCLUDE_PDB_NISHI_
class pdb_nishi{
protected:
        int bfanum,bfrnum;
        float bfcoox,bfcooy,bfcooz,bfoccu,bftemf;
        char bfreco[7],bfatmn[5],bfresn[5],bfchai[2],bfelem[5];
	//char bfreco[9],bfatmn[9],bfresn[9],bfchai[9],bfelem[9];
public:
        vector<int> anum,rnum; // pdb elements
        vector<double> coox,cooy,cooz,occu,temf; // pdb elements
        vector<string> reco,atmn,resn,chai,elem; // pdb elements
	vector<int> resi_mark; // final internal num. of the residue

	string pdb_name;

	unsigned int total_atom;//anum.size()
	unsigned int total_residue;//resi_mark.size()

        pdb_nishi(); // constructer 1 for error
        pdb_nishi(const char *pdbname); // constructer for input of pdb
        int disp_line(int n); // n is internal number
	int write_pdb(const char* filename, char mode); // put output filename with w or a option
	int write_pdb(const char* filename); // put output filename
	int search_n(char a,int aa); // return internal num. of chain ID "a" and residue num. "aa"; if return -1, then could not find it
        int search_n_end(char a,int aa); // search the final internal num. of chain ID "a" and residue num. "aa"
	int fix_step(const char *filename,float fxcell,float fycell,float fzcell);

	//please use the following variables after function center_r()
	vector<double> comx_r,comy_r,comz_r;// center of mass of side-chains
	unsigned int total_com_r;//num. of comz_r
	int center_r();
};
#endif
/* ### MEMO ###
	cout<<"TOTAL ATOM = "<<pdb1.rnum.size()<<endl;
	cout<<"TOTAL RESIDUE = "<<pdb1.resi_mark.size()<<endl;
	cout<<"final internal num. = "<<pdb1.resi_mark[pdb1.resi_mark.size()-1]<<endl;
*/
// #########################################################################
//                             tra_nishi
// #########################################################################
// tranishi.cpp
#ifndef _INCLUDE_TRA_NISHI_
#define _INCLUDE_TRA_NISHI_
class tra_nishi 
{
protected:
//private:
	int bufint;
	float buffloat;
        int bfanum,bfrnum;
        float bfcoox,bfcooy,bfcooz,bfoccu,bftemf;
        //char bfreco[7],bfatmn[6],bfresn[5],bfchai[2],bfelem[5];
	char bfreco[9],bfatmn[9],bfresn[9],/*bfchai[9],*/bfelem[9];
public:
	vector<int> loopnum,num15svw,num15hyd;
	vector<double> sitime,cputime,totalE,kineticE,temp,potent,rmsf,rmsd;
	vector<double> cordx,cordy,cordz,length_x,length_y,length_z;
	unsigned int total_step, total_sel; // total_* >= 1
	pdb_nishi* pdb1;
	string atom_sel;
	string cod_name, pdb_name;

	void constructor(const char *codname, const char *pdbname,int stride,string atomsel);
	tra_nishi(const char *codname, const char *pdbname); // default constructer
	tra_nishi(const char *codname, const char *pdbname,int stride); // constructer for input of cod file
	tra_nishi(const char *codname, const char *pdbname,int stride,string atomsel); // constructer for input of cod file
	tra_nishi(const char *codname, const char *pdbname,string atomsel); // constructer for input of cod file

	int disp_line(int step); // display info. of step without coordinates
	int write_cod(const char* filename, int stride);
	int write_cod(const char* filename);
	int write_step(const char* filename, int n); // write pdb at n step
	int fix_step(const char *filename, int n,float fxcell,float fycell,float fzcell);
	int fix_cod(float fxcell,float fycell,float fzcell);
	int fix_cod_npt();
	~tra_nishi();
};
#endif
///////////////////////////////////////////////////////////////
// #########################################################################
//        inheritance of tra_nishi  ---> Fix_cod
// #########################################################################
#ifndef _INCLUDE_FIX_COD_
#define _INCLUDE_FIX_COD_
class Fix_cod: public tra_nishi // inherit tra_nishi
{
private:
public:
	float fxcell, fycell, fzcell;

        Fix_cod(const char *codname, const char *pdbname,int stride,float a,float b,float c)
        : tra_nishi(codname,pdbname,stride) {
		fxcell=a; fycell=b; fzcell=c;
        	float rdiff;
        	int imove;
        	for(unsigned int n=0;n<total_step;n++){
                	//fprintf(fout,"MODEL %d\n",n+1);
       		for(unsigned int j = n * pdb1->total_atom +1; j < (n+1) * pdb1->total_atom;j++){
                	rdiff = cordx[j]-cordx[j-1];
        	        if(rdiff>=0) imove = int(rdiff/fxcell+0.5);
                	else imove = int(rdiff/fxcell-0.5);
	                cordx[j] = cordx[j] - imove*fxcell;
        	        rdiff = cordy[j]-cordy[j-1];
	                if(rdiff>=0) imove = int(rdiff/fycell+0.5);
        	        else imove = int(rdiff/fycell-0.5);
                	cordy[j] = cordy[j] - imove*fycell;
	                rdiff = cordz[j]-cordz[j-1];
                	if(rdiff>=0) imove = int(rdiff/fzcell+0.5);
	                else imove = int(rdiff/fzcell-0.5);
        	        cordz[j] = cordz[j] - imove*fzcell;
	        }
        	}
        //write_cod("zzz.pdb",1);
	} // constructer
};
#endif
///////////////////////////////////////////////////////////////
// #########################################################################
//        inheritance of tra_nishi  ---> Tra_ana
// #########################################################################
#ifndef _INCLUDE_TRA_ANA_
#define _INCLUDE_TRA_ANA_
class Tra_ana: public tra_nishi // inherit tra_nishi
{
private:
public:
        //struct inp_se{ int sta_rnum; int end_rnum; string sta_chai; string end_chai; }
        //struct inp_se input;
        int intra_sta,intra_end,sta_rnum,end_rnum;
        char sta_chai,end_chai;

        Tra_ana(const char *codname, const char *pdbname,int stride,char a,char b,int aa,int bb)
        : tra_nishi(codname,pdbname,stride) {
                sta_rnum = aa; end_rnum = bb;
                sta_chai = a ; end_chai = b ;
                intra_sta = pdb1->search_n(a,aa);
                intra_end = pdb1->search_n(b,bb)-1;
        } // constructer
};
#endif
// #########################################################################
//                         pdb_nishi ---> centnishi
// #########################################################################
// centnishi.cpp

#ifndef _INCLUDE_AVE_
#define _INCLUDE_AVE_
        inline double ave(float *all, int i){         // to determine averagei
                float total=0;
                int u;
                for(u=0;u<i;u++){
                        total+=all[u];
                }
                return total/i;
        }
#endif

#ifndef _INCLUDE_CENT_NISHI_
#define _INCLUDE_CENT_NISHI_
class Cent_nishi : public pdb_nishi // inherit pdb_nishi 
{
private:
public:
        vector<int> pair1_R,pair2_R;
        vector<double> dist_residue,comx_r,comy_r,comz_r;
        unsigned int total_com_r;
	//char label; // for checking dist mode
/*
	inline double ave(float *all, int i){         // to determine averagei
        	float total=0;
        	int u;
                for(u=0;u<i;u++){
	                total+=all[u];
                }
		return total/i;
	}
*/
        Cent_nishi(const char *pdbname):pdb_nishi(pdbname) // constructor
	{
        float *bufx,*bufy,*bufz;
        bufx=new float[100];bufy=new float[100];bufz=new float[100];
        int jj=0;
        for(int k= 0; k<resi_mark[0]; k++ ){
                //strcpy(bfatmn,atmn[k].c_str());
                if(! (strncmp( atmn[k].c_str(),"N",4 ) == 0  ||
                        strncmp( atmn[k].c_str(),"CA",4 ) == 0 ||
                        strncmp( atmn[k].c_str(),"C",4 ) == 0  ||
                        strncmp( atmn[k].c_str(),"O",4 ) == 0)  ){
                        //cout<<"anum "<<anum[k]<<", atmn "<<atmn[k]<<", resn "<<resn[k] <<endl;
                        //cout <<"rnum "<<rnum[k]<<endl;
                        //cout << "bufx " << bufx[k-resnum[j-1]] << endl;
                        bufx[jj] = coox[k];
                        bufy[jj] = cooy[k];
                        bufz[jj] = cooz[k];
                        //cout<<"jj= "<<jj<<", bufz "<<bufz[jj]<<endl;
                        jj++;
                }
        }
        //cout<<"total jj= "<<jj<<endl;
        comx_r.push_back( ave(bufx,jj) );
        comy_r.push_back( ave(bufy,jj) );
        comz_r.push_back( ave(bufz,jj) );
        //cout<<"comx0= "<<comx_r[0]<<", comy0= "<<comy_r[0]<<", comz0= "<<comz_r[0]<<endl;
        for(unsigned int j=1 ; j<total_residue ; j++ ){
                jj=0;
                for(int k= resi_mark[j-1]+1; k<resi_mark[j]; k++ ){
                        if(! (strncmp( atmn[k].c_str(),"N",4 ) == 0  ||
                        strncmp( atmn[k].c_str(),"CA",4 ) == 0 ||
                        strncmp( atmn[k].c_str(),"C",4 ) == 0  ||
                        strncmp( atmn[k].c_str(),"O",4 ) == 0)  ){
                                //cout<<"anum "<<" "<<anum[k]<<endl;
                                //cout <<"rnum "<<rnum[k]<<endl;
                                //cout << "bufx " << bufx[k-resnum[j-1]] << endl;
                                bufx[jj] = coox[k];
                                bufy[jj] = cooy[k];
                                bufz[jj] = cooz[k];
                                jj++;
                        }
                }
                comx_r.push_back( ave(bufx,jj) );
                comy_r.push_back( ave(bufy,jj) );
                comz_r.push_back( ave(bufz,jj) );
                //cout <<"comx[j]:"<<j<<" " <<comx[j]<<endl;
                //cout <<"res 2 "<< ave(bufx,resnum[j]-resnum[j-1]) << endl;
        }
        total_com_r=comz_r.size();
        //cout<<"total_com_r: "<<total_com_r<<", total_residue: "<<total_residue<<endl;
        if( total_com_r != total_residue){
                printf("ERROR: dist_nishi::dist_resi, disagreement\n");
                cout<<"total_com_r: "<<total_com_r<<", total_residue: "<<total_residue<<endl;
                exit(1);
        }
	} // constructer
        //int dist_resi(const char* filename, float dist);
        //~Cent_nishi();
};
#endif

// #########################################################################
//                         Inp_nishi
// #########################################################################
// inpnishi.cpp
#ifndef _INCLUDE_INPUT_PARAMETERS_
#define _INCLUDE_INPUT_PARAMETERS_
class Inp_nishi{
private:
public:
   string filename;

   Inp_nishi( const char *inputname );
      string read( string );
   };
#endif
			    

// #########################################################################
//                select_atom
// atomsel:
//	all: all atoms
//	protein: without resn WAT, CIP and CIM
//	heavy: protein without hydrogens
//	mainchain: only atom name CA, N, C and O
//	ca: only atom name CA
// #########################################################################
// select_atom() in tranishi.cpp
int select_atom( pdb_nishi &pdb1, vector<double> &vec, string &atomsel, int i );
int select_atom( pdb_nishi &pdb1, double x, double y, double z, vector<double> &vec, string &atomsel, int i );

// search_sel in tranishi.cpp
int search_sel( pdb_nishi &pdb1, string chai, int resn, string atmn, string atomsel);
