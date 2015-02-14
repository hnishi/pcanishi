#include"nlib.h"
/*double ave(float *all, int i){         // to determine average
        float total=0;
        int u;
                for(u=0;u<i;u++){
                total+=all[u];
                }
return total/i;
}*/
// #########################################################################
//                           CLASS: Cent_nishi
// #########################################################################
// include nlib.h
/* version 0.0

class Cent_nishi
{
private:
public:
        vector<int> pair1_R,pair2_R;
        vector<float> dist_residue,comx_r,comy_r,comz_r;
        unsigned int total_dist;
        //char label; // for checking dist mode
        pdb_nishi* pdb1;

        Cent_nishi(const char *pdbname);//:pdb_nishi(pdbname); // constructer for input of cod file
        int dist_resi(const char* filename,unsigned float dist);
        ~Cent_nishi();
};
*/
/*
Cent_nishi::Cent_nishi(const char *pdbname){
	pdb1=new pdb_nishi(pdbname);
	cout<<"successfully open "<<pdbname<<endl;
}*/

/*int Cent_nishi::dist_resi(const char* filename,float dist,char* chai1,char* char2,int rnum1,int rnum2){
	float *bufx,*bufy,*bufz;
	bufx=new float[100];bufy=new float[100];bufz=new float[100];
	int jj=0;
        for(int k= 0; k<pdb1->resi_mark[0]; k++ ){
		//strcpy(bfatmn,pdb1->atmn[k].c_str());
                if(! (strncmp( pdb1->atmn[k].c_str(),"N",4 ) == 0  ||
                        strncmp( pdb1->atmn[k].c_str(),"CA",4 ) == 0 ||
                        strncmp( pdb1->atmn[k].c_str(),"C",4 ) == 0  ||
                        strncmp( pdb1->atmn[k].c_str(),"O",4 ) == 0)  ){
                        cout<<"anum "<<pdb1->anum[k]<<", atmn "<<pdb1->atmn[k]<<", resn "<<pdb1->resn[k] <<endl;
                        //cout <<"rnum "<<pdb1->rnum[k]<<endl;
                        //cout << "bufx " << bufx[k-resnum[j-1]] << endl;
                        bufx[jj] = pdb1->coox[k];
                        bufy[jj] = pdb1->cooy[k];
                        bufz[jj] = pdb1->cooz[k];
			cout<<"jj= "<<jj<<", bufz "<<bufz[jj]<<endl;
        	        jj++;
        	}
        }
	cout<<"total jj= "<<jj<<endl;
        comx_r.push_back( ave(bufx,jj) );
        comy_r.push_back( ave(bufy,jj) );
        comz_r.push_back( ave(bufz,jj) );
	cout<<"comx0= "<<comx_r[0]<<", comy0= "<<comy_r[0]<<", comz0= "<<comz_r[0]<<endl;
	for(unsigned int j=1 ; j<pdb1->total_residue ; j++ ){
		jj=0;
		for(int k= pdb1->resi_mark[j-1]+1; k<pdb1->resi_mark[j]; k++ ){
	                if(! (strncmp( pdb1->atmn[k].c_str(),"N",4 ) == 0  ||
        	        strncmp( pdb1->atmn[k].c_str(),"CA",4 ) == 0 ||
                	strncmp( pdb1->atmn[k].c_str(),"C",4 ) == 0  ||
                	strncmp( pdb1->atmn[k].c_str(),"O",4 ) == 0)  ){
                        	//cout<<"anum "<<" "<<pdb1->anum[k]<<endl;
				//cout <<"rnum "<<pdb1->rnum[k]<<endl;
                        	//cout << "bufx " << bufx[k-resnum[j-1]] << endl;
                                bufx[jj] = pdb1->coox[k];
                                bufy[jj] = pdb1->cooy[k];
                                bufz[jj] = pdb1->cooz[k];			
				jj++;
			}
		}
		comx_r.push_back( ave(bufx,jj) );
                comy_r.push_back( ave(bufy,jj) );
                comz_r.push_back( ave(bufz,jj) );
		//cout <<"comx[j]:"<<j<<" " <<comx[j]<<endl;
		//cout <<"res 2 "<< ave(bufx,resnum[j]-resnum[j-1]) << endl;
	}
	total_distR=comz_r.size();
	cout<<"total_distR: "<<total_distR<<", total_residue: "<<pdb1->total_residue<<endl;
	if( total_distR != pdb1->total_residue){
		printf("ERROR: Cent_nishi::dist_resi, disagreement\n");
		cout<<"total_distR: "<<total_distR<<", total_residue: "<<pdb1->total_residue<<endl;
		exit(1);
	}
	delete [] bufx;delete [] bufy;delete [] bufz;


	return 0;
}*/

/*Cent_nishi::~Cent_nishi(){
        delete pdb1;
}*/
