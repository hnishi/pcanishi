#include"nlib.h"
#include"math_nishi.h"
#include<Eigen/Dense>  // Eigen/Core, Geometry, LU, Cholesky, SVD, QR and Eigenvalues
//#include<sstream>
using namespace Eigen;

// user define
//#define INSTART_CHAIN ' ' 
//#define INSTART_NUM   1       
//#define INEND_CHAIN   ' '  
//#define INEND_NUM     27 
#define INSTART_CHAIN  'B'
#define INEND_CHAIN    'B'  
#define INSTART_NUM   210       
#define INEND_NUM     222 
#define TEMPERATURE 300.0
#define BIN_DE 5 
#define BIN_EMIN1 -200

// constant values
#define GAS_CONST 8.31451  // gas constant, R ( joule/mol*k )
//#define BOLTZMAN_CONST 1.380658e-23  // Boltzman constant, kb ( joule/k )
#define BOLTZMAN_CONST 8.3144621  // Boltzman constant, kb ( joule/k/mol )
// kb = R / Na, where Na is Avogadro constant (6.02214129e-23)
// printf("Boltzman constant kb = %e \n", kb);
#define JOULE_CALORIE 4.184 // joule-calorie conversion unit

#define DEBUG_PCA 1 // for debugging, instead of comment out
//      #if DEBUG_PCA == 1
//            ........................
//      #endif


int pca_prctc(char *codname,char *pdbname,int stride){
// if stride is set to -1, then
if(stride == -1){
// END 
  return -1;
}

// ############# trajectory  ########################################
/* (1) load trajectory by tra_nishi
*/
   cout<<endl<<"REPORT> (1) load trajectory \n";
   cout<<"loading\n"<<codname<<"        "<<pdbname<<endl;
	tra_nishi* tra1;
	tra1 = new tra_nishi(codname, pdbname, stride);
	cout<<"TOTAL FRAME = "<<tra1->total_step<<endl;
	cout<<"TOTAL ATOM = "<<tra1->pdb1->total_atom<<endl;
	unsigned int frame = tra1->total_step;
/* (2) setting of region
  from intra_start to intra_end (internal number of atom)
*/
   cout<<endl<<"REPORT> (2) specify the region \n";
   cout<<"test\nINSTART_CHAIN:\""<<INSTART_CHAIN<<"\", INSTART_NUM:\""<<INSTART_NUM<<"\"\n";
   cout<<"test\nINEND_CHAIN:\""<<INEND_CHAIN<<"\", INEND_NUM:\""<<INEND_NUM<<"\"\n";
	int intra_start, intra_end;
	intra_start = tra1->pdb1->search_n(INSTART_CHAIN,INSTART_NUM);
	//intra_end = tra1->pdb1->search_n(INEND_CHAIN,INEND_NUM);  // for final residue
	if( tra1->pdb1->rnum[tra1->pdb1->total_atom - 1] == INEND_NUM ){
		intra_end = tra1->pdb1->search_n(INEND_CHAIN,INEND_NUM );
	}else{
		intra_end = tra1->pdb1->search_n(INEND_CHAIN,INEND_NUM + 1) -1;
	}
	cout<<"intra_start = "<<intra_start<<endl;
	cout<<"intra_end = "<<intra_end<<endl;

	unsigned int dim_0 = intra_end - intra_start +1;
	//cout<<"!!!!!!!!!!!! dim_0 = "<<dim_0<<endl;

	//tra1.pdb1->disp_line(11860);

	tra1->pdb1->disp_line(intra_start-1);
	tra1->pdb1->disp_line(intra_start);
	cout<<"...\n";
        tra1->pdb1->disp_line(intra_end);
        tra1->pdb1->disp_line(intra_end+1);

        tra1->pdb1->write_pdb("zzz.pdb"); cout<<"output pdb file (zzz.pdb)\n";

/* (3) create vector which includes only atoms within h3-loop to create vecors for variance-covariance matrix
  method : minimum "heavy atom - heavy atom" distance between residues except for ones of |i-j| <= 2 
*/
   cout<<endl<<"REPORT> (3) PCA calculation starts \n";
/* ********
 * 3-1   create vector Q = (q1x,q1y,q1z,q2x,...,q20z) 
 *     
 * */
   cout<<" 3-1  create vector Q \n";
	//int * Q; Q = new double [dim_0 * 3][frame] //
	double Q[dim_0 * 3][frame];
	for(unsigned int n=0;n<frame;n++){
	        for(unsigned int i=0;i<dim_0*3;i++){
			Q[i][n] = 0;
		}
	}
	//cout<<"!!!! frame = "<<frame<<endl;
	unsigned int ii = 0;
	for(unsigned int n=0;n<frame;n++){
		ii = 0;
	        for(int i=intra_start;i<=intra_end;i++){
			//if( tra1->pdb1->elem[i] == "H" )continue;
			if( tra1->pdb1->atmn[i] != "CA" )continue;
			Q[ii][n] = tra1->cordx[n*tra1->pdb1->total_atom+i];
                	Q[ii+1][n] = tra1->cordy[n*tra1->pdb1->total_atom+i];
	                Q[ii+2][n] = tra1->cordz[n*tra1->pdb1->total_atom+i];
        		ii = ii + 3;
		}
	}

	//delete tra1;

	unsigned int dim_Q = ii;
	cout<<"dimensionality of Q is "<<dim_Q<<endl;
	cout<<"total num. of Q is "<<frame<<endl;

/* *********
 * 3-2   create variance-covariance matrix
 *
 * */
   cout<<" 3-2  create variance-covariance matrix \n";
	double average_Q[dim_Q], average_QQ[ dim_Q * dim_Q ];  // initialize 2d-array
        for(unsigned int n=0;n<frame;n++){
                for(unsigned int i=0;i<dim_Q;i++){
                	average_Q[i] = 0;
			for(unsigned int j=0;j<dim_Q;j++){
				average_QQ[ j + i * dim_Q ] = 0;  // 
			}
        	}
	}
	for(unsigned int n=0;n<frame;n++){
                for(unsigned int i=0;i<dim_Q;i++){
                	average_Q[i] += Q[i][n];
			//cerr<<average_Q[i]<<" ";
			//cout<<"average_Q["<<i<<"] = "<<average_Q[i]<<endl;
			for(unsigned int j=0;j<dim_Q;j++){
				average_QQ[ j + i * dim_Q ] += Q[j][n] * Q[i][n];  // 
			}
		}
        }
        for(unsigned int i=0;i<dim_Q;i++){
                average_Q[i] = average_Q[i] / frame;
		//cerr<<"average_Q["<<i<<"] = "<<average_Q[i]<<endl;
        }
        for(unsigned int i=0;i<dim_Q*dim_Q;i++){
                average_QQ[i] = average_QQ[i] / frame;
		//cout<<"average_QQ["<<i<<"] = "<<average_QQ[i]<<endl;
        }
	//delete[] Q;

        MatrixXd VCV(dim_Q,dim_Q);// Eigen/Dense
        for(unsigned int n=0;n<dim_Q;n++){
                for(unsigned int i=0;i<dim_Q;i++){
                	VCV(i,n) = average_QQ[ i + n * dim_Q ] - average_Q[i] * average_Q[n];
                	//cout<<"\nVCV(i,n) = "<<average_QQ[ i + n * dim_Q ]<<" - "<<average_Q[i] * average_Q[n];
			//cout<<" "<<VCV(i,n);
			//cerr<<VCV(i,n)<<" ";
		}
        	//cout<<endl;
	}
	//cerr<<VCV<<endl;
        for(unsigned int n=0;n<dim_Q;n++){
                for(unsigned int i=0;i<dim_Q;i++){
			//cout<<"DEBUG: VCV("<<i<<","<<n<<") = "<<VCV(i,n)<<", VCV("<<n<<","<<i<<") = "<<VCV(n,i)<<endl;
			if( VCV(i,n) != VCV(n,i) )cout<<"ERROR: !!!!!!!!!!!! VCV is not symmetric !!!!!!!!!!!!!!\n";
                }
        }


//#if DEBUG_PCA == 1
//#endif

/* **********
 * 3-3   calcuate eigen-value and eigen-vector
 *
 * */
   cout<<" 3-3  calculate eigen-value \n";
	SelfAdjointEigenSolver<MatrixXd> es(VCV);
	//cout<<"\neigenvalues in ascending order\n"<<es.eigenvalues()<<endl<<endl;
	//cout<<"!!!!\n"<<es.eigenvalues().rows()<<", dim_Q = "<<dim_Q<<endl;
	//cout<<"\neigenvectors\n"<<es.eigenvectors()<<endl<<endl;
	cout<<"The maximum eigenvalue of VCV\n"<<es.eigenvalues()[dim_Q - 1]<<endl;
	//cout<<"The maximum eigenvalue of VCV\n"<<es.eigenvalues().row(dim_Q - 1)<<endl;
	cout<<"The eigenvector of the maximum eigenvalue\n"<<es.eigenvectors().col(dim_Q-1).transpose()<<endl;
	cout<<"The second eigenvalue of VCV\n"<<es.eigenvalues()[dim_Q-2]<<endl;
	cout<<"The eigenvector of the second eigenvalue\n"<<es.eigenvectors().col(dim_Q-2).transpose()<<endl;

        ofstream ofs_eigen;
        ofs_eigen.open("out_eigen.txt");
        ofs_eigen<<"> Eigen values\n"<<es.eigenvalues()<<endl;
 
   for(unsigned int n=0;n<dim_Q;n++){
      ofs_eigen<<"> Eigen vector no"<<n<<"\n"<<es.eigenvectors().col(n).transpose()<<endl;
   }

        ofs_eigen.close(); 
	cout<<"ouput out_eigen.txt\n";


/* **********
 *  3-5  calculate the coordinates of the ensemble along the PCA axes
 *
 * */
   cout<<" 3-5  calculate components along principle axis \n";
	double c1[frame], c2[frame];

	vector<double> buf_vec;

        for(unsigned int i=0;i<dim_Q;i++){
	        buf_vec.push_back(average_Q[i]);
        }
	VectorXd Q2=Map<VectorXd>(&buf_vec[0],buf_vec.size());

        for(unsigned int n=0;n<frame;n++){
		buf_vec.clear();
                for(unsigned int i=0;i<dim_Q;i++){
			buf_vec.push_back(Q[i][n]);
		}
		VectorXd Q1=Map<VectorXd>(&buf_vec[0],buf_vec.size());
        	c1[n] = es.eigenvectors().col(dim_Q-1).transpose() * ( Q1 - Q2 );
		//cout<<"c1["<<n<<"] = "<<c1[n]<<endl;
        	c2[n] = es.eigenvectors().col(dim_Q-2).transpose() * ( Q1 - Q2 );
		//cout<<"c1["<<n<<"] = "<<c1[n]<<endl;
	}	

/* (4) output (write) c1 vs c2 
*/
	ofstream ofs;
	ofs.open("c1vsc2.dat");

        for(unsigned int n=0;n<frame;n++){
        	ofs<<c1[n]<<"   "<<c2[n]<<endl;
	}

	ofs.close();
   cout<<"output c1vsc2.dat (2-D graph) \n";

/* (5) write pdb of max and min 
 *    calculate min max along first axis
 * */
   cout<<endl<<"REPORT> (5) find min max along the axis \n";
   double max_c = -999999, min_c = 999999;
   int max_c_n = 0, min_c_n = 0;
   for(unsigned int n=0;n<frame;n++){
      if( c1[n] >= max_c ){
         max_c = c1[n]; max_c_n = n;
      }
      if( c1[n] <= min_c ){
         min_c = c1[n]; min_c_n = n;
      }
   }
   cout<<"max_c = "<<max_c<<",  max_c_n = "<<max_c_n<<endl;
   cout<<"min_c = "<<min_c<<",  min_c_n = "<<min_c_n<<endl;
   tra1->write_step("max_c.pdb", max_c_n);  //write pdb at step n+1  
   tra1->write_step("min_c.pdb", min_c_n);  //write pdb at step n+1
   //delete tra1;
   cout<<"output min_c.pdb and max_c.pdb \n";

/* (6) calculate contribution ratio of c1 and c2
 *     
 * */
   cout<<endl<<"REPORT> (6) calculate contribution ratio of c1vsc2 \n";
   double sum_lambda = 0;
   for(unsigned int i=0;i<dim_Q;i++){
      sum_lambda += es.eigenvalues()[i];
   }
   double cntr_rt = (es.eigenvalues()[dim_Q - 1] + es.eigenvalues()[dim_Q - 2]) / sum_lambda;
   cout<<"contribution ratio = "<<cntr_rt<<endl;

/* (7) PMF calculation
 *
 * */
   cout<<endl<<"REPORT> (7) PMF calculation \n";
   //int num_bin = 1000;
   //double length_bin = (max_c - min_c) / num_bin;
   double length_bin = BIN_DE;
   int emax, emin;
   emax = (int)max_c + 3*length_bin;
   emin = (int)min_c - 3*length_bin; cout<<"emax, emin = "<<emax<<", "<<emin<<endl;
   cout<<"length_bin = "<<length_bin<<endl;
   int num_bin = ( emax - emin ) / length_bin + 1; cout<<"num_bin = "<<num_bin<<endl;
   double pmf[num_bin][num_bin];
   for(int i=0;i<num_bin;i++){  //initialize array
      for(int j=0;j<num_bin;j++){
         pmf[j][i] = 0;
      }
   }
   int count_pmf = 0;
   for(unsigned int n=0;n<frame;n++){ //count
      for(int i=0;i<num_bin;i++){
         for(int j=0;j<num_bin;j++){
            if(c1[n] > emin + length_bin * j 
 	    && c1[n] <= emin + length_bin * (j + 1) 
	    && c2[n] > emin + length_bin * i
	    && c2[n] <= emin + length_bin * (i + 1)     ){
	       pmf[j][i] ++ ;
	       count_pmf ++ ;
	       //goto NEXT_PMF;
	    }
         }
      }

//NEXT_PMF:
   //continue;
   }
   cout<<"count_pmf / frame = "<<count_pmf<<" / "<<frame<<endl;
   if( count_pmf != (int)frame )cout<<"WARNING: count_pmf != frame; "<<count_pmf<<" != "<<frame<<endl;
   double min_pmf = 999999, max_pmf = -999999; 
   for(int i=0;i<num_bin;i++){ //normalization
      for(int j=0;j<num_bin;j++){
         pmf[j][i] = pmf[j][i] / frame;
         if( pmf[j][i] <= min_pmf ){
            min_pmf = pmf[j][i];
         }
         if( pmf[j][i] >= max_pmf ){
            max_pmf = pmf[j][i];
         }
      }
   }
   cout<<"Maximum Probability = "<<max_pmf<<endl;
   cout<<"Minimum Probability = "<<min_pmf<<endl;

   //double const_boltz = 1.380658e-23, pmf_temp = 300;
   min_pmf = 999999, max_pmf = -999999; 
   int min_pmf_n[2]; double min_pmf_c[2];
   min_pmf_n[0] = 0; min_pmf_n[1] = 0;
   int max_pmf_n[2]; double max_pmf_c[2];
   max_pmf_n[0] = 0; max_pmf_n[1] = 0;
   for(int i=0;i<num_bin;i++){ //PMF calculation
      for(int j=0;j<num_bin;j++){
         if( pmf[j][i] == 0 ){
	    //pmf[j][i] = 0;
	    continue;
	 }
         pmf[j][i] = -1 * BOLTZMAN_CONST * TEMPERATURE * log( pmf[j][i] );
         if( pmf[j][i] < min_pmf ){
            min_pmf = pmf[j][i];
            min_pmf_n[0] = j; min_pmf_n[1] = i;
         }
         if( pmf[j][i] > max_pmf ){
            max_pmf = pmf[j][i];
            max_pmf_n[0] = j; max_pmf_n[1] = i;
            cout<<"DEBUG> max_pmf = "<<max_pmf<<", max_pmf_n[0], max_pmf_n[1] = "<<max_pmf_n[0]<<", "<<max_pmf_n[1]<<endl;
         }
      }
   }
   min_pmf_c[0] = emin + min_pmf_n[0] * length_bin + length_bin / 2;
   min_pmf_c[1] = emin + min_pmf_n[1] * length_bin + length_bin / 2;
   max_pmf_c[0] = emin + max_pmf_n[0] * length_bin + length_bin / 2;
   max_pmf_c[1] = emin + max_pmf_n[1] * length_bin + length_bin / 2;
   
   cout<<"Maximum PMF = "<<max_pmf<<" (J/mol)"<<endl;
   cout<<"Minimum PMF = "<<min_pmf<<" (J/mol)"<<endl;
   cout<<"(c1, c2) of Maximum PMF = ("<<max_pmf_c[0]<<", "<<max_pmf_c[1]<<")\n"; 
   cout<<"(c1, c2) of Minimum PMF = ("<<min_pmf_c[0]<<", "<<min_pmf_c[1]<<")\n"; 
   
/*   ofstream ofs_pmf;
   ofs_pmf.open("out_pmf.dat");
   for(int i=0;i<num_bin;i++){
      for(int j=0;j<num_bin;j++){
         ofs_pmf<<min_c + length_bin * j<<"   "<<min_c + length_bin * i<<"   "<<pmf[j][i]<<endl;
         //ofs_pmf<<format("%8.3f%8.3f%8.3f \n")  %(min_c + length_bin * j), %(min_c + length_bin * i), %pmf[j][i]; //not work
      }
   }
   ofs_pmf.close();*/
   FILE *fout;
   if((fout = fopen("out_pmf.dat","w")) == NULL ){
      printf("cannot open output file: %s\n","out_pmf.dat");
      return 1;
   }
   for(int i=0;i<num_bin;i++){
      for(int j=0;j<num_bin;j++){
         if( pmf[j][i] == 0 ){
         fprintf(fout,"%12.3f%12.3f%12.3f \n", emin + length_bin * j + length_bin / 2, emin + length_bin * i + length_bin / 2, 0.0 );}else{
         fprintf(fout,"%12.3f%12.3f%12.3f \n", emin + length_bin * j + length_bin / 2, emin + length_bin * i + length_bin / 2, (pmf[j][i] - max_pmf) / JOULE_CALORIE /1000 );
         }
      }
   }
   fclose( fout );
   cout<<"output out_pmf.dat"<<endl;
/* (8) Find the conformation which has the lowest energy
 *
 * */
   cout<<endl<<"REPORT> (8) Find the conformation which has the lowest energy \n";
   double min_8_1 = 999999, min_8_2 = 999999, min_8_3 = 999999, min_8_4 = 999999;
   unsigned int min_n_8 = 0, min_n_8_2 = 0;
   for(unsigned int n=0;n<frame;n++){ //count
      if( fabs(c1[n] - min_pmf_c[0]) <= min_8_1 && fabs(c2[n] - min_pmf_c[1]) <= min_8_2 ){
         min_8_1 = fabs(c1[n] - min_pmf_c[0]);
         min_8_2 = fabs(c2[n] - min_pmf_c[1]);
         min_n_8 = n;
      }
      if( fabs(c1[n] - max_pmf_c[0]) <= min_8_3 && fabs(c2[n] - max_pmf_c[1]) <= min_8_4 ){
         min_8_3 = fabs(c1[n] - max_pmf_c[0]);
         min_8_4 = fabs(c2[n] - max_pmf_c[1]);
         min_n_8_2 = n;
      }
   }
   cout<<"min_8_1, min_8_2 = "<<min_8_1<<", "<<min_8_2<<endl;
   cout<<"closest conformation (most stable): step no. "<<min_n_8<<endl;
   cout<<"(c1, c2) = ("<<c1[min_n_8]<<", "<<c2[min_n_8]<<")\n";
   cout<<"min_8_3, min_8_4 = "<<min_8_3<<", "<<min_8_4<<endl;
   cout<<"closest conformation (most unstable): step no. "<<min_n_8_2<<endl;
   cout<<"(c1, c2) = ("<<c1[min_n_8_2]<<", "<<c2[min_n_8_2]<<")\n";

   tra1->write_step("mostStable.pdb", min_n_8);
   cout<<"output mostStable.pdb\n";
   tra1->write_step("mostUnStable.pdb", min_n_8_2);
   cout<<"output mostUnStable.pdb\n";

   vector<int> mem_8;
   char filename[256]={0};   //string filename;     
   system("rm -r ./yyy");
   system("mkdir ./yyy");
   for(unsigned int n=0;n<frame;n++){ //count
      if(c1[n] >= emin + length_bin * min_pmf_n[0] 
            && c1[n] <= emin + length_bin * (min_pmf_n[0] + 1)
            && c2[n] >= emin + length_bin * min_pmf_n[1] 
            && c2[n] <= emin + length_bin * (min_pmf_n[1] + 1)     ){
         mem_8.push_back(n);
         sprintf( filename,"./yyy/%i.pdb",n );   //filename = "./yyy/" + to_string( (long long int)n ) + ".pdb";     
         tra1->write_step( filename , n);   //tra1->write_step( filename.c_str() , n);
      }
   }
   cout<<"output the most stable pdbs in ./yyy/ \n";
   cout<<"the coordinates of most stable bin. (total num = "<<mem_8.size()<<")\n";
/*   cout<<"step num.s\n";
   for(unsigned int i=0;i<mem_8.size();i++){
      cout<<mem_8[i]<<"   (c1, c2) = ("<<c1[mem_8[i]]<<", "<<c2[mem_8[i]]<<")\n";
   }*/

// END
	delete tra1;
	return 0;
}

int pca_prctc(char *codname,char *pdbname){
   return pca_prctc( codname, pdbname, 1);
}


