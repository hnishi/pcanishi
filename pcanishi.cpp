#include"nlib.h"
#include"math_nishi.h"
#include<Eigen/Dense>  // Eigen/Core, Geometry, LU, Cholesky, SVD, QR and Eigenvalues
//#include"inpnishi.h"

//#include<sstream>
using namespace Eigen;

// user define
#define TEMPERATURE 300.0
#define BIN_DE 5 
#define BIN_EMIN1 -200

// constant values
#define GAS_CONST 8.31451  // gas constant, R ( joule/mol*k )
//#define BOLTZMAN_CONST 1.380658e-23  // Boltzman constant, kb ( joule/k )
#define BOLTZMAN_CONST 8.3144621  // Boltzman constant, kb ( joule/k*mol )
// kb = R / Na, where Na is Avogadro constant (6.02214129e-23)
// printf("Boltzman constant kb = %e \n", kb);
#define JOULE_CALORIE 4.184 // joule-calorie conversion unit

#define DEBUG_PCA 1 // for debugging, instead of comment out
//      #if DEBUG_PCA == 1
//            ........................
//      #endif


vector<double> quaternion( vector<double> &vec_ref, vector<double> &vec_tar );  // give two vectors as pointer, so please notice the changes of these vectors in this function will remain out of it.
int transfer_quat( vector<double> &vec, vector<double> &transf );
int rotate_quat( vector<double> &vec, vector<double> &rot_mat );


/* ********************************************
 *    PCA by ATOM COORDINATES
 * ********************************************/
int pcanishi(  Inp_nishi inp1 ){

// ############# trajectory  ########################################
/* (1) load trajectory by tra_nishi
 * 
 */
   //cout<<endl<<"REPORT> (1) load trajectory \n";
   cout<<endl<<"--- TRAJECTORY INFORMATION --- \n";
   string codname = inp1.read("COD1");
   string pdbname = inp1.read("REFPDBNAME");
   string pcaatom = inp1.read("PCAATOM") ;
   int stride = atoi(inp1.read("STRIDE").c_str());
   int startframe = atoi( inp1.read("STARTFRAME").c_str() ) - 1 ;
   string superpbase = inp1.read("SUPERPBASE");   //v.1.1.0

   if(stride <= 0){
      return -1;
   }
   cout<<"loading\n"<<codname<<"        "<<pdbname<<endl;
	tra_nishi* tra1;
	tra1 = new tra_nishi(codname.c_str(), pdbname.c_str(), stride, pcaatom);
	cout<<"TOTAL FRAME = "<<tra1->total_step<<endl;
	cout<<"TOTAL ATOM = "<<tra1->pdb1->total_atom<<endl;
	cout<<"TOTAL SELECTED ATOM = "<<tra1->total_sel<<endl;
	//unsigned int frame = tra1->total_step;
/* (2) setting of region
 *     from intra_start to intra_end (internal number of atom)
 */
   //cout<<endl<<"REPORT> (2) specify the region \n";
   cout<<endl<<"--- RESIDUE RANGE --- \n";
   char startchain = inp1.read("STARTCHAIN").c_str()[0];
   char endchain = inp1.read("ENDCHAIN").c_str()[0];
   int startres = atoi(inp1.read("STARTRES").c_str());
   int endres = atoi(inp1.read("ENDRES").c_str());
   int intra_start, intra_end;
	intra_start = tra1->pdb1->search_n( startchain , startres );
	intra_end = tra1->pdb1->search_n_end( endchain , endres );
	cout<<"intra_start = "<<intra_start<<endl;
	cout<<"intra_end = "<<intra_end<<endl;


	tra1->pdb1->disp_line(intra_start-1);
	tra1->pdb1->disp_line(intra_start);
	cout<<"...\n";
        tra1->pdb1->disp_line(intra_end);
        tra1->pdb1->disp_line(intra_end+1);

        //tra1->pdb1->write_pdb("zzz.pdb"); cout<<"output pdb file (zzz.pdb)\n";

/* (3) ATOM SELECTION
*/
   cout<<"\n--- ATOM SELECTION --- \n";
   int rej_ca=0, rej_h=0, rej_wat=0, rej_cim=0, rej_cip=0, rej_mainchain=0;
   int rtrn_sel, count=0, start_sel = 0, end_sel = 0;
   vector<double> vec_ref;
      for( int i=0;(signed)i<(signed)tra1->pdb1->total_atom;i++){
      //for( int i=intra_start;i<=intra_end;i++){
         rtrn_sel = select_atom( *tra1->pdb1, vec_ref, pcaatom, i );
	 if( i>=intra_start && i<=intra_end ){
	    switch( rtrn_sel ){
	    case 0: 
	       if( start_sel == 0){
	          start_sel = count;
		  tra1->pdb1->disp_line(i);
   	       }
	       end_sel = count ;
	       tra1->pdb1->disp_line(i);
	       break;
	    case 1: rej_ca++; break;
	    case 2: rej_mainchain++; break;
	    case 3: rej_h++; break;
	    case 4: rej_wat++; break;
	    case 5: rej_cim++; break;
	    case 6: rej_cip++; break;
	    default: cout<<"Unknown value of rtrn_sel \n";
	    }
	 }
	 if(rtrn_sel==0)count++;
      }
      
      int cycle_frame = end_sel - start_sel + 1;
      cout<<"!!!! start_sel and end_sel = "<<start_sel<<" and "<<end_sel<<endl;
      cout<<"Atoms in region between STARTRES and ENDRES = "<<intra_end - intra_start + 1<<endl;
      cout<<"Num. of atoms selected for PCA calculation = "<<end_sel - start_sel +1<<" in reference"<<endl;
      cout<<"Rejected atoms are as follows (Selected: "<<pcaatom<<")\n";
      cout<<"!CA atoms = "<<rej_ca<<endl;
      cout<<"!CA & !N & !C & !O atoms = "<<rej_mainchain<<endl;
      cout<<"element H (Hydrogens) = "<<rej_h<<endl;
      cout<<"residue name WAT (Water molecules) = "<<rej_wat<<endl;
      cout<<"residue name CIM (minus ions) = "<<rej_cim<<endl;
      cout<<"residue name CIP (plus ions) = "<<rej_cip<<endl;
   
   //vec.clear();
   int cod_num_i=2;
   cout<<"!!!!! cod_num = COD1"<<endl;

   vector<double> vec, vec_buf;
   if( superpbase == "YES" ){
      for(int i=0;i<start_sel*3;i++){
        vec_buf.push_back(vec_ref[i]);
      }
      for(unsigned int i=end_sel*3+3;i<vec_ref.size();i++){
        vec_buf.push_back(vec_ref[i]);
      }
   }
   else{
      for(int i=start_sel*3;i<=end_sel*3+2;i++){
        vec_buf.push_back(vec_ref[i]);
      }
   }
   //cout<<"!!!!!!! vec_buf.size() = "<<vec_buf.size()<<endl;
   //cout<<"!!!!!!! vec_ref.size() = "<<vec_ref.size()<<endl;
   //cout<<"!!! vec_ref[0] = "<<vec_buf[0]<<", vec_ref[n] = "<<vec_buf[vec_buf.size() -1]<<endl;

flag100:
   cout<<"TOTAL FRAME = "<<tra1->total_step<<endl;
   cout<<"TOTAL ATOM = "<<tra1->pdb1->total_atom<<endl;
   cout<<"TOTAL SELECTED ATOM = "<<tra1->total_sel<<endl;
   vector<double> vec_tar;
   //cout<<"!!!!! tra1->cordx.size() = "<<tra1->cordx.size()<<endl;
   for(unsigned int n=startframe;n<tra1->total_step;n++){
      vec_tar.clear();
      if( superpbase == "YES" ){
         for(int i=0;i<start_sel;i++){
            vec_tar.push_back(tra1->cordx[n*tra1->total_sel+i]);
            vec_tar.push_back(tra1->cordy[n*tra1->total_sel+i]);
            vec_tar.push_back(tra1->cordz[n*tra1->total_sel+i]);
         }
         for(unsigned int i=end_sel+1;i<tra1->total_sel;i++){
            vec_tar.push_back(tra1->cordx[n*tra1->total_sel+i]);
            vec_tar.push_back(tra1->cordy[n*tra1->total_sel+i]);
            vec_tar.push_back(tra1->cordz[n*tra1->total_sel+i]);
         }
         //cout<<"!!!!!!! vec_tar.size() = "<<vec_tar.size()<<endl;
         //cout<<"!!!!!!! total_sel*3 = "<<tra1->total_sel*3<<endl;
         vector<double> rot_mat = quaternion( vec_buf, vec_tar );  // give two vectors as pointer, so please notice the changes of these vectors in this function will remain out of it.
	 vec_tar.clear();
         for(int i=start_sel;i<=end_sel;i++){
            vec_tar.push_back(tra1->cordx[n*tra1->total_sel+i]);
            vec_tar.push_back(tra1->cordy[n*tra1->total_sel+i]);
            vec_tar.push_back(tra1->cordz[n*tra1->total_sel+i]);
         }
         vector<double> transf;  // transfer target
         transf.push_back( rot_mat[12] );
         transf.push_back( rot_mat[13] );
         transf.push_back( rot_mat[14] );
         transfer_quat( vec_tar, transf );

         /*transf.clear();  //no need to transfer reference because do not calculate RMSD
         transf.push_back( rot_mat[ 9] );
         transf.push_back( rot_mat[10] );
         transf.push_back( rot_mat[11] );
         transfer_quat( pdb_ref->coox, pdb_ref->cooy, pdb_ref->cooz, transf );
*/
         rotate_quat( vec_tar, rot_mat );  // rotate target
       
      }
      else {
         for(int i=start_sel;i<=end_sel;i++){
            vec_tar.push_back(tra1->cordx[n*tra1->total_sel+i]);
            vec_tar.push_back(tra1->cordy[n*tra1->total_sel+i]);
            vec_tar.push_back(tra1->cordz[n*tra1->total_sel+i]);
         }
         ///cout<<n<<" rmsd before = "<<rmsd2(vec_tar,vec_buf)<<endl;
         quaternion( vec_buf, vec_tar );  // give two vectors as pointer, so please notice the changes of these vectors in this function will remain out of it.
      }
      //cout<<"!!! vec_tar[0] = "<<vec_tar[0]<<", vec_tar[n] = "<<vec_tar[vec_tar.size() -1]<<endl;
      //cout<<n<<" rmsd after = "<<rmsd2(vec_tar,vec_buf)<<endl;
      //cout<<"!!! vec_ref[0] = "<<vec_buf[0]<<", vec_ref[n] = "<<vec_buf[vec_buf.size() -1]<<endl;
      //cout<<"!!! vec_tar[0] = "<<vec_tar[0]<<", vec_tar[n] = "<<vec_tar[vec_tar.size() -1]<<endl;
   
      //cout<<"!!!!! size of cycle_frame and vec_tar = "<<cycle_frame<<" and "<<vec_tar.size()<<endl;
      for(unsigned int i=0;i<vec_tar.size();i++){
         vec.push_back( vec_tar[i] );
      }
   }
   //cout<<"!!!!! size of vec_ref and vec_tar = "<<vec_buf.size()<<" and "<<vec_tar.size()<<endl;
   //tra1->write_step("zzz.pdb",tra1->total_step -1);

   string cod_num;  char buf[32];
   sprintf(buf,"COD%d",cod_num_i);  //itoa()
   cod_num = buf;
   cout<<"!!!!! cod_num = "<<cod_num<<endl;
   cod_num_i++;
   codname = inp1.read(cod_num);
   if( codname == "nothing" ){
      cout<<"Section of reading trajectories was ended normaly\n";
      delete tra1;
   }else{
      delete tra1;
      tra1 = new tra_nishi(codname.c_str(), pdbname.c_str(), stride, pcaatom);
      goto flag100;
   }
   //unsigned int dim_0 = end_sel - start_sel +1;

   cout<<"\n--- PCA CALCULATION --- \n";
   //cout<<endl<<"REPORT> (3) PCA calculation starts \n";
/* ********
 * 3-1   create vector Q = (q1x,q1y,q1z,q2x,...,q20z) 
 *     
 * */
   //cout<<" 3-1  create vector Q \n";
   unsigned int dim_Q = vec_tar.size();
   cout<<"dimensionality of Q (components per one structure) is "<<dim_Q<<endl;
   cout<<"cycle_frame = "<<cycle_frame<<endl;
   unsigned int frame = vec.size() / dim_Q;
   cout<<"Thr number of structures is "<<frame<<endl;

/* *********
 * 3-2   create variance-covariance matrix
 *
 * */
   cout<<" 3-2  create variance-covariance matrix \n";
   double *average_Q; double *average_QQ;
   average_Q = new double [dim_Q];
   average_QQ = new double [dim_Q*dim_Q];
   // initialize 2d-array
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
         average_Q[i] += vec[i+n*dim_Q];
			//cerr<<average_Q[i]<<" ";
			//cout<<"average_Q["<<i<<"] = "<<average_Q[i]<<endl;
	 for(unsigned int j=0;j<dim_Q;j++){
	    average_QQ[ j + i * dim_Q ] += vec[j+n*dim_Q] * vec[i+n*dim_Q];  // 
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
	//delete average_Q;
	delete[] average_QQ;
	//cerr<<VCV<<endl;
        for(unsigned int n=0;n<dim_Q;n++){
                for(unsigned int i=0;i<dim_Q;i++){
			//cout<<"DEBUG: VCV("<<i<<","<<n<<") = "<<VCV(i,n)<<", VCV("<<n<<","<<i<<") = "<<VCV(n,i)<<endl;
			if( VCV(i,n) != VCV(n,i) )cout<<"ERROR: !!!!!!!!!!!! VCV is not symmetric !!!!!!!!!!!!!!\n";
                }
        }
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

	string outeigen = inp1.read("OUTEIGEN");
        ofstream ofs_eigen;
        ofs_eigen.open( outeigen.c_str() );
        ofs_eigen<<"> Eigen values\n"<<es.eigenvalues()<<endl;
 
   for(unsigned int n=0;n<dim_Q;n++){
      ofs_eigen<<"> Eigen vector no"<<n+1<<"\n"<<es.eigenvectors().col(n).transpose()<<endl;
   }

        ofs_eigen.close(); 
	cout<<"ouput "<<outeigen<<"\n";

/* **********
 *  3-5  calculate the coordinates of the ensemble along the PCA axes
 *
 * */
   cout<<" 3-5  calculate components along principle axis \n";
	double* c1; double* c2;
	double* c3; double* c4;
	c1 = new double [frame]; c2 = new double [frame];
	c3 = new double [frame]; c4 = new double [frame];

	vector<double> buf_vec;

        for(unsigned int i=0;i<dim_Q;i++){
	        buf_vec.push_back(average_Q[i]);
        }
	VectorXd Q2=Map<VectorXd>(&buf_vec[0],buf_vec.size());
	delete[] average_Q;
        for(unsigned int n=0;n<frame;n++){
		buf_vec.clear();
                for(unsigned int i=0;i<dim_Q;i++){
			buf_vec.push_back(vec[i+n*dim_Q]);
		}
		VectorXd Q1=Map<VectorXd>(&buf_vec[0],buf_vec.size());
        	c1[n] = es.eigenvectors().col(dim_Q-1).transpose() * ( Q1 - Q2 );
		//cout<<"c1["<<n<<"] = "<<c1[n]<<endl;
        	c2[n] = es.eigenvectors().col(dim_Q-2).transpose() * ( Q1 - Q2 );
		//cout<<"c1["<<n<<"] = "<<c1[n]<<endl;
        	c3[n] = es.eigenvectors().col(dim_Q-3).transpose() * ( Q1 - Q2 );
        	c4[n] = es.eigenvectors().col(dim_Q-4).transpose() * ( Q1 - Q2 );
	}	

/* (4) output (write) c1 vs c2 
*/
	string c1c2out;
	c1c2out = inp1.read("C1C2OUT").c_str();
	if( c1c2out != "NO" ){
		ofstream ofs;
		ofs.open( c1c2out.c_str() );
	        for(unsigned int n=0;n<frame;n++){
        		ofs<<c1[n]<<"   "<<c2[n]<<"   "<<c3[n]<<"   "<<c4[n]<<"   "<<n+1<<endl;
		}
		ofs.close();
   		cout<<"output c1vsc2.dat (2-D graph) \n";
   	}
	else{
		cout<<"do not output C1C2OUT \n";
	}
	delete[] c1; delete[] c2; 
	delete[] c3; delete[] c4; 
/* (6) calculate contribution ratio of c1 and c2
 *     
 * */
   cout<<endl<<"REPORT> (6) calculate contribution ratio of c1 and c2 \n";
   double sum_lambda = 0;
   for(unsigned int i=0;i<dim_Q;i++){
      sum_lambda += es.eigenvalues()[i];
   }
   double cntr_rt = (es.eigenvalues()[dim_Q - 1]) / sum_lambda;
   cout<<"contribution ratio (PC1) = "<<cntr_rt<<endl;
   cntr_rt = (es.eigenvalues()[dim_Q - 1] + es.eigenvalues()[dim_Q - 2]) / sum_lambda;
   cout<<"contribution ratio (PC1 + PC2) = "<<cntr_rt<<endl;
   cntr_rt = (es.eigenvalues()[dim_Q - 1] + es.eigenvalues()[dim_Q - 2] + es.eigenvalues()[dim_Q - 3]) / sum_lambda;
   cout<<"contribution ratio (PC1 + PC2 + PC3) = "<<cntr_rt<<endl;
   cntr_rt = (es.eigenvalues()[dim_Q - 1] + es.eigenvalues()[dim_Q - 2] + es.eigenvalues()[dim_Q - 3] + es.eigenvalues()[dim_Q - 4]) / sum_lambda;
   cout<<"contribution ratio (PC1 + PC2 + PC3 + PC4) = "<<cntr_rt<<endl;


// END
        return 0;
}
