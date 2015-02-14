#include"math_nishi.h"
//#include"inpnishi.h"
//#include<Eigen/Eigen>
#include<Eigen/Dense>
//#include<Eigen/Eigenvalues>
using Eigen::MatrixXd;
//using Eigen::EigenSolver;
using Eigen::SelfAdjointEigenSolver;


struct inp_se{ 
   vector<double> vec_tar;
   vector<double> vec_ref; 
};
//struct inp_se vec_quat;

/* **********************
 *   proto-type delaration
 * *********************/
vector<double> quaternion( vector<double> &vec_ref, vector<double> &vec_tar );  // give two vectors as pointer, so please notice the changes of these vectors in this function will remain out of it.
int transfer_quat( vector<double> &x, vector<double> &y,  vector<double> &z, vector<double> &transf );
int transfer_quat( vector<double> &vec, vector<double> &transf );
int rotate_quat( vector<double> &x, vector<double> &y,  vector<double> &z, vector<double> &rot_mat );
int rotate_quat( vector<double> &vec, vector<double> &rot_mat );
int select_quat( pdb_nishi &pdb1, vector<double> &vec, string &rmsdatom, int i );
int select_quat( pdb_nishi &pdb1, vector<double> &x, vector<double> &y,vector<double> &z, vector<double> &vec, string &rmsdatom, int i );

/* ************************************* 
 *   quatnishi 
 *
 * *************************************/
//int quatnishi(char *codname,char *pdbname,int stride){
int quatnishi( Inp_nishi inp1 ){
   // (1) input depending on RMSDMODE and RMSDATOM
   cout<<"\n------ MODE INPUT ------\n";
   string rmsdmode = inp1.read("RMSDMODE") ;  //char rmsdmode[5];
   string rmsdatom = inp1.read("RMSDATOM") ;  
   string inversermsd = inp1.read("INVERSERMSD") ;  
   
   // (2) read pdb and memorize the range of residues for rmsd calculation
   pdb_nishi* pdb_tmp;
   //pdb_tmp = new pdb_nishi( inp1.read("REFPDBNAME").c_str() );
   char refpdbname[100];  // = inp1.read("REFPDBNAME").c_str();
   strcpy(refpdbname,inp1.read("REFPDBNAME").c_str() );
   pdb_tmp = new pdb_nishi( refpdbname );
   cout<<"\n------ REFERENCE PDB INFORMATION ------\n";
   cout<<"Total_atom of pdb = "<<pdb_tmp->total_atom<<endl;

   cout<<"\n------ RANGE AND ATOM SELECTION ------\n";
   char startchain = inp1.read("STARTCHAIN").c_str()[0];
   char endchain = inp1.read("ENDCHAIN").c_str()[0];
   int startres = atoi(inp1.read("STARTRES").c_str());
   int endres = atoi(inp1.read("ENDRES").c_str());
   int intra_start, intra_end;
   //intra_start = pdb_tmp->search_n( 'B', atoi(inp1.read("STARTRES").c_str() ));
   //cout<<"DEBUG "<<inp1.read("STARTCHAIN").c_str()[0]<<endl;
   //intra_start = pdb_tmp->search_n( inp1.read("STARTCHAIN").c_str()[0], atoi(inp1.read("STARTRES").c_str() ));
   intra_start = pdb_tmp->search_n( startchain, startres );
   if( pdb_tmp->rnum[pdb_tmp->total_atom - 1] == endres ){
      intra_end = pdb_tmp->total_atom -1;
   }else{
      intra_end = pdb_tmp->search_n(endchain,endres + 1) -1;
   }
   cout<<"\nplease confirm the boader line of residue-range for rmsd calculation \n";
   cout<<"intra_start = "<<intra_start<<endl;
   cout<<"intra_end = "<<intra_end<<endl;

   pdb_tmp->disp_line(intra_start-1);
   pdb_tmp->disp_line(intra_start);
   cout<<"...\n";
   pdb_tmp->disp_line(intra_end);
   pdb_tmp->disp_line(intra_end+1);

   int intra_end2=intra_end, intra_start2=intra_start, flag=0;
   if( inversermsd == "YES" ){
      intra_start = 0;
      intra_end = intra_start2 - 1;
   }

      int rej_ca=0, rej_h=0, rej_wat=0, rej_cim=0, rej_cip=0, rej_mainchain=0;
      int rtrn_sel;
      vector<double> vec_ref;
flag100:
      for(int i=intra_start;i<=intra_end;i++){
         rtrn_sel = select_quat( *pdb_tmp, vec_ref, rmsdatom, i );
	 switch( rtrn_sel ){
	 case 0: break;
	 case 1: rej_ca++; break;
	 case 2: rej_mainchain++; break;
	 case 3: rej_h++; break;
	 case 4: rej_wat++; break;
	 case 5: rej_cim++; break;
	 case 6: rej_cip++; break;
	 default: cout<<"Unknown value of rtrn_sel \n";
	 }
      }
      if( inversermsd == "YES" && flag == 0 ){
         intra_start = intra_end2 + 1;
	 intra_end = pdb_tmp->total_atom - 1;
	 flag = 100;
         goto flag100;
      }
      if( flag == 100 ){
         intra_start = 0;
	 intra_end = intra_start2 - 1;
      }

      cout<<"Atoms in region between STARTRES and ENDRES = "<<intra_end - intra_start + 1<<endl;
      cout<<"Num. of atoms selected for rmsd calculation = "<<vec_ref.size()/3<<" in reference"<<endl;
      cout<<"Rejected atoms are as follows (Selected: "<<rmsdatom<<")\n";
      cout<<"!CA atoms = "<<rej_ca<<endl;
      cout<<"!CA & !N & !C & !O atoms = "<<rej_mainchain<<endl;
      cout<<"element H (Hydrogens) = "<<rej_h<<endl;
      cout<<"residue name WAT (Water molecules) = "<<rej_wat<<endl;
      cout<<"residue name CIM (minus ions) = "<<rej_cim<<endl;
      cout<<"residue name CIP (plus ions) = "<<rej_cip<<endl;

   delete pdb_tmp;//free dynamic memory

   // **********************************************************************************
   // (3-1) superpose pdb to pdb and calculate RMSD
   // **********************************************************************************
   // 3-1-1 create vectors which include only atoms for rmsd calculation
   if( rmsdmode == "pdb" ){
      cout<<"\n------ PDB ------\n";
      cout<<"begin to calculate RMSD between pdb and pdb \n";
      pdb_nishi* pdb_ref;
      pdb_ref = new pdb_nishi( refpdbname );
      pdb_nishi* pdb_tar;
      pdb_tar = new pdb_nishi( inp1.read("PDBNAME").c_str() );
      /*if( pdb_ref->total_atom != pdb_tar->total_atom ){
         cerr<<"ERROR: Total num.s of target pdb and reference pdb are different\n";
	 cerr<<"exit this program by quatnishi(3-1-1) \n";
         exit(1);
      }*/

      vector<double> vec_tar;
flag200:
      for(int i=intra_start;i<=intra_end;i++){
         select_quat( *pdb_tar, vec_tar, rmsdatom, i );
      }

      if( inversermsd == "YES" && flag == 100 ){
         intra_start = intra_end2 + 1;
	 intra_end = pdb_tar->total_atom - 1;
	 flag = 200;
         goto flag200;
      }
      if( flag == 200 ){
         intra_start = 0;
	 intra_end = intra_start2 - 1;
      }

      cout<<"Num. of atoms selected for rmsd calculation = "<<vec_tar.size()/3<<" in target"<<endl;
      if( vec_ref.size() != vec_tar.size() ){
         cerr<<"ERROR: num. of selected atoms is different of reference and target \n";
	 cerr<<"reference atoms: "<<vec_ref.size()/3<<", target atoms: "<<vec_tar.size()/3<<endl;
	 exit(1);
      }


      // 3-1-2  get rotation matrix and transform coordinates
      vector<double> rot_mat;
      rot_mat = quaternion( vec_ref, vec_tar );  // get rotation matrix

      vector<double> transf;  // transfer target
      transf.push_back( rot_mat[12] );
      transf.push_back( rot_mat[13] );
      transf.push_back( rot_mat[14] );
      transfer_quat( pdb_tar->coox, pdb_tar->cooy, pdb_tar->cooz, transf );

      transf.clear();  // transfer reference
      transf.push_back( rot_mat[ 9] );
      transf.push_back( rot_mat[10] );
      transf.push_back( rot_mat[11] );
      transfer_quat( pdb_ref->coox, pdb_ref->cooy, pdb_ref->cooz, transf );

      rotate_quat( pdb_tar->coox, pdb_tar->cooy, pdb_tar->cooz, rot_mat );  // rotate target
       
      vector<double> ax,ay,az,bx,by,bz;  // format
      for(unsigned int i=0;i<vec_tar.size();i=i+3){ //for rmsd()
         ax.push_back( vec_tar[i] );
         ay.push_back( vec_tar[1+i] );
         az.push_back( vec_tar[2+i] );
         bx.push_back( vec_ref[i] );
         by.push_back( vec_ref[1+i] );
         bz.push_back( vec_ref[2+i] );
     } 
       	//FILE *fout2;  char filename2[] = "superpose.pdb";  if((fout2 = fopen(filename2,"w")) == NULL ){  printf("cannot open output file: %s\n",filename2);  exit(1); }  for(unsigned int n=0;n<ax.size();n++){  fprintf( fout2, "%4s%2s%5d%1s%-4s%8s%6s%8.3f%8.3f%8.3f%22s%-2s\n","ATOM"," ",n+1," ","O","MOL E 1"," ",ax[n],ay[n],az[n]," ","O" );  }  for(unsigned int n=0;n<ax.size();n++){  fprintf( fout2, "%4s%2s%5d%1s%-4s%8s%6s%8.3f%8.3f%8.3f%22s%-2s\n","ATOM"," ",n+1," ","C","MOL F 1"," ",bx[n],by[n],bz[n]," ","C" );   	}  fclose( fout2 );

      // 3-1-3  output results
      double rmsdq = rmsd(pdb_tar->coox,pdb_tar->cooy,pdb_tar->cooz,pdb_ref->coox,pdb_ref->cooy,pdb_ref->cooz);
      cout<<"RMSD of all atoms in pdb (including water) = "<<rmsdq<<" A \n";
      double rmsd_sel = rmsd(ax,ay,az,bx,by,bz);
      cout<<"RMSD of selected atoms = "<<rmsd_sel<<" A, selection = "<<rmsdatom<<endl;;

      char superpdb[100];
      strcpy( superpdb, inp1.read("SUPERPDB").c_str() );
      FILE *fout;
      if((fout = fopen(superpdb,"w")) == NULL ){
		printf("cannot open output file: %s\n",superpdb);
		exit(1);
      }
      fprintf(fout,"REMARK  RMSD of all atoms in pdb (including water) = %f A\n",rmsdq);
      fprintf(fout,"REMARK  RMSD of selected region = %f A, ATOM SELECTION = %s\n",rmsd_sel,rmsdatom.c_str());
      fprintf(fout,"REMARK  STARTCHAIN: %c, STARTRES: %i, ENDCHAIN: %c, ENDRES: %i\n",startchain,startres,endchain,endres);
      fclose(fout);
      pdb_tar->write_pdb(superpdb,'a');
      pdb_ref->write_pdb(superpdb,'a');

      delete pdb_tar;//free dynamic memory
      delete pdb_ref;//free dynamic memory
   }

   else if( rmsdmode == "cod" || rmsdmode == "crd"){
   // ********************************************************************
   // (3-2) calculate rmsd about trajectory to reference
   // ********************************************************************
   // 3-2-1 create vectors which include only atoms for rmsd calculation
      cout<<"\n------ COD ------\n";
      cout<<"begin to calculate RMSD of trajectory \n";
      tra_nishi* tra1;
      tra1 = new tra_nishi( inp1.read("CODNAME").c_str() ,refpdbname );
      cout<<"TOTAL FRAME = "<<tra1->total_step<<endl;

      vector<double> vec_tar;
      vector<double> rmsd_tra;
      vector<double> buf_x, buf_y, buf_z;
      for(unsigned int n=0;n<tra1->total_step;n++){
         for(unsigned int ii=0;ii<tra1->pdb1->total_atom;ii++){
	    buf_x.push_back( tra1->cordx[n*tra1->pdb1->total_atom+ii] );
	    buf_y.push_back( tra1->cordy[n*tra1->pdb1->total_atom+ii] );
	    buf_z.push_back( tra1->cordz[n*tra1->pdb1->total_atom+ii] );
	 }
flag300:
         for(int i=intra_start;i<=intra_end;i++){
            select_quat( *tra1->pdb1, buf_x, buf_y, buf_z, vec_tar, rmsdatom, i );
	 }

         if( inversermsd == "YES" && flag == 100 ){
            intra_start = intra_end2 + 1;
       	    intra_end = tra1->pdb1->total_atom - 1;
	    flag = 300;
            goto flag300;
         }
         if( flag == 300 ){
            intra_start = 0;
	    intra_end = intra_start2 - 1;
	    flag = 100;  // 
         } 

         // 3-2-2  get rotation matrix and transform coordinates
         quaternion( vec_ref, vec_tar );  // get rotation matrix
         vector<double> ax,ay,az,bx,by,bz;  // format
         for(unsigned int i=0;i<vec_tar.size();i=i+3){ //for rmsd()
            ax.push_back( vec_tar[i] );
            ay.push_back( vec_tar[1+i] );
            az.push_back( vec_tar[2+i] );
            bx.push_back( vec_ref[i] );
            by.push_back( vec_ref[1+i] );
            bz.push_back( vec_ref[2+i] );
         }  
         rmsd_tra.push_back( rmsd(ax,ay,az,bx,by,bz) );
	 buf_x.clear(); buf_y.clear(); buf_z.clear(); vec_tar.clear();
      }
      string rmsdfile = inp1.read("RMSDFILE") ;  
      FILE *fout;
      if((fout = fopen( rmsdfile.c_str() ,"w")) == NULL ){
		printf("cannot open output file: %s\n", rmsdfile.c_str() );
		exit(1);
      }
      for(unsigned int i=0;i<tra1->total_step;i++){
         fprintf( fout,"%12i %12.5f \n", i+1, rmsd_tra[i] );
         //cout<<"Frame "<<i+1<<",  RMSD = "<<rmsd_tra[i]<<endl;
      }
      fclose(fout);

      double tot=0;
      for(unsigned int i=0;i<rmsd_tra.size();i++){
         tot = tot + rmsd_tra[i];
      }
      tot = tot / rmsd_tra.size();
      cout<<"Average RMSD of trajectory = "<<tot<<" A \n";
   }
   // *************************************************************
   //  else
   // *************************************************************
   else{
      cerr<<"ERROR: cannot recognaize RMSDMODE; "<<rmsdmode<<endl;
      cerr<<"exit this program by quatnishi() \n";
      exit(1);
   }

   // END of  quatnishi()
   return 0;
}

/* ************************************* 
 *  quaternion() in quatnishi.cpp
 *
 * *************************************/
vector<double> quaternion( vector<double> &vec_ref, vector<double> &vec_tar ){
   vector<double> rot_mat;

   if( vec_ref.size() != vec_tar.size() ){
      printf("ERROR: in quaternion(), the size of two vectors are different\n");
      cerr<<"this program was ended by quaternion() in quatnishi.cpp\n";
      exit(1);
   }

   // (1) calculate center of mass
   double com_vecs[6]; for(int i=0;i<6;i++){com_vecs[i]=0;}
   for( unsigned int i=0;i<vec_ref.size();i = i + 3 ){
      com_vecs[0] = com_vecs[0] + vec_ref[i];
      com_vecs[1] = com_vecs[1] + vec_ref[i+1];
      com_vecs[2] = com_vecs[2] + vec_ref[i+2];
      com_vecs[3] = com_vecs[3] + vec_tar[i];
      com_vecs[4] = com_vecs[4] + vec_tar[i+1];
      com_vecs[5] = com_vecs[5] + vec_tar[i+2];
   }
   for(int i=0;i<6;i++){
      com_vecs[i] = com_vecs[i] / ( vec_ref.size() / 3 );
      //cout<<"center of mass of vec-num. "<<i<<" = "<<com_vecs[i]<<endl;
   }

   // (2) transfer c.o.m. to origin
   for( unsigned int i=0;i<vec_ref.size();i = i + 3 ){
      vec_ref[i] = vec_ref[i] - com_vecs[0];
      vec_ref[i+1] = vec_ref[i+1] - com_vecs[1];
      vec_ref[i+2] = vec_ref[i+2] - com_vecs[2];
      vec_tar[i] = vec_tar[i] - com_vecs[3];
      vec_tar[i+1] = vec_tar[i+1] - com_vecs[4];
      vec_tar[i+2] = vec_tar[i+2] - com_vecs[5];
   }
   
   // (3) for skew matrix Ak, vec_ref + vec_tar and vec_ref - vec_tar
   double ak_plus[vec_ref.size()], ak_minus[vec_ref.size()];
   for( unsigned int i=0;i<vec_ref.size();i++ ){
      ak_plus[i] = vec_ref[i] + vec_tar[i];
      ak_minus[i] = vec_ref[i] - vec_tar[i];
   }
   
   // (4) Symmetric matrix B = Akt*Ak
   double B11=0,B12=0,B13=0,B14=0,B21,B22=0,B23=0,B24=0,B31,B32,B33=0,B34=0,B41,B42,B43,B44=0;
   for(unsigned int i=0;i<vec_ref.size();i=i+3){ // Symmetric matrix B = Akt*Ak
      B11 +=  ak_minus[i]*ak_minus[i] + ak_minus[i+1]*ak_minus[i+1] + ak_minus[i+2]*ak_minus[i+2];
      B12 +=  ak_plus[i+2]*ak_minus[i+1] - ak_minus[i+2]*ak_plus[i+1];
      B13 += -ak_plus[i+2]*ak_minus[i] + ak_plus[i]*ak_minus[i+2];
      B14 +=  ak_plus[i+1]*ak_minus[i] - ak_plus[i]*ak_minus[i+1];
      B22 +=  ak_minus[i]*ak_minus[i] + ak_plus[i+1]*ak_plus[i+1] + ak_plus[i+2]*ak_plus[i+2];
      B23 +=  ak_minus[i]*ak_minus[i+1] - ak_plus[i]*ak_plus[i+1];
      B24 +=  ak_minus[i]*ak_minus[i+2] - ak_plus[i+2]*ak_plus[i];
      B33 +=  ak_plus[i]*ak_plus[i] + ak_minus[i+1]*ak_minus[i+1] + ak_plus[i+2]*ak_plus[i+2];
      B34 +=  ak_minus[i+1]*ak_minus[i+2] - ak_plus[i+1]*ak_plus[i+2];
      B44 +=  ak_plus[i]*ak_plus[i] + ak_plus[i+1]*ak_plus[i+1] + ak_minus[i+2]*ak_minus[i+2];
      }
	unsigned int n = vec_ref.size();
	B11/=n;B12/=n;B13/=n;B14/=n;B22/=n;B23/=n;B24/=n;B33/=n;B34/=n;B44/=n;

	B21 = B12;
	B31 = B13;
	B41 = B14;
	B32 = B23;
	B42 = B24;
	B43 = B34;
	//puts("B matrix");
	//printf("%f %f %f %f \n",B11,B12,B13,B14);
	//printf("%f %f %f %f \n",B21,B22,B23,B24);
	//printf("%f %f %f %f \n",B31,B32,B33,B34);
	//printf("%f %f %f %f \n",B41,B42,B43,B44);
	
	MatrixXd m(4,4);// Eigen/Dense
	m(0,0)=B11;m(0,1)=B12;m(0,2)=B13;m(0,3)=B14;
	m(1,0)=B21;m(1,1)=B22;m(1,2)=B23;m(1,3)=B24;
	m(2,0)=B31;m(2,1)=B32;m(2,2)=B33;m(2,3)=B34;
	m(3,0)=B41;m(3,1)=B42;m(3,2)=B43;m(3,3)=B44;
	//cout<<"symmetric matrix B\n"<<m<<endl;

   // (5) eigenvector
	SelfAdjointEigenSolver<MatrixXd> es(m);
/*	cout<<"The eigenvalues of B\n"<<es.eigenvalues()<<endl;
	cout<<"The matrix of eigenvectors V\n"<<es.eigenvectors()<<endl;

	cout<<"The first eigenvalue, lambda = "<<es.eigenvalues()[0]<<endl;
	cout<<"The first eigenvector, v = \n"<<es.eigenvectors().col(0)<<endl;
	//VectorXcd v = es.eigenvectors().col(0);
	cout<<es.eigenvectors().col(0)[3]<<endl;
	//cout<<"lambds * v\n"<<es.eigenvectors().col(0)*es.eigenvalues()[0]<<endl;
*/
	double q0,q1,q2,q3;  //,q11,q12,q13,q21,q22,q23,q31,q32,q33;
	q0 = -1 * es.eigenvectors().col(0)[0];
	q1 = -1 * es.eigenvectors().col(0)[1];
	q2 = -1 * es.eigenvectors().col(0)[2];
	q3 = -1 * es.eigenvectors().col(0)[3];
	//printf("v\n%f %f %f %f \n",q0,q1,q2,q3);
	rot_mat.push_back( 2*q0*q0 + 2*q1*q1 - 1 );
	rot_mat.push_back( 2*q1*q2 - 2*q0*q3     );
	rot_mat.push_back( 2*q1*q3 + 2*q0*q2     );
	rot_mat.push_back( 2*q1*q2 + 2*q0*q3     );
	rot_mat.push_back( 2*q0*q0 + 2*q2*q2 - 1 );
	rot_mat.push_back( 2*q2*q3 - 2*q0*q1     );
	rot_mat.push_back( 2*q1*q3 - 2*q0*q2     );
	rot_mat.push_back( 2*q2*q3 + 2*q0*q1     );
	rot_mat.push_back( 2*q0*q0 + 2*q3*q3 - 1 );

	for(int i=0;i<6;i++){
		rot_mat.push_back( com_vecs[i] );
	}

	/*puts("rotation matrix q");
	for(int i=0;i<9;i=i+3){
	   printf("%10f %10f %10f\n",rot_mat[i],rot_mat[i+1],rot_mat[i+2]);
	}*/

   double buf[3];
   for(unsigned int i=0;i<vec_tar.size();i=i+3){ // Rotate by R(q)*b
      buf[0] = rot_mat[0]*vec_tar[i]+rot_mat[1]*vec_tar[1+i]+rot_mat[2]*vec_tar[2+i];
      buf[1] = rot_mat[3]*vec_tar[i]+rot_mat[4]*vec_tar[1+i]+rot_mat[5]*vec_tar[2+i];
      buf[2] = rot_mat[6]*vec_tar[i]+rot_mat[7]*vec_tar[1+i]+rot_mat[8]*vec_tar[2+i];
      vec_tar[i] = buf[0];
      vec_tar[1+i] = buf[1];
      vec_tar[2+i] = buf[2];
   }

   return rot_mat;
}

int transfer_quat( vector<double> &x, vector<double> &y,  vector<double> &z, vector<double> &transf ){
      for(unsigned int i=0;i<x.size();i++){ // Transfer
         x[i] = x[i] - transf[0];
         y[i] = y[i] - transf[1];
         z[i] = z[i] - transf[2];
      }
   return 0;
}
int transfer_quat( vector<double> &vec, vector<double> &transf ){
   for(unsigned int i=0;i<vec.size();i=i+3){ // Transfer
         vec[i] = vec[i] - transf[0];
         vec[i+1] = vec[i+1] - transf[1];
         vec[i+2] = vec[i+2] - transf[2];
   }
   return 0; //transfer_quat();
}
int rotate_quat( vector<double> &x, vector<double> &y,  vector<double> &z, vector<double> &rot_mat ){
      double buf[3];
      for(unsigned int i=0;i<x.size();i++){ // Rotate by R(q)*b
         buf[0] = rot_mat[0]*x[i]+rot_mat[1]*y[i]+rot_mat[2]*z[i];
         buf[1] = rot_mat[3]*x[i]+rot_mat[4]*y[i]+rot_mat[5]*z[i];
         buf[2] = rot_mat[6]*x[i]+rot_mat[7]*y[i]+rot_mat[8]*z[i];
	 x[i] = buf[0];
	 y[i] = buf[1];
	 z[i] = buf[2];
      }
      return 0;
}
int rotate_quat( vector<double> &vec, vector<double> &rot_mat ){
      double buf[3];
      for(unsigned int i=0;i<vec.size();i=i+3){ // Rotate by R(q)*b
         buf[0] = rot_mat[0]*vec[i]+rot_mat[1]*vec[i+1]+rot_mat[2]*vec[i+2];
         buf[1] = rot_mat[3]*vec[i]+rot_mat[4]*vec[i+1]+rot_mat[5]*vec[i+2];
         buf[2] = rot_mat[6]*vec[i]+rot_mat[7]*vec[i+1]+rot_mat[8]*vec[i+2];
	 vec[i] = buf[0];
	 vec[i+1] = buf[1];
	 vec[i+2] = buf[2];
      }
      return 0;
}
int select_quat( pdb_nishi &pdb1, vector<double> &x, vector<double> &y,vector<double> &z, vector<double> &vec, string &rmsdatom, int i ){
	 if( pdb1.atmn[i] != "CA" && rmsdatom == "ca" )return 1;
	 if( pdb1.atmn[i] != "CA" && pdb1.atmn[i] != "N" && pdb1.atmn[i] != "C" && pdb1.atmn[i] != "O" && rmsdatom == "mainchain" )return 2;
         if( pdb1.elem[i] == "H" && rmsdatom != "all" )return 3;
         if( pdb1.resn[i] == "WAT" )return 4;
	 if( pdb1.resn[i] == "CIM" )return 5;
	 if( pdb1.resn[i] == "CIP" )return 6;
         vec.push_back(x[i]);
         vec.push_back(y[i]);
         vec.push_back(z[i]);
   return 0;
}
int select_quat( pdb_nishi &pdb1, vector<double> &vec, string &rmsdatom, int i ){
   return select_quat( pdb1, pdb1.coox, pdb1.cooy, pdb1.cooz, vec, rmsdatom, i );
}
