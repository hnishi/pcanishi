/*  CLASS tra_nishi
*/

#include"nlib.h"

void tra_nishi::constructor( const char *codname, const char *pdbname, int n, string atomsel )
{
   cod_name = codname;
   pdb_name = pdbname;
	atom_sel = atomsel ;
	ifstream ifs( codname, ios::in | ios::binary);//open file, codname

	if(ifs.fail()){//error handling
		cerr<<"cannot open file; "<<codname<<endl;
		exit(1);
	}

	pdb1=new pdb_nishi(pdbname);//call class pdb_nishi from nlib.h
	//cout<<"TOTAL NUM. OF ATOMS = "<<pdb1.anum.size()<<endl;
	//cout<<"TOTAL NUM. OF ATOMS = "<<pdb1.total_atom<<endl;
	while(!ifs.eof()){// to read binary
                ifs.ignore(4);
                ifs.read( (char *) &bufint, sizeof( int ) );
                loopnum.push_back (bufint);
                ifs.read( (char *) &buffloat, 4 );
                sitime.push_back (buffloat);
                ifs.read( (char *) &buffloat, 4 );
                cputime.push_back (buffloat);
                ifs.read( (char *) &buffloat, 4 );
                totalE.push_back (buffloat);
                ifs.read( (char *) &buffloat, 4 );
                kineticE.push_back (buffloat);
                ifs.read( (char *) &buffloat, 4 );
                temp.push_back (buffloat);
                ifs.read( (char *) &buffloat, 4 );
                potent.push_back (buffloat);
                ifs.read( (char *) &buffloat, 4 );
                rmsf.push_back (buffloat);
                ifs.read( (char *) &bufint, 4 );
                num15svw.push_back (bufint);
                ifs.read( (char *) &bufint, 4 );
                num15hyd.push_back (bufint);
                ifs.read( (char *) &buffloat, 4 );
                rmsd.push_back (buffloat);
		ifs.ignore(8);

		double minx=99999,maxx=-99999,miny=99999,maxy=-99999,minz=99999,maxz=-99999;
		float bufx, bufy, bufz;   int rtrn, rtrn_c=0;   vector<double> vec;
                for(unsigned int w=0; w < pdb1->total_atom; w++){
                        ifs.read( (char *) &bufx, 4 );
			if(minx>buffloat&&pdb1->resn[w]!="WAT")minx=bufx;
			if(maxx<buffloat&&pdb1->resn[w]!="WAT")maxx=bufx;
                        ifs.read( (char *) &bufy, 4 );
			if(miny>buffloat)miny=bufy;
			if(maxy<buffloat)maxy=bufy;
                        ifs.read( (char *) &bufz, 4 );
			if(minz>buffloat)minz=bufz;
			if(maxz<buffloat)maxz=bufz;
			
                        //ifs.read( (char *) &bufx, 4 );
                        //ifs.read( (char *) &bufy, 4 );
                        //ifs.read( (char *) &bufz, 4 );
			rtrn = select_atom( *pdb1, bufx, bufy, bufz, vec, atomsel, w );
			if( rtrn != 0 ) rtrn_c++;
		}			
		length_x.push_back(maxx - minx);
		length_y.push_back(maxy - miny);
		length_z.push_back(maxz - minz);

                ifs.ignore(4);  // totally 4*11+4*4+4*3*total_atom=60+12*total_atom
                ifs.ignore(((60+12*pdb1->total_atom)*(n-1)));   //ignore for stride
        	
                for(unsigned int w=0; w < vec.size(); w=w+3){
			//cout<<"!!! vec = "<<vec[w]<<endl;
			cordx.push_back( vec[w] );
			cordy.push_back( vec[w+1] );
			cordz.push_back( vec[w+2] );
		}
		total_sel = vec.size() /3 ;
		//cout<<"!!!! rtrn_c = "<<rtrn_c<<endl;
        }
	//cout<<"!!!! total_sel = "<<total_sel<<endl;
        //total_step = loopnum.size();
        total_step = rmsd.size();
        int loopdist=loopnum[1]-loopnum[0];
        //cout<<"loopdist= "<<loopdist<<endl;
        while( loopnum[total_step-1] != ( loopnum[total_step-2] + loopdist ) ){
                total_step--;   // check whether total_step is wrong or correct
		cout<<"WARNING: total_step was decreased because the final loop-number of a trajectory was not correctly read\n";
        }
	//cout<<"!!!! total_step = "<<total_step<<endl;
        ifs.close();
}
tra_nishi::tra_nishi(const char *codname,const char *pdbname){//default constructor
	int n=1;
	constructor(codname,pdbname,n, "all");
}
tra_nishi::tra_nishi(const char *codname,const char *pdbname, int n){
	constructor(codname,pdbname,n, "all");
}
tra_nishi::tra_nishi(const char *codname,const char *pdbname, int n,string atomsel){
	constructor(codname,pdbname,n, atomsel);
}
tra_nishi::tra_nishi(const char *codname,const char *pdbname,string atomsel){
	constructor(codname,pdbname, 1 , atomsel);
}


int tra_nishi::disp_line(int n){//display one step data at internal num. of n
	printf("loopnum: %12d, sitime: %12f, cputime: %12f\n",
	loopnum[n],sitime[n],cputime[n]);
	printf("totalE: %10f, kineticE: %10f, temperature: %10f, potential: %10f\n",
	totalE[n],kineticE[n],temp[n],potent[n]);
	printf("rmsf: %12f, num15svw: %6d, num15hyd: %6d, rmsd: %12f\n",    
        //printf("rmsf: %12f, rmsd: %12f %12f %12f %12f %12f %12f %12f\n",
	rmsf[n],num15svw[n],num15hyd[n],rmsd[n]/*,zzz1[n],zzz2[n],zzz3[n],z1[n],z2[n],z3[n]*/);
	return 0;
}


int tra_nishi::write_step(const char* filename,int n){//output pdb at n step
   n = n-1;
	FILE *fout;
        if((fout = fopen(filename,"w")) == NULL ){//error handling
                printf("cannot open output file: %s\n",filename);
                exit(1);
        }
   if( (unsigned) n >= total_step ){ //added 2015.1.30
      printf("cannot write_step because (step num.) %i >= (total_step) %i \n",n,total_step);
      exit(1);
   }
   fprintf(fout,"REMARK   this file was output by tra_nishi::write_step at frame %i/%i\n", n+1,total_step);
   fprintf(fout,"REMARK   original crd file is %s\n", cod_name.c_str() );
   fprintf(fout,"REMARK   original pdb file is %s\n", pdb_name.c_str() );
        //int bfanum,bfrnum;
        //float bfcoox,bfcooy,bfcooz,bfoccu,bftemf;
        //char bfreco[6],bfatmn[6],bfresn[5],bfchai[2],bfelem[5];
	char bfch;   int rtrn_sel ,ccc=0;   vector<double> vec;
	//cout<<"!!! atomsel = "<<atom_sel<<endl;
        for( unsigned int j=0; j < pdb1->total_atom; j++ ){
		rtrn_sel = select_atom( *pdb1, vec, atom_sel, j );
		if( rtrn_sel != 0 ) continue;
                
                strcpy(bfreco,pdb1->reco[j].c_str());
                bfanum=pdb1->anum[j];
                strcpy(bfatmn,pdb1->atmn[j].c_str());
                strcpy(bfresn,pdb1->resn[j].c_str());
                //strcpy(bfchai,pdb1->chai[j].c_str());
                bfch=pdb1->chai[j].c_str()[0];
		bfrnum=pdb1->rnum[j];
                bfcoox=cordx[n * total_sel + ccc];
                bfcooy=cordy[n * total_sel + ccc];
                bfcooz=cordz[n * total_sel + ccc];
                bfoccu=pdb1->occu[j];
                bftemf=pdb1->temf[j];
                strcpy(bfelem,pdb1->elem[j].c_str());
                if( strlen(bfatmn)==4 ){
                //fprintf(fout,"%-6s%5d%5s%4s%c%4d%12.3f%8.3f%8.3f%6.2f%8.4f%12s\n",
                //fprintf(fout,"%-6s%5d%5s %-4s%5d%12.3f%8.3f%8.3f%6.2f%6.2f\n",
                fprintf(fout,WRITE_FORMAT_1,
                bfreco,bfanum,bfatmn,bfresn,bfch,bfrnum,
                bfcoox,bfcooy,bfcooz,bfoccu,bftemf,bfelem);
                }
                else{
                //fprintf(fout,"%-6s%5d  %-3s%4s%c%4d%12.3f%8.3f%8.3f%6.2f%8.4f%12s\n",
                //fprintf(fout,"%-6s%5d  %-3s %-4s%5d%12.3f%8.3f%8.3f%6.2f%6.2f\n",
                fprintf(fout,WRITE_FORMAT_2,
                bfreco,bfanum,bfatmn,bfresn,bfch,bfrnum,
                bfcoox,bfcooy,bfcooz,bfoccu,bftemf,bfelem);
                }
		ccc++;
        }
        fclose(fout);
	return 0;
}

int tra_nishi::write_cod(const char* filename,int stride){//output trajectory in ASCII code
        FILE *fout;
        if((fout = fopen(filename,"w")) == NULL ){
                printf("cannot open output file: %s\n",filename);
                exit(1);
        }

        //int bfanum,bfrnum;
        //float bfcoox,bfcooy,bfcooz,bfoccu,bftemf;
        //char bfreco[6],bfatmn[6],bfresn[5],bfchai[2],bfelem[5];
        for(unsigned int n=0;n<total_step;n=n+stride){
	char bfch;   int rtrn_sel , ccc=0 ;   vector<double> vec;
        //fprintf(fout,"MODEL %d\n",n+1);
        fprintf(fout,"MODEL %d\n",loopnum[n]);
        for( unsigned int j=0; j < pdb1->total_atom; j++ ){
		rtrn_sel = select_atom( *pdb1, vec, atom_sel, j );
		if( rtrn_sel != 0 ) continue;

                strcpy(bfreco,pdb1->reco[j].c_str());
                bfanum=pdb1->anum[j];
                strcpy(bfatmn,pdb1->atmn[j].c_str());
                strcpy(bfresn,pdb1->resn[j].c_str());
                //strcpy(bfchai,pdb1->chai[j].c_str());
                bfch=pdb1->chai[j].c_str()[0];
		bfrnum=pdb1->rnum[j];
                bfcoox=cordx[n * total_sel+ ccc];
                bfcooy=cordy[n * total_sel + ccc];
                bfcooz=cordz[n * total_sel + ccc];
                bfoccu=pdb1->occu[j];
                bftemf=pdb1->temf[j];
                strcpy(bfelem,pdb1->elem[j].c_str());
                if( strlen(bfatmn)==4 ){
                //fprintf(fout,"%-6s%5d%5s%4s%2s%4d%12.3f%8.3f%8.3f%6.2f%8.4f%12s\n",
                //fprintf(fout,"%-6s%5d%5s %-4s%5d%12.3f%8.3f%8.3f%6.2f%6.2f\n",
                fprintf(fout,WRITE_FORMAT_1,
                bfreco,bfanum,bfatmn,bfresn,bfch,bfrnum,
                bfcoox,bfcooy,bfcooz,bfoccu,bftemf,bfelem);
                }
                else{
                //fprintf(fout,"%-6s%5d  %-3s%4s%2s%4d%12.3f%8.3f%8.3f%6.2f%8.4f%12s\n",
                //fprintf(fout,"%-6s%5d  %-3s %-4s%5d%12.3f%8.3f%8.3f%6.2f%6.2f\n",
                fprintf(fout,WRITE_FORMAT_2,
                bfreco,bfanum,bfatmn,bfresn,bfch,bfrnum,
                bfcoox,bfcooy,bfcooz,bfoccu,bftemf,bfelem);
                }
		ccc++;
        }
	fprintf(fout,"ENDMDL\n");
        }
        fclose(fout);
        return 0;
}
int tra_nishi::write_cod(const char* filename){
  return write_cod(filename, 1);
}

int tra_nishi::fix_step(const char *filename, int n,float fxcell,float fycell,float fzcell){//not yet implemented
	
	return 0;
}


int tra_nishi::fix_cod(float fxcell,float fycell,float fzcell){//replicate atoms because of periodic boundary condition
	float rdiff;
	int imove;

        for(unsigned int n=0;n<total_step;n++){
                //fprintf(fout,"MODEL %d\n",n+1);
                //printf("%i: %f %f %f \n",n,fxcell,fycell,fzcell);
        int qqq=0, ccc=0, rtrn_sel;
	vector<double> vec;
        for(unsigned int j = n * pdb1->total_atom +1; j < (n+1) * pdb1->total_atom;j++){
                qqq++;

		rtrn_sel = select_atom( *pdb1, vec, atom_sel, qqq );
		if( rtrn_sel != 0 ) continue;
		ccc++;

                if( pdb1->resn[qqq] == "WAT" )continue;
                if( pdb1->resn[qqq] == "CIM" )continue;
                if( pdb1->resn[qqq] == "CIP" )continue;
                //printf("%s %s \n",pdb1->chai[qqq].c_str(),pdb1->chai[qqq-1].c_str());
                //cout<<pdb1->chai[qqq].c_str()<<" "<<pdb1->chai[qqq-1].c_str()<<endl;
                if( pdb1->chai[qqq] != pdb1->chai[qqq-1])continue;
                rdiff = cordx[n*total_sel + ccc]-cordx[n*total_sel + ccc-1];
                //cout<<"frame "<<n<<", atom "<<qqq<<", rdiff of x = "<<rdiff<<endl;
                //rdiff = cordx[j]-cordx[j-pdb1->total_atom];
                if(rdiff>=0){ imove = int(rdiff/fxcell+0.5); }
                else{ imove = int(rdiff/fxcell-0.5); }
                //cout<<"imove of x = "<<imove<<endl;
                //cout<<"cordx before = "<<cordx[[n*total_sel + ccc]]<<endl;
                cordx[n*total_sel + ccc] = cordx[n*total_sel + ccc] - imove*fxcell;
                //cout<<"cordx after = "<<cordx[n*total_sel + ccc]<<endl;
                //if(rdiff<0&&rdiff>fxcell/2)cordx[n*total_sel + ccc] = cordx[n*total_sel + ccc] + fxcell;
                rdiff = cordy[n*total_sel + ccc]-cordy[n*total_sel + ccc-1];
                //rdiff = cordy[n*total_sel + ccc]-cordy[n*total_sel + ccc-pdb1->total_atom];
                if(rdiff>=0) imove = int(rdiff/fycell+0.5);
                else imove = int(rdiff/fycell-0.5);
                cordy[n*total_sel + ccc] = cordy[n*total_sel + ccc] - imove*fycell;
                //if(rdiff<0&&rdiff>fycell/2)cordy[n*total_sel + ccc] = cordy[n*total_sel + ccc] + fycell;
                rdiff = cordz[n*total_sel + ccc]-cordz[n*total_sel + ccc-1];
                //rdiff = cordz[n*total_sel + ccc]-cordz[n*total_sel + ccc-pdb1->total_atom];
                if(rdiff>=0) imove = int(rdiff/fzcell+0.5);
                else imove = int(rdiff/fzcell-0.5);
                cordz[n*total_sel + ccc] = cordz[n*total_sel + ccc] - imove*fzcell;
                //if(rdiff<0&&rdiff>fzcell/2)cordz[n*total_sel + ccc] = cordz[n*total_sel + ccc] + fzcell;
                //cout<<j<<" ";
        }
        }
	//write_cod(filename,stride);
//### the following comment is a version of not changing cords ###
        /*FILE *fout;
        if((fout = fopen(filename,"w")) == NULL ){
                printf("cannot open output file: %s\n",filename);
                exit(1);
        }
	float bfx,bfy,bfz;

        for(unsigned int n=0;n<total_step;n=n+stride){
        //fprintf(fout,"MODEL %d\n",n+1);
        fprintf(fout,"MODEL %d\n",loopnum[n]);
        for(unsigned int j = n * pdb1->total_atom +1; j < (n+1) * pdb1->total_atom;j++){
        }
        for(unsigned int j=0;j<pdb1->total_atom;j++){
                rdiff = cordx[n*pdb1->total_atom+1+j]-cordx[j-1];
                //rdiff = cordx[n*pdb1->total_atom+1+j]-cordx[j-pdb1->total_atom];
                if(rdiff>=0) imove = int(rdiff/fxcell+0.5);
                else imove = int(rdiff/fxcell-0.5);
                bfx = cordx[n*pdb1->total_atom+1+j] - imove*fxcell;
                rdiff = cordy[n*pdb1->total_atom+1+j]-cordy[j-1];
                //rdiff = cordy[n*pdb1->total_atom+1+j]-cordy[j-pdb1->total_atom];
                if(rdiff>=0) imove = int(rdiff/fycell+0.5);
                else imove = int(rdiff/fycell-0.5);
                bfy = cordy[n*pdb1->total_atom+1+j] - imove*fycell;
                rdiff = cordz[n*pdb1->total_atom+1+j]-cordz[j-1];
                //rdiff = cordz[n*pdb1->total_atom+1+j]-cordz[j-pdb1->total_atom];
                if(rdiff>=0) imove = int(rdiff/fzcell+0.5);
                else imove = int(rdiff/fzcell-0.5);
                bfz = cordz[n*pdb1->total_atom+1+j] - imove*fzcell;

                strcpy(bfreco,pdb1->reco[j].c_str());
                bfanum=pdb1->anum[j];
                strcpy(bfatmn,pdb1->atmn[j].c_str());
                strcpy(bfresn,pdb1->resn[j].c_str());
                strcpy(bfchai,pdb1->chai[j].c_str());
                bfrnum=pdb1->rnum[j];
                bfcoox=bfx;
                bfcooy=bfy;
                bfcooz=bfz;
                bfoccu=pdb1->occu[j];
                bftemf=pdb1->temf[j];
                strcpy(bfelem,pdb1->elem[j].c_str());
                if( strlen(bfatmn)==4 ){
                //fprintf(fout,"%-6s%5d%5s%4s%2s%4d%12.3f%8.3f%8.3f%6.2f%8.4f%12s\n",
                //fprintf(fout,"%-6s%5d%5s %-4s%5d%12.3f%8.3f%8.3f%6.2f%6.2f\n",
                fprintf(fout,WRITE_FORMAT_1,
                bfreco,bfanum,bfatmn,bfresn,bfchai,bfrnum,
                bfcoox,bfcooy,bfcooz,bfoccu,bftemf,bfelem);
                }
                else{
                //fprintf(fout,"%-6s%5d  %-3s%4s%2s%4d%12.3f%8.3f%8.3f%6.2f%8.4f%12s\n",
                //fprintf(fout,"%-6s%5d  %-3s %-4s%5d%12.3f%8.3f%8.3f%6.2f%6.2f\n",
                fprintf(fout,WRITE_FORMAT_2,
                bfreco,bfanum,bfatmn,bfresn,bfchai,bfrnum,
                bfcoox,bfcooy,bfcooz,bfoccu,bftemf,bfelem);
                }
        }
        fprintf(fout,"ENDMDL\n");
        }
        fclose(fout);*/
	return 0;
}

int tra_nishi::fix_cod_npt(){
        float rdiff;
        int imove;
	float fxcell,fycell,fzcell;
        for(unsigned int n=0;n<total_step;n++){
                //fprintf(fout,"MODEL %d\n",n+1);
		fxcell = length_x[n];
		fycell = length_y[n];
		fzcell = length_z[n];
		//printf("%i: %f %f %f \n",n,fxcell,fycell,fzcell);
	int qqq=0;
        int ccc=0, rtrn_sel;  vector<double> vec;
        for(unsigned int j = n * pdb1->total_atom +1; j < (n+1) * pdb1->total_atom;j++){
		qqq++;
		rtrn_sel = select_atom( *pdb1, vec, atom_sel, qqq );
		if( rtrn_sel != 0 ) continue;
		ccc++;
                if( pdb1->resn[qqq] == "WAT" )continue;
                if( pdb1->resn[qqq] == "CIM" )continue;
                if( pdb1->resn[qqq] == "CIP" )continue;
		//printf("%s %s \n",pdb1->chai[qqq].c_str(),pdb1->chai[qqq-1].c_str());
		//cout<<pdb1->chai[qqq].c_str()<<" "<<pdb1->chai[qqq-1].c_str()<<endl;	
		if( pdb1->chai[qqq] != pdb1->chai[qqq-1])continue;
                rdiff = cordx[n*total_sel + ccc]-cordx[n*total_sel + ccc-1];
		//cout<<"rdiff of x = "<<rdiff<<endl;
                //rdiff = cordx[n*total_sel + ccc]-cordx[n*total_sel + ccc-pdb1->total_atom];
                if(rdiff>=0){ imove = int(rdiff/fxcell+0.5); }
                else{ imove = int(rdiff/fxcell-0.5); }
		//cout<<"imove of x = "<<imove<<endl;
                //cout<<"cordx before = "<<cordx[n*total_sel + ccc]<<endl;
		cordx[n*total_sel + ccc] = cordx[n*total_sel + ccc] - imove*fxcell;
		//cout<<"cordx after = "<<cordx[n*total_sel + ccc]<<endl;
                //if(rdiff<0&&rdiff>fxcell/2)cordx[n*total_sel + ccc] = cordx[n*total_sel + ccc] + fxcell;
		rdiff = cordy[n*total_sel + ccc]-cordy[n*total_sel + ccc-1];
                //rdiff = cordy[n*total_sel + ccc]-cordy[n*total_sel + ccc-pdb1->total_atom];
                if(rdiff>=0) imove = int(rdiff/fycell+0.5);
                else imove = int(rdiff/fycell-0.5);
                cordy[n*total_sel + ccc] = cordy[n*total_sel + ccc] - imove*fycell;
                //if(rdiff<0&&rdiff>fycell/2)cordy[n*total_sel + ccc] = cordy[n*total_sel + ccc] + fycell;
                rdiff = cordz[n*total_sel + ccc]-cordz[n*total_sel + ccc-1];
                //rdiff = cordz[n*total_sel + ccc]-cordz[n*total_sel + ccc-pdb1->total_atom];
                if(rdiff>=0) imove = int(rdiff/fzcell+0.5);
                else imove = int(rdiff/fzcell-0.5);
                cordz[n*total_sel + ccc] = cordz[n*total_sel + ccc] - imove*fzcell;
                //if(rdiff<0&&rdiff>fzcell/2)cordz[n*total_sel + ccc] = cordz[n*total_sel + ccc] + fzcell;
                //cout<<j<<" ";
        }
        }
	return 0;
}

tra_nishi::~tra_nishi(){//deconstructor
	delete pdb1;//free dynamic memory
}


int select_atom( pdb_nishi &pdb1, double x, double y, double z, vector<double> &vec, string &atomsel, int i ){
	 if( pdb1.atmn[i] != "CA" && atomsel == "ca" )return 1;
	 if( pdb1.atmn[i] != "CA" && pdb1.atmn[i] != "N" && pdb1.atmn[i] != "C" && pdb1.atmn[i] != "O" && atomsel == "mainchain" )return 2;
         if( pdb1.elem[i] == "H" && atomsel == "heavy" )return 3;
         if( pdb1.resn[i] == "WAT" && atomsel != "all" )return 4;
	 if( pdb1.resn[i] == "CIM" && atomsel != "all" )return 5;
	 if( pdb1.resn[i] == "CIP" && atomsel != "all" )return 6;
         vec.push_back(x);
         vec.push_back(y);
         vec.push_back(z);
   return 0;
}
int select_atom( pdb_nishi &pdb1, vector<double> &vec, string &atomsel, int i ){
   return select_atom( pdb1, pdb1.coox[i], pdb1.cooy[i], pdb1.cooz[i], vec, atomsel, i );
}

int search_sel( pdb_nishi &pdb1, string chai, int resn, string atmn, string atomsel){
  int total_sel = 0, intra_num = -1, check1=0;
  vector<double> vec;
  for(unsigned int w=0; w < pdb1.total_atom; w++){
    int rtrn = select_atom( pdb1, vec, atomsel, w );
    if( rtrn == 0 ){
      //cout<<pdb1.chai[w]<<" "<<chai<<" "<<pdb1.rnum[w]<<" "<<resn<<" "<<pdb1.atmn[w]<<" "<<atmn<<endl;
      if( pdb1.chai[w] == chai && pdb1.rnum[w] == resn && pdb1.atmn[w] == atmn ){
        intra_num = total_sel;
        check1++;
      }
      total_sel++; 
    }
  }
  if(check1 > 1){
    cout<<"WARNING in tranishi.cpp: search_sel() detected more than 1 selected atom"<<endl; 
  }
  //cout<<"DEBUG: total_sel = "<<total_sel<<endl;
  return intra_num;
}
