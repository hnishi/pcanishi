#include"nlib.h"

//#define SCAN_FORMAT "%s%d%s%s%s%d%f%f%f%f%f%s" // no missing
//#define WRITE_FORMAT_1 "%-6s%5d%5s %-4s%s%4d%12.3f%8.3f%8.3f%6.2f%6.2f%12s\n"
//#define WRITE_FORMAT_2 "%-6s%5d  %-3s %-4s%s%4d%12.3f%8.3f%8.3f%6.2f%6.2f%12s\n"
//#define SCAN_FORMAT "%s%d%s%s%d%f%f%f%f%f" // missing chain element
//#define WRITE_FORMAT_1 "%-6s%5d%5s %-4s%5d%12.3f%8.3f%8.3f%6.2f%6.2f\n"
//#define WRITE_FORMAT_2 "%-6s%5d  %-3s %-4s%5d%12.3f%8.3f%8.3f%6.2f%6.2f\n"


// #########################################################################
//                           CLASS: pdb_nishi
// #########################################################################
// include nlib.h

/*
class pdb_nishi{
private:
public:
	vector<int> anum,rnum;
	vector<double> coox,cooy,cooz,occu,temf;
	vector<string> reco,atmn,resn,chai,elem;

	pdb_nishi(); // constructer 1 for error
	pdb_nishi(char *pdbname); // constructer for input of pdb
	void disp_line(int n);
};
*/

//        int bfanum,bfrnum;
//        float bfcoox,bfcooy,bfcooz,bfoccu,bftemf;
//        char bfreco[8],bfatmn[5],bfresn[5],bfchai[5],bfelem[5];


// #########################################################################
//                        CLASS FUNCTION constructor
// #########################################################################
pdb_nishi::pdb_nishi(){ //for error
	cout<<"cannot read pdb file\n";
	exit(1);
}
pdb_nishi::pdb_nishi(const char *pdbname){
	pdb_name = pdbname;  //public variable of class pdb_nishi
	char buf1[100];
	FILE *fin;

	if((fin = fopen(pdbname,"r")) == NULL ){
		printf("cannot open input pdb: %s\n",pdbname);
		exit(1);
	}
	//printf("successfully open input pdb: %s\n",pdbname);

	int i=0;
	char buf2[9];
	while( fgets(buf1,sizeof(buf1),fin)!=NULL ){
		if( !((strncmp(buf1,"ATOM",4)==0) || (strncmp(buf1,"HETATM",6)==0) )){
			continue;
		}
		//sscanf(buf1,"%s%d%s%s%d%f%f%f%f%f",
		//sscanf(buf1,SCAN_FORMAT,
		//bfreco,&bfanum,bfatmn,bfresn,bfchai,&bfrnum,&bfcoox,&bfcooy,&bfcooz,&bfoccu,&bftemf,bfelem);
		//strncpy(bfreco,&buf1[0],6);
		strncpy(buf2,&buf1[0],6); buf2[6]='\0';
		//cout<<"buf2: "<<buf2<<endl;
		//strcpy(bfreco,buf2);
		//cout<<"bfreco: "<<bfreco<<endl;
		sscanf(buf2,"%s", bfreco);
		strncpy(buf2,&buf1[6],5); buf2[5]='\0';
		sscanf(buf2,"%d",&bfanum);
		strncpy(buf2, &buf1[12], 4); buf2[4] = '\0';
		sscanf(buf2, "%s", bfatmn);
		//strncpy(atom->altLoc, &textbuf[16], 1);
		//strncpy(bfresn, &buf1[17], 3);
		strncpy(buf2, &buf1[17], 4); buf2[4]='\0';
		sscanf(buf2,"%s",bfresn);
		//strncpy(bfchai, &buf1[21], 1);
		//if(bfchai) sscanf(buf2,"%s",bfchai);
		//else bfchai[0] = ' ';
		strncpy(buf2, &buf1[21], 1); buf2[1]='\0';
		//if(*buf2) cout<<buf2<<"ARUYO!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
		sscanf(buf2,"%c",bfchai);
		//if(bfchai[0]=='\0') cout<<"NULL\n";
		//else cout<<"bfchai= "<<bfchai<<"  ,endl\n";
		//bfchai[1]=' ';
		strncpy(buf2, &buf1[22], 4); buf2[4] = '\0';
		sscanf(buf2, "%d", &bfrnum);
		//strncpy(atom->iCode, &textbuf[26], 1);
		strncpy(buf2, &buf1[30], 8); buf2[8] = '\0';
		sscanf(buf2, "%f", &bfcoox);
		strncpy(buf2, &buf1[38], 8); buf2[8] = '\0';
		sscanf(buf2, "%f", &bfcooy);
		strncpy(buf2, &buf1[46], 8); buf2[8] = '\0';
		sscanf(buf2, "%f", &bfcooz);
                strncpy(buf2, &buf1[54], 6); buf2[6] = '\0';
                sscanf(buf2, "%f", &bfoccu);
                strncpy(buf2, &buf1[60], 6); buf2[6] = '\0';
                sscanf(buf2, "%f", &bftemf);

                //cout<<"rtrn_elem = "<<strncpy(buf2, &buf1[76], 3)<<endl; buf2[3] = '\0';  //could not see what return-value
		string str_buf = buf1;  //v1.1.0
		//cout<<str_buf<<", size = "<<str_buf.size()<<endl;
		//cout<<"size of buf1 = "<<sizeof(buf1)<<endl;
		//cout<<"buf1[77] = "<<buf1[77]<<endl;
                if( str_buf.size() < 74 ){  //bfelem[0] = ' ';
		  if( strcmp( bfresn,"CIP" ) == 0 ){
		    strcpy( bfelem, "NA" );
		  }
		  else if( strcmp( bfresn,"CIM" ) == 0 ){
		    strcpy( bfelem, "CL" );
		  }
		  else if( strcmp( bfatmn,"FE" ) == 0 ){
		    strcpy( bfelem, "FE" );
		  }
		  else if( strcmp( bfatmn,"MG" ) == 0 ){
		    strcpy( bfelem, "MG" );
		  }
		  else if( strcmp( bfatmn,"CU" ) == 0 ){
		    strcpy( bfelem, "CU" );
		  }
		  else if( strcmp( bfatmn,"AG" ) == 0 ){
		    strcpy( bfelem, "AG" );
		  }
		  else if( strcmp( bfatmn,"ZN" ) == 0 ){
		    strcpy( bfelem, "ZN" );
		  }
		  else if( strlen(bfatmn)==4 ){  //stable expression for elemnt
                    strncpy(buf2, &buf1[12], 1); buf2[1]='\0';
		    sscanf(buf2, "%s", bfelem);
		  }
                  else{
		    strncpy(buf2, &buf1[13], 1); buf2[1]='\0';
		    sscanf(buf2, "%s", bfelem);
		  }
		}
		else{ 
                  strncpy(buf2, &buf1[76], 3); buf2[3] = '\0';
		  buf2[3] = '\0';
                  sscanf(buf2, "%s", bfelem);
                }
		/*if( strlen(bfatmn)==4 ){  //stable expression for elemnt
                strncpy(buf2, &buf1[12], 1); buf2[1]='\0';
		sscanf(buf2, "%s", bfelem);
		}
                else{
		strncpy(buf2, &buf1[13], 1); buf2[1]='\0';
		sscanf(buf2, "%s", bfelem);
		}*/
		//cout<<bfanum<<" chai: "<<bfchai<<endl;
                //cout<<bfanum<<" reco: "<<bfreco<<endl;
		/*printf("%-6s%5d%5s %-4s%s%4d%12.3f%8.3f%8.3f%6.2f%6.2f%12s\n",
                bfreco,bfanum,bfatmn,bfresn,bfchai,bfrnum,
                bfcoox,bfcooy,bfcooz,bfoccu,bftemf,bfelem);
		*/
		reco.push_back(bfreco);
                anum.push_back(bfanum);
                atmn.push_back(bfatmn);
                resn.push_back(bfresn);
                chai.push_back(bfchai);
		rnum.push_back(bfrnum);
                coox.push_back(bfcoox);
                cooy.push_back(bfcooy);
                cooz.push_back(bfcooz);
                occu.push_back(bfoccu);
                temf.push_back(bftemf);
		elem.push_back(bfelem);

		if((i!=0) && (rnum[i-1]!=rnum[i])){
			resi_mark.push_back(i-1);
		}
		i++;
	}
	resi_mark.push_back(i-1);//final internal number
	//cout<<"TOTAL NUM. OF ATOM: "<<anum.size()<<endl;

	total_atom = anum.size();
	total_residue = resi_mark.size();	
	fclose(fin);
}

// #########################################################################
//                           CLASS FUNCTION disp_line
// #########################################################################
int pdb_nishi::disp_line(int j){  // "int j" should be internal number
   /* 2015.01.08 added 
   ERROR CHECK        */
   if( j >= (signed)total_atom || j<0 ){
      cout<<"WARNING: pdb_nishi::disp_line ; j >= total_atom or j < 0\n";
      cout<<" intra num. j = "<<j<<",   total_atom = "<<total_atom<<endl;
      return -1;
   }
	char bfch;
	/*char *cstr1 = new char[reco[j].length() + 1];
	strcpy(cstr1, reco[j].c_str());
	delete [] cstr1;
        */strcpy(bfreco,reco[j].c_str());
        bfanum=anum[j];
        /*char *cstr2 = new char[atmn[j].length() + 1];
        strcpy(cstr2, atmn[j].c_str());
        delete [] cstr2;
        */strcpy(bfatmn,atmn[j].c_str());
        /*char *cstr3 = new char[resn[j].length() + 1];
        strcpy(cstr3, resn[j].c_str());
        delete [] cstr3;
        */strcpy(bfresn,resn[j].c_str());
        /*char *cstr4 = new char[chai[j].length() + 1];
        strcpy(cstr4, chai[j].c_str());
        delete [] cstr4;
	*///strcpy(bfchai,chai[j].c_str());
        bfch=chai[j].c_str()[0];
	bfrnum=rnum[j];
        bfcoox=coox[j];
        bfcooy=cooy[j];
        bfcooz=cooz[j];
        bfoccu=occu[j];
        bftemf=temf[j];
        /*char *cstr5 = new char[elem[j].length() + 1];
        strcpy(cstr5, elem[j].c_str());
        delete [] cstr5;
	*/strcpy(bfelem,elem[j].c_str());
	if( strlen(bfatmn)==4 ){
                //fprintf(fout,"%-6s%5d%5s%4s%2s%4d%12.3f%8.3f%8.3f%6.2f%8.4f%12s\n",
                //printf("%-6s%5d%5s %-4s%5d%12.3f%8.3f%8.3f%6.2f%6.2f\n",
                printf(WRITE_FORMAT_1,
		bfreco,bfanum,bfatmn,bfresn,bfch,bfrnum,
                bfcoox,bfcooy,bfcooz,bfoccu,bftemf,bfelem);
                }
        else{
                //fprintf(fout,"%-6s%5d  %-3s%4s%2s%4d%12.3f%8.3f%8.3f%6.2f%8.4f%12s\n",
                //printf("%-6s%5d  %-3s %-4s%5d%12.3f%8.3f%8.3f%6.2f%6.2f\n",
                printf(WRITE_FORMAT_2,
		bfreco,bfanum,bfatmn,bfresn,bfch,bfrnum,
                bfcoox,bfcooy,bfcooz,bfoccu,bftemf,bfelem);
        }
	return 0;
}

// #########################################################################
//                           CLASS FUNCTION write_pdb
// #########################################################################
int pdb_nishi::write_pdb(const char* filename, char mode){
	FILE *fout;
   if( mode == 'w'){
	if((fout = fopen(filename,"w")) == NULL ){
		printf("cannot open output file: %s\n",filename);
		exit(1);
	}
   }
   else if( mode == 'a' ){
	if((fout = fopen(filename,"a")) == NULL ){
		printf("cannot open output file: %s\n",filename);
		exit(1);
	}
   }
   else{
      cerr<<"in pdb_nishi::write_pdb, failed to read mode \n";
      exit(1);
   }

   fprintf(fout,"REMARK  original pdb is %s \n",pdb_name.c_str() );
	char bfch;
	for(unsigned int j=0;j<anum.size();j++){
		strcpy(bfreco,reco[j].c_str());
                bfanum=anum[j];
                strcpy(bfatmn,atmn[j].c_str());
                strcpy(bfresn,resn[j].c_str());
		//strcpy(bfchai,chai[j].c_str());
                bfch=chai[j].c_str()[0];
		bfrnum=rnum[j];
                bfcoox=coox[j];
                bfcooy=cooy[j];
                bfcooz=cooz[j];
                bfoccu=occu[j];
                bftemf=temf[j];
		strcpy(bfelem,elem[j].c_str());
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
	}
        fprintf(fout,"END\n");
	fclose(fout);
	return 0;
}
int pdb_nishi::write_pdb(const char* filename){
   return write_pdb( filename, 'w' );
}

// #########################################################################
//                           CLASS FUNCTION fix_step
// #########################################################################
int pdb_nishi::fix_step(const char *filename,float fxcell,float fycell,float fzcell){
	float rdiff;
	int imove;
        for(unsigned int j = 1; j < total_atom;j++){
                rdiff = coox[j]-coox[j-1];
                if(rdiff>=0)
			imove = (signed int)(rdiff/fxcell+0.5);
                else
                        imove = (signed int)(rdiff/fxcell-0.5);
		coox[j] = coox[j] - imove*fxcell;
                rdiff = cooy[j]-cooy[j-1];
		if(rdiff>=0)
                	imove = (signed int)(rdiff/fycell+0.5);
		else
			imove = (signed int)(rdiff/fycell-0.5);
                //cout<<"rdiff/fycell: "<<rdiff/fycell;
		//cout<<", rdiff_y: "<<rdiff<<", y: "<<cooy[j]<<", imove: "<<imove<<endl;
                cooy[j] = cooy[j] - imove*fycell;
                rdiff = cooz[j]-cooz[j-1];
                if(rdiff>=0)
                	imove = int(rdiff/fzcell+0.5);
                else
                        imove = (signed int)(rdiff/fzcell-0.5);
                cooz[j] = cooz[j] - imove*fzcell;
                //cout<<j<<" ";
        }
	write_pdb(filename);
	return 0;
}

// #########################################################################
//                           CLASS FUNCTION search_n
// #########################################################################
int pdb_nishi::search_n(char a,int aa) // search internal num. of chain ID "a" and residue num. "aa"
{
        for(unsigned int j=0;j<total_atom;j++){
		if( chai[j].c_str()[0] == a )
		{
			if( rnum[j] == aa ){
				return j;
			}
		}
        }
	cout<<"ERROR: CHAIN ID "<<a<<" AND RESIDUE NUM. "<<aa<<" DOES NOT EXIST.\n";
	return -1;
}
// #########################################################################
//                           CLASS FUNCTION search_n_end
// #########################################################################
int pdb_nishi::search_n_end(char a,int aa) // search internal num. of chain ID "a" and residue num. "aa"
{
        for(unsigned int j=0;j<total_atom;j++){
		if( chai[j].c_str()[0] == a )
		{
			if( rnum[j] == aa && rnum[j+1] != aa){
				return j;
			}
		}
        }
	cout<<"WARNING: pdb_nishi::search_n_end( char CHAIN-ID, int RESIDUE-NUMBER ) \n";
	cout<<"WARNING: CHAIN ID "<<a<<" AND RESIDUE NUM. "<<aa<<" DOES NOT EXIST.\n";
	return -1;
}

// #########################################################################
//                     CLASS FUNCTION center_r and ave()
// #########################################################################
/*inline double ave(float *all, int i){         // to determine average
        float total=0;
        int u;
                for(u=0;u<i;u++){
                total+=all[u];
                }
return total/i;
}*/

int pdb_nishi::center_r(){
        float *bufx,*bufy,*bufz;
        bufx=new float[100];bufy=new float[100];bufz=new float[100];
        int jj=0;
        for(int k= 0; k<resi_mark[0]; k++ ){
                //strcpy(bfatmn,atmn[k].c_str());
                if(! (strncmp( atmn[k].c_str(),"N",4 ) == 0  ||
                        strncmp( atmn[k].c_str(),"CA",4 ) == 0 ||
                        strncmp( atmn[k].c_str(),"C",4 ) == 0  ||
                        strncmp( atmn[k].c_str(),"O",4 ) == 0)  ){
                        cout<<"anum "<<anum[k]<<", atmn "<<atmn[k]<<", resn "<<resn[k] <<endl;
                        //cout <<"rnum "<<rnum[k]<<endl;
                        //cout << "bufx " << bufx[k-resnum[j-1]] << endl;
                        bufx[jj] = coox[k];
                        bufy[jj] = cooy[k];
                        bufz[jj] = cooz[k];
                        cout<<"jj= "<<jj<<", bufz "<<bufz[jj]<<endl;
                        jj++;
                }
        }
        cout<<"total jj= "<<jj<<endl;
        comx_r.push_back( ave(bufx,jj) );
        comy_r.push_back( ave(bufy,jj) );
        comz_r.push_back( ave(bufz,jj) );
        cout<<"comx0= "<<comx_r[0]<<", comy0= "<<comy_r[0]<<", comz0= "<<comz_r[0]<<endl;
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
        cout<<"total_com_r: "<<total_com_r<<", total_residue: "<<total_residue<<endl;
        if( total_com_r != total_residue){
                printf("ERROR: dist_nishi::dist_resi, disagreement\n");
                cout<<"total_com_r: "<<total_com_r<<", total_residue: "<<total_residue<<endl;
                exit(1);
        }
        delete [] bufx;delete [] bufy;delete [] bufz;

	return 0;
}

