#include"nlib.h"

int calc(){
	cout<<"Program calc() starts!!!\n";
	
// OUTPUT TRAJECTORIES IN MD0
	/*string codname1, pdbname1, outname1;
	char bfcod[200], bfpdb[200], bfout[200];
	
	char buf[3];
	for(int i=1;i<9;i++){
	pdbname1="/Users/nishigami/AbModeling/02mcmd_ab10/cargo_md0/10_npt_43000.pdb";
	codname1="/Users/nishigami/AbModeling/02mcmd_ab10/cargo_md0/no/md0.crd";
	outname1="/Users/nishigami/AbModeling/02mcmd_ab10/cargo_md0/cod150_no.pdb";
	//outname1="cod150_no.pdb";
	//cout<<codname1<<endl;
	//itoa(i,buf,10);
	sprintf(buf,"%d",i);
	codname1.insert(52,buf); //insert "1" behind charcter at 52 in codname1; that is "o" 
	outname1.insert(59,buf);
	//outname1.insert(9,"1");
	
	cout<<codname1<<endl;
	cout<<outname1<<endl;
	strcpy(bfpdb,pdbname1.c_str());
	strcpy(bfcod,codname1.c_str());
	strcpy(bfout,outname1.c_str());
	//cout<<bfpdb<<" "<<bfcod<<endl;
	//tra_nishi tra1(codname1,pdbname1);
	tra_nishi tra1(bfcod,bfpdb);
	tra1.write_cod(bfout,10);
	//tra1.write_cod("zzz.pdb",10);
	}*/
	return 0;
}
