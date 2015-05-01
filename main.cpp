#include"nlib.h"

using namespace std;

int calc(Inp_nishi);  // calc() is empty func. (written in calc.cpp)
int quatnishi(Inp_nishi);
int pcanishi(Inp_nishi);

int main(int argc, char *argv[]){
  cout<<"Version info. pcanishi v1.2.4 \n";
// ##################### ARGUMENT HANDLING ##########################
// argv[1]: input parameter file
  if( argv[1]==NULL ){
    puts("No ARGUMEMTS");
    puts("USAGE: ./a.out (argv[1]: input parameter file)" );
    return 1;
  }
  cout<<"Your input-parameter file: "<<argv[1]<<endl;

// INPUT_PARAMETERS
   //Inp_nishi inp1( "parameter.txt" );
   Inp_nishi inp1( argv[1] );
   
// DO quatnishi
   //quatnishi(inp1);

// DO pcanishi
   int rtrn = pcanishi( inp1 );
   if(rtrn==0)cout<<"\nend of pcanishi\n";


// END
	cout<<"\n\nit took "<<(float)clock()/CLOCKS_PER_SEC<<" sec of CPU to execute this program"<<endl;
	return 0;
}
