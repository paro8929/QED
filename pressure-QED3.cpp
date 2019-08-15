#include <QED3.h>


using namespace std;


//Run parameters which you might want to change
const int N0=400;
const int N1=202;
const int write_to_file=0; //whether or not to (over-)write the results file
const double IR=0.000001; //IR cutoff for spatial integration
const double alpha_default=1.0; //default value of interaction if not writing to file
//

int ABORT=0;

double ahat;  // is dimensionless coupling ahat=e^2 Nf/(4 Pi T)

double wn[N0];
double wnweights[N0];
double p0vals[N1];
double p0weights[N1];
double modp0weights[N1];
double s2p0vals[N1],c2p0vals[N1],tp0vals[N1];



std::complex<double> I (0,1);


void load_quadrature_data()
{
  

  char fname[255];
  sprintf(fname,"LegendreData/LegQuad-N%i.dat",N1-1);
  
  printf("===> Info: Loading Quadrature Data from File %s ...\n",fname);

  fstream parameters;
  parameters.open(fname, ios::in);
  if(parameters.is_open())
    {
      
      for (int i=0;i<N1;i++)
	{
	  parameters >> p0vals[i];
	  parameters >> p0weights[i];
	  //printf("i=%i (%.16f,%.16f)\n",i,p0vals[i],p0weights[i]);
	}

      printf("\t\t\t\t\t\t\t\t\t\t ...finished!\n");
    }
  else
    {
      printf("\t\t\t\t\t\t\t\t\t\t ...FAILED\n");
    }
  parameters.close();
}


void derived_data()
{

  for (int i=0;i<N1;i++)
    {
      s2p0vals[i]=sin(p0vals[i]*M_PI*0.5)*sin(p0vals[i]*M_PI*0.5);
      c2p0vals[i]=cos(p0vals[i]*M_PI*0.5)*cos(p0vals[i]*M_PI*0.5);
      tp0vals[i]=tan(p0vals[i]*M_PI*0.5);
      modp0weights[i]=p0weights[i]/c2p0vals[i];

    }
  modp0weights[0]=0.5*p0weights[0]/c2p0vals[0];

  for (int i=0;i<N0;i++)
    {
      wn[i]=2*M_PI*i;
      wnweights[i]=1;
      //printf("hello i=%i %f\n",i,wn[i]);
    }
  wnweights[0]=0.5;

}


void allocate_memory()
{
  
}

void free_memory()
{
  
  
}

double get_PiB(int n,double p1)
{

  double res=0;

  if (n==0)
    {

      double cut=atan(p1*0.5);

      for (int j=0;j<N1;j++)
	{
	  //double k=p0vals[j];
	  //res-=4*ahat*p1*p0weights[j]/(1+exp(k*p1*0.5))*(sqrt(1-k*k));

	  //double k=tp0vals[j];
	   //if (k<p1*0.5)
	   // res-=modp0weights[j]/(1+exp(k))*sqrt(1-4*k*k/(p1*p1));


	   double q=tan(p0vals[j]*cut);
	   double weight=modp0weights[j]*c2p0vals[j]/(cos(p0vals[j]*cut)*cos(p0vals[j]*cut));
	   res-=weight/(1+exp(q))*sqrt(1-4*q*q/(p1*p1));
	  
	  //printf("k=%f add=%f\n",k,4*ahat*p1*p0weights[j]/(1+exp(k*p1*0.5))*(sqrt(1-k*k)));
	  //printf("j=%i k=%f old %f q=%f new %f\n",j,k,M_PI*0.5*modp0weights[j]/(1+exp(k))*sqrt(1-4*k*k/(p1*p1)), q,cut*weight/(1+exp(q))*sqrt(1-4*q*q/(p1*p1)));
	  //printf("where %f %f\n",M_PI*0.5*modp0weights[j],cut*weight);
	}
      //res*=M_PI*4*ahat;
      res*=8*ahat*cut;

      res+=log(2)*8*ahat;
    }
  else
    {
      double p0=wn[n];
      double P2=p0*p0+p1*p1;
      double SP2=sqrt(P2);
      
      for (int j=0;j<N1;j++)
	{
	  double k=tp0vals[j];
	  std::complex<double> a1=P2-4*k*k-4.*I*p0*k;      
	  res+=modp0weights[j]/(1+exp(k))*(1-sqrt(a1).real()/SP2);
	}
      res*=0.5*M_PI;


      res*=8*ahat*P2/(p1*p1);
    }
  //res*=pow(ahat,-0.1);

 
  return res;
  
}

double get_PiA(int n,double p1)
{

  double res=0;

  if (n==0)
    {
      double cut=atan(p1*0.5);
      for (int j=0;j<N1;j++)
	{
	  //double k=p0vals[j];
	  //res-=4*ahat*p1*p0weights[j]/(1+exp(k*p1*0.5))*k*k/sqrt(1-k*k);

	  //double k=tp0vals[j];
	  //if (k<p1*0.5)
	  //res+=modp0weights[j]/(1+exp(k))*(4*k*k/sqrt(1-4*k*k/(p1*p1)));

	  double q=tan(p0vals[j]*cut);
	  double weight=modp0weights[j]*c2p0vals[j]/(cos(p0vals[j]*cut)*cos(p0vals[j]*cut));
	  res+=weight/(1+exp(q))*(4*q*q/sqrt(1-4*q*q/(p1*p1)));

	  
	}
      //res*=-4*ahat*M_PI/(p1*p1);
      res*=-8*ahat*cut/(p1*p1);
    }
  else
    {
      double p0=wn[n];
      double P2=p0*p0+p1*p1;
      double SP2=sqrt(P2);
      
      for (int j=0;j<N1;j++)
	{
	  double k=tp0vals[j];
	  std::complex<double> a1=P2-4*k*k-4.*I*p0*k;
	  std::complex<double> a2=(I*p0+2*k);
	  
      res+=modp0weights[j]/(1+exp(k))*(p0*p0+SP2*(a2*a2/sqrt(a1)).real());
	}
      res*=0.5*M_PI;

      
      res*=-8*ahat/(p1*p1);
    }
  //res*=pow(ahat,-0.1);
  
  return res;
  
}

double vacuum()
{
  double res=0;
  double e2=ahat*M_PI*0.5;
  for (int j=0;j<N1;j++)
    {
      double k=tp0vals[j]+IR;
      
      double myarg=1/(exp(k)-1)*((k*k+e2*e2)*atan(k/e2)-e2*k);
      res+=myarg*modp0weights[j];
    }
  res/=8*M_PI;
  return -(0.0956566+res);
}

double partitionB()
{
  double res=0;
  for (int i=0;i<N0;i++)
    for (int j=0;j<N1;j++)
      {
	double mypi=get_PiB(i,sqrt(tp0vals[j]+IR));
	double p2=wn[i]*wn[i]+tp0vals[j]+IR;
	double myarg=log(1+mypi/(p2+ahat*M_PI*0.5*sqrt(p2)));
	res+=myarg*modp0weights[j]*wnweights[i];
	//if (i==0)
	//  printf("i=%i j=%i p2=%f res=%f mypi*p2=%f added=%f\n",i,j,p2,res,mypi*p2,myarg*modp0weights[j]*wnweights[i]);
      }
  
  return res/(8.);
}

double partitionA()
{
  double res=0;
  for (int i=0;i<N0;i++)
    for (int j=0;j<N1;j++)
      {
	double mypi=get_PiA(i,sqrt(tp0vals[j]+IR));
	double p2=wn[i]*wn[i]+tp0vals[j]+IR;
	double myarg=log(1+mypi/(p2+ahat*M_PI*0.5*sqrt(p2)));
	res+=myarg*modp0weights[j]*wnweights[i];
	//if (i==0)
	//  printf("i=%i j=%i p2=%f res=%f mypi*p2=%f added=%f\n",i,j,p2,res,mypi*p2,myarg*modp0weights[j]*wnweights[i]);
      }
  
  return res/(8.);
}

double piece1()
{
  double res=0;
  for (int i=0;i<N0;i++)
    for (int j=0;j<N1;j++)
      {
	double mypi=get_PiB(i,sqrt(tp0vals[j]+IR));
	double p2=wn[i]*wn[i]+tp0vals[j]+IR;
	double myarg=mypi/(p2+ahat*M_PI*0.5*sqrt(p2));
	res+=myarg*modp0weights[j]*wnweights[i];
      }
  
  return res/(8.);
}

double piece2()
{
  double res=0;
  for (int i=0;i<N0;i++)
    for (int j=0;j<N1;j++)
      {
	double mypi=get_PiB(i,sqrt(tp0vals[j]+IR));
	double p2=wn[i]*wn[i]+tp0vals[j]+IR;
	double myarg=(mypi/(p2+ahat*M_PI*0.5*sqrt(p2)))*(mypi/(p2+ahat*M_PI*0.5*sqrt(p2)));
	res+=myarg*modp0weights[j]*wnweights[i];
      }
  
  return res/(-16.);
}

double check(double m2)
{
  double res=0;
  for (int i=0;i<N0;i++)
    for (int j=0;j<N1;j++)
      {
	double p2=wn[i]*wn[i]+tp0vals[j]+IR;
	double myarg=1./(p2+m2)/(p2+m2);
	res+=myarg*modp0weights[j]*wnweights[i];
      }
  res*=sqrt(m2);
  
  return res/(8.);
}


int main(int argc, char* argv[])
{

  load_quadrature_data();
  derived_data();
  allocate_memory();

  char buffer[255];
  fstream out;
  sprintf(buffer,"results-NT%i-N%i-IRcut%f.dat",N0,N1,IR);
  if (write_to_file)
    out.open(buffer,ios::out);

  if (write_to_file)
    printf("===> Info: Will store results in file %s\n",buffer);
  else
    printf("===> Info: Not storing results in file\n");

  if (write_to_file)
    {
      sprintf(buffer,"#alpha/T\t\t f_A\t f_B\t f_V\t sum\n");
      out << buffer;
      out.flush();
    }
  
  
  printf("===> Info: working in units where T=1\n");

  //printf("test pib=%f\n",get_PiB(0.1,0.1));
  //printf("test pib=%f\n",get_PiB(0.1,0.001));
  //printf("test pib=%f\n",get_PiB(0.01,0.0001));
  //printf("test partition %f\n",partitionB());
  //printf("test piece1 %f\n",piece1());
  //printf("test piece1 %f\n",piece2());
  //printf("check(1)=%f check(100)=%f  check(10000)=%f\n",check(1.),check(100.),check(10000));



    
  double inc=0.1;
    
  if (write_to_file)
    {
      /*
	for (int i=0;i<301;i++)
	{
	ahat=inc*i;
	double fa=partitionA();
	double fb=partitionB();
	double v=vacuum();
	
	printf("alpha/T=%f fA=%f fB=%f fV=%f sum=%f \n",ahat,fa,fb,v,fa+fb+2*v);
	
	sprintf(buffer,"%f\t\t %f\t %f\t %f\t %f\n",ahat,fa,fb,v,fa+fb+2*v);
	out << buffer;
	out.flush();
	
	}*/
      
      inc=1.1;
      ahat=0.001;
      for (int i=0;i<201;i++)
	{
	  
	  double fa=partitionA();
	  double fb=partitionB();
	  double v=vacuum();
	  
	  printf("alpha/T=%f fA=%f fB=%f fV=%f sum=%f (sum/dof=%f) \n",ahat,fa,fb,v,fa+fb+2*v,(fa+fb+2*v)/0.191313);
	  
	  sprintf(buffer,"%f\t\t %f\t %f\t %f\t %f\n",ahat,fa,fb,v,fa+fb+2*v);
	  out << buffer;
	  out.flush();
	  ahat*=inc;
	}
    }
  
  if (!write_to_file)
    {
      
      ahat=alpha_default;
      double fa=partitionA();
      double fb=partitionB();
      double v=vacuum();
      
      printf("alpha/T=%f fA=%f fB=%f fV=%f sum=%f (sum/dof=%f) \n",ahat,fa,fb,v,fa+fb+2*v,(fa+fb+2*v)/0.191313);
    }

  //printf("echk %f %f\n",get_PiA(0,100),get_PiB(0,100.));
  
  free_memory();
  if (write_to_file)
    out.close();
}
