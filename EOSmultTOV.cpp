//Generator of a number of random EOS
//check the speed of sound
//check for mass limit

using namespace std;
#include <cmath>
#define  G 6.6732*TMath::Power(10,-8)
#define  c 2.9979*TMath::Power(10,10)
#define  Mo 1.987*TMath::Power(10,33)
#define  eo  TMath::Power(10,15)
#define  Mb  1.66*TMath::Power(10,-24)
#define  kappa c*c/(G*eo)

  
class findParameter {
public:
  double p1;
  double Gamma1;
  double Gamma2;
  double Gamma3;
  double K1;
  double K2;
  double K3;
  double p2;
  double rho1;
  double rho2;
  double rho0;
  void setParameters(double pa1,double pa2,double pa3,double pa4,double pa5, double par6, double par7){
    p1=TMath::Power(10,pa1)*kappa*G/(c*c*c*c);
    Gamma1=pa2;
    Gamma2=pa3;
    Gamma3=pa4;
    rho1= par6*2.28*TMath::Power(10,14)*kappa*G/(c*c);
    rho2= par7*2.28*TMath::Power(10,14)*kappa*G/(c*c);
    K1=p1/TMath::Power(rho1,Gamma1);
    K2=p1/TMath::Power(rho1,Gamma2);
    K3=K2*TMath::Power(rho2,Gamma2-Gamma3);
    p2=K2*TMath::Power(rho2,Gamma2);
    rho0=pa5*kappa*G/(c*c);
  }




  //needed for the crust
  int crust(double p,vector<double> vp,vector<double> vd){//crust EOS
    int pos;
    for(int i=0;i<vp.size()-1;i++){
      if(vp.at(i)<=p && vp.at(i+1)>=p) return i;
    }
    
  }
  double dmdr(double r,double m,double p,vector<double> vp,vector<double> ve){//energy density
    int pos;
    for(int i=0;i<vp.size()-1;i++){
      if(vp.at(i)<=p && vp.at(i+1)>=p) pos= i;
    }


    double a=(TMath::Log10(ve.at(pos))-TMath::Log10(ve.at(pos+1)))/(TMath::Log10(vp.at(pos))-TMath::Log10(vp.at(pos+1)));
    double b=TMath::Power(10,TMath::Log10(ve.at(pos+1))-(TMath::Log10(ve.at(pos))-TMath::Log10(ve.at(pos+1)))/(TMath::Log10(vp.at(pos))-TMath::Log10(vp.at(pos+1)))*TMath::Log10(vp.at(pos+1)));
    double en=b*TMath::Power(p,a);
    return 4.*TMath::Pi()*r*r*en;
  }


  double dpdr(double r,double m,double p,vector<double> vp,vector<double> ve){
    int pos;
    for(int i=0;i<vp.size()-1;i++){
      if(vp.at(i)<=p && vp.at(i+1)>=p) pos= i;
    }

    double a=(TMath::Log10(ve.at(pos))-TMath::Log10(ve.at(pos+1)))/(TMath::Log10(vp.at(pos))-TMath::Log10(vp.at(pos+1)));
    double b=TMath::Power(10,TMath::Log10(ve.at(pos+1))-(TMath::Log10(ve.at(pos))-TMath::Log10(ve.at(pos+1)))/(TMath::Log10(vp.at(pos))-TMath::Log10(vp.at(pos+1)))*TMath::Log10(vp.at(pos+1)));
    double en=b*TMath::Power(p,a);
  
    return -(en+p)*(m+4.*TMath::Pi()*r*r*r*p)/(r*(r-2.*m));

  }
  //need for crust
  

  
  double pressure(double d){//EOS
    if(d<=rho1){
      return K1*TMath::Power(d,Gamma1);
    }
    else if(d>rho1 && d<rho2){
      return K2*TMath::Power(d,Gamma2);
    }
    else if(d>=rho2){
      return K3*TMath::Power(d,Gamma3);
    }
  }



  double h(double d){
    //y=2.56x-3.37 FPS
    //y=1.5x+11.6 NV
    double a0FPS,a0NV,a1,a2,a3,K0NV,K0FPS;
    if (rho0<TMath::Power(10,14.35)*kappa*G/(c*c)){//NV and FPS
      K0NV=2.122*TMath::Power(10,30)/(c*c*c*c/(G*kappa))/TMath::Power(1.192786036339965*1.66*TMath::Power(10,14)/(c*c/(G*kappa)),1.5);
      a0NV=-1-K0NV*TMath::Power(TMath::Power(10,12)*kappa*G/(c*c),1.5-1)/(1.5-1)+1;

      a1=-1-K1*TMath::Power(rho0,Gamma1-1)/(Gamma1-1)+(1+a0NV)+K0NV*TMath::Power(rho0,1.5-1)/(1.5-1);
    }
    else{//FPS
      K0NV=1.257*TMath::Power(10,30)/((c*c*c*c/(G*kappa))*TMath::Power(5.96927*1.66*TMath::Power(10,13)/(c*c/(G*kappa)),1.5));
      a0NV=-1-K0NV*TMath::Power(TMath::Power(10,12)*kappa*G/(c*c),1.5-1)/(1.5-1)+1;

      K0FPS=5.19662*TMath::Power(10,35)/(c*c*c*c/(G*kappa))/TMath::Power(9.989852383181389*1.66*TMath::Power(10,14)/(c*c/(G*kappa)),2.56);
      a0FPS=-1-K0FPS*TMath::Power(TMath::Power(10,14.35)*kappa*G/(c*c),2.56-1)/(2.56-1)+(1+a0NV)+K0NV*TMath::Power(TMath::Power(10,14.35)*kappa*G/(c*c),2.56-1)/(2.56-1);
      
      a1=-1-K1*TMath::Power(rho0,Gamma1-1)/(Gamma1-1)+(1+a0FPS)+K0FPS*TMath::Power(rho0,2.56-1)/(2.56-1);
    }
    
    
    a2=-1-K2*TMath::Power(rho1,Gamma2-1)/(Gamma2-1)+(1+a1)+K1*TMath::Power(rho1,Gamma1-1)/(Gamma1-1);
    //a1=-1-K1*TMath::Power(rho0,Gamma1-1)/(Gamma1-1)+1;
    //a2=-1-K2*TMath::Power(rho1,Gamma2-1)/(Gamma2-1)+(1+a1)+K1*TMath::Power(rho1,Gamma1-1)/(Gamma1-1);
    a3=-1-K3*TMath::Power(rho2,Gamma3-1)/(Gamma3-1)+(1+a2)+K2*TMath::Power(rho2,Gamma2-1)/(Gamma2-1);
    
    if(d<=rho1){
      return 1.+a1+K1*Gamma1*TMath::Power(d,Gamma1-1.)/(Gamma1-1.);
    }
    else if(d>rho1&&d<rho2){
      return 1.+a2+K2*Gamma2*TMath::Power(d,Gamma2-1.)/(Gamma2-1.);
    }
    else if(d>=rho2){
      return 1.+a3+K3*Gamma3*TMath::Power(d,Gamma3-1.)/(Gamma3-1.);
    }
  
 
  }
  double e(double d){//energy density
    //y=2.56x-3.37 FPS
    //y=1.5x+11.6 NV
    double a0FPS,a0NV,a1,a2,a3,K0NV,K0FPS;
    if (rho0<TMath::Power(10,14.35)*kappa*G/(c*c)){//NV and FPS
      K0NV=2.122*TMath::Power(10,30)/(c*c*c*c/(G*kappa))/TMath::Power(1.192786036339965*1.66*TMath::Power(10,14)/(c*c/(G*kappa)),1.5);
      a0NV=-1-K0NV*TMath::Power(TMath::Power(10,12)*kappa*G/(c*c),1.5-1)/(1.5-1)+1;

      a1=-1-K1*TMath::Power(rho0,Gamma1-1)/(Gamma1-1)+(1+a0NV)+K0NV*TMath::Power(rho0,1.5-1)/(1.5-1);
    }
    else{//FPS
      K0NV=1.257*TMath::Power(10,30)/((c*c*c*c/(G*kappa))*TMath::Power(5.96927*1.66*TMath::Power(10,13)/(c*c/(G*kappa)),1.5));
      a0NV=-1-K0NV*TMath::Power(TMath::Power(10,12)*kappa*G/(c*c),1.5-1)/(1.5-1)+1;

      K0FPS=5.19662*TMath::Power(10,35)/(c*c*c*c/(G*kappa))/TMath::Power(9.989852383181389*1.66*TMath::Power(10,14)/(c*c/(G*kappa)),2.56);
      a0FPS=-1-K0FPS*TMath::Power(TMath::Power(10,14.35)*kappa*G/(c*c),2.56-1)/(2.56-1)+(1+a0NV)+K0NV*TMath::Power(TMath::Power(10,14.35)*kappa*G/(c*c),2.56-1)/(2.56-1);
      
      a1=-1-K1*TMath::Power(rho0,Gamma1-1)/(Gamma1-1)+(1+a0FPS)+K0FPS*TMath::Power(rho0,2.56-1)/(2.56-1);
    }
    
    
    a2=-1-K2*TMath::Power(rho1,Gamma2-1)/(Gamma2-1)+(1+a1)+K1*TMath::Power(rho1,Gamma1-1)/(Gamma1-1);
    //a1=-1-K1*TMath::Power(rho0,Gamma1-1)/(Gamma1-1)+1;
    //a2=-1-K2*TMath::Power(rho1,Gamma2-1)/(Gamma2-1)+(1+a1)+K1*TMath::Power(rho1,Gamma1-1)/(Gamma1-1);
    a3=-1-K3*TMath::Power(rho2,Gamma3-1)/(Gamma3-1)+(1+a2)+K2*TMath::Power(rho2,Gamma2-1)/(Gamma2-1);
    
    if(d<=rho1){
      return (1+a1)*d+K1*TMath::Power(d,Gamma1)/(Gamma1-1.);
    }
    else if(d>rho1&&d<rho2){
      return (1+a2)*d+K2*TMath::Power(d,Gamma2)/(Gamma2-1.);
    }
    else if(d>=rho2){
      return (1+a3)*d+K3*TMath::Power(d,Gamma3)/(Gamma3-1.);
    }  
  }




  double f(double p,double en){
  
    return 1./(en*(kappa*G)/(c*c)+(p*(G*kappa)/(c*c*c*c)));
  }


  double enthalpy(vector<double> p,vector<double> en,double h0,int num,int data){

    double dx=(p.at(num)-p.at(data))/((num-data)*1.);
    double integr=0;
    double press=p.at(0);
    
    integr=f(p.at(data),en.at(data))+f(p.at(num),en.at(num));

    
    for(int i=1;i<(num-data)*1-1;i++){
      press+=dx;
      int pos=crust(press,p,en);

      double a=(TMath::Log10(en.at(pos))-TMath::Log10(en.at(pos+1)))/(TMath::Log10(p.at(pos))-TMath::Log10(p.at(pos+1)));
      double b=TMath::Power(10,TMath::Log10(en.at(pos+1))-(TMath::Log10(en.at(pos))-TMath::Log10(en.at(pos+1)))/(TMath::Log10(p.at(pos))-TMath::Log10(p.at(pos+1)))*TMath::Log10(p.at(pos+1)));
      double ener=b*TMath::Power(press,a);
      if(i%2==0){
	integr+=2.*f(press,ener);
      }
      else{
	integr+=4.*f(press,ener);
      }
    }

    integr=h0+integr*dx*(kappa*G)/(c*c*3.);
 
     
    return integr;
  }




  bool soundSpeed(double p,double d){

    if(d<=rho1){
      if(TMath::Sqrt(Gamma1*p/(e(d)+p))>1){
	return false;
      }
      else{
	return true;
      }
    }
    else if(d>rho1&&d<rho2){
      if(TMath::Sqrt(Gamma2*p/(e(d)+p))>1){
	return false;
      }
      else{
	return true;
      }
    }
    else if(d>=rho2){
      if(TMath::Sqrt(Gamma3*p/(e(d)+p))>1){
	return false;
      }
      else{
	return true;
      }
    } 

  }



  double maxMassOrRadius(vector<double> v_d,vector<double> v_p,vector<double> v_e,int m_or_r){

    int count=0;
    double stepsize=100./(TMath::Power(kappa,0.5));
    double r0=TMath::Power(10,-15)/(TMath::Power(kappa,0.5));
    double m0=TMath::Power(10,-15)*G/(c*c*TMath::Power(kappa,0.5));
    double d0=v_d.at(v_d.size()-1);
    double p0=v_p.at(v_p.size()-1);

    
    double r,m,d,p,k0,k1,k2,k3,l0,l1,l2,l3;

    d=d0;
    p=p0;
    m=m0;
  
    //RK4 Begining
    while(d>TMath::Power(10,0.9)*kappa*G/(c*c)){   
      if(count==0){
	/*
	  r=r0;//radius
	  m=m0;//mass
	  d=d0;//density
	  p=p0;//pressure
	*/
	r=r0;
	m=4*TMath::Pi()*d0*r0*r0*r0/3-TMath::Abs((TMath::Log10(d0)-TMath::Log10(v_d.size()-2))/(TMath::Log10(p0)-TMath::Log10(v_p.size()-2)))*8/15*3.14*3.14*d0/p0*(d0+p0)*(d0+3*p0)*TMath::Power(r0,5);
	d=d0-2/3*3.14*d0/p0*(d0+p0)*(d0+3*p0)*r0*r0;
	p=p0-TMath::Abs((TMath::Log10(d0)-TMath::Log10(v_d.at(v_p.size()-2)))/(TMath::Log10(p0)-TMath::Log10(v_p.at(v_p.size()-2))))*2*3.14*(d0+p0)*(d0+3*p0)*r0*r0/3;
      }
 
      else{
	
	double mold=m;
	double rold=r;
	double dold=d;
	double pold=p;
	double rnew,mnew,dnew,pnew;
	if(d<TMath::Power(10,14.5)*kappa*G/(c*c)) stepsize=10./(TMath::Power(kappa,0.5));
	if(d<TMath::Power(10,5.5)*kappa*G/(c*c)) stepsize=1./(TMath::Power(kappa,0.5));
	if(d<TMath::Power(10,3.5)*kappa*G/(c*c)) stepsize=0.01/(TMath::Power(kappa,0.5));
	if(d<TMath::Power(10,1.5)*kappa*G/(c*c)) stepsize=0.0001/(TMath::Power(kappa,0.5));
	
	k0=dmdr(rold,mold,pold,v_p,v_e);
	l0=dpdr(rold,mold,pold,v_p,v_e);

	if(pold+stepsize*0.5*l0<v_p.at(0))break;
	k1= dmdr(rold+stepsize/2.,mold+stepsize*0.5*k0,pold+stepsize*0.5*l0,v_p,v_e);
	l1= dpdr(rold+stepsize/2.,mold+stepsize*0.5*k0,pold+stepsize*0.5*l0,v_p,v_e);

	if(pold+stepsize*0.5*l1<v_p.at(0))break;
	k2 = dmdr(rold+stepsize/2.,mold+stepsize*0.5*k1,pold+stepsize*0.5*l1,v_p,v_e);
	l2 = dpdr(rold+stepsize/2.,mold+stepsize*0.5*k1,pold+stepsize*0.5*l1,v_p,v_e);

	if(pold+stepsize*l2<v_p.at(0))break;
	k3 = dmdr(rold+stepsize,mold+stepsize*k2,pold+stepsize*l2,v_p,v_e);
	l3 = dpdr(rold+stepsize,mold+stepsize*k2,pold+stepsize*l2,v_p,v_e);

	rnew=rold+stepsize;
	mnew=mold+(k0+2.*k1+2.*k2+k3)*stepsize/6.;
	pnew=pold+(l0+2.*l1+2.*l2+l3)*stepsize/6.;

	if(pnew<v_p.at(0))break;
	int pos=crust(pnew,v_p,v_d);

	double a=(TMath::Log10(v_d.at(pos))-TMath::Log10(v_d.at(pos+1)))/(TMath::Log10(v_p.at(pos))-TMath::Log10(v_p.at(pos+1)));
	double b=TMath::Power(10,TMath::Log10(v_d.at(pos+1))-(TMath::Log10(v_d.at(pos))-TMath::Log10(v_d.at(pos+1)))/(TMath::Log10(v_p.at(pos))-TMath::Log10(v_p.at(pos+1)))*TMath::Log10(v_p.at(pos+1)));
	dnew=b*TMath::Power(pnew,a);

 
      
	m=mnew;
	r=rnew;
	p=pnew;
	d=dnew;
	

      }
      count++;
    }

    if(m_or_r==0){
      return (m*c*c*TMath::Power(kappa,0.5))/(Mo*G);
    }else if(m_or_r==1){
      return (r*TMath::Power(kappa,0.5))/100000.;//Radius
    }
  }
  
  
};



void EOSmultTOV(){
  int num;
  cout<<"Give the number of EOS"<<endl;
  cin>>num;

  ofstream eosTable;
  eosTable.open("EosTable",ios::out);
  eosTable<<"EOS  log(p1)  Gamma1  Gamma2  Gamma3"<<endl;
  bool Continue= false;

  //choose surface EOS
  vector<string> surfEOS;
  surfEOS.push_back("eosNV");
  surfEOS.push_back("eosFPS");
  int sEOSnum=0;//0 = NV   1 = FPS
  int fileNum=0;


  double d0=1.58489*TMath::Power(10,15)*kappa*G/(c*c);//central density (rho not e)

  //Random EOS Generation Begins
  while(fileNum<num){
    


    findParameter fp;
    double par1,par2,par3,par4,par5,par6,par7;
    TRandom3* rand =new TRandom3();
    rand->SetSeed(0);
    par1=34.331;//rand->Uniform(34,35);//log(p1)
    par2=3.418;//rand->Uniform(1.5,4.5);//Gamma1
    par3=2.835;//rand->Uniform(3,6);//Gamma2
    par4=2.832;//rand->Uniform(2.5,3);//Gamma3
    par6=2.198189621;//rand->Uniform(1.5,3);//rho1*2.28 10^14 g/cm^3
    par7=4.385964912;//rand->Uniform(TMath::Max(2.00,par6),4);//rho2*2.28 10^14 g/cm^3




    //finds the core-crust density "border"
    if(sEOSnum==0){//y=1.5x+11.6
      par5=TMath::Power(10,(11.6-TMath::Log10(TMath::Power(10,par1)/(TMath::Power(par6*2.28*TMath::Power(10,14),par2))))/(par2-1.5));//rho0
      if((par5<=TMath::Power(10,14.35)) && (par5>=TMath::Power(10,12.2967))){
	Continue=true;
      }
      else{
	Continue=false;
      }
    }
    else if(sEOSnum==1){//y=2.56x-3.37
      par5=TMath::Power(10,(-3.37-TMath::Log10(TMath::Power(10,par1)/(TMath::Power(par6*2.28*TMath::Power(10,14),par2))))/(par2-2.56));//rho0
      if((par5>=TMath::Power(10,14.35)) && (par5<=TMath::Power(10,16.22))){
	Continue=true;
      }
      else if(par5<TMath::Power(10,14.35)){
	par5=TMath::Power(10,(11.6-TMath::Log10(TMath::Power(10,par1)/(TMath::Power(par6*2.28*TMath::Power(10,14),par2))))/(par2-1.5));//rho0
	if((par5<=TMath::Power(10,14.35)) && (par5>=TMath::Power(10,12.2967))){
	  Continue=true;
	}
	else{
	  Continue=false;
	}
      }
      else{ 
	Continue=false;
      }

    }

    
    if(!Continue) continue;


    //call class
    fp.setParameters(par1,par2,par3,par4,par5,par6,par7);




    //Take the Surface EOS - Begining
    bool BoolSpeed=true;
    int count=0;
    vector<double> v_m,v_r,v_p,v_d,v_e,v_h;
    vector<double> vEnergyDensity,vEnthalpy,vDensity,vPressure; //non dimensional
    ifstream myinputfile;
    myinputfile.open(surfEOS.at(sEOSnum),ios::in);
    string line;
    string a1,a2,a3,a4,a5;
    int data;
    while( getline(myinputfile,line) ){
      if(count==0){
	stringstream liness = stringstream(line);
	getline(liness,a1,' ');
	data=stod(a1);
      }
      else{
	stringstream liness = stringstream(line);
	getline(liness,a1,' ');
	getline(liness,a2,' ');
	getline(liness,a3,' ');
	getline(liness,a4,' ');
	if(sEOSnum==0){getline(liness,a5,' ');}
	if((stod(a4)*1.66*TMath::Power(10,-24))>(fp.rho0*c*c/(kappa*G))) break;
	v_e.push_back(stod(a1));
	v_p.push_back(stod(a2));
	v_h.push_back(stod(a3));
	v_d.push_back(stod(a4)*Mb);
	vPressure.push_back(stod(a2)*kappa*G/(c*c*c*c));
	vEnergyDensity.push_back(stod(a1)*kappa*G/(c*c));
	vEnthalpy.push_back(stod(a3)/(c*c));
	vDensity.push_back(stod(a4)*Mb*kappa*G/(c*c));
      }
      count++;
    }
    count=0;
    //Take the Surface EOS -END
 






    
    //EOS Begining
    double stepsize=(d0-fp.rho0)/22.;
    double p0,d,p;
    d=fp.rho0;
    p=fp.pressure(fp.rho0);
 
    while(d<=d0){
 
	v_d.push_back(d*c*c/(G*kappa));
	v_p.push_back(p*c*c*c*c/(G*kappa));
	v_e.push_back(fp.e(d)*c*c/(G*kappa));

	vPressure.push_back(p);
	vEnergyDensity.push_back(fp.e(d));
	vDensity.push_back(d);	
	//Finds N by using the politropic equation
	//v_h.push_back(TMath::Log(fp.h(d))*(c*c));
	//vEnthalpy.push_back(TMath::Log(fp.h(d)));
 
      BoolSpeed=fp.soundSpeed(p,d);
      if(!BoolSpeed)break;
      d=d+stepsize;
      p=fp.pressure(d);
      count++;
    }
    count=0;
    //EOS end


    

    //Check the speed of sound
    if(!BoolSpeed) continue;



    //Finds the N by integration
    /*
    for(int i=data;i<v_d.size();i++){
      double tempH=fp.enthalpy(v_p,v_e,v_h.at(data-2),i,data-1);
      v_h.push_back(tempH);
      vEnthalpy.push_back(tempH/(c*c));
    }
    */
    //Get e,p,h,d ENDS
   
    

    //check the Maximum Mass
    if(fp.maxMassOrRadius(vDensity,vPressure,vEnergyDensity,0)>1.97) continue;

    
    //Write EOS parameters Begining
    string title="";
    title=title+"eosSA"+fileNum;
    fileNum++;
    eosTable<<title<<" "<<par1<<" "<<par2<<" "<<par3<<" "<<par4<<endl;
    ofstream output;
    output.open(title,ios::out);
    //Write EOS parameters End

  
    //output Data
    output<<v_d.size()/*-data*/<<endl;//remove data when suurface included
    for(int i=0;i<v_d.size();i++){
      output<<v_e.at(i)<<" "<<v_p.at(i)<<" "<<v_h.at(i)<<" "<<v_d.at(i)/(Mb)<<" "<<endl;
    }


    //reset the vectors
    output.close();
    vEnergyDensity.clear();
    vEnthalpy.clear();
    vDensity.clear();
    vPressure.clear();
  }
  //Random EOS Generation Begins
}
