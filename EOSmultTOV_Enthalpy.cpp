//Generator of a number of random EOS
//check the speed of sound
//check for mass limit
//is based on the enthalpy

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
  double a1;
  double a2;
  double a3;
  
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

  void setAlphas(double ec, double dc,double pc){

    a1=-1-K1*TMath::Power(dc,Gamma1-1)/(Gamma1-1)+ec/dc;
    a2=-1-K2*TMath::Power(rho1,Gamma2-1)/(Gamma2-1)+(1+a1)+K1*TMath::Power(rho1,Gamma1-1)/(Gamma1-1);
    a3=-1-K3*TMath::Power(rho2,Gamma3-1)/(Gamma3-1)+(1+a2)+K2*TMath::Power(rho2,Gamma2-1)/(Gamma2-1);


  }


  //finds the location
int loc(double p,vector<double> vp){//crust EOS
  int pos;
  for(int i=vp.size()-1;i>0;i--){
    if(vp.at(i)>=p && vp.at(i-1)<=p){
      return i-1;
    }
    else if(vp.at(vp.size()-1)<p){
      return vp.size()-2;
    }
    else if(vp.at(0)>p){
      return 0;
    }
  }
    
}

  //need for crust


 double h(double d){

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

double rho(double eta){//energy density
 
  if(eta <= (h(rho1)-1) ){
    return TMath::Power((eta-a1)/(K1*(1/(Gamma1-1)+1)),1/(Gamma1-1));
  }
  else if(eta > (h(rho1)-1) && eta < (h(rho2)-1)){
    return TMath::Power((eta-a2)/(K2*(1/(Gamma2-1)+1)),1/(Gamma2-1));
  }
  else if(eta >= (h(rho2)-1)){
    return TMath::Power((eta-a3)/(K3*(1/(Gamma3-1)+1)),1/(Gamma3-1));
  }
  
}

  double pressure(double eta){
  
    if(eta <= (h(rho1)-1) ){
      return K1*TMath::Power((eta-a1)/(K1*(1/(Gamma1-1)+1)),1/(Gamma1-1)+1);
    }
    else if(eta > (h(rho1)-1) && eta < (h(rho2)-1)){
      return K2*TMath::Power((eta-a2)/(K2*(1/(Gamma2-1)+1)),1/(Gamma2-1)+1);
    }
    else if(eta >= (h(rho2)-1)){
      return K3*TMath::Power((eta-a3)/(K3*(1/(Gamma3-1)+1)),1/(Gamma3-1)+1);
    }
 
  }



  double e(double eta){//energy density
   
    if(eta <= (h(rho1)-1) ){
      return rho(eta)*(1.+(a1+eta/(Gamma1-1.))/(1.+1./(Gamma1-1.)));
    }
    else if(eta > (h(rho1)-1) && eta < (h(rho2)-1)){
      return rho(eta)*(1.+(a2+eta/(Gamma2-1.))/(1.+1./(Gamma2-1.)));
    }
    else if(eta >= (h(rho2)-1)){
      return rho(eta)*(1.+(a3+eta/(Gamma3-1.))/(1.+1./(Gamma3-1.)));;
    }
  
 
  }

double drhodh(double eta){
	
	double result;

	if(eta <= (h(rho1)-1) ){
	  result = TMath::Power((eta-a1)/(K1*(1/(Gamma1-1)+1)),1/(Gamma1-1))*(1./(Gamma1-1.)/(1.+1./(Gamma1-1.)))* TMath::Power((eta-a1)/(K1*(1/(Gamma1-1)+1)),1/(Gamma1-1)-1)*1/(Gamma1-1)/(K1*(1/(Gamma1-1)+1))*(1+(a1+eta/(Gamma1-1))/(1+1/(Gamma1-1)));
	}
	else if(eta > (h(rho1)-1) && eta < (h(rho2)-1)){
	  result= TMath::Power((eta-a2)/(K2*(1/(Gamma2-1)+1)),1/(Gamma2-1))*(1./(Gamma2-1.)/(1.+1./(Gamma2-1.)))* TMath::Power((eta-a2)/(K2*(1/(Gamma2-1)+1)),1/(Gamma2-1)-1)*1/(Gamma2-1)/(K2*(1/(Gamma2-1)+1))*(1+(a2+eta/(Gamma2-1))/(1+1/(Gamma2-1)));
	}
	else if(eta >= (h(rho2)-1)){
	  result = TMath::Power((eta-a3)/(K3*(1/(Gamma3-1)+1)),1/(Gamma3-1))*(1./(Gamma3-1.)/(1.+1./(Gamma3-1.)))* TMath::Power((eta-a3)/(K3*(1/(Gamma3-1)+1)),1/(Gamma3-1)-1)*1/(Gamma3-1)/(K3*(1/(Gamma3-1)+1))*(1+(a3+eta/(Gamma3-1))/(1+1/(Gamma3-1)));
	}
	

	return result;
}

double drdh(double eta,double r,double m,vector<double> veta,vector<double> vp,vector<double> ve){//energy density

  double p;

    int pos=loc(eta,veta);
    double a=(TMath::Log10(vp.at(pos))-TMath::Log10(vp.at(pos+1)))/(TMath::Log10(veta.at(pos))-TMath::Log10(veta.at(pos+1)));
    double b=TMath::Power(10,TMath::Log10(vp.at(pos+1))-(TMath::Log10(vp.at(pos))-TMath::Log10(vp.at(pos+1)))/(TMath::Log10(veta.at(pos))-TMath::Log10(veta.at(pos+1)))*TMath::Log10(veta.at(pos+1)));
    p=b*TMath::Power(eta,a);
  
    return -(r*(r-2.*m))/((m+4.*TMath::Pi()*r*r*r*p)*(eta+1.));
}


double dmdh(double eta,double r,double m,vector<double> veta,vector<double> vp,vector<double> ve){

  double en;
  
  int pos=loc(eta,veta);

  double a=(TMath::Log10(ve.at(pos))-TMath::Log10(ve.at(pos+1)))/(TMath::Log10(veta.at(pos))-TMath::Log10(veta.at(pos+1)));
  double b=TMath::Power(10,TMath::Log10(ve.at(pos+1))-(TMath::Log10(ve.at(pos))-TMath::Log10(ve.at(pos+1)))/(TMath::Log10(veta.at(pos))-TMath::Log10(veta.at(pos+1)))*TMath::Log10(veta.at(pos+1)));
  en=b*TMath::Power(eta,a);

  return 4.*TMath::Pi()*r*r*en*drdh(eta,r,m,veta,vp,ve);

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


  
  double maxMassOrRadius(vector<double> v_eta,vector<double> v_p,vector<double> v_e,int m_or_r){

    int count=0;
    double eta0=v_eta.at(v_eta.size()-1);
    double stepsize=-0.0001;
    
    double p0=v_p.at(v_eta.size()-1);
    double eta=eta0;
    
    double r,m,d,p,k0,k1,k2,k3,l0,l1,l2,l3;

    //RK4 Begining
    while(eta>v_eta.at(0)){  
      if(eta > (eta0-0.0001)){

	//this is for small distances


	r=TMath::Sqrt(3*(eta0-eta)/(2*TMath::Pi()*(e(eta0)+3.*p0)))*(1.-1./4.*(e(eta0)-3.*p0-3.*drhodh(eta0)/5.)*(eta0-eta)/(e(eta0)+3.*p0));

	m=4.*TMath::Pi()*e(eta0)*r*r*r*(1-3.*(eta0-eta)*drhodh(eta0)/(e(eta0)*5.))/3.;

	eta-=0.000001;

      }
 
      else{
	if (eta<v_eta.at(25)) stepsize=-0.00001;
	if (eta<v_eta.at(15)) stepsize=-0.0000001;
	if (eta<v_eta.at(10)) stepsize=-0.000000001;
	if (eta<v_eta.at(5)) stepsize=-0.00000000001;

	double mold=m;
	double rold=r;
	double etaold=eta;
	double rnew,mnew,etanew;

	if(eta+stepsize<v_eta.at(0))break;

	k0=drdh(etaold,rold,mold,v_eta,v_p,v_e);
	l0=dmdh(etaold,rold,mold,v_eta,v_p,v_e);
	
	k1= drdh(etaold+stepsize/2.,rold+stepsize*0.5*k0,mold+stepsize*0.5*l0,v_eta,v_p,v_e);
	l1= dmdh(etaold+stepsize/2.,rold+stepsize*0.5*k0,mold+stepsize*0.5*l0,v_eta,v_p,v_e);


	k2 = drdh(etaold+stepsize/2.,rold+stepsize*0.5*k1,mold+stepsize*0.5*l1,v_eta,v_p,v_e);
	l2 = dmdh(etaold+stepsize/2.,rold+stepsize*0.5*k1,mold+stepsize*0.5*l1,v_eta,v_p,v_e);

	k3 = drdh(etaold+stepsize,rold+stepsize*k2,mold+stepsize*l2,v_eta,v_p,v_e);
	l3 = dmdh(etaold+stepsize,rold+stepsize*k2,mold+stepsize*l2,v_eta,v_p,v_e);
	

	etanew=etaold+stepsize;

	rnew=rold+(k0+2.*k1+2.*k2+k3)*stepsize/6.;
	mnew=mold+(l0+2.*l1+2.*l2+l3)*stepsize/6.;
	

	m=mnew;
	r=rnew;
	eta=etanew;

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



void EOSmultTOV_Enthalpy(){
  int num;
  cout<<"Give the number of EOS"<<endl;
  cin>>num;

  ofstream eosTable;
  eosTable.open("EosTable",ios::out);
  eosTable<<"EOS  log(p1)  Gamma1  Gamma2  Gamma3 rho1 (g/cm^3) rho2 (g/cm^3)  Surface EOS"<<endl;
  bool Continue= false;

  //choose surface EOS
  vector<string> surfEOS;
  surfEOS.push_back("eosNV");
  surfEOS.push_back("eosFPS");
  int sEOSnum=1;//0 = NV   1 = FPS
  int fileNum=0;


  double d0=1.58489*TMath::Power(10,15)*kappa*G/(c*c);//central density (rho not e)

  //Random EOS Generation Begins
  while(fileNum<num){
    


    findParameter fp;
    double par1,par2,par3,par4,par5,par6,par7;
    TRandom3* rand =new TRandom3();
    rand->SetSeed(0);
    par1=/*34.331;*/rand->Uniform(34,35);//log(p1)
    par2=/*3.418;*/rand->Uniform(1.5,4.5);//Gamma1
    par3=/*2.835;*/rand->Uniform(3,6);//Gamma2
    par4=/*2.832;*/rand->Uniform(2.5,3);//Gamma3
    par6=/*2.198189621;*/rand->Uniform(1.5,3);//rho1*2.28 10^14 g/cm^3
    par7=/*4.385964912;*/rand->Uniform(TMath::Max(2.00,par6),4);//rho2*2.28 10^14 g/cm^3




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
      if((par5>=TMath::Power(10,14.35)) && (par5<=par6*2.28*TMath::Power(10,14)/*TMath::Power(10,16.22)*/)){
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
    vector<double> v_m,v_r,v_p,v_d,v_e,v_enth;
    vector<double> vEnergyDensity,vEnthalpy,vDensity,vPressure,vEta; //non dimensional
    ifstream myinputfile;
    myinputfile.open(surfEOS.at(sEOSnum),ios::in);
    string line;
    string read1,read2,read3,read4,read5;
    int data;
    while( getline(myinputfile,line) ){
      if(count==0){
	stringstream liness = stringstream(line);
	getline(liness,read1,' ');
	data=stod(read1);
      }
      else{
	stringstream liness = stringstream(line);
	getline(liness,read1,' ');
	getline(liness,read2,' ');
	getline(liness,read3,' ');
	getline(liness,read4,' ');
	if(sEOSnum==0){getline(liness,read5,' ');}
	if((stod(read4)*1.66*TMath::Power(10,-24))>(fp.rho0*c*c/(kappa*G))) break;
	v_e.push_back(stod(read1));
	v_p.push_back(stod(read2));
	v_enth.push_back((stod(read3)));
	v_d.push_back(stod(read4)*Mb);
	
	vPressure.push_back(stod(read2)*kappa*G/(c*c*c*c));
	vEnergyDensity.push_back(stod(read1)*kappa*G/(c*c));
	vEnthalpy.push_back((stod(read3))/(c*c));
	vDensity.push_back(stod(read4)*Mb*kappa*G/(c*c));
	vEta.push_back(TMath::Exp(stod(read3)/(c*c))-1);
      }
      count++;
    }
    data=v_e.size();
    count=0;
    //Take the Surface EOS -END
 
    fp.setAlphas(vEnergyDensity.at(data-1),vDensity.at(data-1),vPressure.at(data-1));





    
    //EOS Begining
    double stepsize=(fp.h(d0)-fp.h(fp.rho0))/22.;
    double p0,d,p,eta;
    bool CheckCont=true;
    eta=fp.h(fp.rho0)-1;
 
    while(eta<=(fp.h(d0)-1)){
 
	v_enth.push_back(TMath::Log(eta+1)*c*c);
	v_p.push_back(fp.pressure(eta)/(kappa*G/(c*c*c*c)));
	v_d.push_back(fp.rho(eta)/(kappa*G/(c*c)));
	v_e.push_back(fp.e(eta)/(kappa*G/(c*c)));

	
	vPressure.push_back(fp.pressure(eta));
	vEnergyDensity.push_back(fp.e(eta));
	vDensity.push_back(fp.rho(eta));	
	vEnthalpy.push_back(TMath::Log(eta+1));
	vEta.push_back(eta);

	if( (v_enth.at(v_enth.size()-1)<v_enth.at(v_enth.size()-2))  || (v_p.at(v_p.size()-1)<v_p.at(v_p.size()-2)) || (v_d.at(v_d.size()-1)<v_d.at(v_d.size()-2)) || (v_e.at(v_e.size()-1)<v_e.at(v_e.size()-2))){
	  CheckCont=false;
	  break;
	}
	
      BoolSpeed=fp.soundSpeed(fp.pressure(eta),fp.rho(eta));
      if(!BoolSpeed)break;
      eta+=stepsize;
      count++;
    }
    count=0;
    //EOS end


    

    //Check the speed of sound
    if(!BoolSpeed) continue;
    if(!CheckCont) continue;

    

    //check the Maximum Mass
    if(fp.maxMassOrRadius(vEta,vPressure,vEnergyDensity,0)>1.97) continue;

    //Write EOS parameters Begining
    string title="";
    title=title+"eosSA"+fileNum;
    fileNum++;
    eosTable<<title<<" "<<par1<<" "<<par2<<" "<<par3<<" "<<par4<<" "<<fp.rho1*c*c/(kappa*G)<<" "<<fp.rho2*c*c/(kappa*G)<<" "<<surfEOS.at(sEOSnum)<<endl;
    ofstream output;
    output.open(title,ios::out);
    //Write EOS parameters End

  
    //output Data
    output<<v_d.size()/*-data*/<<endl;//remove data when suurface included
    for(int i=0;i<v_d.size();i++){
      output<<v_e.at(i)<<" "<<v_p.at(i)<<" "<<v_enth.at(i)<<" "<<v_d.at(i)/(Mb)<<" "<<endl;
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
