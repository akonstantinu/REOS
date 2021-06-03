//finds M-R relation including the crust and Max mass and Radius by reading an EOS file
//but use the enthalpy derivatives
using namespace std;
#include <cmath>
#define  G 6.6732*TMath::Power(10,-8)
#define  c 2.9979*TMath::Power(10,10)
#define  Mo 1.987*TMath::Power(10,33)
#define  eo  TMath::Power(10,15)
#define  kappa c*c/(G*eo)
#define  rho0 TMath::Power(10,14.35)*kappa*G/(c*c)
#define  rho1 TMath::Power(10,14.7)*kappa*G/(c*c)
#define  rho2 TMath::Power(10,15)*kappa*G/(c*c)
#define  p1 TMath::Power(10,34.331)*kappa*G/(c*c*c*c)
#define  Gamma1 3.418
#define  Gamma2 2.835
#define  Gamma3 2.832
#define  K1 p1/TMath::Power(rho1,Gamma1)
#define  K2 p1/TMath::Power(rho1,Gamma2)
#define  K3 K2*TMath::Power(rho2,Gamma2-Gamma3)
#define  p2 K2*TMath::Power(rho2,Gamma2)






//needed for the crust
int crust(double p,vector<double> vp,vector<double> vd){//crust EOS
  int pos;
  for(int i=0;i<vp.size()-1;i++){
    if(vp.at(i)<=p && vp.at(i+1)>=p) return i;
  }
    
}

double dmdr1(double r,double m,double p,vector<double> vp,vector<double> ve){//energy density
  int pos;
  for(int i=0;i<vp.size()-1;i++){
    if(vp.at(i)<=p && vp.at(i+1)>=p) pos= i;
  }


  double a=(TMath::Log10(ve.at(pos))-TMath::Log10(ve.at(pos+1)))/(TMath::Log10(vp.at(pos))-TMath::Log10(vp.at(pos+1)));
  double b=TMath::Power(10,TMath::Log10(ve.at(pos+1))-(TMath::Log10(ve.at(pos))-TMath::Log10(ve.at(pos+1)))/(TMath::Log10(vp.at(pos))-TMath::Log10(vp.at(pos+1)))*TMath::Log10(vp.at(pos+1)));
  double en=b*TMath::Power(p,a);
  return 4.*TMath::Pi()*r*r*en;
}


double dpdr1(double r,double m,double p,vector<double> vp,vector<double> ve){
  int pos;
  for(int i=0;i<vp.size()-1;i++){
    if(vp.at(i)<=p && vp.at(i+1)>=p) pos= i;
  }

   double a=(TMath::Log10(ve.at(pos))-TMath::Log10(ve.at(pos+1)))/(TMath::Log10(vp.at(pos))-TMath::Log10(vp.at(pos+1)));
  double b=TMath::Power(10,TMath::Log10(ve.at(pos+1))-(TMath::Log10(ve.at(pos))-TMath::Log10(ve.at(pos+1)))/(TMath::Log10(vp.at(pos))-TMath::Log10(vp.at(pos+1)))*TMath::Log10(vp.at(pos+1)));
  double en=b*TMath::Power(p,a);
  
  return -(en+p)*(m+4.*TMath::Pi()*r*r*r*p)/(r*(r-2.*m));

}


///needed for the crust





double h(double d){
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







double rho(double eta){//energy density
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
      
      a1=-1-K1*TMath::Power(rho0,Gamma1-1)/(Gamma1-1)+(1+a0FPS)+K0FPS*TMath::Power(rho0,2.5-1)/(2.56-1);
    }
    
    
    a2=-1-K2*TMath::Power(rho1,Gamma2-1)/(Gamma2-1)+(1+a1)+K1*TMath::Power(rho1,Gamma1-1)/(Gamma1-1);
    //a1=-1-K1*TMath::Power(rho0,Gamma1-1)/(Gamma1-1)+1;
    //a2=-1-K2*TMath::Power(rho1,Gamma2-1)/(Gamma2-1)+(1+a1)+K1*TMath::Power(rho1,Gamma1-1)/(Gamma1-1);
    a3=-1-K3*TMath::Power(rho2,Gamma3-1)/(Gamma3-1)+(1+a2)+K2*TMath::Power(rho2,Gamma2-1)/(Gamma2-1);

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
      
      a1=-1-K1*TMath::Power(rho0,Gamma1-1)/(Gamma1-1)+(1+a0FPS)+K0FPS*TMath::Power(rho0,2.5-1)/(2.56-1);
    }
    
    
    a2=-1-K2*TMath::Power(rho1,Gamma2-1)/(Gamma2-1)+(1+a1)+K1*TMath::Power(rho1,Gamma1-1)/(Gamma1-1);
    //a1=-1-K1*TMath::Power(rho0,Gamma1-1)/(Gamma1-1)+1;
    //a2=-1-K2*TMath::Power(rho1,Gamma2-1)/(Gamma2-1)+(1+a1)+K1*TMath::Power(rho1,Gamma1-1)/(Gamma1-1);
    a3=-1-K3*TMath::Power(rho2,Gamma3-1)/(Gamma3-1)+(1+a2)+K2*TMath::Power(rho2,Gamma2-1)/(Gamma2-1);
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



double drdh(double eta,double r,double m){//energy density

  return -(r*(r-2.*m))/((m+4.*TMath::Pi()*r*r*r*pressure(eta))*(eta+1.));
}


double dmdh(double eta,double r,double m){
  
  
  return 4.*TMath::Pi()*r*r*e(eta)*drdh(eta,r,m);

}//dr/dh



void TOV_enthalpy(){


  
  int count=0;

  //define the crust
  vector<string> surfEOS;
  surfEOS.push_back("eosNV");
  surfEOS.push_back("eosFPS");
  int sEOSnum=0;//0 = NV   1 = FPS

  
  
  vector<double> v_p,v_d,v_e,v_eta,v_m,v_r,maxm,maxr,v_d1,v_p1,cr_e,cr_p,cr_d;

  //Maximum masses and Radii from rns
  vector<double> maxm1={1.8504,1.81534,1.7729,1.72214,1.66212,1.59184
		        ,1.5104,1.41704,1.31115,1.19264};
 
  vector<double> maxr1= {10.5166,10.6373,10.7572,10.875,10.989,11.0973
			 ,11.1975,11.2872,11.3636,11.424};


  //Ask for the Surface EOS - Begining
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
      if((stod(a4)*1.66*TMath::Power(10,-24))>(rho0*c*c/(kappa*G))) break;
      cr_e.push_back(stod(a1)*kappa*G/(c*c));
      cr_p.push_back(stod(a2)*kappa*G/(c*c*c*c));
      cr_d.push_back(stod(a4)*1.66*TMath::Power(10,-24)*kappa*G/(c*c));

      v_e.push_back(stod(a1)*kappa*G/(c*c));
      v_p.push_back(stod(a2)*kappa*G/(c*c*c*c));
      v_eta.push_back(stod(a3)/(c*c));
      v_d.push_back(stod(a4)*1.66*TMath::Power(10,-24)*kappa*G/(c*c));

      
    }
    count++;
  }
  cr_d.push_back(rho0);
  cr_p.push_back(pressure(rho0));
  cr_e.push_back(e(rho0));
  
  //Ask for the Surface EOS -END
  count=0;
  for(int dD=0;dD<15;dD++){//the loop is for different central densities 
    double stepsize=-0.001;
    double eta0=h((9.54753*1.66*TMath::Power(10,14)-dD*TMath::Power(10,14)*0.5)*kappa*G/(c*c))-1;
    double d0=(9.54753*1.66*TMath::Power(10,14)-dD*TMath::Power(10,14)*0.05)*kappa*G/(c*c);
    double r0=TMath::Power(10,-7)/(TMath::Power(kappa,0.5));
    double m0=TMath::Power(10,-15)*G/(TMath::Power(kappa,0.5)*c*c);
    double p0=pressure(eta0);
    //cout<<d0/(kappa*G/(c*c))<<" "<<rho(eta0)/(kappa*G/(c*c))<<" "<<p0/(kappa*G/(c*c*c*c))<<endl;
    double r,m,eta,p,d,k0,k1,k2,k3,l0,l1,l2,l3;

    eta=eta0;
    d=d0;
    //RK4 Begining
    while(d>TMath::Power(10,1)*kappa*G/(c*c)/*rho0*/){   
      if(eta > (eta0-0.001)){
	//this is for small distances
	double a1,a2,a3;
	a1=-1-K1*TMath::Power(rho0,Gamma1-1)/(Gamma1-1)+1;
	a2=-1-K2*TMath::Power(rho1,Gamma2-1)/(Gamma2-1)+(1+a1)+K1*TMath::Power(rho1,Gamma1-1)/(Gamma1-1);
	a3=-1-K3*TMath::Power(rho2,Gamma3-1)/(Gamma3-1)+(1+a2)+K2*TMath::Power(rho2,Gamma2-1)/(Gamma2-1);
	double drhodh= TMath::Power((eta0-a3)/(K3*(1/(Gamma3-1)+1)),1/(Gamma3-1))*(1./(Gamma3-1.)/(1.+1./(Gamma3-1.)))* TMath::Power((eta0-a3)/(K3*(1/(Gamma3-1)+1)),1/(Gamma3-1)-1)*1/(Gamma3-1)/(K3*(1/(Gamma3-1)+1))*(1+(a3+eta0/(Gamma3-1))/(1+1/(Gamma3-1)));

	r=TMath::Sqrt(3*(eta0-eta)/(2*TMath::Pi()*(e(eta0)+3.*p0)))*(1.-1./4.*(e(eta0)-3.*p0-3.*drhodh/5.)*(eta0-eta)/(e(eta0)+3.*p0));

	m=4.*TMath::Pi()*e(eta0)*r*r*r*(1-3.*(eta0-eta)*drhodh/(e(eta0)*5.))/3.;

	eta=eta-0.00001;

	p=pressure(eta0);
	d=rho(eta0);
	
      } 
      else{
	if(d>rho0){
	  //This is until the plotropic EOS stops
	  double mold=m;
	  double rold=r;
	  double pold=p;
	  double etaold=eta;
	  double rnew,mnew,etanew,pnew,dnew;

	  k0=drdh(etaold,rold,mold);
	  l0=dmdh(etaold,rold,mold);
	
	  k1= drdh(etaold+stepsize/2.,rold+stepsize*0.5*k0,mold+stepsize*0.5*l0);
	  l1= dmdh(etaold+stepsize/2.,rold+stepsize*0.5*k0,mold+stepsize*0.5*l0);


	  k2 = drdh(etaold+stepsize/2.,rold+stepsize*0.5*k1,mold+stepsize*0.5*l1);
	  l2 = dmdh(etaold+stepsize/2.,rold+stepsize*0.5*k1,mold+stepsize*0.5*l1);

	  k3 = drdh(etaold+stepsize,rold+stepsize*k2,mold+stepsize*l2);
	  l3 = dmdh(etaold+stepsize,rold+stepsize*k2,mold+stepsize*l2);
	
	  etanew=etaold+stepsize;
	
	  if(etanew<0)break;

	  rnew=rold+(k0+2.*k1+2.*k2+k3)*stepsize/6.;
	  mnew=mold+(l0+2.*l1+2.*l2+l3)*stepsize/6.;
	  pnew=pressure(eta);
	  dnew=rho(eta);
	  if(pnew<cr_p.at(0))break;

	
	  m=mnew;
	  r=rnew;
	  p=pnew;
	  eta=etanew;
	  d=dnew;
	}
	else{//this is for the crust
	  double mold=m;
	  double rold=r;
	  double dold=d;
	  double pold=p;
	  
	  double rnew,mnew,dnew,pnew;
	  if(d<rho0) stepsize=10./(TMath::Power(kappa,0.5));
	  if(d<TMath::Power(10,5.5)*kappa*G/(c*c)) stepsize=1./(TMath::Power(kappa,0.5));
	  if(d<TMath::Power(10,3.5)*kappa*G/(c*c)) stepsize=0.01/(TMath::Power(kappa,0.5));
	  if(d<TMath::Power(10,1.5)*kappa*G/(c*c)) stepsize=0.0001/(TMath::Power(kappa,0.5));
	
	  k0=dmdr1(rold,mold,pold,cr_p,cr_e);
	  l0=dpdr1(rold,mold,pold,cr_p,cr_e);
	  
	  if(pold+stepsize*0.5*l0<cr_p.at(0))break;
	  k1= dmdr1(rold+stepsize/2.,mold+stepsize*0.5*k0,pold+stepsize*0.5*l0,cr_p,cr_e);
	  l1= dpdr1(rold+stepsize/2.,mold+stepsize*0.5*k0,pold+stepsize*0.5*l0,cr_p,cr_e);

	  if(pold+stepsize*0.5*l1<cr_p.at(0))break;
	  k2 = dmdr1(rold+stepsize/2.,mold+stepsize*0.5*k1,pold+stepsize*0.5*l1,cr_p,cr_e);
	  l2 = dpdr1(rold+stepsize/2.,mold+stepsize*0.5*k1,pold+stepsize*0.5*l1,cr_p,cr_e);

	  if(pold+stepsize*l2<cr_p.at(0))break;
	  k3 = dmdr1(rold+stepsize,mold+stepsize*k2,pold+stepsize*l2,cr_p,cr_e);
	  l3 = dpdr1(rold+stepsize,mold+stepsize*k2,pold+stepsize*l2,cr_p,cr_e);
	  
	  rnew=rold+stepsize;
	  mnew=mold+(k0+2.*k1+2.*k2+k3)*stepsize/6.;
	  pnew=pold+(l0+2.*l1+2.*l2+l3)*stepsize/6.;
	  
	  if(pnew<cr_p.at(0))break;
	  int pos=crust(pnew,cr_p,cr_d);

	  double a=(TMath::Log10(cr_d.at(pos))-TMath::Log10(cr_d.at(pos+1)))/(TMath::Log10(cr_p.at(pos))-TMath::Log10(cr_p.at(pos+1)));
	  double b=TMath::Power(10,TMath::Log10(cr_d.at(pos+1))-(TMath::Log10(cr_d.at(pos))-TMath::Log10(cr_d.at(pos+1)))/(TMath::Log10(cr_p.at(pos))-TMath::Log10(cr_p.at(pos+1)))*TMath::Log10(cr_p.at(pos+1)));
	  dnew=b*TMath::Power(pnew,a);


	  m=mnew;
	  r=rnew;
	  p=pnew;
	  d=dnew;

	  
	}
	
      }
      if(count!=0 && count%2==0){
      	v_m.push_back((m*c*c*TMath::Power(kappa,0.5))/(Mo*G));
	v_r.push_back(r*TMath::Power(kappa,0.5)/100000.);
      }
      //cout<<(m*c*c*TMath::Power(kappa,0.5))/(Mo*G)<<" "<<r*TMath::Power(kappa,0.5)/100000.<<endl;
      count++;
    }
    //RK4 end
    count=0;
    maxm.push_back(v_m.at(v_m.size()-1));
    maxr.push_back(v_r.at(v_r.size()-1));
    //plot M-R relation

    v_m.clear();
    v_r.clear();
  }

  
  TGraph* m_r =new TGraph(maxm.size(),&maxr.at(0),&maxm.at(0));
  m_r->SetTitle("Max Mass vs Max Radius");
  m_r->GetYaxis()->SetTitle("Mass (Mo)");
  m_r->GetXaxis()->SetTitle("Radius (Km)");
  m_r->SetLineColor(4);
  m_r->SetMarkerColor(4);
  m_r->SetMarkerStyle(20);
  m_r->SetMarkerSize(0.5);
  TGraph* m_r1 =new TGraph(maxr1.size(),&maxr1.at(0),&maxm1.at(0));
  m_r1->SetLineColor(3);
  m_r1->SetMarkerColor(3);
  m_r1->SetMarkerStyle(20);
  m_r1->SetMarkerSize(0.5);

  /*TGraph* p_d =new TGraph(v_r.size(),&v_r.at(0),&v_m.at(0));
  p_d->SetTitle(" Mass vs Radius");
  p_d->GetYaxis()->SetTitle("Mass (Mo)");
  p_d->GetXaxis()->SetTitle("Radius (Km)");
  p_d->SetLineColor(4);
  p_d->SetMarkerColor(4);
  p_d->SetMarkerStyle(20);
  p_d->SetMarkerSize(0.5);
  TGraph* p_d1 =new TGraph(v_r1.size(),&v_r1.at(0),&v_m1.at(0));
  p_d1->SetLineColor(3);
  p_d1->SetMarkerColor(3);
  p_d1->SetMarkerStyle(20);
  p_d1->SetMarkerSize(0.5);
  */
  TCanvas* canvas =new TCanvas();
  canvas->cd();
  TImage *img = TImage::Create();

  //img->FromPad(5, 10, 10, 300, 200);
  img->FromPad(canvas);
  //p_d1->Draw("APL");
  //p_d->Draw("PL SAME");
  //p_d->Draw("AP");
  m_r->Draw("APL");
  m_r1->Draw("PL SAME");
  canvas->SaveAs("Max_Mass_Radius(N_D).png");
  /*canvas->SetLogx();
  canvas->SetLogy();
  p_d1->Draw("APL");
  //p_d->Draw("PL SAME");
  canvas->SaveAs("EOS2(N_D).png");
  */
}



