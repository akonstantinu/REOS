#define n0 0.16*TMath::Power(10,39)//1/cm^3
#define mN 1.675*TMath::Power(10,-24)//gr
#define G 6.6732*TMath::Power(10,-8)
#define c 2.9979*TMath::Power(10,10)
#define Mo 1.987*TMath::Power(10,33)
#define eo  TMath::Power(10,15)
#define Mb  1.66*TMath::Power(10,-24)
#define kappa c*c/(G*eo)
#define hbar 1.054571628*TMath::Power(10,-34)



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
  


double f(double e,double a1,double a2,double a3,double a4,double a5,double a6){
  
  double speed;
  if(e/(mN*n0)>1.5){
    speed=a1*TMath::Exp(-0.5*TMath::Power(e/(mN*n0)-a2,2)/(a3*a3))+a6+(1./3.-a6)/(1.+TMath::Exp(-a5*((e/(mN*n0)-a4))));
  }
  else{
    //speed=(1-0.52)*(hbar*hbar)*TMath::Power(3.*TMath::Power(TMath::Pi(),2)*e/(mN),2./3.)/(0.92*3.*mN*mN*c*c);
    speed=(1-0.52)*TMath::Power(3.*TMath::Power(TMath::Pi(),2)*e/(mN),2./3.)/(0.92*3.*2.2676*TMath::Power(10,27));
  }


  if (speed<0){
    return 0;
  }
  else{
    return speed;
  }
}

double pressure(double e,double e0,double p0,double a1,double a2,double a3,double a4,double a5,double a6){
  double dx=(e-e0)/1000.;
  double integ=0;
  double ener=e0;
  integ=f(e0,a1,a2,a3,a4,a5,a6)+f(e,a1,a2,a3,a4,a5,a6);

    
  for(int i=1;i<1000-1;i++){
    ener+=dx;
    if(i%2==0){
      integ+=2.*f(ener,a1,a2,a3,a4,a5,a6);
    }
    else{
      integ+=4.*f(ener,a1,a2,a3,a4,a5,a6);
    }
  }

  return p0+integ*dx*c*c/3.;
}

double g(double e,double e0, double p0,double a1,double a2,double a3,double a4,double a5,double a6){
  
  return 1./(e+pressure(e,e0,p0,a1,a2,a3,a4,a5,a6)/(c*c));
}


double enthalpy(double e,double e0,double p0,double h0,double a1,double a2,double a3,double a4,double a5,double a6){
  double dx=(e-e0)/100.;
  double integ=0;
  double ener=e0;
  integ=g(e0,e0,p0,a1,a2,a3,a4,a5,a6)+g(e,e0,p0,a1,a2,a3,a4,a5,a6);

    
  for(int i=1;i<100-1;i++){
    ener+=dx;
    if(i%2==0){
      integ+=2.*g(ener,e0,p0,a1,a2,a3,a4,a5,a6);
    }
    else{
      integ+=4.*g(ener,e0,p0,a1,a2,a3,a4,a5,a6);
    }
  }
  return h0+integ*(pressure(e,e0,p0,a1,a2,a3,a4,a5,a6)-pressure(e0,e0,p0,a1,a2,a3,a4,a5,a6))/(3*100.);
}




void grief(){

  //asks for EOS num
  int num;
  cout<<"Give the number of EOS"<<endl;
  cin>>num;

  int interv=40;
  double e0=3.5*TMath::Power(10,15); /* kappa*G/(c*c) *///central energy density 
  double a1,a2,a3,a4,a5,a6;

  ofstream eosTable;
  eosTable.open("EosTable",ios::out);
  eosTable<<"EOS  a1  a2  a3  a4  a5  a6  central_energy"<<endl;
  
  TRandom3* rand =new TRandom3();
  rand->SetSeed(0);
  
  int fileNum=0;
  
  //Loop for EOS
  while(fileNum<num){
    
    vector<double> v_p,v_d,v_e,v_h;
    vector<double> vEnergyDensity,vEnthalpy,vDensity,vPressure; //non dimensional
    
    //choose surface EOS
    vector<string> surfEOS;
    surfEOS.push_back("eosNV");
    surfEOS.push_back("eosFPS");
    surfEOS.push_back("eosBBB2");
    int sEOSnum=2;//0 = NV   1 = FPS 2 = BBB2


    ifstream myinputfile;
    myinputfile.open(surfEOS.at(sEOSnum),ios::in);
    string line;
    string read1,read2,read3,read4,read5;
    int data;
    int count=0;
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
	if(stod(read4)>1.1*n0) break;
	v_e.push_back(stod(read1));
	v_p.push_back(stod(read2));
	v_h.push_back(stod(read3));
	v_d.push_back(stod(read4)*Mb);
	vPressure.push_back(stod(read2)*kappa*G/(c*c*c*c));
	vEnergyDensity.push_back(stod(read1)*kappa*G/(c*c));
	vEnthalpy.push_back(stod(read3)/(c*c));
	vDensity.push_back(stod(read4)*Mb*kappa*G/(c*c));

      }
      count++;
    }


     
    double e=v_e.at(v_e.size()-1);
    double de=(e0-e)/double(interv);      
    bool BoolSpeed=true;

    a1=/*0.820174;*/rand->Uniform(0.1,1.5);
    a2=/*6.93428;*/rand->Uniform(1.5,12);
    a3=/*2.63873;*/rand->Uniform(0.05*a2,2*a2);
    a4=/*24.8329;*/rand->Uniform(1.5,37);
    a5=/*0.894644;*/rand->Uniform(0.1,1);

    double u =TMath::Sqrt((v_p.at(v_p.size()-1)-v_p.at(v_p.size()-2))/(v_e.at(v_e.size()-1)-v_e.at(v_e.size()-2)))/(c);

    a6=/*-0.0395129;*/(u*u-a1*TMath::Exp(-0.5*TMath::Power(v_e.at(v_e.size()-1)/(mN*n0)-a2,2)/(a3*a3))-1./(3.*(1.+TMath::Exp(-a5*((v_e.at(v_e.size()-1)/(mN*n0)-a4))))))/(1.-1./(1.+TMath::Exp(-a5*((v_e.at(v_e.size()-1)/(mN*n0)-a4)))));
    //a6=/*-0.0395129;*/(0.05-a1*TMath::Exp(-0.5*TMath::Power(TMath::Power(10,14)/(mN*n0)-a2,2)/(a3*a3))-1./(3.*(1.+TMath::Exp(-a5*(3.*TMath::Power(10,14)/(mN*n0)-a4)))))/(1.-1./(1.+TMath::Exp(-a5*((TMath::Power(10,14)/(mN*n0)-a4)))));
    
    //cout<<u<<endl;
    if(f(TMath::Power(10,16),a1,a2,a3,a4,a5,a6)>(1./3.+0.001) || f(TMath::Power(10,16),a1,a2,a3,a4,a5,a6)<(1./3.-0.001)) continue;
 
    for(int i=0;i<interv;i++){
      e+=de;
      double pr=pressure(e,v_e.at(v_e.size()-1),v_p.at(v_p.size()-1),a1,a2,a3,a4,a5,a6);
      double h=enthalpy(e,v_e.at(v_e.size()-1),v_p.at(v_p.size()-1),v_h.at(v_h.size()-1),a1,a2,a3,a4,a5,a6);
      double n= (e+pr/(c*c))*exp(-h/(c*c))/(Mb);
      //cout<<e<<" "<<pr<<" "<<f(e,a1,a2,a3,a4,a5,a6)<<" "<<TMath::Sqrt((pr-v_p.at(v_p.size()-1))/(e-v_e.at(v_e.size()-1)))/(c)<<endl;
      v_p.push_back(pr);
      v_e.push_back(e);
      v_d.push_back(n*Mb);
      v_h.push_back(h);
      vPressure.push_back(pr*kappa*G/(c*c*c*c));
      vEnergyDensity.push_back(e*kappa*G/(c*c));
      vDensity.push_back(n*Mb*kappa*G/(c*c));
      vEnthalpy.push_back(h/(c*c));
      
      if(f(e,a1,a2,a3,a4,a5,a6)>1 || f(e,a1,a2,a3,a4,a5,a6)<=0){
	BoolSpeed=false;
	break;
      }
      if((e/(mN*n0)<1.5) && f(e,a1,a2,a3,a4,a5,a6)>0.163){
	BoolSpeed=false;
	break;
      }

    }
    
    if(!BoolSpeed) continue;
    //if(maxMassOrRadius(vDensity,vPressure,vEnergyDensity,0)<1.97) continue;
    //Write EOS parameters Begining
    string title="";
    title=title+"eosSA"+fileNum;
    eosTable<<title<<" "<<a1<<" "<<a2<<" "<<a3<<" "<<a4<<" "<<a5<<" "<<a6<<" "<<e0<<" "<<endl;
    
    ofstream output;
    output.open(title,ios::out);
    //Write EOS parameters End
    //output Data
    output<<v_d.size()<<endl;
    for(int i=0;i<v_d.size();i++){
      output<<v_e.at(i)<<" "<<v_p.at(i)<<" "<<v_h.at(i)<<" "<<v_d.at(i)/(Mb)<<endl;
    }
    output.close();
    v_e.clear();
    v_p.clear();
    v_d.clear();
    v_h.clear();
    fileNum++;
  }
}
