//finds M-R relation including the crust and Max mass and Radius by reading an EOS file
using namespace std;
#include <cmath>
#define  G 6.6732*TMath::Power(10,-8)
#define  c 2.9979*TMath::Power(10,10)
#define  Mo 1.987*TMath::Power(10,33)
#define  eo  TMath::Power(10,15)
#define  kappa c*c/(G*eo)
#define  rho0 TMath::Power(10,14.3445)*kappa*G/(c*c)
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


///needed for the crust




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

}//drho/dr



void TOV_reading_my(){



  int num;
  cout<<"Give the number of EOS"<<endl;
  cin>>num;
  TCanvas* canvas =new TCanvas();
  canvas->cd();
  TImage *img = TImage::Create();
  //img->FromPad(5, 10, 10, 300, 200);
  img->FromPad(canvas);

  TRandom3* rand =new TRandom3();
  rand->SetSeed(0);


  //Plot limits
  vector<double> v2,v4;
  v4.push_back(6);
  v2.push_back(0);
  v4.push_back(6);
  v2.push_back(3);
  v4.push_back(18);
  v2.push_back(0);
  v4.push_back(18);
  v2.push_back(3);

  
  TGraph* lim =new TGraph(v2.size(),&v4.at(0),&v2.at(0));
  lim->SetTitle("Mass vs Radius Curves");
  lim->GetYaxis()->SetTitle("Mass (Mo)");
  lim->GetXaxis()->SetTitle("Radius (Km)");
  lim->SetLineColor(0);
  lim->SetMarkerColor(0);
  lim->SetMarkerStyle(20);
  lim->SetMarkerSize(1.5);
  
  lim->Draw("AP");
  
  int fileNum=0;

  for(int k=0;k<num;k++){
    int count=0;
    vector<double> rns_r,rns_m,rns_d,rns_p;
    ifstream input;
    input.open("nrs",ios::in);
    string lines;
  
    vector<double> v_p,v_d,v_e,v_m,v_r,maxm,maxr,v_d1,v_p1;
  
    string title="";
    title=title+"eosSA"+fileNum+"table";
    fileNum++;
    ofstream output;
    output.open(title,ios::out);
    output<<"Mass Radius MaxSoundSpeed CentralEnergyDensity"<<endl;

    
    //Ask for the Surface EOS - Begining
    ifstream myinputfile;
    myinputfile.open("eosPol"+std::to_string(k),ios::in);
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
	v_e.push_back(stod(a1)*kappa*G/(c*c));
	v_p.push_back(stod(a2)*kappa*G/(c*c*c*c));
	v_d.push_back(stod(a4)*1.66*TMath::Power(10,-24)*kappa*G/(c*c));
      }
      count++;
    }
    //Ask for the Surface EOS -END
    count=0;

    for(int dD=0;dD<int(data)-70;dD++){
      double stepsize=100./(TMath::Power(kappa,0.5));
      double r0=TMath::Power(10,-15)/(TMath::Power(kappa,0.5));
      double m0=TMath::Power(10,-15)*G/(c*c*TMath::Power(kappa,0.5));
      double d0=v_d.at(v_d.size()-1-dD);//(TMath::Power(10,14.9)+dD*TMath::Power(10,15)*0.1)*kappa*G/(c*c);
      double p0=v_p.at(v_p.size()-1-dD);//pressure(d0);
      double maxSoundSpeed=0;

      
      //cout<<d0/(1.66*TMath::Power(10,-24)*kappa*G/(c*c))<<endl;

    
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
	  m=4*TMath::Pi()*d0*r0*r0*r0/3-TMath::Abs((TMath::Log10(d0)-TMath::Log10(v_d.size()-1-dD-1))/(TMath::Log10(p0)-TMath::Log10(v_p.size()-1-dD-1)))*8/15*3.14*3.14*d0/p0*(d0+p0)*(d0+3*p0)*TMath::Power(r0,5);
	  d=d0-2/3*3.14*d0/p0*(d0+p0)*(d0+3*p0)*r0*r0;
	  p=p0-TMath::Abs((TMath::Log10(d0)-TMath::Log10(v_d.at(v_p.size()-1-dD-1)))/(TMath::Log10(p0)-TMath::Log10(v_p.at(v_p.size()-1-dD-1))))*2*3.14*(d0+p0)*(d0+3*p0)*r0*r0/3;
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

	  v_m.push_back((m*c*c*TMath::Power(kappa,0.5))/(Mo*G));
	  v_r.push_back(r*TMath::Power(kappa,0.5)/100000.);
	  v_p1.push_back(p/(kappa*G)*(c*c*c*c));
	  v_d1.push_back(d/(kappa*G)*(c*c));
	}
	count++;
      }
      //RK4 end
      count=0;
      //cout<<(m*c*c*TMath::Power(kappa,0.5))/(Mo*G)<<" "<<r*TMath::Power(kappa,0.5)/100000.<<endl;
      maxm.push_back(v_m.at(v_m.size()-1));
      maxr.push_back(v_r.at(v_r.size()-1));

      for (int i=1;i<(v_e.size()-1-dD);i++){
	double us=TMath::Sqrt((TMath::Log10(v_p.at(i))-TMath::Log10(v_p.at(i-1)))*v_p.at(i)/(TMath::Log10(v_d.at(i))-TMath::Log10(v_d.at(i-1)))/(v_e.at(i)+v_p.at(i)));
	if (us>maxSoundSpeed) maxSoundSpeed=us;
      }

      output<<v_m.at(v_m.size()-1)<<" "<<v_r.at(v_r.size()-1)<<" "<<maxSoundSpeed<<" "<<v_e.at(v_e.size()-1-dD)*(c*c)/(kappa*G)<<endl;
      
      //plot M-R relation
 
      v_m.clear();
      v_r.clear();
      v_p1.clear();
      v_d1.clear();
    }

 
    TGraph* m_r =new TGraph(maxm.size(),&maxr.at(0),&maxm.at(0));
    m_r->SetLineWidth(3.5);
    //m_r->SetLineColor(int(rand->Uniform(0.5,9.5)));
    m_r->SetLineColor(k+1);
    m_r->SetMarkerColor(k+1);
    m_r->SetMarkerStyle(20);
    m_r->SetMarkerSize(0.5);


    m_r->Draw("L SAME");
  canvas->SaveAs("Max_Mass_Radius(N_D).png");
  }


  
}



