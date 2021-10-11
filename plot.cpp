using namespace std;


void plot(){
    int num;
  cout<<"Give the number of EOS"<<endl;
  cin>>num;
  int count=0;
  vector<double> v1,v2,v3,v4;
  TCanvas* canvas =new TCanvas();
  canvas->cd();
  TImage *img = TImage::Create();

  //img->FromPad(5, 10, 10, 300, 200);
  img->FromPad(canvas);
  canvas->SetLogx();
  canvas->SetLogy();


  //Ask for the Surface EOS - Begining
  ifstream myinputfile;
  string line;
  string a1,a2,a3,a4,a5;
  int data;

  
  //Plot limits 
  v4.push_back(TMath::Power(10,14));
  v2.push_back(TMath::Power(10,32));
  v4.push_back(TMath::Power(10,14));
  v2.push_back(TMath::Power(10,36.4));
  v4.push_back(TMath::Power(10,15.7));
  v2.push_back(TMath::Power(10,32));
  v4.push_back(TMath::Power(10,15.7));
  v2.push_back(TMath::Power(10,36.4));
  
  
  TGraph* eos =new TGraph(v2.size(),&v4.at(0),&v2.at(0));
  eos->SetTitle("Eos");
  eos->GetYaxis()->SetTitle("Pressure (dynes/cm^2)");
  eos->GetXaxis()->SetTitle("Energy density (gr/cm^3)");
  eos->SetLineColor(0);
  eos->SetMarkerColor(0);
  eos->SetMarkerStyle(20);
  eos->SetMarkerSize(1.5);

  eos->Draw("AP");

  v1.clear();
  v2.clear();
  v3.clear();
  v4.clear();
  
  //rest
  TRandom3* rand =new TRandom3();
  rand->SetSeed(0);


  
  for(int k=0;k<num;k++){
    
    string title="";
    title=title+"eosPol"+k;
    ifstream myinput;
    
    myinput.open(title,ios::in);
 
    while( getline(myinput,line) ){
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
	v1.push_back(stod(a1));
	v2.push_back(stod(a2));
	v3.push_back(stod(a3));
	v4.push_back(stod(a4)*1.66*TMath::Power(10,-24));
	
      }
      count++;
    }
    //Ask for the Surface EOS -END
    count=0;
  
 
    TGraph* eos1 =new TGraph(v2.size(),&v1.at(0),&v2.at(0));
    //eos1->SetLineColor(int(rand->Uniform(0.5,9.5)));
    eos1->SetLineColor(k+1);
    eos1->SetMarkerColor(k+1);
    eos1->SetMarkerStyle(20);
    eos1->SetMarkerSize(1.5);
    eos1->SetLineWidth(2);
    eos1->Draw("L SAME");
    v1.clear();
    v2.clear();
    v3.clear();
    v4.clear();
  }

  


  canvas->SaveAs("plotEOSHNG.png");
}
