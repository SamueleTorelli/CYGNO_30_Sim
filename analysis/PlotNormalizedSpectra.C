/*
These are the real maps

std::map<std::string,double> ElementMass ={ {"Cathodes",0.145746}, {"GEMs",193.536}, {"Lens",0.4995}, {"Rings",1114.74}, {"Sensors",0.1392}, {"Vessel",1102.24}  }; // in kilograms



std::map<std::string, std::map<std::string, double >> Contaminant = {  
  {"GEMs",          { {"U238",1.32e-2}, {"Th232",5.45e-3},                    {"U235",2.8e-2},  {"K40",6.31e-2},  {"Co60",2.34e-3},  {"Cs137",1.56e-3}   } },
  {"Lens",          { {"U238",1.23e-4}, {"Th232",4.07e-5},                                      {"K40",3.1e-4}                                           } },
  {"Vessel",        { {"U238",2.96e-4}, {"Th232",5.69e-5},                                      {"K40",7.12e-5}                                          } },
  {"Rings",         { {"U238",1.2e-5 }, {"Th232",4.1e-6},                                       {"K40",6.1e-5},   {"Co60",2.4e-4},   {"Cs137",2.9e-4}    } },
  {"Cathodes",      { {"U238",9.01e-1},                     {"U234",4.07e-5}                                                                             } },
  {"Sensors",       { {"U238",6.8e-3},  {"Th232",5.20e-3},                    {"U235",9.1e-4},  {"K40",3.5},                         {"Cs137",4.2e-4}    } }
};  // call Contaminant["Cat"]["U238"]



std::map<std::string, std::map<std::string, double >> NEvents = {  
  {"GEMs",      { {"U238",1e6}, {"Th232",1e6},               {"U235",1e6},   {"K40",1e7},  {"Co60",1e7},  {"Cs137",1e7}   } },
  {"Lens",      { {"U238",1e6}, {"Th232",1e6},                               {"K40",1e6}                                  } },
  {"Vessel",    { {"U238",1e7}, {"Th232",1e7},                               {"K40",1e7}                                  } },
  {"Rings",     { {"U238",1e6}, {"Th232",1e6},                               {"K40",1e7},  {"Co60",1e7},  {"Cs137",1e7}   } },
  {"Cathodes",  { {"U238",1e6},                 {"U234",1e6}                                                              } },
  {"Sensors",   { {"U238",1e7}, {"Th232",1e7},               {"U235",1e7},   {"K40",1e7},                 {"Cs137",1e7}   } }
};  // call Contaminant["Cat"]["U238"]
*/

//creating now fake testing maps

std::map<std::string,double> ElementMass ={ {"Cathodes",0.145746}, {"GEMs",193.536}, {"Rings",1114.74},  };

std::map<std::string, std::map<std::string, double >> Contaminant = {  
  {"GEMs",          { {"Co60",2.34e-3},  {"Cs137",1.56e-3}   } },
  {"Rings",         { {"Th232",4.1e-6}                       } },
  {"Cathodes",      { {"U238",9.01e-1}                       } }
  
};  // call Contaminant["Cat"]["U238"]

std::map<std::string, std::map<std::string, double >> NEvents = {  
  {"GEMs",          { {"Co60",1e7},  {"Cs137",1e7}   } },
  {"Rings",         { {"Th232",1e6}                  } },
  {"Cathodes",      { {"U238",1e6}                   } }
  
};  // call Contaminant["Cat"]["U238"]



void BuildDetectorMap(std::map<Int_t,TVector3>& aMap);
bool isWithin(const std::map<Int_t,TVector3>& aMap,const Double_t x,const Double_t y,const Double_t z,const Int_t Volnum);
  
void GenNormalizedSpectra(){

  //Build centroid detector map for fiducialization
  std::map<Int_t,TVector3> VolumeMap;
  BuildDetectorMap(VolumeMap);

  //declare tree variables
  Int_t evNumber;
  std::string* PartName = nullptr;
  std::vector<double>* EDep=nullptr;
  std::vector<double>* VolNum=nullptr;
  std::string* Nucleus= nullptr;
  std::vector<double>* X_Vertex=nullptr;
  std::vector<double>* Y_Vertex=nullptr;
  std::vector<double>* Z_Vertex=nullptr;
  
  
  TFile* f;
  TTree* tree;

  TH1D* temp_histo;
  TH1D* temp_histocut;
  
  std::vector<TH1D*> Histo;
  std::vector<TH1D*> HistoCut;
  
  for(auto& component: ElementMass){ //for components in the map Element Mass
    for(auto& val : Contaminant[component.first]){// for val (that ia map) that get for each detector component the map with <nuclide,contamination>  

      std::cout << "Mass of " << component.first << " is " << component.second << " Kg with contamination of \t" << val.first
		<<" equal to \t" << val.second <<" Bq/Kg nEvents: \t" << NEvents[component.first][val.first] <<"\n";
      
      f= TFile::Open(Form("elab_%s_%s_t0.root",(component.first).c_str(),(val.first).c_str()),"r");
      std::cout << "file " << Form("elab_%s_%s_t0.root", (component.first).c_str(), (val.first).c_str()) << "\n";

      tree = (TTree*)f->Get("elabHits");

      tree->SetBranchAddress("evNumber",&evNumber);
      tree->SetBranchAddress("PartName",&PartName);
      tree->SetBranchAddress("EDep_Out",&EDep);
      tree->SetBranchAddress("VolNnum_Out",&VolNum);
      tree->SetBranchAddress("Nucleus",&Nucleus);
      tree->SetBranchAddress("X_Vertex",&X_Vertex);
      tree->SetBranchAddress("Y_Vertex",&Y_Vertex);
      tree->SetBranchAddress("Z_Vertex",&Z_Vertex);

      temp_histo = new TH1D(Form("%s_%s",(component.first).c_str(), (val.first).c_str() ), Form("%s_%s",(component.first).c_str(), (val.first).c_str() ),900,0,2000);
      temp_histocut = new TH1D(Form("%s_%s_cut",(component.first).c_str(), (val.first).c_str() ), Form("%s_%s_cut",(component.first).c_str(), (val.first).c_str() ),900,0,2000);
      
      for(int i=0;i< tree->GetEntries();i++){
	tree->GetEntry(i);
	
	if(i%40000 == 0) std::cout << i << "/" << tree->GetEntries()<< std::endl;
	
	if(std::strcmp((*PartName).c_str(),"e+")==0 || std::strcmp((*PartName).c_str(),"e-")==0 ){
	  for(int j=0; j<(*EDep).size();j++ ){
	    
	    if(isWithin( VolumeMap,(*X_Vertex)[j],(*Y_Vertex)[j],(*Z_Vertex)[j],(*VolNum)[j] )){
	      temp_histocut->Fill( (*EDep)[j]*1000 );
	    }//chiudo if within the volume
	    
	    temp_histo->Fill( (*EDep)[j]*1000 ); 
	    
	  }//chiudo for on vector
	}//chiudo if on e-
	
	
      }//tree over entries

      temp_histo->Sumw2();
      temp_histocut->Sumw2();
      
      temp_histo->Scale( val.second*component.second*60*60*24*365/NEvents[component.first][val.first] );
      temp_histocut->Scale( val.second*component.second*60*60*24*365/NEvents[component.first][val.first] );

      Histo.push_back(temp_histo);
      HistoCut.push_back(temp_histocut);
      
                  
    }// for on contaminant map with access to the element
  }//for on mass element

  THStack* Hstack("Hstack","Hstack");
  THStack* Hstack_cut("Hstack_cut","Hstack_cut");
  
  TFile* outDef = new TFile("NormalizedHisto.root","recreate");
  outDef->cd();

  for(int i=0;i<Histo.size();i++){
    Histo[i]->SetMarkerColor(30+i);
    Histo[i]->SetLineColor(30+i);
    Histo[i]->Write();
    Hstack->Add(Histo[i]);

    HistoCut[i]->SetMarkerColor(30+i);
    HistoCut[i]->SetLineColor(30+i);
    HistoCut[i]->Write();
    Hstack_cut->Add(HistoCut[i]);
  }

  Hstack->Write();
  Hstack_cut->Write();
  
  outDef->Save();
  outDef->Close();
  
  
}















void BuildDetectorMap(std::map<Int_t,TVector3>& aMap){

  Double_t VolumeSize_x=500;
  Double_t VolumeSize_y=800;
  Double_t VolumeSize_z = 500;

  Double_t DetectorSpace=3;
  
  int counter=0;
  
  for(int i=-12;i<13;i++){
    for(int j=-1;j<2;j++){
      aMap[counter] = TVector3( i*(VolumeSize_x+DetectorSpace) , j*(VolumeSize_y+DetectorSpace) , VolumeSize_z/2  );
      counter++;
    }
  }
  
  
  for(int i=-12;i<13;i++){
    for(int j=-1;j<2;j++){
      aMap[counter] = TVector3( i*(VolumeSize_x+DetectorSpace) , j*(VolumeSize_y+DetectorSpace) , -VolumeSize_z/2  );
      counter++;
    }
  }

}


bool isWithin(const std::map<Int_t,TVector3>& aMap,const Double_t x,const Double_t y,const Double_t z,const Int_t Volnum){

  Double_t VolumeSize_x=500;
  Double_t VolumeSize_y=800;
  Double_t VolumeSize_z = 500;

  Double_t FiducialCut_xy = 20;
  Double_t FiducialCut_z = 20; 

  Double_t dx = abs(x-aMap.at(Volnum).X() );
  Double_t dy = abs(y-aMap.at(Volnum).Y() );
  Double_t dz = abs(z-aMap.at(Volnum).Z() );

  bool cond = (dx <= (VolumeSize_x/2-FiducialCut_xy) && dy <= (VolumeSize_y/2-FiducialCut_xy) && dz <= (VolumeSize_z/2-FiducialCut_z) );
    
  //std::cout << "VolNum: " << Volnum << "  Center:  " << aMap.at(Volnum).X()<<"\t" << aMap.at(Volnum).Y()<<"\t" << aMap.at(Volnum).Z() << "  dist: " << dx<<"\t" << dy<<"\t" << dz << "\t" << cond << "\n";
  
  return cond;
}
