std::map<std::string,double> ElementMass ={ {"Cat",0.145746}, {"GEM",193.536}, {"Lens",0.4995}, {"Rings",1114.74}, {"Sens",0.1392}, {"Vessel",1102.24}  }; // in kilograms

std::map<std::string, std::map<std::string, double >> Contaminant = {  
  {"GEM",      { {"U238",1.32e-2}, {"Th232",5.45e-3},  {"U235",2.8e-2},  {"K40",6.31e-2},  {"Co60",2.34e-3},  {"Cs137",1.56e-3}   } },
  {"Lens",     { {"U238",1.23e-4}, {"Th232",4.07e-5},                    {"K40",3.1e-4}                                           } },
  {"Vessel",   { {"U238",2.96e-4}, {"Th232",5.69e-5},                    {"K40",7.12e-5}                                          } },
  {"Rings",    { {"U238",1.2e-5 }, {"Th232",4.1e-6},                     {"K40",6.1e-5},   {"Co60",2.4e-4},   {"Cs137",2.9e-4}    } },
  {"Cat",      { {"U238",9.01e-1}, {"U234",4.07e-5}                                                                               } },
  {"Sens",     { {"U238",6.8e-3},  {"Th232",5.20e-3},  {"U235",9.1e-4},  {"K40",3.5},                         {"Cs137",4.2e-4}    } }
};  // call Contaminant["Cat"]["U238"]



void GenNormalizedSpectra(){

  for(auto& component: ElementMass){
    for(auto& val : Contaminant[component.first]){

      std::cout << "Mass of " << component.first << " is " << component.second << " with contamination of " << val.first <<" equal to  " << val.second << "\n";
      
    }// for on contaminant map with access to the element
  }//for on mass element

  
  
}
