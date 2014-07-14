{     
  // include path for RooFit
  TString rfitpath("$ROOFITSYS/include");
  TString path = gSystem->GetIncludePath();
  path += "-I. -I$ROOTSYS/src -I";
  path += rfitpath;
  gSystem->SetIncludePath(path.Data());

  gROOT->Macro("CPlot.cc+");
  gROOT->Macro("MitStyleRemix.cc+");
  	    
  // Show which process needs debugging
  gInterpreter->ProcessLine(".! ps |grep root.exe");  
}
