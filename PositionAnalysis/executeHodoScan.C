{
  gROOT->ProcessLine(".L hodoScan.C");
  TTree* tree=_file0->Get("tree_passedEvents");
  hodoScan t(tree);
  t.Loop();
}
