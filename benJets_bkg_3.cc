#include "Pythia8/Pythia.h"
#include "Pythia8/HeavyIons.h"
#include "Pythia8Plugins/FastJet3.h"

#include "TFile.h"
#include "TH1D.h"
#include "TH2D.H"
#include "TString.h"
#include "TTree.h"

#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
//#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/Selector.hh"
#include "fastjet/tools/GridMedianBackgroundEstimator.hh"
//#include "fastjet/tools/JetMedianBackgroundEstimator.hh"
//#include "fastjet/tools/Subtractor.hh"
#include "fastjet/contrib/IterativeConstituentSubtractor.hh"
#include "fastjet/contrib/SoftDrop.hh"



using namespace Pythia8;
using namespace fastjet;

int main(){
  //  int nEvents = 1;

  bool print_debug = false;

  double pTmin_jet = 20;
  double pTmin_hadron = 2;

  double R = 0.4;
  double etaMax_hadron = 2.0;
  double etaMax_jet = 1.0;
  
  Pythia pythia;
  Event& event = pythia.event;


  bool AAMode = false;
  
  pythia.readString("HardQCD:all = on");
  pythia.readString("PhaseSpace:pTHatMin = 20.");

  if(AAMode){
    pythia.readString("Beams:idA = 1000791970");
    pythia.readString("Beams:idB = 1000791970");
  }
  pythia.readString("Beams:eCM = 200.");

  //  pythia.readString("HeavyIon:SigFitErr = "
  //		    "0.02,0.02,0.1,0.05,0.05,0.0,0.1,0.0");
  //  pythia.readString("HeavyIon:SigFitDefPar = "
  //		    "17.24,2.15,0.33,0.0,0.0,0.0,0.0,0.0");
  //  pythia.readString("HeavyIon:SigFitNGen = 20");
  
  
  pythia.init();

  JetDefinition jetDef( antikt_algorithm, R);
  //JetDefinition jetDefCA( kt_algorithm, 0.4);
  //  AreaDefinition areaDef(active_area,GhostedAreaSpec(etaMax+1));
  //  AreaDefinition areaDef_bkg(active_area_explicit_ghosts,GhostedAreaSpec(etaMax+1));
  contrib::SoftDrop sd(0.0,0.1);
  contrib::SoftDrop sd2(0.0,0.1);

  GridMedianBackgroundEstimator bge(etaMax_hadron,0.5);

  contrib::IterativeConstituentSubtractor sub;
  sub.set_distance_type(contrib::ConstituentSubtractor::deltaR);
  vector<double> max_distance;
  max_distance.push_back(0.2);
  max_distance.push_back(0.1);
  vector<double> alpha;
  alpha.push_back(1.0);
  alpha.push_back(1.0);
  sub.set_parameters(max_distance,alpha);
  sub.set_ghost_removal(true);
  sub.set_ghost_area(0.004);
  sub.set_max_eta(etaMax_hadron);
  sub.set_background_estimator(&bge);

  Selector sel_max_pt = SelectorPtMax(10);
  sub.set_particle_selector(&sel_max_pt);

  sub.initialize();

  cout << sub.description() << endl;
  
  TFile *outfile;
  if(AAMode) outfile =  new TFile("v3_bkg_benJets_AuAu.root","RECREATE");
  else   outfile =  new TFile("v3_bkg_benJets.root","RECREATE");

  TH1D *jet_pT = new TH1D("jet_pT","Jet p_{T};p_{T} [GeV/c];1/N_{jet} dN_{jet}/dp_{T} [GeV/c]^{-1}",500,0.0,50.0);
  TH1D *jet_zg = new TH1D("jet_zg","Jet z_{g};z_{g};1/N_{jet} dN_{jet}/dz_{g}",20,0.0,0.5);
  TH1D *jet_rg = new TH1D("jet_rg","Jet R_{g};R_{g};1/N_{jet} dN_{jet}/dR_{g}",20,0.0,0.5);
  TH1D *jet_time = new TH1D("jet_time","Jet Formation Time;ln(#tau_{f}) [fm/c];1/N_{jet} dN_{jet}/dln(#tau_{f}) [fm/c]^{-1}",20,-2,4);

  TH1D *jet_pT_sub = new TH1D("jet_pT_sub","BG Subtracted Jet p_{T};p_{T} [GeV/c];1/N_{jet} dN_{jet}/dp_{T} [GeV/c]^{-1}",500,0.0,50.0);
  TH1D *jet_zg_sub = new TH1D("jet_zg_sub","BG Subtracted Jet z_{g};z_{g};1/N_{jet} dN_{jet}/dz_{g}",20,0.0,0.5);
  TH1D *jet_rg_sub = new TH1D("jet_rg_sub","BG Subtracted Jet R_{g};R_{g};1/N_{jet} dN_{jet}/dR_{g}",20,0.0,0.5);
  TH1D *jet_time_sub = new TH1D("jet_time_sub","BG Subtracted Jet Formation Time;ln(#tau_{f}) [fm/c];1/N_{jet} dN_{jet}/dln(#tau_{f}) [fm/c]^{-1}",20,-2,4);

  TH1D *jet_pT_pT = new TH1D("jet_pT_pT","Jet p_{T} hadron p_{T}>2 GeV/c;p_{T} [GeV/c];1/N_{jet} dN_{jet}/dp_{T} [GeV/c]^{-1}",500,0.0,50.0);
  TH1D *jet_zg_pT = new TH1D("jet_zg_pT","Jet z_{g} hadron p_{T}>2 GeV/c;z_{g};1/N_{jet} dN_{jet}/dz_{g}",20,0.0,0.5);
  TH1D *jet_rg_pT = new TH1D("jet_rg_pT","Jet R_{g} hadron p_{T}>2 GeV/c;R_{g};1/N_{jet} dN_{jet}/dR_{g}",20,0.0,0.5);
  TH1D *jet_time_pT = new TH1D("jet_time_pT","Jet Formation Time hadron p_{T}>2 GeV/c;ln(#tau_{f}) [fm/c];1/N_{jet} dN_{jet}/dln(#tau_{f}) [fm/c]^{-1}",20,-2,4);

  TH2D *jet_zg_correlation = new TH2D("jet_zg_correlation","Correlation between z_{g} with and without background subtraction;z_{g} bkg subtracted; z_{g} no subtraction",200,0.0,0.5,200,0.0,0.5);
  TH2D *jet_zg_correlation_hadCuts = new TH2D("jet_zg_correlation_hadCuts","Correlation between z_{g} with and without hadron p_{T} cuts;z_{g} p_{T}^{hadron}>2 GeV/c; z_{g} p_{T}^{hadron}>0 GeV/c",200,0.0,0.5,200,0.0,0.5);


  TH2D *jet_pT_correlation = new TH2D("jet_pT_correlation","Correlation between p_{T} with and without background subtraction;p_{T} bkg subtracted; p_{T} no subtraction",500,0.0,50.0,500,0.0,50.0);
  TH2D *jet_pT_correlation_noCut = new TH2D("jet_pT_correlation_noCut","Correlation between p_{T} with and without background subtraction;p_{T} bkg subtracted; p_{T} no subtraction",500,0.0,50.0,500,0.0,50.0);
  TH2D *jet_pT_correlation_hadCuts = new TH2D("jet_pT_correlation_hadCuts","Correlation between p_{T} with and without hadron p_{T} cuts;p_{T} p_{T}^{hadron}>2 GeV/c; p_{T} p_{T}^{hadron}>0 GeV/c",500,0.0,50.0,500,0.0,50.0);

  TTree *tree = new TTree("tree","");
  double data[24];
  //,"pT:phi:eta:pTg:phig:etag:zg:rg:pTs:phis:etas:pTsg:phisg:etasg:zgs:rgs");
  tree->Branch("pT",&data[0]);
  tree->Branch("phi",&data[1]);
  tree->Branch("eta",&data[2]);
  tree->Branch("pTg",&data[3]);
  tree->Branch("phig",&data[4]);
  tree->Branch("etag",&data[5]);
  tree->Branch("zg",&data[6]);
  tree->Branch("rg",&data[7]);
  tree->Branch("pTs",&data[8]);
  tree->Branch("phis",&data[9]);
  tree->Branch("etas",&data[10]);
  tree->Branch("pTsg",&data[11]);
  tree->Branch("phisg",&data[12]);
  tree->Branch("etasg",&data[13]);
  tree->Branch("zgs",&data[14]);
  tree->Branch("rgs",&data[15]);
  tree->Branch("pTpT",&data[16]);
  tree->Branch("phipT",&data[17]);
  tree->Branch("etapT",&data[18]);
  tree->Branch("pTpTs",&data[19]);
  tree->Branch("phipTs",&data[20]);
  tree->Branch("etapTs",&data[21]);
  tree->Branch("zgpT",&data[22]);
  tree->Branch("rgpT",&data[23]);
  
  
  int iEvent = 0;
  int nJets_gt_20 = 0;
  int nJets_sub_gt_20 = 0;
  int nWrong = 0;
  while (nJets_sub_gt_20 < 10000){
  //  for(int iEvent =0; iEvent < nEvents; iEvent++){
    if(!pythia.next()) continue;
    cout << "working on event " << iEvent << endl;
    int nCharged = 0;
    int nJets_20_30 = 0;
    vector<PseudoJet> particles;
    vector<PseudoJet> particles_no_low_pT;
    for(int i=0; i < event.size(); ++i){
      if(event[i].isFinal() && event[i].isCharged()) nCharged++;
      if(event[i].isFinal() && event[i].isVisible() && abs(event[i].eta()) < etaMax_hadron) particles.push_back( PseudoJet(event[i].px(), event[i].py(), event[i].pz(), event[i].e() ));
      if(event[i].isFinal() && event[i].isVisible() && abs(event[i].eta()) < etaMax_hadron && event[i].pT() > pTmin_hadron) particles_no_low_pT.push_back( PseudoJet(event[i].px(), event[i].py(), event[i].pz(), event[i].e() ));
    }//end loop over tracks

    
    ClusterSequence cs(particles, jetDef);
    vector<PseudoJet> sortedJets = sorted_by_pt( cs.inclusive_jets() );

    ClusterSequence cspT(particles_no_low_pT, jetDef);
    vector<PseudoJet> pTJets = sorted_by_pt( cspT.inclusive_jets() );
    
    bge.set_particles(particles);

    vector<PseudoJet> corrEvent = sub.subtract_event(particles);

    ClusterSequence csCorr(corrEvent, jetDef);
    vector<PseudoJet> corrJets = sorted_by_pt( csCorr.inclusive_jets() );
    
    for (unsigned int i=0; i<sortedJets.size(); i++){
      //      cout << "working on base jet " << i << "/" << sortedJets.size() << endl;
      PseudoJet sdJet = sd(sortedJets[i]);

      double zg = sdJet.structure_of<contrib::SoftDrop>().symmetry();
      double rg = sdJet.structure_of<contrib::SoftDrop>().delta_R();
      double energy =sdJet.E();
      double t_form = 1.0/(zg*(1.0-zg)*rg*rg*energy);
      double log_t_form = log(t_form);

      jet_pT->Fill(sdJet.pt());
      
      if(sdJet.pt() >= pTmin_jet && sdJet.pt() <= 30 && fabs(sdJet.eta()) < 1){
	jet_zg->Fill(zg);
	jet_rg->Fill(rg);
	jet_time->Fill(log_t_form);
	nJets_gt_20++;
      }
    }

    for( unsigned int i=0; i<pTJets.size(); i++){
      //      cout << "working on pT jet " << i << "/" << pTJets.size() << endl;
      PseudoJet sdJet_pT = sd(pTJets[i]);

      double minDeltaR = 100000;
      int minK = -1;
      for(unsigned int k=0; k<sortedJets.size(); k++){
	if(pTJets[i].delta_R(sortedJets[k]) < minDeltaR){
	  minDeltaR = pTJets[i].delta_R(sortedJets[k]);
	  minK = k;
	}
      }



      
      double zg_pT = sdJet_pT.structure_of<contrib::SoftDrop>().symmetry();
      double rg_pT = sdJet_pT.structure_of<contrib::SoftDrop>().delta_R();
      double energy_pT =sdJet_pT.E();
      double t_form_pT = 1.0/(zg_pT*(1.0-zg_pT)*rg_pT*rg_pT*energy_pT);
      double log_t_form_pT = log(t_form_pT);

      PseudoJet sdJet = sd(sortedJets[minK]);
      //PseudoJet sdJet = sortedJets[i];

      double zg = sdJet.structure_of<contrib::SoftDrop>().symmetry();
      double rg = sdJet.structure_of<contrib::SoftDrop>().delta_R();
      double energy =sdJet.E();
      double t_form = 1.0/(zg*(1.0-zg)*rg*rg*energy);
      double log_t_form = log(t_form);

      
      jet_pT_pT->Fill(sdJet_pT.pt());
      
      if(sdJet_pT.pt() >= pTmin_jet && sdJet_pT.pt() <= 30 && fabs(sdJet_pT.eta()) < 1){
	jet_zg_pT->Fill(zg_pT);
	jet_rg_pT->Fill(rg_pT);
	jet_time_pT->Fill(log_t_form_pT);

	jet_zg_correlation_hadCuts->Fill(zg_pT,zg);

	jet_pT_correlation_hadCuts->Fill(sdJet_pT.pt(),sdJet.pt());
	

      }

    }
    
    for(unsigned int i=0; i<corrJets.size(); i++){
      if(corrJets[i] == 0) continue;
      //      cout << "working on subtracted jet " << i << "/" << corrJets.size() << endl;
      PseudoJet sdJet_sub = sd2(corrJets[i]);
      
      
      double zg_sub = sdJet_sub.structure_of<contrib::SoftDrop>().symmetry();
      double rg_sub = sdJet_sub.structure_of<contrib::SoftDrop>().delta_R();
      double energy_sub =sdJet_sub.E();
      double t_form_sub = 1.0/(zg_sub*(1.0-zg_sub)*rg_sub*rg_sub*energy_sub);
      double log_t_form_sub = log(t_form_sub);

      double minDeltaR = 100000;
      int minJ = -1;
      for(unsigned int j=0; j<sortedJets.size(); j++){
	if(corrJets[i].delta_R(sortedJets[j]) < minDeltaR){
	  minDeltaR = corrJets[i].delta_R(sortedJets[j]);
	  minJ = j;
	}
      }

      double pT;
      double eta;
      double phi;
      double zg;
      double rg;
      double energy;
      double t_form;
      double log_t_form;
      
      if( minJ != -1){
      
	PseudoJet sdJet = sd(sortedJets[minJ]);
	//PseudoJet sdJet = sortedJets[i];

	pT = sdJet.pt();
	eta = sdJet.eta();
	phi = sdJet.phi();
	zg = sdJet.structure_of<contrib::SoftDrop>().symmetry();
	rg = sdJet.structure_of<contrib::SoftDrop>().delta_R();
	energy =sdJet.E();
	t_form = 1.0/(zg*(1.0-zg)*rg*rg*energy);
	log_t_form = log(t_form);

      }
      
      double minDeltaRK = 100000;
      int minK = -1;
      for(unsigned int k=0; k<pTJets.size(); k++){
	if(sortedJets[minJ].delta_R(pTJets[k]) < minDeltaRK){
	  minDeltaRK = sortedJets[minJ].delta_R(pTJets[k]);
	  minK = k;
	}
      }

      double pT_pT = -999;
      double eta_pT = -999;
      double phi_pT = -999;
      double zg_pT = -999;
      double rg_pT = -999;
      double energy_pT = -999;
      double t_form_pT = -999;
      double log_t_form_pT = -999;

      //      cout << "looped over pT Jets" << endl;
      
      if( minK!= -1){
      
	PseudoJet sdJet_pT = sd(pTJets[minK]);
	
	pT_pT = sdJet_pT.pt();
	eta_pT = sdJet_pT.eta();
	phi_pT = sdJet_pT.phi();
	zg_pT = sdJet_pT.structure_of<contrib::SoftDrop>().symmetry();
	rg_pT = sdJet_pT.structure_of<contrib::SoftDrop>().delta_R();
	energy_pT =sdJet_pT.E();
	t_form_pT = 1.0/(zg_pT*(1.0-zg_pT)*rg_pT*rg_pT*energy_pT);
	log_t_form_pT = log(t_form_pT);
      }
      //      cout << "set pT jet stuff" << endl;
      /*
      cout << "jet constituents: " << sortedJets[minJ].constituents().size() << endl;
      cout << "groomed jet constituents: " << sdJet.constituents().size() << endl;
      cout << "zg: " << zg << endl;
      cout << "bg sub jet constituents: " << corrJets[i].constituents().size() << endl;
      cout << "bg sub groomed jet constituents: " << sdJet_sub.constituents().size() << endl;
      cout << "bg sub zg: " << zg_sub << endl;
      */

      jet_pT_sub->Fill(sdJet_sub.pt());

      //      cout << "filled pT hist" << endl;
      
      if(sdJet_sub.pt() > corrJets[i].pt()) nWrong++;

      //      cout << "nWrong" << endl;
      
      if(sdJet_sub.pt() >= pTmin_jet && sdJet_sub.pt() <= 30 && fabs(sdJet_sub.eta()) < 1){
	//	cout << "filling hists" << endl;
	jet_zg_sub->Fill(zg_sub);
	jet_rg_sub->Fill(rg_sub);
	jet_time_sub->Fill(log_t_form_sub);
	
	jet_zg_correlation->Fill(zg_sub,zg);

	jet_pT_correlation->Fill(sdJet_sub.pt(),pT);
	
	nJets_sub_gt_20++;
	//	cout << "done with hists" << endl;
      }
      //      cout << "out of hists" << endl;
      jet_pT_correlation_noCut->Fill(sdJet_sub.pt(),pT);

      //      cout << "making data" << endl;

      double pTJet_pT = -999;
      double pTJet_phi = -999;
      double pTJet_eta = -999;
      if(pTJets.size() != 0){
	pTJet_pT = pTJets[minK].pt();
	pTJet_phi = pTJets[minK].phi();
	pTJet_eta = pTJets[minK].eta();
      }


      double Jet_pT = -999;
      double Jet_phi = -999;
      double Jet_eta = -999;
      if(sortedJets.size() != 0){
	Jet_pT = sortedJets[minJ].pt();
	Jet_phi = sortedJets[minJ].phi();
	Jet_eta = sortedJets[minJ].eta();
      }

      double data2[] = {Jet_pT,Jet_phi,Jet_eta,pT,phi,eta,zg,rg,corrJets[i].pt(),corrJets[i].phi(),corrJets[i].eta(),sdJet_sub.pt(),sdJet_sub.phi(),sdJet_sub.eta(),zg_sub,rg_sub,pTJet_pT,pTJet_phi,pTJet_eta,pT_pT,phi_pT,eta_pT,zg_pT,rg_pT};

      //      cout << "setting data" << endl;
      
      for(int k=0; k<24 ;k++){
	data[k] = data2[k];
	cout << "data[" << k << "]: " << data[k] << endl;

      }
      //      cout << "filling tree" << endl;
      tree->Fill();
      //      cout << "finished jet " << i << "/" << sortedJets.size() << endl;
    }//end loop over jets

    iEvent++;
    cout << "Number of jets after subtraction with pT>20 GeV/c: " << nJets_sub_gt_20 << endl;
    //    cout << "Number of wrong groomed jets (groomed pT > ungroomed pT): " << nWrong << endl;
  }//end while over events

  outfile->cd();

  jet_pT->Scale(1.0,"width");
  jet_pT->Scale(1.0/jet_pT->GetEntries());

  jet_zg->Scale(1.0,"width");
  jet_zg->Scale(1.0/jet_zg->GetEntries());

  jet_rg->Scale(1.0,"width");
  jet_rg->Scale(1.0/jet_rg->GetEntries());

  jet_time->Scale(1.0,"width");
  jet_time->Scale(1.0/jet_time->GetEntries());


  jet_pT_pT->Scale(1.0,"width");
  jet_pT_pT->Scale(1.0/jet_pT_pT->GetEntries());

  jet_zg_pT->Scale(1.0,"width");
  jet_zg_pT->Scale(1.0/jet_zg_pT->GetEntries());

  jet_rg_pT->Scale(1.0,"width");
  jet_rg_pT->Scale(1.0/jet_rg_pT->GetEntries());

  jet_time_pT->Scale(1.0,"width");
  jet_time_pT->Scale(1.0/jet_time_pT->GetEntries());


  jet_pT_sub->Scale(1.0,"width");
  jet_pT_sub->Scale(1.0/jet_pT_sub->GetEntries());

  jet_zg_sub->Scale(1.0,"width");
  jet_zg_sub->Scale(1.0/jet_zg_sub->GetEntries());

  jet_rg_sub->Scale(1.0,"width");
  jet_rg_sub->Scale(1.0/jet_rg_sub->GetEntries());

  jet_time_sub->Scale(1.0,"width");
  jet_time_sub->Scale(1.0/jet_time_sub->GetEntries());

  jet_pT_correlation->Scale(1.0,"width");
  jet_pT_correlation->Scale(1.0/jet_pT_correlation->GetEntries());

  jet_pT_correlation_noCut->Scale(1.0,"width");
  jet_pT_correlation_noCut->Scale(1.0/jet_pT_correlation_noCut->GetEntries());

  jet_pT_correlation_hadCuts->Scale(1.0,"width");
  jet_pT_correlation_hadCuts->Scale(1.0/jet_pT_correlation_hadCuts->GetEntries());

  jet_zg_correlation->Scale(1.0,"width");
  jet_zg_correlation->Scale(1.0/jet_zg_correlation->GetEntries());

  jet_zg_correlation_hadCuts->Scale(1.0,"width");
  jet_zg_correlation_hadCuts->Scale(1.0/jet_zg_correlation_hadCuts->GetEntries());



  tree->Write();
  
  jet_pT->Write();
  jet_zg->Write();
  jet_rg->Write();
  jet_time->Write();

  jet_pT_pT->Write();
  jet_zg_pT->Write();
  jet_rg_pT->Write();
  jet_time_pT->Write();

  jet_pT_sub->Write();
  jet_zg_sub->Write();
  jet_rg_sub->Write();
  jet_time_sub->Write();
  
  jet_pT_correlation->Write();
  jet_pT_correlation_noCut->Write();
  jet_pT_correlation_hadCuts->Write();
  
  jet_zg_correlation->Write();
  jet_zg_correlation_hadCuts->Write();

  outfile->Close();
  
  return 0;
}//end main
