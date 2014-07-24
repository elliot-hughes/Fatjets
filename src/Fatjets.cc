/* C++
Path: Analyzers/Fatjets/src/Fatjets.cc
Package: Fatjets
Class: Fatjets

Author: Elliot Hughes
Last Updated: 140718
Purpose: To make and ntuplize fatjets.

I'm going to find 5 subjets for the leading and second-leading fatjets in a number of ways:
* Recluster with KT and R, then pull back 5 times.
* Recluster with CA and R, then pull back 5 times.
* Filter.
* Trim or prune.

I will store everything in an ntuple and use a python script to make plots later.  For each event, I will store pT, m, nsubjettiness variables, and pT and m of each of the 5 subjets under each method above, all for leading and second-leading fatjets.

*/

// INCLUDES
// system include files
#include <iostream>

// user include files
	// basic includes
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/ESHandle.h"

	// class includes
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/JetReco/interface/BasicJet.h"
#include "DataFormats/JetReco/interface/BasicJetCollection.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"		// JetDefinition, ClusterSequence, etc.
#include "fastjet-contrib/Nsubjettiness/Nsubjettiness.hh"
#include "fastjet/tools/Pruner.hh"


#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

	// ROOT includes
#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
//#include "TFileDirectory.h"
#include "TTree.h"
#include "TBranch.h"
#include "TNtuple.h"
// \INCLUDES
// NAMESPACES
using namespace std;
using namespace reco;
using namespace edm;
using namespace fastjet;
using namespace fastjet::contrib;
// \NAMESPACES

// CLASS DEFINITION
class Fatjets : public edm::EDAnalyzer {
	public:
		explicit Fatjets(const edm::ParameterSet&);	// Set the class argument to be (a reference to) a parameter set (?)
		~Fatjets();	// Create the destructor.
		static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

	private:
		virtual void beginJob();
		virtual void analyze(const edm::Event&, const edm::EventSetup&);
		virtual void endJob();
		virtual void beginRun(const edm::Run&, const edm::EventSetup&);
		virtual void endRun(const edm::Run&, const edm::EventSetup&);
		virtual void beginLuminosityBlock(const edm::LuminosityBlock&, const edm::EventSetup&);
		virtual void endLuminosityBlock(const edm::LuminosityBlock&, const edm::EventSetup&);

	// member data
	double R_;
	bool make_gen_;
	bool make_pf_;
	bool pruning;
	vector<string> jet_strings;
	// basic fatjet variables
	int njets;
	vector<double> vec_nsub1;
	vector<double> vec_nsub2;
	vector<double> vec_nsub3;
	vector<double> phi;
	vector<double> y;
	vector<double> pt;
	vector<double> e;
	vector<double> m;
	int njets_gen;
	vector<double> vec_nsub1_gen;
	vector<double> vec_nsub2_gen;
	vector<double> vec_nsub3_gen;
	vector<double> phi_gen;
	vector<double> y_gen;
	vector<double> pt_gen;
	vector<double> e_gen;
	vector<double> m_gen;
	// variables that involve subjets
	double dpt0;
	vector<double> pt0_sub;
	vector<double> e0_sub;
	vector<double> m0_sub;
	double dpt1;
	vector<double> pt1_sub;
	vector<double> e1_sub;
	vector<double> m1_sub;
	
	// CONSTRUCT ROOT OBJECT CONTAINERS
	map<string, TNtuple*> ntuple;
	map<string, TTree*> trees;
	map<string, TBranch*> branches;
	map<string, TH1F*> h1;
	map<string, TH2F*> h2;

	// Algorithm Variables
};

//
// constants, enums and typedefs
//
	// DEFINE CUTS
//	float cut_er = 5;

//
// static data member definitions
//


// constructors and destructors
//
Fatjets::Fatjets(const edm::ParameterSet& iConfig) :
	R_(iConfig.getParameter<double>("R")),
	make_gen_(iConfig.getParameter<bool>("make_gen")),
	make_pf_(iConfig.getParameter<bool>("make_pf"))
{
//do what ever initialization is needed
	jet_strings.push_back("ca");
	jet_strings.push_back("kt");
	jet_strings.push_back("ktr");
	// Set up output root file
	edm::Service<TFileService> fs;
	trees["fat"] = fs->make<TTree>("ttree_fat", "");
	trees["fat_gen"] = fs->make<TTree>("ttree_fat_gen", "");
	trees["ca"] = fs->make<TTree>("ttree_ca", "");
	trees["kt"] = fs->make<TTree>("ttree_kt", "");
	trees["ktr"] = fs->make<TTree>("ttree_ktr", "");
	h2["eta_phi"] =  fs->make<TH2F>("eta_phi", "Geometry Of Event", 40, -5, 5, 30, -3.2, 3.2);
	h2["y_phi"] =  fs->make<TH2F>("y_phi", "Geometry Of Fatjets", 40, -5, 5, 30, 0, 6.4);
	// Debug
	cout << "R = " << R_ << ", make_gen is " << make_gen_ << ", make_pf is " << make_pf_ << endl;
	//Nsubjettiness
}

// DEFINE THE DESTRUCTOR
Fatjets::~Fatjets()
{
}


// CLASS METHODS ("method" = "member function")

// ------------ called for each event  ------------
void Fatjets::analyze(
	const edm::Event& iEvent,
	const edm::EventSetup& iSetup
){
	pruning = true;
	map<string, JetDefinition> jet_defs;
	jet_defs[jet_strings[0]] = JetDefinition(cambridge_algorithm, R_);
	jet_defs[jet_strings[1]] = JetDefinition(kt_algorithm, R_);
	jet_defs[jet_strings[2]] = JetDefinition(kt_algorithm, 0.7);
	
	if (make_gen_) {		// This whole section is completely experimental.  Right now, it clusters every single gen particle (like, even the starting protons...).
		Handle<GenParticleCollection> objects_gen;
		iEvent.getByLabel("genParticles", objects_gen);
	
		cout << "The number of gen objects is " << objects_gen->size() << endl;
		
		vector<PseudoJet> jobjects_gen;
		int n = -1;
		for (GenParticleCollection::const_iterator gen = objects_gen->begin(); gen != objects_gen->end(); ++ gen) {
			if (gen->status() == 1) {
				n ++;
				if (n < 10){
					cout << n << ": pt = " << gen->pt() << ", y = " << gen->rapidity() << ", pdgid = " << gen->pdgId() << endl;
				}
				jobjects_gen.push_back( PseudoJet(gen->px(), gen->py(), gen->pz(), gen->energy()) );
			}
		}
		for (unsigned j = 0; j < jobjects_gen.size(); j++) {
			h2["y_phi"]->Fill(jobjects_gen[j].rapidity(), jobjects_gen[j].phi());
		}
		cout << "The number of gen objects with Status 1 is " << n << endl;
		ClusterSequence cs_gen(jobjects_gen, jet_defs[jet_strings[0]]);
		vector<PseudoJet> fatjets_gen = sorted_by_pt(cs_gen.inclusive_jets());
		
		
	    double beta = 1.0; // power for angular dependence, e.g. beta = 1 --> linear k-means, beta = 2 --> quadratic/classic k-means
	    double R0 = 1.2; // Characteristic jet radius for normalization	      
	    double Rcut = 1.2; // maximum R particles can be from axis to be included in jet
	    Nsubjettiness nsub1(1, Njettiness::onepass_kt_axes, beta, R0, Rcut);
	    Nsubjettiness nsub2(2, Njettiness::onepass_kt_axes, beta, R0, Rcut);
	    Nsubjettiness nsub3(3, Njettiness::onepass_kt_axes, beta, R0, Rcut);
	    Nsubjettiness nsub4(4, Njettiness::onepass_kt_axes, beta, R0, Rcut);
	    Nsubjettiness nsub5(5, Njettiness::onepass_kt_axes, beta, R0, Rcut);
	    Nsubjettiness nsub6(6, Njettiness::onepass_kt_axes, beta, R0, Rcut);
		
		vec_nsub1_gen.clear();
		vec_nsub2_gen.clear();
		vec_nsub3_gen.clear();
		phi_gen.clear();
		y_gen.clear();
		pt_gen.clear();
		e_gen.clear();
		m_gen.clear();
		njets_gen = fatjets_gen.size();
		int njets_gen_cut = -1;
		for (unsigned i = 0; i < fatjets_gen.size(); i++) {
			if (abs(fatjets_gen[i].rapidity()) < 6) {
				njets_gen_cut ++;
			}
			phi_gen.push_back( fatjets_gen[i].phi_std() );		// http://fastjet.fr/repo/doxygen-3.0.2/classfastjet_1_1PseudoJet.html
			y_gen.push_back( fatjets_gen[i].rapidity() );
			pt_gen.push_back( fatjets_gen[i].pt() );
			e_gen.push_back( fatjets_gen[i].e() );
			m_gen.push_back( fatjets_gen[i].m() );
			vec_nsub1_gen.push_back( nsub1(fatjets_gen[i]) );
			vec_nsub2_gen.push_back( nsub2(fatjets_gen[i]) );
			vec_nsub3_gen.push_back( nsub3(fatjets_gen[i]) );
		}
		cout << "-----: " << njets_gen_cut << endl;
		cout << "-----: " << njets_gen << endl;
		branches["njets"] = trees["fat_gen"]->Branch("njets", &njets_gen, "njets/I");
		branches["phi"] = trees["fat_gen"]->Branch("phi", &phi_gen, 64000, 0);
		branches["y"] = trees["fat_gen"]->Branch("y", &y_gen, 64000, 0);
		branches["pt"] = trees["fat_gen"]->Branch("pt", &pt_gen, 64000, 0);
		branches["e"] = trees["fat_gen"]->Branch("e", &e_gen, 64000, 0);
		branches["m"] = trees["fat_gen"]->Branch("m", &m_gen, 64000, 0);
		branches["nsub1"] = trees["fat_gen"]->Branch("nsub1", &vec_nsub1_gen, 64000, 0);
		branches["nsub2"] = trees["fat_gen"]->Branch("nsub2", &vec_nsub2_gen, 64000, 0);
		branches["nsub3"] = trees["fat_gen"]->Branch("nsub3", &vec_nsub3_gen, 64000, 0);
		trees["fat_gen"]->Fill();
		for (unsigned i = 0; i < 3; i++) {
			double phi1 = fatjets_gen[i].phi();
			double y1 = fatjets_gen[i].rapidity();
			for (int k = 0; k < 50; k++) {
				h2["y_phi"]->Fill(y1, phi1);
			}
			for (GenParticleCollection::const_iterator j = objects_gen->begin(); j != objects_gen->end(); ++ j) {
				double phi2 = j->phi();
				double y2 = j->rapidity();
				if (j->mass() > 200 && j->status() == 22) {
					double dphi = abs(phi1-phi2);
					if (dphi > M_PI) dphi -= 2 * M_PI;
					double dR2 = (y1-y2)*(y1-y2) + dphi*dphi;
					cout << i << ": pt_jet = " << fatjets_gen[i].pt() << ", m_jet = " << fatjets_gen[i].m() << " phi,y " << phi1 << "," << y1 << ", m_part = " << j->mass() << ", status = " << j->status() << ", dR2 = " << dR2 << endl;
				}
			}
		}
	}
	if (make_pf_) {
		Handle<PFCandidateCollection> objects_pf;
		iEvent.getByLabel("particleFlow", objects_pf);
	
		cout << "The number of PF objects is " << objects_pf->size() << endl;
		
		vector<PseudoJet> jobjects_pf;
		int n = -1;
		for (PFCandidateCollection::const_iterator pf = objects_pf->begin(); pf != objects_pf->end(); ++ pf) {
			n ++;
			if (n < 5){
				cout << n << ": px = " << pf->px() << ", py = " << pf->py() << ", pz = " << pf->pz() << ", E = " << pf->energy() << ": pdgId = " << pf->pdgId() << endl;
			}
			h2["eta_phi"]->Fill(pf->eta(), pf->phi());
			jobjects_pf.push_back( PseudoJet(pf->px(), pf->py(), pf->pz(), pf->energy()) );
		}
		JetDefinition thing(cambridge_algorithm, R_);
		ClusterSequence cs(jobjects_pf, thing);
		cout << "The number of objects to be used in clustering is " << cs.n_particles() << endl;
		vector<PseudoJet> fatjets_pf = sorted_by_pt(cs.inclusive_jets());

	    double beta = 1.0; // power for angular dependence, e.g. beta = 1 --> linear k-means, beta = 2 --> quadratic/classic k-means
	    double R0 = 1.2; // Characteristic jet radius for normalization	      
	    double Rcut = 1.2; // maximum R particles can be from axis to be included in jet
	    Nsubjettiness nsub1(1, Njettiness::onepass_kt_axes, beta, R0, Rcut);
	    Nsubjettiness nsub2(2, Njettiness::onepass_kt_axes, beta, R0, Rcut);
	    Nsubjettiness nsub3(3, Njettiness::onepass_kt_axes, beta, R0, Rcut);
	    Nsubjettiness nsub4(4, Njettiness::onepass_kt_axes, beta, R0, Rcut);
	    Nsubjettiness nsub5(5, Njettiness::onepass_kt_axes, beta, R0, Rcut);
	    Nsubjettiness nsub6(6, Njettiness::onepass_kt_axes, beta, R0, Rcut);
		
		vec_nsub1.clear();
		vec_nsub2.clear();
		vec_nsub3.clear();
		phi.clear();
		y.clear();
		pt.clear();
		e.clear();
		m.clear();
		njets = fatjets_pf.size();
		if (pruning != true) {
			// Loop over fatjets:
			for (unsigned i = 0; i < fatjets_pf.size(); i++) {
				phi.push_back( fatjets_pf[i].phi_std() );		// http://fastjet.fr/repo/doxygen-3.0.2/classfastjet_1_1PseudoJet.html
				y.push_back( fatjets_pf[i].rapidity() );
				pt.push_back( fatjets_pf[i].pt() );
				e.push_back( fatjets_pf[i].e() );
				m.push_back( fatjets_pf[i].m() );
				vec_nsub1.push_back( nsub1(fatjets_pf[i]) );
				vec_nsub2.push_back( nsub2(fatjets_pf[i]) );
				vec_nsub3.push_back( nsub3(fatjets_pf[i]) );
			}
			branches["njets"] = trees["fat"]->Branch("njets", &njets, "njets/I");
			branches["phi"] = trees["fat"]->Branch("phi", &phi, 64000, 0);
			branches["y"] = trees["fat"]->Branch("y", &y, 64000, 0);
			branches["pt"] = trees["fat"]->Branch("pt", &pt, 64000, 0);
			branches["e"] = trees["fat"]->Branch("e", &e, 64000, 0);
			branches["m"] = trees["fat"]->Branch("m", &m, 64000, 0);
			branches["nsub1"] = trees["fat"]->Branch("nsub1", &vec_nsub1, 64000, 0);
			branches["nsub2"] = trees["fat"]->Branch("nsub2", &vec_nsub2, 64000, 0);
			branches["nsub3"] = trees["fat"]->Branch("nsub3", &vec_nsub3, 64000, 0);
			trees["fat"]->Fill();
		}
		
		int n_subjets = 4;		// The number of subjets you want to find.
		// Loop over reclustering methods:
		for (unsigned j = 0; j < jet_strings.size(); j++) {
			string jet_string = jet_strings[j];
			cout << endl;
			cout << ">> Performing " << jet_string << " reclustering:" << endl;
			// Loop over the fatjets:
			for (unsigned i = 0; i < fatjets_pf.size(); i++) {
				vector<PseudoJet> jobjects;
				// Recluster the fatjet, so its cluster sequence is specific to the jet in question.
				double e_original = fatjets_pf[i].e();	//debug variable
				vector<PseudoJet> subjets;
				if (pruning) {
					Pruner pruner(jet_defs[jet_strings[1]], 1, 0.5);
					PseudoJet fatjet = pruner(fatjets_pf[i]);
					jobjects = fatjet.constituents();
					if (j == 0) {
						phi.push_back( fatjet.phi_std() );		// http://fastjet.fr/repo/doxygen-3.0.2/classfastjet_1_1PseudoJet.html
						y.push_back( fatjet.rapidity() );
						pt.push_back( fatjet.pt() );
						e.push_back( fatjet.e() );
						m.push_back( fatjet.m() );
						vec_nsub1.push_back( nsub1(fatjet) );
						vec_nsub2.push_back( nsub2(fatjet) );
						vec_nsub3.push_back( nsub3(fatjet) );
					}
				}
				else {
					jobjects = fatjets_pf[i].constituents();		// Get all of the original objects that formed into this jet.
				}
				ClusterSequence cs_jet(jobjects, jet_defs[jet_strings[j]]);		// Cluster these objects again (back into the same jet).
				vector<PseudoJet> reclustered_fatjets = sorted_by_pt(cs_jet.inclusive_jets());		// This is saved for debugging; read down a few lines before asking, "what?".
				if (j == 0 || j == 1) {
					subjets = sorted_by_pt(cs_jet.exclusive_jets_up_to(n_subjets));		// This basically does what my entire algorithm does...
				}
				else if (j == 2){
					subjets = reclustered_fatjets;
				}
				if (reclustered_fatjets.size() != 1 && j != 2) {
					cout << "WARNING (clustering): Reclustering yielded more than one fatjet." << endl;
					cout << "WARNING cont: energy of first jet = " << reclustered_fatjets[0].e() << " GeV, e2 = " << reclustered_fatjets[1].e() << " GeV" << endl;
				}
				if (i == 0) {
					pt0_sub.clear();
					e0_sub.clear();
					m0_sub.clear();
					double e_total = 0;
					for (unsigned k = 0; k < subjets.size(); k++) {
						pt0_sub.push_back( subjets[k].pt() );
						e0_sub.push_back( subjets[k].e() );
						m0_sub.push_back( subjets[k].m() );
						cout << "    Subjet " << k << ": pt = " << subjets[k].pt() << " GeV,\te = " << subjets[k].e() << " GeV" << endl;
						e_total += subjets[k].e();
					}
	//				cout << "    max_index = " << max_index << "  target = " << target << endl;'
	//				double e_range;
					dpt0 = subjets[0].pt() - subjets[subjets.size()-1].pt();
					cout << "    pt_range = " << dpt0 << endl;
				}
				if (i == 1) {
					pt1_sub.clear();
					e1_sub.clear();
					m1_sub.clear();
					double e_total = 0;
					for (unsigned k = 0; k < subjets.size(); k++) {
						pt1_sub.push_back( subjets[k].pt() );
						e1_sub.push_back( subjets[k].e() );
						m1_sub.push_back( subjets[k].m() );
						cout << "    Subjet " << k << ": pt = " << subjets[k].pt() << " GeV,\te = " << subjets[k].e() << " GeV" << endl;
						e_total += subjets[k].e();
					}
	//				cout << "    max_index = " << max_index << "  target = " << target << endl;'
	//				double e_range;
					dpt1 = subjets[0].pt() - subjets[subjets.size()-1].pt();
					cout << "    pt_range = " << dpt1 << endl;
				}
				
//				// BEGIN SUBJET ALGORITHM (Find the n subjets by unfolding fatjet back n-1 times.)
//				int max_index = fatjet.cluster_hist_index();
//				vector<PseudoJet> pjets = cs_jet.jets();		// The history as seen through jets. Last one is the last one formed, 5th one is 5th jet formed.
//			
//				int target = max_index - n + 1;		// Algorith variable.
//				vector<int> subjet_indices;		// Vector of indices of "partners" <- my definition.
//				subjet_indices.push_back(max_index);
//				sort(subjet_indices.begin(), subjet_indices.end());		// Sorts vector: lowest -> highest.
//				int max_partner = subjet_indices.back();
//			
//				while (max_partner > target) {
//					subjet_indices.erase( remove(subjet_indices.begin(), subjet_indices.end(), max_partner), subjet_indices.end() ); // Remove max_partner (because it will be replaced by the indices of its parents in the following code).
//					PseudoJet s1;
//					PseudoJet s2;
//					pjets[max_partner].has_parents(s1, s2);
//					subjet_indices.push_back(s1.cluster_hist_index());
//					subjet_indices.push_back(s2.cluster_hist_index());
//					sort(subjet_indices.begin(), subjet_indices.end());		// Sorts vector: lowest -> highest.
//					max_partner = subjet_indices.back();
//				}
//				vector<PseudoJet> subjets;
//				for (unsigned j = 0; j < subjet_indices.size(); j++) {
//					subjets.push_back(pjets[subjet_indices[j]]);
//				} 
//				// END SUBJET ALGORITHM (The subjet indices are stored in subjet_indices and the subjets are stored in subjets.
			}
			if (pruning) {
				branches["njets"] = trees["fat"]->Branch("njets", &njets, "njets/I");
				branches["phi"] = trees["fat"]->Branch("phi", &phi, 64000, 0);
				branches["y"] = trees["fat"]->Branch("y", &y, 64000, 0);
				branches["pt"] = trees["fat"]->Branch("pt", &pt, 64000, 0);
				branches["e"] = trees["fat"]->Branch("e", &e, 64000, 0);
				branches["m"] = trees["fat"]->Branch("m", &m, 64000, 0);
				branches["nsub1"] = trees["fat"]->Branch("nsub1", &vec_nsub1, 64000, 0);
				branches["nsub2"] = trees["fat"]->Branch("nsub2", &vec_nsub2, 64000, 0);
				branches["nsub3"] = trees["fat"]->Branch("nsub3", &vec_nsub3, 64000, 0);
				trees["fat"]->Fill();
			}
			branches["dpt0"] = trees[jet_string]->Branch("dpt0", &dpt0, "dpt0/D");
			branches["pt0_sub"] = trees[jet_string]->Branch("pt0_sub", &pt0_sub, 64000, 0);
			branches["e0_sub"] = trees[jet_string]->Branch("e0_sub", &e0_sub, 64000, 0);
			branches["m0_sub"] = trees[jet_string]->Branch("m0_sub", &m0_sub, 64000, 0);
			branches["dpt1"] = trees[jet_string]->Branch("dpt1", &dpt1, "dpt1/D");
			branches["pt1_sub"] = trees[jet_string]->Branch("pt1_sub", &pt1_sub, 64000, 0);
			branches["e1_sub"] = trees[jet_string]->Branch("e1_sub", &e1_sub, 64000, 0);
			branches["m1_sub"] = trees[jet_string]->Branch("m1_sub", &m1_sub, 64000, 0);
			trees[jet_string]->Fill();
		}
	}
}


// ------------ called once each job just before starting event loop  ------------
void Fatjets::beginJob()
{
}

// ------------  called once each job just after ending the event loop  ------------
void Fatjets::endJob()
{
}

// ------------ method called when starting to processes a run  ------------
void 
Fatjets::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
Fatjets::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
Fatjets::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
Fatjets::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
Fatjets::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(Fatjets);
