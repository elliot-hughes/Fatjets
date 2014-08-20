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
#include <typeinfo>

// user include files
//// basic includes
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/ESHandle.h"

//// important class includes
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/JetReco/interface/BasicJet.h"
#include "DataFormats/JetReco/interface/BasicJetCollection.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"		// JetDefinition, ClusterSequence, etc.
#include "fastjet-contrib/Nsubjettiness/Nsubjettiness.hh"
#include "fastjet/tools/Pruner.hh"
#include "fastjet/Selector.hh"
#include "fastjet/tools/Filter.hh"


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

// OTHER
struct sort_by_m {
	bool operator() (PseudoJet jet1, PseudoJet jet2) {
		return (jet1.m() > jet2.m());
	}
};





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
	// These variables are filled by setting the python configuration file:
	bool v_;		// Verbose control
	double R_;		// Fatjet distance parameter
	bool make_gen_, make_pf_;		// Controls to make gen fatjets or pf fatjets
	bool do_true_;		// Control to look only at squark products (DOESN'T WORK)
	bool do_prune_, do_trim_;		// Controls for fatjet transformers
	double trim_R_, trim_ptf_;		// Trimmer transformer parameters
	vector<string> jet_strings;
	// basic fatjet variables
	vector<string> fat_vars;
	map<string, vector<double>> fat_gen;		// This is a container for the branches I eventually save as gen output.
	map<string, vector<double>> fat_pf;
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
	int n_event, counter, n_error_g, n_error_sq, n_error_m, n_error_sort;
	map<string, JetDefinition> jet_defs;
	vector<int> outliers;
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
	v_(iConfig.getParameter<bool>("v")),
	R_(iConfig.getParameter<double>("R")),
	make_gen_(iConfig.getParameter<bool>("make_gen")),
	make_pf_(iConfig.getParameter<bool>("make_pf")),
	do_true_(iConfig.getParameter<bool>("do_true")),
	do_prune_(iConfig.getParameter<bool>("do_prune")),
	do_trim_(iConfig.getParameter<bool>("do_trim")),
	trim_R_(iConfig.getParameter<double>("trim_R")),
	trim_ptf_(iConfig.getParameter<double>("trim_ptf"))
{
//do what ever initialization is needed
	// Simple algorithm variables
	n_event = 0;
	n_error_g = 0;
	n_error_sq = 0;
	n_error_m = 0;
	n_error_sort = 0;
	counter = 0;
	// Jet algorithm variables
	jet_strings.push_back("ca");
	jet_strings.push_back("kt");
	jet_strings.push_back("ktr");
	jet_defs[jet_strings[0]] = JetDefinition(cambridge_algorithm, R_);
	jet_defs[jet_strings[1]] = JetDefinition(kt_algorithm, R_);
	jet_defs[jet_strings[2]] = JetDefinition(kt_algorithm, 0.7);
	// Fatjet variables
	// (When you add a new one here, make sure you fill it!)
	fat_vars.push_back("phi");
	fat_vars.push_back("y");
	fat_vars.push_back("pt");
	fat_vars.push_back("e");
	fat_vars.push_back("m");
	fat_vars.push_back("tau1");
	fat_vars.push_back("tau2");
	fat_vars.push_back("tau3");
	fat_vars.push_back("dr_fj");
	fat_vars.push_back("dr_sq");
	// Set up output root file
	edm::Service<TFileService> fs;
	trees["fat_pf"] = fs->make<TTree>("fat_pf", "");		// This will contain branches for each pf fatjet variable.
	trees["fat_gen"] = fs->make<TTree>("fat_gen", "");		// This will contain branches for each gen fatjet variable.
	trees["ca"] = fs->make<TTree>("sub_ca", "");
	trees["kt"] = fs->make<TTree>("sub_kt", "");
	trees["ktr"] = fs->make<TTree>("sub_ktr", "");
	h1["dR"] =  fs->make<TH1F>("dR", "#Delta R", 40, 0, 2);
	h2["y_phi"] =  fs->make<TH2F>("y_phi", "Geometry Of Fatjets", 40, -5, 5, 30, 0, 6.4);
	// Debug
	cout << endl;
	cout << "Starting the Fatjet analyzer..." << endl;
	cout << "R = " << R_ << ", make_gen is " << make_gen_ << ", make_pf is " << make_pf_ << endl;
}

// DEFINE THE DESTRUCTOR
Fatjets::~Fatjets()
{
}

// ------------ called once each job just before starting event loop  ------------
void Fatjets::beginJob()
{
}

// CLASS METHODS ("method" = "member function")

// ------------ called for each event  ------------
void Fatjets::analyze(
	const edm::Event& iEvent,
	const edm::EventSetup& iSetup
){
	n_event ++;
	// Declare fatjet transformers:
	Pruner pruner(jet_defs[jet_strings[1]], 0.3, 1.0);		// jobject is discarded if both are true: zcut: pt1 or pt2 < zcut*pt(1+2), Rcut_factor: dR12 > Rcut_factor*2m/pt (m, pt are the fatjet's).
	Filter trimmer(trim_R_, SelectorPtFractionMin(trim_ptf_));		// Recluster FJ using CA with R = trim_R_, then discard any subjets that have pt less than trim_ptf_*pt_fatjet.
	if (make_gen_) {
		// Declare some gen-specific variables:
		vector<Candidate*> goddesses;		// This is a confused idea that I need to rethink.
		vector<GenParticle*> squarks;		// Containers for the initial squarks
		double dR_sq;		// Distance between the two initial squarks
		
		// Get gen particles from event:
		Handle<GenParticleCollection> objects_gen;
		iEvent.getByLabel("genParticles", objects_gen);
		
		if (v_) cout << ">> Starting to make gen fatjets:" << endl;
		if (v_) cout << ">> The number of gen objects in the event is " << objects_gen->size() << "." << endl;
		
		// Cluster gen particles:
		vector<PseudoJet> jobjects_gen;		// Container for gen particles that will be clustered into a fatjet
		vector<PseudoJet> jobjects_gen_true;		// Container for specialized gen particles, if a more exotic selection is needed.
		int n_object = -1;		// Counter for gen particles
		int n_object_1 = -1;		// Counter for gen particles with Status 1 (stable)
		//// Loop over gen particles:
		for (GenParticleCollection::const_iterator object = objects_gen->begin(); object != objects_gen->end(); ++ object) {
			n_object ++;
//			cout << "n_gen_particle = " << n_object << endl;
			int status = object->status();		// Get the gen particle status.
			int pdgid = object->pdgId();		// Get the gen particle PDG ID.
			int n_mothers = object->numberOfMothers();		// Get the number of mothers the gen particle has.
			double mass = object->mass(), pt = object->pt(), y = object->rapidity(), phi = object->phi();		// Get kinematic variables.
			bool is_goddess = false;
			Candidate* goddess;
			
			//// Get the squarks:
			if (status == 22 && abs(pdgid) == 1000005) {
				squarks.push_back(object->clone());		// Remember, clone makes a copy of the object and returns a pointer to it.
			}
			
			//// Play around with gen particles:
//			if (status == 22 || status == 23 || status == 4) {
//				if (n_mothers > 0) {
//					cout << pdgid << "   " << status << "   " << y << "   " << phi << "   " << n_mothers << "   " << object->mother()->pdgId() << endl;
//				}
//				else {
//					cout << pdgid << "   " << status << "   " << y << "   " << phi << "   " << n_mothers << endl;
//				}
//			}
			// Find object's goddess (supreme mother):
//			GenParticle object_temp = *object;
//			const Candidate* goddess_temp = (object->mother()) ? object->mother() : &object_temp;
//			map<int, Candidate> lineage;
			vector<int> lineage;
			lineage.push_back(status);
			if (object->mother() && object->numberOfMothers() < 2){		// THIS IS A PROBLEM (it only looks one step back...)
				if (object->numberOfMothers() > 1) {
					n_error_m ++;
				}
				const Candidate* goddess_temp = object->mother();
				while ( !( (goddess_temp->status() == 22 && abs(goddess_temp->pdgId()) < 1000021) || goddess_temp->status() == 4) ) {
					lineage.push_back(goddess_temp->status());
					goddess_temp = goddess_temp->mother();
				}
				lineage.push_back(goddess_temp->status());
				goddess = goddess_temp->clone();
			}
			else {
				is_goddess = true;
			}
			bool add = true;
			if ( !(is_goddess) ) {
//				if (n_object < 100) {
//					cout << goddess->pdgId() << "  " << goddess->status() << "  " << goddess->rapidity() << "  " << goddess->phi() << endl;
//				}
				for (auto g : goddesses) {
					if ( (goddess->pdgId() == g->pdgId() && goddess->rapidity() == g->rapidity() && goddess->phi() == g->phi()) ) {
//						cout << goddess->pdgId() << "  " << g->pdgId() << endl;
						add = false;
					}
				}
				if (add == true) {
					goddesses.push_back(goddess);
				}
			}
			
			//// Fill jobjects with the particles that stick around:
			if (status == 1) {
				n_object_1 ++;
				jobjects_gen.push_back( PseudoJet(object->px(), object->py(), object->pz(), object->energy()) );
//				if (goddess->pdgId() == 1000005 && n_mothers < 2) {		// I make a collection of objects that only came from the squarks. (DOESN'T WORK)
//					jobjects_gen_true.push_back( PseudoJet(object->px(), object->py(), object->pz(), object->energy()) );
//				}
//				double dR = sqrt( pow(goddess->rapidity()-object->rapidity(), 2) + pow(deltaPhi(*goddess, *object), 2) );
//				h1["dR"]->Fill(dR);
				
				
//				if (n_object_1 < 100) {
//					cout << n_object_1 << ": pdgid = " << pdgid << ", pt = " << pt << ", lineage = ";
//					for (auto i : lineage) {
//						cout << i << " ";
//					}
//					cout << endl;
//				}
			}
//			if (object->pdgId() > 1000021) {
//				cout << n_object << ": m = " << object->mass() << ", phi = " << object->phi() << ", y = " << object->rapidity() << ", pdgid = " << object->pdgId() << ", status = " << object->status() << ", goddess = " << goddess->pdgId() << "  " << goddess->status() << endl;
//			}
		}		// End of gen particle loop.
//		cout << goddesses.size() << endl;
//		if (goddesses.size() != 4) {
//			n_error_g ++;
//		}
		// Do we only see 2 squarks in the event (as we expected to)?
		if (squarks.size() != 2) {
			n_error_sq ++;
			cout << "Error: " << squarks.size() << " were found!" << endl;
		}
		else {
			// Calculate dR between the squarks:
			double y0 = squarks[0]->rapidity();
			double y1 = squarks[1]->rapidity();
			double phi0 = squarks[0]->phi();
			double phi1 = squarks[1]->phi();
			dR_sq = sqrt( pow( (y0 - y1), 2 ) + pow( deltaPhi(phi0, phi1), 2 ) );
//			cout << dR_sq << endl;
		}
		
		// More stuff I need to fix up.
		int first_squark = 0;
//		Candidate* first_squark;
//		if (n_event == 404) {
		bool found = false;
		int n_goddess = -1;
		for (auto g : goddesses) {
			n_goddess ++;
			if (g->pdgId() == 1000005 && found == false) {
				first_squark = n_goddess;
//				first_squark = goddesses[n_goddess];
				found = true;
//				cout << "HERE" << endl;
//				cout << first_squark << endl;
////				cout << first_squark->pdgId() << endl;
			}
//			cout << g->pdgId() << endl; // "  " << g->rapidity() << "  " << g->phi() << endl;
		}
//		if (found == true) {
//			cout << first_squark << endl;
//			cout << goddesses[first_squark]->pdgId() << endl; // << "  " << goddessfirst_squark->rapidity() << "  " << first_squark->phi() << endl;
//		}
//		else {
//			cout << "Could not find first squark." << endl;
//		}
		
		// Make the fatjets!:
		ClusterSequence cs_gen = (do_true_) ? ClusterSequence(jobjects_gen_true, jet_defs[jet_strings[0]]) : ClusterSequence(jobjects_gen, jet_defs[jet_strings[0]]);		// Perform the clustering with FastJet.
		vector<PseudoJet> fatjets_gen = cs_gen.inclusive_jets();
//		vector<PseudoJet> fatjets_unsorted = sorted_by_pt(cs_gen.inclusive_jets());
//		vector<PseudoJet> fatjets_gen = fatjets_unsorted;
		sort(fatjets_gen.begin(), fatjets_gen.end(), sort_by_m());		// Sort fatjets by invariant mass.
		
		if (v_) cout << ">> Number of gen fatjets: " << fatjets_gen.size() << endl;
		
//		vector<double> m(fatjets_unsorted.size());
//		for (size_t i = 0; i < fatjets_unsorted.size(); i++) {m[i] = fatjets_unsorted[i].m();}
//		vector<PseudoJet> fatjets_gen = sorted_by_pt(fatjets_unsorted);
//		vector<PseudoJet> fatjets_gen = objects_sorted_by_values(fatjets_unsorted, m);		// Make a vector of fatjets, sorted by m.
		
//		// get a vector of indices
//		vector<int> indices(m.size());
//		for (size_t i = 0; i < indices.size(); i++) {indices[i] = i;}
//		// sort the indices
//		sort_indices(indices, m);
//		// copy the objects 
//		vector<PseudoJet> fatjets_gen(fatjets_unsorted.size());
//		// place the objects in the correct order
//		for (size_t i = 0; i < indices.size(); i++) {
//			fatjets_gen[i] = fatjets_unsorted[indices[i]];
//		}

		// Check fatjet sorting (delete eventually):
		if (fatjets_gen.size() > 1){
//			cout << "HERE 275 " << fatjets_gen[0].m() << "  " << fatjets_gen[1].m() << "  " << fatjets_gen[2].m() << "  " << fatjets_gen[3].m() << "  " << fatjets_gen[4].m() << "   " << fatjets_gen.size() << endl;
			if (fatjets_gen[0].m() < fatjets_gen[1].m() || fatjets_gen[0].m() < fatjets_gen[2].m()) {
				n_error_sort ++;
			}
		}

		// ADD: FIX NSUBJETTINESS WITH DINKO'S ALGO:
		double beta = 1.0; // power for angular dependence, e.g. beta = 1 --> linear k-means, beta = 2 --> quadratic/classic k-means
		double R0 = 1.2; // Characteristic jet radius for normalization	      
		double Rcut = 1.2; // maximum R particles can be from axis to be included in jet
		Nsubjettiness nsub1(1, Njettiness::onepass_kt_axes, beta, R0, Rcut);
		Nsubjettiness nsub2(2, Njettiness::onepass_kt_axes, beta, R0, Rcut);
		Nsubjettiness nsub3(3, Njettiness::onepass_kt_axes, beta, R0, Rcut);
		
		// Clear fatjet variable containers:
		for (vector<string>::iterator var = fat_vars.begin(); var != fat_vars.end(); ++ var){
			fat_gen[*var].clear();
		}
		
		// Loop over fatjets:
		int n_fj = -1;
		for (vector<PseudoJet>::iterator fj = fatjets_gen.begin(); fj != fatjets_gen.end(); ++fj) {
			n_fj ++;
			
			// Declare variables specific to fatjet loop:
			vector<double> dR_fj_sq;		// Distances between fajet and squarks
			double dR_min = 100;
			int sq_close = -1;
		
			// Apply transformers:
			PseudoJet fatjet;
			if (do_prune_) {		// Apply pruning if need be.
				fatjet = pruner(*fj);
			}
			else if (do_trim_){
				fatjet = trimmer(*fj);
			}
			else {
				fatjet = *fj;
			}
			
			if (n_fj < 2){
				if (v_) cout << ">> FJ " << n_fj << ": m = " << fatjet.m() << " GeV" << endl;
			}
			
			// Find squark closest to fatjet:
			//// Loop over squarks:
			for (unsigned sq = 0; sq < squarks.size(); sq ++) {
				double phi_sq = squarks[sq]->phi();
				double dR_temp = sqrt( pow(squarks[sq]->rapidity()-fatjet.rapidity(), 2) + pow(deltaPhi(phi_sq, fatjet.phi()), 2) );
				dR_fj_sq.push_back(dR_temp);
				if (dR_temp < dR_min) {
					dR_min = dR_temp;
					sq_close = sq;
				}
			}
//			cout << "The closest squark is Squark " << sq_close << " with dR = " << dR_min << ",   ";
//			for (auto r : dR_fj_sq) {
//				cout << r << "  ";
//			}
//			cout << endl;
			
			// Old stuff:
//			if (fatjet.m() > 600) {
//				outliers.push_back(n_event);
//				vector<PseudoJet> pieces = fatjet.constituents();
//				sort(pieces.begin(), pieces.end(), sort_by_m());
//				cout << "Big fatjet: m = " << fatjet.m() << ", n_constituents = " << pieces.size() << endl;
//				int n_piece = -1;
//				for (vector<PseudoJet>::iterator piece = pieces.begin(); piece != pieces.end(); ++ piece) {
//					n_piece ++;
//					if (n_piece < 10) {
//						cout << "	" << n_piece << ": m = " << piece->m() << ": y = " << piece->rapidity() << ": phi = " << piece->phi() << endl;
//					}
//				}
//			}
			
			// Fill the vectors that I want to store as branches.
			fat_gen["phi"].push_back( fatjet.phi() );
			fat_gen["y"].push_back( fatjet.rapidity() );
			fat_gen["pt"].push_back( fatjet.pt() );
			fat_gen["e"].push_back( fatjet.e() );
			fat_gen["m"].push_back( fatjet.m() );
			fat_gen["tau1"].push_back( nsub1(fatjet) );
			fat_gen["tau2"].push_back( nsub2(fatjet) );
			fat_gen["tau3"].push_back( nsub3(fatjet) );
			fat_gen["dr_fj"].push_back( dR_min );
		}
		fat_gen["dr_sq"].push_back( dR_sq );
		
		// Make branches from the fatjet variable vectors:
		for (vector<string>::iterator var = fat_vars.begin(); var != fat_vars.end(); ++ var){
			branches[*var] = trees["fat_gen"]->Branch(var->c_str(), &(fat_gen[*var]), 64000, 0);
		}
		trees["fat_gen"]->Fill();		// Save branches to the tree.
		
//		// Gen particle exploration:
//		for (unsigned i = 0; i < 3; i++) {
//			double phi1 = fatjets_gen[i].phi();
//			double y1 = fatjets_gen[i].rapidity();
//			for (int k = 0; k < 50; k++) {
//				h2["y_phi"]->Fill(y1, phi1);
//			}
//			for (GenParticleCollection::const_iterator j = objects_gen->begin(); j != objects_gen->end(); ++ j) {
//				double phi2 = j->phi();
//				double y2 = j->rapidity();
//				if (j->mass() > 200 && j->status() == 22) {
//					double dphi = abs(phi1-phi2);
//					if (dphi > M_PI) dphi -= 2 * M_PI;
//					double dR2 = (y1-y2)*(y1-y2) + dphi*dphi;
//					cout << i << ": pt_jet = " << fatjets_gen[i].pt() << ", m_jet = " << fatjets_gen[i].m() << " phi,y " << phi1 << "," << y1 << ", m_part = " << j->mass() << ", status = " << j->status() << ", dR2 = " << dR2 << endl;
//				}
//			}
//		}
	}
	
	// REVISIONS STOP HERE...
	if (make_pf_) {
		
		// Get PF objects from the event:
		Handle<PFCandidateCollection> objects_pf;
		iEvent.getByLabel("particleFlow", objects_pf);
	
		if (v_) cout << ">> Starting to make PF fatjets:" << endl;
		if (v_) cout << ">> The number of PF objects is " << objects_pf->size() << "." << endl;
		
		// Cluster PF objects:
		vector<PseudoJet> jobjects_pf;		// Container for PF particles that will be clustered into a fatjet
		int n_object = -1;		// Counter for PF objects
		//// Loop over PF objects:
		for (PFCandidateCollection::const_iterator object = objects_pf->begin(); object != objects_pf->end(); ++ object) {
			n_object ++;
			
			// ADD: define some variables to store object values like pt
			
//			if (n_object < 5){
//				cout << "PF object Number " << n_object << ": px = " << object->px() << ", py = " << object->py() << ", pz = " << object->pz() << ", E = " << object->energy() << ": pdgId = " << object->pdgId() << endl;
//			}
//			h2["eta_phi"]->Fill(object->eta(), object->phi());
			
			jobjects_pf.push_back( PseudoJet(object->px(), object->py(), object->pz(), object->energy()) );
		}
		//// Cluster:
		ClusterSequence cs_pf(jobjects_pf, jet_defs[jet_strings[0]]);		// Perform the clustering with FastJet.
		vector<PseudoJet> fatjets_pf = cs_pf.inclusive_jets();
		sort(fatjets_pf.begin(), fatjets_pf.end(), sort_by_m());		// Sort fatjets by invariant mass.
		
		if (v_) cout << ">> Number of PF fatjets: " << fatjets_pf.size() << endl;
		
		
		// Check fatjet sorting (delete eventually):
		if (fatjets_pf.size() >= 1){
			if (fatjets_pf[0].m() < fatjets_pf[1].m()) {
				cout << ">> ERROR: sorting by m" << endl;
				n_error_sort ++;
			}
		}
		
		
		// ADD: FIX NSUBJETTINESS WITH DINKO'S ALGO:
		double beta = 1.0; // power for angular dependence, e.g. beta = 1 --> linear k-means, beta = 2 --> quadratic/classic k-means
		double R0 = 1.2; // Characteristic jet radius for normalization	      
		double Rcut = 1.2; // maximum R particles can be from axis to be included in jet
		Nsubjettiness nsub1(1, Njettiness::onepass_kt_axes, beta, R0, Rcut);
		Nsubjettiness nsub2(2, Njettiness::onepass_kt_axes, beta, R0, Rcut);
		Nsubjettiness nsub3(3, Njettiness::onepass_kt_axes, beta, R0, Rcut);
		
		// Clear fatjet variable containers:
		for (vector<string>::iterator var = fat_vars.begin(); var != fat_vars.end(); ++ var){
			fat_pf[*var].clear();
		}
		
		// Loop over fatjets:
		int n_fj = -1;
		for (vector<PseudoJet>::iterator fj = fatjets_pf.begin(); fj != fatjets_pf.end(); ++fj) {
			n_fj ++;
			// Declare variables specific to PF fatjet loop:
		
			// Apply transformers:
			PseudoJet fatjet;
			if (do_prune_) {		// Apply pruning if need be.
				fatjet = pruner(*fj);
			}
			else if (do_trim_){		// Apply trimming if need be.
				fatjet = trimmer(*fj);
			}
			else {
				fatjet = *fj;
			}
			
			if (n_fj < 2){
				if (v_) cout << ">> FJ " << n_fj << ": m = " << fatjet.m() << " GeV" << endl;
			}
			
			// Fill the vectors that I want to store as branches.
			fat_pf["phi"].push_back( fatjet.phi() );
			fat_pf["y"].push_back( fatjet.rapidity() );
			fat_pf["pt"].push_back( fatjet.pt() );
			fat_pf["e"].push_back( fatjet.e() );
			fat_pf["m"].push_back( fatjet.m() );
			fat_pf["tau1"].push_back( nsub1(fatjet) );
			fat_pf["tau2"].push_back( nsub2(fatjet) );
			fat_pf["tau3"].push_back( nsub3(fatjet) );
		}
		
		// Make branches from the fatjet variable vectors:
		for (vector<string>::iterator var = fat_vars.begin(); var != fat_vars.end(); ++ var){
			branches[*var] = trees["fat_pf"]->Branch(var->c_str(), &(fat_pf[*var]), 64000, 0);
		}
		trees["fat_pf"]->Fill();		// Save branches to the tree.

//		// OLD SUBJET THINGS
//		int n_subjets = 4;		// The number of subjets you want to find.
//		// Loop over reclustering methods:
//		for (unsigned j = 0; j < jet_strings.size(); j++) {
//			string jet_string = jet_strings[j];
//			cout << endl;
//			cout << ">> Performing " << jet_string << " reclustering:" << endl;
//			// Loop over the fatjets:
//			for (unsigned i = 0; i < fatjets_pf.size(); i++) {
//				vector<PseudoJet> jobjects;
//				// Recluster the fatjet, so its cluster sequence is specific to the jet in question.
////				double e_original = fatjets_pf[i].e();	//debug variable
//				vector<PseudoJet> subjets;
//				if (do_prune_) {
//					Pruner pruner(jet_defs[jet_strings[1]], 1, 0.5);
//					PseudoJet fatjet = pruner(fatjets_pf[i]);
//					jobjects = fatjet.constituents();
//					if (j == 0) {
//						phi.push_back( fatjet.phi_std() );		// http://fastjet.fr/repo/doxygen-3.0.2/classfastjet_1_1PseudoJet.html
//						y.push_back( fatjet.rapidity() );
//						pt.push_back( fatjet.pt() );
//						e.push_back( fatjet.e() );
//						m.push_back( fatjet.m() );
//						vec_nsub1.push_back( nsub1(fatjet) );
//						vec_nsub2.push_back( nsub2(fatjet) );
//						vec_nsub3.push_back( nsub3(fatjet) );
//					}
//				}
//				else {
//					jobjects = fatjets_pf[i].constituents();		// Get all of the original objects that formed into this jet.
//				}
//				ClusterSequence cs_jet(jobjects, jet_defs[jet_strings[j]]);		// Cluster these objects again (back into the same jet).
//				vector<PseudoJet> reclustered_fatjets = sorted_by_pt(cs_jet.inclusive_jets());		// This is saved for debugging; read down a few lines before asking, "what?".
//				if (j == 0 || j == 1) {
//					subjets = sorted_by_pt(cs_jet.exclusive_jets_up_to(n_subjets));		// This basically does what my entire algorithm does...
//				}
//				else if (j == 2){
//					subjets = reclustered_fatjets;
//				}
//				if (reclustered_fatjets.size() != 1 && j != 2) {
//					cout << "WARNING (clustering): Reclustering yielded more than one fatjet." << endl;
//					cout << "WARNING cont: energy of first jet = " << reclustered_fatjets[0].e() << " GeV, e2 = " << reclustered_fatjets[1].e() << " GeV" << endl;
//				}
//				if (i == 0) {
//					pt0_sub.clear();
//					e0_sub.clear();
//					m0_sub.clear();
//					double e_total = 0;
//					for (unsigned k = 0; k < subjets.size(); k++) {
//						pt0_sub.push_back( subjets[k].pt() );
//						e0_sub.push_back( subjets[k].e() );
//						m0_sub.push_back( subjets[k].m() );
//						cout << "    Subjet " << k << ": pt = " << subjets[k].pt() << " GeV,\te = " << subjets[k].e() << " GeV" << endl;
//						e_total += subjets[k].e();
//					}
//	//				cout << "    max_index = " << max_index << "  target = " << target << endl;'
//	//				double e_range;
//					dpt0 = subjets[0].pt() - subjets[subjets.size()-1].pt();
//					cout << "    pt_range = " << dpt0 << endl;
//				}
//				if (i == 1) {
//					pt1_sub.clear();
//					e1_sub.clear();
//					m1_sub.clear();
//					double e_total = 0;
//					for (unsigned k = 0; k < subjets.size(); k++) {
//						pt1_sub.push_back( subjets[k].pt() );
//						e1_sub.push_back( subjets[k].e() );
//						m1_sub.push_back( subjets[k].m() );
//						cout << "    Subjet " << k << ": pt = " << subjets[k].pt() << " GeV,\te = " << subjets[k].e() << " GeV" << endl;
//						e_total += subjets[k].e();
//					}
//	//				cout << "    max_index = " << max_index << "  target = " << target << endl;'
//	//				double e_range;
//					dpt1 = subjets[0].pt() - subjets[subjets.size()-1].pt();
//					cout << "    pt_range = " << dpt1 << endl;
//				}
				
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
//			}
//			if (do_prune_) {
//				branches["njets"] = trees["fat"]->Branch("njets", &njets, "njets/I");
//				branches["phi"] = trees["fat"]->Branch("phi", &phi, 64000, 0);
//				branches["y"] = trees["fat"]->Branch("y", &y, 64000, 0);
//				branches["pt"] = trees["fat"]->Branch("pt", &pt, 64000, 0);
//				branches["e"] = trees["fat"]->Branch("e", &e, 64000, 0);
//				branches["m"] = trees["fat"]->Branch("m", &m, 64000, 0);
//				branches["nsub1"] = trees["fat"]->Branch("nsub1", &vec_nsub1, 64000, 0);
//				branches["nsub2"] = trees["fat"]->Branch("nsub2", &vec_nsub2, 64000, 0);
//				branches["nsub3"] = trees["fat"]->Branch("nsub3", &vec_nsub3, 64000, 0);
//				trees["fat"]->Fill();
//			}
//			branches["dpt0"] = trees[jet_string]->Branch("dpt0", &dpt0, "dpt0/D");
//			branches["pt0_sub"] = trees[jet_string]->Branch("pt0_sub", &pt0_sub, 64000, 0);
//			branches["e0_sub"] = trees[jet_string]->Branch("e0_sub", &e0_sub, 64000, 0);
//			branches["m0_sub"] = trees[jet_string]->Branch("m0_sub", &m0_sub, 64000, 0);
//			branches["dpt1"] = trees[jet_string]->Branch("dpt1", &dpt1, "dpt1/D");
//			branches["pt1_sub"] = trees[jet_string]->Branch("pt1_sub", &pt1_sub, 64000, 0);
//			branches["e1_sub"] = trees[jet_string]->Branch("e1_sub", &e1_sub, 64000, 0);
//			branches["m1_sub"] = trees[jet_string]->Branch("m1_sub", &m1_sub, 64000, 0);
//			trees[jet_string]->Fill();
//		}
	}
}

// ------------  called once each job just after ending the event loop  ------------
void Fatjets::endJob()
{
	cout << "Error records for gen clustering:" << endl;
//	cout << n_error_g << endl;
	cout << "* Squark number errors: " << n_error_sq << endl;
	cout << "* Sorting errors: " << n_error_sort << endl;
//	cout << outliers.size() << endl;
//	for (auto o : outliers) {
//		cout << o << "  ";
//	}
//	cout << endl;
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
