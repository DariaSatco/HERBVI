// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableFinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/MissingMomentum.hh"

namespace Rivet {


  /// Generic analysis looking at various distributions of final state particles
  /// 
  class MY_ANALYSIS_NUMBER : public Analysis {
  public:

    /// Constructor
    MY_ANALYSIS_NUMBER()
      : Analysis("MY_ANALYSIS_NUMBER")
    {    }


  public:

    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Projections
      const FinalState cnfs(Cuts::abseta < 5.0 && Cuts::pT > 500*MeV);
      declare(cnfs, "FS");
      declare(UnstableFinalState(Cuts::abseta < 5.0 && Cuts::pT > 500*MeV), "UFS");
      declare(FastJets(cnfs, FastJets::ANTIKT, 0.4), "Jets");
      declare(MissingMomentum(cnfs), "InclMET");

      // Histograms
      // @todo Choose E/pT ranged based on input energies... can't do anything about kin. cuts, though

      _histStablePIDs  = bookHisto1D("MultsStablePIDs", 3335, -0.5, 3334.5);
      //_histDecayedPIDs = bookHisto1D("MultsDecayedPIDs", 3335, -0.5, 3334.5);
      _histAllPIDs     = bookHisto1D("MultsAllPIDs", 3335, -0.5, 3334.5);

      _histZ       = bookHisto1D("Z_bosons", 50, -0.5, 49.5);
      _histW       = bookHisto1D("W_bosons", 100, -0.5, 70.5);
      _histH       = bookHisto1D("Higgs", 15, -0.5, 14.5);

      _histELECTRON = bookHisto1D("ELECTRONs", 50, -0.5, 49.5);
      _histMUON     = bookHisto1D("MUONs", 30, -0.5, 29.5);
      _histTAU      = bookHisto1D("TAUs", 25, -0.5, 24.5);
      _histNEUTRINO = bookHisto1D("NEUTRINOs", 50, -0.5, 49.5);
      _histLEPTON   = bookHisto1D("LEPTONs", 50, -0.5, 49.5);

      _histTOP      = bookHisto1D("TOP_quarks", 50, -0.5, 49.5);
      _histBOTTOM   = bookHisto1D("BOTTOM_quarks", 50, -0.5, 49.5);
      _hist_n_jet   = bookHisto1D("n-jet", 100, -0.5, 150.5);

      _h_met_incl  = bookHisto1D("met_incl", logspace(30, 1, 6000));
      _hist_ht_jet = bookHisto1D("HT_jets", 50, -0.5, 40000);
      _hist_m_eff  = bookHisto1D("M_eff_jets", 50, -0.5, 40000);

      _hist_lept_jet = bookHisto2D("leptons_jets_2D", 160, -0.5, 159.5, 60, -0.5, 59.5); 
      _hist_lept_and_jet = bookHisto1D("mult_leptons_jets", 150, -0.5, 199.5);
    }



    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const double weight = event.weight();

      // Unphysical (debug) plotting of all PIDs in the event, physical or otherwise
      foreach (const GenParticle* gp, particles(event.genEvent())) {
        _histAllPIDs->fill(abs(gp->pdg_id()), weight);
      }

      // Charged + neutral final state PIDs
      const FinalState& cnfs = apply<FinalState>(event, "FS");
      foreach (const Particle& p, cnfs.particles()) {
        _histStablePIDs->fill(p.abspid(), weight);
      }

      // Number of Jets, transverse momentum
      const FastJets& jetpro = apply<FastJets>(event, "Jets");
      const Jets jets = jetpro.jetsByPt();
      MSG_DEBUG("Jet multiplicity = " << jets.size());
      _hist_n_jet->fill(jets.size(), weight);
      double HT=0.0;
      int njet=0;
      foreach (const Jet& j, jets) {
      const FourMomentum& pj = j.momentum();
      ++njet;
      HT+=pj.pT()/GeV;
      }

      _hist_ht_jet->fill(HT, weight);
	

      //Missing momentum
      const MissingMomentum& mmincl = apply<MissingMomentum>(event, "InclMET");
      _h_met_incl->fill(mmincl.met()/GeV, weight);
      
      //M_eff
      _hist_m_eff->fill(HT+mmincl.met()/GeV, weight);

      // Multiplicity distribution for different final state PIDs
      const UnstableFinalState& ufs = apply<UnstableFinalState>(event, "UFS");

      int nZ=0;
      int nW=0;
      int nH=0;
      int nELECTRON=0;
      int nNEUTRINO=0;
      int nMUON=0;
      int nTAU=0;
      int nBOTTOM=0;
      int nTOP=0;
      

      foreach (const GenParticle* gp, particles(event.genEvent())) {
	//        _histDecayedPIDs->fill(p.pid(), weight);
	//        const double eta_abs = p.abseta();
        const PdgId pid = abs(gp->pdg_id()); //if (PID::isMeson(pid) && PID::hasStrange()) {
        if (pid == 23) ++nZ;
        else if (pid == 24 ) ++nW;
        else if (pid == 25) ++nH;

	else if (pid == 11) ++nELECTRON;
        else if (pid == 12 || pid == 14 || pid == 16) ++nNEUTRINO;
        else if (pid == 13) ++nMUON;

	else if (pid == 15) ++nTAU;
	else if (pid == 5) ++nBOTTOM;
        else if (pid == 6) ++nTOP;

      }
      
      _histZ->fill(nZ, weight);
      _histW->fill(nW, weight);
      _histH->fill(nH, weight);
      _histELECTRON->fill(nELECTRON, weight);
      _histNEUTRINO->fill(nNEUTRINO, weight);
      _histMUON->fill(nMUON, weight);
      _histTAU->fill(nTAU, weight);
      _histLEPTON->fill(nELECTRON+nMUON, weight);
      _histBOTTOM->fill(nBOTTOM, weight);
      _histTOP->fill(nTOP, weight);
      _hist_lept_jet->fill(njet, nELECTRON+nMUON, weight);
      _hist_lept_and_jet->fill(njet+nELECTRON+nMUON, weight);

    }



    /// Finalize
    void finalize() {
      scale(_histStablePIDs, 1/sumOfWeights());
      //      scale(_histDecayedPIDs, 1/sumOfWeights());
      scale(_histAllPIDs, 1/sumOfWeights());
     
      scale(_histZ, 1/sumOfWeights());
      scale(_histW,  1/sumOfWeights());
      scale(_histH,  1/sumOfWeights());
      scale(_histELECTRON,  1/sumOfWeights());
      scale(_histMUON,  1/sumOfWeights());
      scale(_histNEUTRINO,  1/sumOfWeights());
      scale(_histTAU,  1/sumOfWeights());
      scale(_histLEPTON, 1/sumOfWeights());
      scale(_histTOP,  1/sumOfWeights());
      scale(_histBOTTOM,  1/sumOfWeights());
      scale(_hist_n_jet, 1/sumOfWeights());
      scale(_hist_ht_jet, 1/sumOfWeights());
      scale(_hist_m_eff, 1/sumOfWeights());
      scale(_hist_lept_jet, 1/sumOfWeights());
      scale(_hist_lept_and_jet, 1/sumOfWeights());
    }

    //@}


  private:

    /// @name Histograms
    //@{
    Histo1DPtr _histStablePIDs, _histAllPIDs;
    Histo1DPtr _histZ,  _histW,  _histH,  _histELECTRON,  _histMUON,  _histNEUTRINO, _histLEPTON;
    Histo1DPtr _histBOTTOM,  _histTOP,  _histTAU;
    Histo1DPtr _hist_n_jet;
    Histo1DPtr _h_met_incl;
    Histo1DPtr _hist_ht_jet, _hist_m_eff;
    Histo2DPtr _hist_lept_jet;
    Histo1DPtr _hist_lept_and_jet;
    //@}

  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(MY_ANALYSIS_NUMBER);

}
