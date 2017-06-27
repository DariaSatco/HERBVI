// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {


  /// @brief Add a short analysis description here
  class MC_PRINT : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(MC_PRINT);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(FinalState(), "FS");


    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      /// @todo Do the event by event analysis here
      foreach(const Particle & particle,applyProjection<FinalState>(event,"FS").particles()) {
	cerr << particle << "\n";
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {


    }

    //@}


    /// @name Histograms
    //@{
    Profile1DPtr _h_XXXX;
    Histo1DPtr _h_YYYY;
    CounterPtr _h_ZZZZ;
    //@}


  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(MC_PRINT);


}
