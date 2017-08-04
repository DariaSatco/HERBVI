extern"C" {
  void herwiginit_(int *,int *,double *);
  void herwiggen_();
  void herwigend_();
  void crosssection_(double *);
}

#include "Rivet/Rivet.hh"
#include "Rivet/AnalysisHandler.hh"
#include "HepMC/IO_HERWIG.h"
int main() {

  HepMC::HEPEVT_Wrapper::set_max_number_entries(15000);
  HepMC::HEPEVT_Wrapper::set_sizeof_real(8);
  HepMC::HEPEVT_Wrapper::set_sizeof_int(4);


  // read the ibeam iproc and energy
  int ibeam;
  std::cin >> ibeam;
  int iproc;
  std::cin >> iproc;
  double energy;
  std::cin >> energy;
  int nevent;
  std::cin >> nevent;
  std::string name_of_file;
  std::cin >> name_of_file;
  std::string path_of_file;
  std::cin >> path_of_file;
  std::ofstream out(path_of_file+name_of_file+".hepmc"); 

  // create Rivet analysis handler
  Rivet::AnalysisHandler *_rivet = new Rivet::AnalysisHandler("");
  _rivet->setIgnoreBeams();
  std::string name;  
  while(true) {
  std::cin  >> name;
  if(name[0]=='#') break;
    _rivet->addAnalysis(name);
  };
  herwiginit_(&ibeam,&iproc,&energy);
  // create HepMC converter
  HepMC::IO_HERWIG _converter;
  _converter.set_print_inconsistency_errors(false);
  // main event loop
  std::cout <<"generating nevents="<< nevent<<std::endl;
  for(unsigned int ix=0;ix<nevent;++ix) {
    herwiggen_();
    HepMC::GenEvent * event = new HepMC::GenEvent();
    event->use_units(HepMC::Units::GEV, HepMC::Units::MM);
    _converter.fill_next_event(event);
    //    HepMC::HEPEVT_Wrapper::print_hepevt();
    
    event->write( out );
    _rivet->analyze(*event);
    delete event;
  }
  herwigend_();
  double cross;
  crosssection_(&cross);
  _rivet->setCrossSection(1000.*cross);
  _rivet->finalize();
  _rivet->writeData(path_of_file+name_of_file+".yoda");
  delete _rivet;
  return 0;
}
