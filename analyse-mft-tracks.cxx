#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"

#include "TDatabasePDG.h"

using namespace o2;
using namespace o2::framework;

using Particles = aod::McParticles;



//First approach to analysing an AO2D.root file
//written thanks to https://aliceo2group.github.io/analysis-framework/docs/tutorials/analysistask.html




struct analyseMFTTracks
{
  int icoll = 0;

  HistogramRegistry registry
  {
     "registry",
     {
       {"TracksPhiEta", "; #varphi; #eta; tracks", {HistType::kTH2F, {{600, 0, 2 * M_PI}, {41, -5., -1.}}}},            //
       {"TracksPhiEtaGen", "; #varphi; #eta; tracks", {HistType::kTH2F, {{600, 0, 2 * M_PI}, {21, -2.1, 2.1}}}},
       {"Multiplicity", "alibi_nightlies/O2DPG_pp_minbias_testbeam.sh/15-11-2021-18:00 - tf13; collisionID; N_{trk}^{MFT}", {HistType::kTH1F, {{101, 0., 100.}}}, true}         //
     }                                                                                                                   //
   };

  void process(o2::aod::Collision const& collision, o2::aod::MFTTracks const& tracks)
    {
      //using o2::aod::MFTTracks;
      //using namespace o2::aod::fwdtrack;
      int n = 0;

        for (auto& track : tracks)
        {
          registry.fill(HIST("TracksPhiEta"), track.phi(), track.eta());
          registry.fill(HIST("Multiplicity"), icoll);
          n++;
        }
      icoll++;

    }
    //end of process


};
//end of MyTask


// struct analyseGenTracks
// {
//   void processGen(Particles const& particles)
//     {
//         for (auto& particle : particles)
//         {
//           registry.fill(HIST("TracksPhiEtaGen"), particle.phi(), particle.eta());
//         }
//
//     }
//     //end of processGen
// };

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<analyseMFTTracks>(cfgc),
    //adaptAnalysisTask<analyseGenTracks>(cfgc),
  };
}
