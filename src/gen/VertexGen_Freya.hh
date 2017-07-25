////////////////////////////////////////////////////////////////////////
/// @class RAT::VertexGen_Freya
///
/// @brief Generates spontaneous fission events using FREYA
///
/// @author Steven Gardiner <sjgardiner@ucdavis.edu>
///
/// @details Allows the user to generate spontaneous fission events for
///          a limited number of isotopes. For supported isotopes, FREYA
///          provides realistic neutron and gamma multiplicities, energies,
///          directional correlations, etc. One may add this generator
///          in RAT-PAC using the macro command
///          @code /generator/add combo freya[:POSITION][:TIME] @endcode
///
////////////////////////////////////////////////////////////////////////

#ifndef __RAT_VertexGen_Freya__
#define __RAT_VertexGen_Freya__

// Geant4 includes
#include "G4PrimaryParticle.hh"

// RAT includes
#include "RAT/GLG4VertexGen.hh"

// LLNL fission library includes
#include "fissionEvent.h"

namespace RAT {

  class VertexGen_Freya : public GLG4VertexGen {
    public:

      VertexGen_Freya(const char* arg_dbname = "freya");

      virtual void GeneratePrimaryVertex(G4Event* event,
        G4ThreeVector& dx, G4double dt) override;

      // The state setting should contain a single ZA value indicating
      // which isotope should be simulated
      virtual void SetState( G4String newValues ) override;
      virtual G4String GetState() override;

    protected:

      // Helper functions for creating Geant4 primary particles from the FREYA
      // results
      G4PrimaryParticle* create_primary(
        const G4ParticleDefinition& particle_definition,
        double kinetic_energy, double dir_x, double dir_y, double dir_z) const;

      G4PrimaryParticle* create_primary_neutron(fissionEvent& fission_event,
        int index) const;

      G4PrimaryParticle* create_primary_gamma(fissionEvent& fission_event,
        int index) const;

      // @brief Integer code representing the fissioning isotope.
      // @details For an atomic number Z and mass number A, the code
      // is given by Z*1000 + A.
      int fIsotopeZA;
  };

} // namespace RAT

#endif
