#include <memory>
#include <stdexcept>
#include <string>

#include "RAT/VertexGen_Freya.hh"
#include "RAT/Log.hh"
#include "RAT/DB.hh"
#include "Randomize.hh" // G4UniformRand()

#include "G4Event.hh"
#include "G4PrimaryVertex.hh"
#include "G4DynamicParticle.hh"
#include "G4ThreeVector.hh"
#include "G4SystemOfUnits.hh"
#include "G4Neutron.hh"
#include "G4Gamma.hh"

// Anonymous namespace for definitions local to this source file
namespace {

  constexpr double DUMMY_AVERAGE_NEUTRON_YIELD = -1.;
  constexpr double DUMMY_INCIDENT_NEUTRON_ENERGY = 0.;

  // FREYA spontaneous fission flag
  constexpr int SPONTANEOUS_FISSION = 0;

  // Flag for LLNL fission library to use FREYA to generate events with full
  // particle correlations
  constexpr int LLNL_USE_FREYA = 3;

  // Maximum string length to reserve for error message strings retrieved from
  // FREYA
  constexpr int FREYA_MESSAGE_MAX_LENGTH = 1000;
}


RAT::VertexGen_Freya::VertexGen_Freya(const char* arg_dbname)
  : GLG4VertexGen(arg_dbname)
{

  // Use the same random number generator as Geant4 in FREYA to create
  // fission events
  fissionEvent::setRNGd( []() -> double { return G4UniformRand(); } );

  // Instruct the LLNL library to use FREYA to generate fission events
  fissionEvent::setCorrelationOption(LLNL_USE_FREYA);

  //SetState("1 400");
}

void RAT::VertexGen_Freya::GeneratePrimaryVertex(G4Event* event,
  G4ThreeVector& dx, G4double dt)
{
  // Use FREYA to make a fission event using the current isotope at the
  // specified time.
  std::unique_ptr<fissionEvent> fe = std::make_unique<fissionEvent>(
    fIsotopeZA, dt / second, DUMMY_AVERAGE_NEUTRON_YIELD,
    DUMMY_INCIDENT_NEUTRON_ENERGY, SPONTANEOUS_FISSION);

  if (fe->getCorrelationOption() != LLNL_USE_FREYA)
    throw std::runtime_error("Error in VertexGen_Freya::"
      "GeneratePrimaryVertex(). FREYA failed to load correctly.");
  else {
    // Check for error complaints from FREYA itself
    int message_length = FREYA_MESSAGE_MAX_LENGTH;
    std::string freya_message(message_length, ' ');
    fe->getFREYAerrors( &message_length, &(freya_message.front()) );
    if (message_length > 1) {
      freya_message.resize(message_length);
      // If FREYA encountered an error, then it will return a string
      // with length > 1
      throw std::runtime_error("FREYA error \"" + freya_message
        + "\" encountered in VertexGen_Freya::GeneratePrimaryVertex().");

    }
  }

  // create a new vertex
  G4PrimaryVertex* vertex = new G4PrimaryVertex(dx, dt);

  double momentum, px, py, pz, kE;

  // Create the primary neutrons
  int num_neutrons = fe->getNeutronNu();
  for (int i = 0; i < num_neutrons; ++i) {

    G4PrimaryParticle* particle = create_primary_neutron(*fe, i);

    vertex->SetPrimary(particle);
  }

  // Create the primary gammas
  int num_gammas = fe->getPhotonNu();
  for (int i = 0; i < num_gammas; ++i) {

    G4PrimaryParticle* particle = create_primary_gamma(*fe, i);

    vertex->SetPrimary(particle);
  }

  event->AddPrimaryVertex(vertex);
}

void RAT::VertexGen_Freya::SetState(G4String newValues)
{
  if (newValues.length() == 0) {
    // print help and current state
    G4cout << "Current state of this VertexGen_Freya:\n"
           << " \"" << GetState() << "\"\n" << G4endl;
    G4cout << "Format of argument to VertexGen_Freya::SetState: \n"
      " \"ZA\"\n" << G4endl << "where Z is the isotope atomic number,"
      " A is the mass number, and ZA = 1000*Z + A\n" << G4endl;
    return;
  }

  std::istringstream is(newValues.c_str());
  int ZA;

  bool ok = static_cast<bool>(is >> ZA);

  if (ok) fIsotopeZA = ZA;
  else Log::Die("VertexGen_Freya: Incorrect vertex setting " + newValues);
}

G4String RAT::VertexGen_Freya::GetState()
{
  return std::to_string(fIsotopeZA);
}

// Creates and returns a pointer to a new G4Primary particle with a given
// particle definition, kinetic energy, and momentum direction unit vector.
G4PrimaryParticle* RAT::VertexGen_Freya::create_primary(
  const G4ParticleDefinition& particle_definition, double kinetic_energy,
  double dir_x, double dir_y, double dir_z) const
{
  G4ThreeVector dir_vec( dir_x, dir_y, dir_z );

  if (dir_vec.mag2() <= 0.) {
    throw std::runtime_error("Invalid"
      " direction vector (" + std::to_string(dir_x)
      + ", " + std::to_string(dir_y) + ", " + std::to_string(dir_z)
      + ") passed to RAT::VertexGen_Freya::create_primary()");
  }

  // Use a temporary G4DynamicParticle so that we can convert
  // the kinetic energy and direction given by FREYA to a momentum
  // vector. This ensures we use the same particle mass as Geant4.
  G4DynamicParticle temp_particle(&particle_definition, dir_vec.unit(),
    kinetic_energy);

  G4ThreeVector momentum = temp_particle.GetMomentum();

  // Create the Geant4 primary particle
  G4PrimaryParticle* primary_particle = new G4PrimaryParticle(
    &particle_definition, momentum.x(), momentum.y(), momentum.z());

  return primary_particle;
}

// Creates a Geant4 primary particle for the index-th neutron from a
// fissionEvent object.
G4PrimaryParticle* RAT::VertexGen_Freya::create_primary_neutron(
  fissionEvent& fission_event, int index) const
{
  if (index >= fission_event.getNeutronNu()) throw std::runtime_error("Index"
    " out of range in RAT::VertexGen_Freya::create_primary_neutron()");

  double kinetic_energy = fission_event.getNeutronEnergy(index) / MeV;

  // Momentum direction unit vector components
  double dir_x = fission_event.getNeutronDircosu(index);
  double dir_y = fission_event.getNeutronDircosv(index);
  double dir_z = fission_event.getNeutronDircosw(index);

  return create_primary( *G4Neutron::Definition(), kinetic_energy, dir_x,
    dir_y, dir_z );
}

// Creates a Geant4 primary particle for the index-th gamma-ray from a
// fissionEvent object.
G4PrimaryParticle* RAT::VertexGen_Freya::create_primary_gamma(
  fissionEvent& fission_event, int index) const
{
  if (index >= fission_event.getPhotonNu()) throw std::runtime_error("Index"
    " out of range in RAT::VertexGen_Freya::create_primary_gamma()");

  double kinetic_energy = fission_event.getPhotonEnergy(index) / MeV;

  // Momentum direction unit vector components
  double dir_x = fission_event.getPhotonDircosu(index);
  double dir_y = fission_event.getPhotonDircosv(index);
  double dir_z = fission_event.getPhotonDircosw(index);

  return create_primary( *G4Gamma::Definition(), kinetic_energy, dir_x,
    dir_y, dir_z );
}
