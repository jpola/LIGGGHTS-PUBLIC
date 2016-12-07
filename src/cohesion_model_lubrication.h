/* ----------------------------------------------------------------------
    This is the

    ██╗     ██╗ ██████╗  ██████╗  ██████╗ ██╗  ██╗████████╗███████╗
    ██║     ██║██╔════╝ ██╔════╝ ██╔════╝ ██║  ██║╚══██╔══╝██╔════╝
    ██║     ██║██║  ███╗██║  ███╗██║  ███╗███████║   ██║   ███████╗
    ██║     ██║██║   ██║██║   ██║██║   ██║██╔══██║   ██║   ╚════██║
    ███████╗██║╚██████╔╝╚██████╔╝╚██████╔╝██║  ██║   ██║   ███████║
    ╚══════╝╚═╝ ╚═════╝  ╚═════╝  ╚═════╝ ╚═╝  ╚═╝   ╚═╝   ╚══════╝®

    DEM simulation engine, released by
    DCS Computing Gmbh, Linz, Austria
    http://www.dcs-computing.com, office@dcs-computing.com

    LIGGGHTS® is part of CFDEM®project:
    http://www.liggghts.com | http://www.cfdem.com

    Core developer and main author:
    Christoph Kloss, christoph.kloss@dcs-computing.com

    LIGGGHTS® is open-source, distributed under the terms of the GNU Public
    License, version 2 or later. It is distributed in the hope that it will
    be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
    of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. You should have
    received a copy of the GNU General Public License along with LIGGGHTS®.
    If not, see http://www.gnu.org/licenses . See also top-level README
    and LICENSE files.

    LIGGGHTS® and CFDEM® are registered trade marks of DCS Computing GmbH,
    the producer of the LIGGGHTS® software and the CFDEM®coupling software
    See http://www.cfdem.com/terms-trademark-policy for details.

-------------------------------------------------------------------------
    Contributing author and copyright for this file:

    Jakub Pola (UWr Wroclaw, Poland)
------------------------------------------------------------------------- */

#ifdef COHESION_MODEL
COHESION_MODEL(COHESION_LUBRICATION,lubrication,16)
#else

#ifndef COHESION_MODEL_LUBRICATION_H
#define COHESION_MODEL_LUBRICATION_H

#include "contact_models.h"
#include "math.h"
#include "math_extra_liggghts.h"
#include "float.h"

#include "global_properties.h"

// Code for managing the model parameters
namespace MODEL_PARAMS
{
inline static ScalarProperty* createFluidDynamicViscosity(PropertyRegistry &registry, const char *caller, bool sanity_checks)
{
  //I use sanity checks to be sure that mu will be positive
  ScalarProperty* mu = MODEL_PARAMS::createScalarProperty(registry, "mu", caller, sanity_checks, 0);
  return mu;
}

inline static ScalarProperty* createCutoffDistance(PropertyRegistry &registry, const char *caller, bool sanity_checks)
{
  ScalarProperty* cutoffDistance = MODEL_PARAMS::createScalarProperty(registry, "cutoff", caller, sanity_checks, DBL_MIN, DBL_MAX);
  return cutoffDistance;
}

inline static ScalarProperty* createEffectiveRoughness(PropertyRegistry &registry, const char *caller, bool sanity_checks)
{
  //JPA: Joseph et al. 2001 eta/R = etaEff = [10e-6, 10e-3]
  ScalarProperty* effectiveRoughness = MODEL_PARAMS::createScalarProperty(registry, "etaEff", caller, sanity_checks, 1e-6, 1e-3);
  return effectiveRoughness;
}

}

namespace LIGGGHTS {
namespace ContactModels {

using namespace std;
using namespace LAMMPS_NS;

template<>
class CohesionModel<COHESION_LUBRICATION> : protected Pointers
{
public:
  static const int MASK = CM_CONNECT_TO_PROPERTIES | CM_SURFACES_INTERSECT;

  int bond_history_offset() { return -1; }

  CohesionModel(LAMMPS *lmp, IContactHistorySetup*, class ContactModelBase *)
    : Pointers(lmp), mu(1.0), cutoff_distance(DBL_MIN)
  { }

  void registerSettings(Settings &settings)
  {

  }

  void connectToProperties(PropertyRegistry &registry)
  {

    //Property registry is dedicated to read from fix/property options;
    registry.registerProperty("mu", &MODEL_PARAMS::createFluidDynamicViscosity);
    registry.registerProperty("cutoff", &MODEL_PARAMS::createCutoffDistance);
    registry.registerProperty("etaEff", &MODEL_PARAMS::createEffectiveRoughness);

    registry.connect("mu", mu, "cohesion_model lubrication");
    registry.connect("cutoff", cutoff_distance, "cohesion_model lubrication");
    registry.connect("etaEff", eta_eff, "cohesion_model lubrication");


    cout << "LUBE: connectToProperties" << endl;
    cout << "LUBE: the mu : " << mu << endl;
    cout << "LUBE: the cutoff_distance : " << cutoff_distance << endl;
    cout << "LUBE: the effective roughness : " << eta_eff << endl;

  }

  void surfacesIntersect(SurfacesIntersectData &sidata, ForceData &i_forces, ForceData &j_forces)
  {
    cout << "Lubrication intersection: does nothing " << endl;
  }

  void beginPass(SurfacesIntersectData&, ForceData&, ForceData&){}

  void endPass(SurfacesIntersectData&, ForceData&, ForceData&){}

  void wallParticleClose(SurfacesCloseData& scdata, ForceData& i_forces, ForceData& j_forces)
  {
    double **v = atom->v;

    const int i = scdata.i;
    const int j = scdata.j;

    const double radi = scdata.radi;

    // r is our delta
    const double r = sqrt(scdata.rsq);
    // dist is our deltan distance between intersection points - Bonometti treats this as dist between particles
    const double dist =  r - radi;

    double R = radi;

    if (dist >= R / 2.0)
    {
      cout << "Particle is too far from wall dist = " << dist << " R/2 = " <<  R/2.0 << " we dont want to calculate" << endl;
      return;
    }

    cout << "Particles approaching to wall distance is ok dist = " << dist << " R/2 = " <<  R/2.0 << endl;

    const double rEff = radi;

    // calculate vn and vt since not in struct
    const double rinv = 1.0 / r;

    const double dx = scdata.delta[0];
    const double dy = scdata.delta[1];
    const double dz = scdata.delta[2];

    const double enx = dx * rinv;
    const double eny = dy * rinv;
    const double enz = dz * rinv;

    // relative translational velocity.
    // surfaceIntersectData have information about velocities of components they are v_i and v_j
    // but surfaceCloseData which is base of it does not have it.
    // I assumed that wall is not moving;
    const double vr1 = v[i][0];
    const double vr2 = v[i][1];
    const double vr3 = v[i][2];

    // normal component
    const double vn = vr1 * enx + vr2 * eny + vr3 * enz;
    const double vn1 = vn * enx;
    const double vn2 = vn * eny;
    const double vn3 = vn * enz;

    //For the future ...
    // tangential component
    const double vt1 = vr1 - vn1;
    const double vt2 = vr2 - vn2;
    const double vt3 = vr3 - vn3;

    // relative rotational velocity
    double const *omega_i = atom->omega[i];

    double wr1 = radi * omega_i[0] * rinv;
    double wr2 = radi * omega_i[1] * rinv;
    double wr3 = radi * omega_i[2] * rinv;

    // relative velocities
    const double vtr1 = vt1 - (dz * wr2 - dy * wr3);
    const double vtr2 = vt2 - (dx * wr3 - dz * wr1);
    const double vtr3 = vt3 - (dy * wr1 - dx * wr2);

    const double stokesPreFactor = -6.0*M_PI*mu*rEff*rEff;


    // JPA: here we should add eta even if the particles are not intersecting?
    //const double deltaijInv = 1.0 / MathExtraLiggghts::max(dist, cutoff_distance) ;
    const double deltaijInv = 1.0 / (dist + (eta_eff*R));


    const double Fn = stokesPreFactor * vn * deltaijInv;

    cout << "FN: " << Fn  << " at distance " << (dist + (eta_eff*R))<< endl;

    const double fx = Fn * enx;
    const double fy = Fn * eny;
    const double fz = Fn * enz;

    scdata.has_force_update = true;

    // return resulting forces
    i_forces.delta_F[0] += fx;
    i_forces.delta_F[1] += fy;
    i_forces.delta_F[2] += fz;
//  i_forces.delta_torque[0] = -radi * tor1; // using radius here, not contact radius
//  i_forces.delta_torque[1] = -radi * tor2;
//  i_forces.delta_torque[2] = -radi * tor3;

    j_forces.delta_F[0] -= fx;
    j_forces.delta_F[1] -= fy;
    j_forces.delta_F[2] -= fz;
//  j_forces.delta_torque[0] = -radj * tor1; // using radius here, not contact radius
//  j_forces.delta_torque[1] = -radj * tor2;
//  j_forces.delta_torque[2] = -radj * tor3;
  }

  void particlesClose(SurfacesCloseData& scdata, ForceData& i_forces, ForceData& j_forces)
  {
    double **v = atom->v;

    const int i = scdata.i;
    const int j = scdata.j;
    const int itype = atom->type[i];
    const int jtype = atom->type[j];
    const double radi = scdata.radi;
    const double radj = scdata.radj;
    // r is our delta
    const double r = sqrt(scdata.rsq);
    // dist is our deltan distance between intersection points - Bonometti treats this as dist between particles
    const double dist =  r - (radi + radj);

    double R = 0.0;
    if (radi < radj || radi > radj )
    {
      cout << "Particles are not equal" << endl;
      R = (radi + radj) / 2.0;
    }
    else
    {
      R = radi;
    }
    if (dist >= R / 2.0)
    {
      cout << "Particles are too far dist = " << dist << " R/2 = " <<  R/2.0 << " we dont want to calculate" << endl;
      return;
    }

    cout << "LUBE: Particles are approaching and distance is ok dist = " << dist << " R/2 = " <<  R/2.0 << endl;


    const double rEff = radi*radj / (radi+radj);

    // calculate vn and vt since not in struct
    const double rinv = 1.0 / r;

    const double dx = scdata.delta[0];
    const double dy = scdata.delta[1];
    const double dz = scdata.delta[2];

    const double enx = dx * rinv;
    const double eny = dy * rinv;
    const double enz = dz * rinv;

    // relative translational velocity
    const double vr1 = v[i][0] - v[j][0];
    const double vr2 = v[i][1] - v[j][1];
    const double vr3 = v[i][2] - v[j][2];

    // normal component
    const double vn = vr1 * enx + vr2 * eny + vr3 * enz;
    const double vn1 = vn * enx;
    const double vn2 = vn * eny;
    const double vn3 = vn * enz;

    //For the future ...
    // tangential component
    const double vt1 = vr1 - vn1;
    const double vt2 = vr2 - vn2;
    const double vt3 = vr3 - vn3;

    // relative rotational velocity
    double wr1, wr2, wr3;
    double const *omega_i = atom->omega[i];
    double const *omega_j = atom->omega[j];

    if(scdata.is_wall) {
      wr1 = radi * omega_i[0] * rinv;
      wr2 = radi * omega_i[1] * rinv;
      wr3 = radi * omega_i[2] * rinv;
    } else {
      wr1 = (radi * omega_i[0] + radj * omega_j[0]) * rinv;
      wr2 = (radi * omega_i[1] + radj * omega_j[1]) * rinv;
      wr3 = (radi * omega_i[2] + radj * omega_j[2]) * rinv;
    }

    // relative velocities
    const double vtr1 = vt1 - (dz * wr2 - dy * wr3);
    const double vtr2 = vt2 - (dx * wr3 - dz * wr1);
    const double vtr3 = vt3 - (dy * wr1 - dx * wr2);


    const double stokesPreFactor = -6.0*M_PI*mu*rEff*rEff;


    // JPA: here we should add eta even if the particles are not intersecting?
    //const double deltaijInv = 1.0 / MathExtraLiggghts::max(dist, cutoff_distance) ;
    const double deltaijInv = 1.0 / (dist + (eta_eff*R));


    const double Fn = stokesPreFactor * vn * deltaijInv;

    cout << "FN: " << Fn  << " at distance " << (dist + (eta_eff*R))<< endl;

    const double fx = Fn * enx;
    const double fy = Fn * eny;
    const double fz = Fn * enz;

    scdata.has_force_update = true;

    // return resulting forces
    if(scdata.is_wall) {
      cout << "WALL?"<< endl;
    }
    else
    {
      i_forces.delta_F[0] += fx;
      i_forces.delta_F[1] += fy;
      i_forces.delta_F[2] += fz;
      //      i_forces.delta_torque[0] = -radi * tor1; // using radius here, not contact radius
      //      i_forces.delta_torque[1] = -radi * tor2;
      //      i_forces.delta_torque[2] = -radi * tor3;

      j_forces.delta_F[0] -= fx;
      j_forces.delta_F[1] -= fy;
      j_forces.delta_F[2] -= fz;
      //      j_forces.delta_torque[0] = -radj * tor1; // using radius here, not contact radius
      //      j_forces.delta_torque[1] = -radj * tor2;
      //      j_forces.delta_torque[2] = -radj * tor3;
    }

  }

  void surfacesClose(SurfacesCloseData& scdata, ForceData& i_forces, ForceData& j_forces)
  {

    if (scdata.is_wall)
    {
      cout << "Particle approaching to wall" << endl;
      wallParticleClose(scdata, i_forces, j_forces);
    }
    else
    {
      cout << "Particle approacing to other particle" << endl;
      particlesClose(scdata, i_forces, j_forces);
    }
  }

private:
  //fluid dynamic viscosity
  double mu;
  double cutoff_distance;
  double eta_eff;

};

}
} //LIGGGHTS


#endif // COHESION_MODEL_LUBRICATION_H
#endif // COHESION_MODEL
