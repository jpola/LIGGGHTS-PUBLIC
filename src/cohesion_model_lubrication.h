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
    {
        cout << "LUBE: Ctor lubrication" << endl;
    }

    void registerSettings(Settings &settings)
    {

        //define my settings like nu, cuttoff
        cout << "LUBE: registerSettings" << endl;
        cout << "\n LUBE: Warning! Be sure that mu is set up correctly!" << endl;

        //TODO:: Define params here? Or in the fix/property which is more restrictive

        settings.registerDoubleSetting("mu", mu, 2.0);
        settings.registerDoubleSetting("cutoff", cutoff_distance, 1.e-2);

    }

    void connectToProperties(PropertyRegistry &registry)
    {
        cout << "LUBE: connectToProperties" << endl;
        cout << "LUBE: the mu : " << mu << endl;
        cout << "LUBE: the cutoff_distance : " << cutoff_distance << endl;


        /**
        //Property registry is dedicated to read from fix/property options;
        registry.registerProperty("mu", &MODEL_PARAMS::createFluidDynamicViscosity);
        registry.registerProperty("cutoff", &MODEL_PARAMS::createFluidDynamicViscosity);

        registry.connect("mu", mu, "cohesion_model lubrication");
        registry.connect("cutoff", cutoff_distance, "cohesion_model lubrication");
        */
    }

    void surfacesIntersect(SurfacesIntersectData &sidata, ForceData &i_forces, ForceData &j_forces)
    {
        cout << "LUBE: surfaceIntersect: ("
             << sidata.delta[0] << ", "
             << sidata.delta[1] << ", "
             << sidata.delta[2] << ") "
             << "\t("
             << sidata.en[0] << ", "
             << sidata.en[1] << ", "
             << sidata.en[2] << ") "<< endl;

    }

    void beginPass(SurfacesIntersectData&, ForceData&, ForceData&){}

    void endPass(SurfacesIntersectData&, ForceData&, ForceData&){}

    void surfacesClose(SurfacesCloseData& scdata, ForceData& i_forces, ForceData& j_forces)
    {

        if (scdata.is_wall)
        {
            cout << "LUBE:  we are close to wall " << scdata.i << ", " << scdata.j
                 << "\t("
                 << scdata.delta[0] << ", "
                 << scdata.delta[1] << ", "
                 << scdata.delta[2] << ") " << endl;

        }
        else
        {
            double **v = atom->v;

            const int i = scdata.i;
            const int j = scdata.j;
            const int itype = atom->type[i];
            const int jtype = atom->type[j];
            const double radi = scdata.radi;
            const double radj = scdata.radj;
            const double r = sqrt(scdata.rsq);
            const double dist =  r - (radi + radj);
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


            const double deltaijInv = 1.0 / MathExtraLiggghts::max(dist, cutoff_distance);


            const double Fn = stokesPreFactor * vn * deltaijInv;

            cout << "FN: " << Fn  << " at distance " << MathExtraLiggghts::max(dist, cutoff_distance) << endl;

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
//              i_forces.delta_torque[0] = -radi * tor1; // using radius here, not contact radius
//              i_forces.delta_torque[1] = -radi * tor2;
//              i_forces.delta_torque[2] = -radi * tor3;

              j_forces.delta_F[0] -= fx;
              j_forces.delta_F[1] -= fy;
              j_forces.delta_F[2] -= fz;
//              j_forces.delta_torque[0] = -radj * tor1; // using radius here, not contact radius
//              j_forces.delta_torque[1] = -radj * tor2;
//              j_forces.delta_torque[2] = -radj * tor3;
            }



//            cout << "LUBE:  surface close " << scdata.i << ", " << scdata.j
//                 <<  "\t("<< dx << ", " << dy << ", " << dz <<")" << endl;


        }


      //  if(scdata.contact_flags) *scdata.contact_flags &= ~CONTACT_COHESION_MODEL;

    }

private:
    //fluid dynamic viscosity
    double mu;
    double cutoff_distance;

};

}
} //LIGGGHTS


#endif // COHESION_MODEL_LUBRICATION_H
#endif // COHESION_MODEL
