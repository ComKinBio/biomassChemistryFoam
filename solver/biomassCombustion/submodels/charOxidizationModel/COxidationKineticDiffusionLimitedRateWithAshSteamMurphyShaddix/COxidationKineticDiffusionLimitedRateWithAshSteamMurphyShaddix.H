/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::COxidationKineticDiffusionLimitedRateWithAsh

Description
    Kinetic/diffusion limited rate surface reaction model for coal parcels.
    Limited to:

        Omiga*C(s) + O2 -> 2*(Omiga-1)CO + (2-Omiga)CO2

    where Omiga is a correlation given by Thumman, H (2002) Combustion of wood particle - a  particle model for Eulerian calculations.

\*---------------------------------------------------------------------------*/

#ifndef COxidationKineticDiffusionLimitedRateWithAshSteamMurphyShaddix_H
#define COxidationKineticDiffusionLimitedRateWithAshSteamMurphyShaddix_H

#include "CharOxidizationModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward class declarations
template<class CloudType>
class COxidationKineticDiffusionLimitedRateWithAshSteamMurphyShaddix;

/*---------------------------------------------------------------------------*\
            Class COxidationKineticDiffusionLimitedRate Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class COxidationKineticDiffusionLimitedRateWithAshSteamMurphyShaddix
:
    public CharOxidizationModel<CloudType>
{
    // Private data

        // Model constants

            //- Omiga correlation constant
            const scalar COmega_;

            //- porosity of the ash layer
            const scalar ep3_;

            //- Mass diffusion limited rate constant, C1 constant in gas phase diffusion coefficient
            const scalar C1_O2_;
            const scalar C1_H2O_;
            const scalar C1_CO2_;
            const scalar C1_H2_;

            //- Kinetics limited rate pre-exponential constant, C2
            const scalar C2_1_;
            const scalar C2_2_;
            const scalar C2_3_;
            const scalar C2_4_;

            //- Kinetics limited rate activation energy RR together
            const scalar E_1_;
            const scalar E_2_;
            const scalar E_3_;
            const scalar E_4_;

        // Static constants

            //- Maximum number of iterations
            static label maxIters_;

            //- Tolerance used in inner iterations
            static scalar tolerance_;
            
        // Addressing

            //- Cs positions in global/local lists
            label CsLocalId_;

            //- O2 position in global list
            label O2GlobalId_;
            
            label CH4GlobalId_;
            
            label H2OGlobalId_;

            //- CO2 positions in global list
            label CO2GlobalId_;

            //- CO positions in global list
            label COGlobalId_;
            
            label H2GlobalId_;
            


        // Local copies of thermo properties

            //- Molecular weight of C [kg/kmol]
            scalar WC_;
            
            scalar WH2O_;
            
            scalar WCO2_;
            
            scalar WH2_;

            //- Molecular weight of O2 [kg/kmol]
            scalar WO2_;

            //- Formation enthalpy for CO2 [J/kg]
            scalar HcCO2_;

            //- Formation enthalpy for CO [J/kg]
            scalar HcCO_;


public:

    //- Runtime type information
    TypeName("COxidationKineticDiffusionLimitedRateWithAshSteamMurphyShaddix");


    // Constructors

        //- Construct from dictionary
        COxidationKineticDiffusionLimitedRateWithAshSteamMurphyShaddix
        (
            const dictionary& dict,
            CloudType& owner
        );

        //- Construct copy
        COxidationKineticDiffusionLimitedRateWithAshSteamMurphyShaddix
        (
            const COxidationKineticDiffusionLimitedRateWithAshSteamMurphyShaddix<CloudType>& srm
        );

        //- Construct and return a clone
        virtual autoPtr<CharOxidizationModel<CloudType>> clone() const
        {
            return autoPtr<CharOxidizationModel<CloudType>>
            (
                new COxidationKineticDiffusionLimitedRateWithAshSteamMurphyShaddix<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~COxidationKineticDiffusionLimitedRateWithAshSteamMurphyShaddix();


    // Member Functions

        //- Update surface reactions
        virtual scalar calculate
        (
            const scalar dt,
            const label celli,
            const scalar d, //outer layer radius
            const scalar T,
            const scalar Tc,
            const scalar pc, 
            const scalar rhoc,
            const scalar mass,
            const scalarField& YGas,
            const scalarField& YLiquid,
            const scalarField& YSolid,
            const scalarField& YMixture,
            const scalar N,  
            scalarField& dMassGas,
            scalarField& dMassLiquid,
            scalarField& dMassSolid,
            scalarField& dMassSRCarrier,
            const scalar di, //inner layer radius
            const scalar muc,//here used for viscosity
            const scalar Rec, //here used for Re number
            const scalar Xi,
            scalarField& dHeatSRCarrier,
            const label particleShape,
            const scalarField& epsilon, 
            const scalar deq
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "COxidationKineticDiffusionLimitedRateWithAshSteamMurphyShaddix.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
