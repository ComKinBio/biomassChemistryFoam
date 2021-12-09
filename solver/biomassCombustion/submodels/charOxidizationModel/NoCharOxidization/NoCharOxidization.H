/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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
    Foam::NoCharOxidization

Description
    Dummy surface reaction model for 'none'

\*---------------------------------------------------------------------------*/

#ifndef NoCharOxidization_H
#define NoCharOxidization_H

#include "CharOxidizationModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                     Class NoCharOxidization Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class NoCharOxidization
:
    public CharOxidizationModel<CloudType>
{
public:

    //- Runtime type information
    TypeName("none");


    // Constructors

        //- Construct from dictionary
        NoCharOxidization(const dictionary& dict, CloudType& owner);

        //- Construct copy
        NoCharOxidization(const NoCharOxidization<CloudType>& srm);

        //- Construct and return a clone
        virtual autoPtr<CharOxidizationModel<CloudType>> clone() const
        {
            return autoPtr<CharOxidizationModel<CloudType>>
            (
                new NoCharOxidization<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~NoCharOxidization();


    // Member Functions

        //- Flag to indicate whether model activates devolatisation model
        virtual bool active() const;

        //- Update surface reactions
        virtual scalar calculate
        (
            const scalar dt,
            const label celli,
            const scalar d,
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
            const scalar di, 
            const scalar muc,
            const scalar Rec,
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
    #include "NoCharOxidization.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
