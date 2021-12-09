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
    Foam::CharOxidizationModel

Description
    Templated surface reaction model class

SourceFiles
    CharOxidizationModel.C
    CharOxidizationModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef CharOxidizationModel_H
#define CharOxidizationModel_H

#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "CloudSubModelBase.H"
#include "scalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class CharOxidizationModel Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class CharOxidizationModel
:
    public CloudSubModelBase<CloudType>
{
protected:

    // Protected data

        //- Mass of lagrangian phase converted
        scalar dMass_;


public:

    //-Runtime type information
    TypeName("charOxidizationModel");


    //- Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        CharOxidizationModel,
        dictionary,
        (
            const dictionary& dict,
            CloudType& cloud
        ),
        (dict, cloud)
    );


    // Constructors

        //- Construct null from owner
        CharOxidizationModel(CloudType& owner);

        //- Construct from dictionary
        CharOxidizationModel
        (
            const dictionary& dict,
            CloudType& cloud,
            const word& type
        );

        //- Construct copy
        CharOxidizationModel(const CharOxidizationModel<CloudType>& srm);

        //- Construct and return a clone
        virtual autoPtr<CharOxidizationModel<CloudType>> clone() const = 0;


    //- Destructor
    virtual ~CharOxidizationModel();


    //- Selector
    static autoPtr<CharOxidizationModel<CloudType>> New
    (
        const dictionary& dict,
        CloudType& cloud
    );


    // Member Functions

        //- Update surface reactions
        //  Returns the heat of reaction
        virtual scalar calculate
        (
            const scalar dt,
            const label celli,
            const scalar d, //outter layer radius
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
            const scalar Xi, //cylinder geometry coefficient
            scalarField& dHeatSRCarrier, //heat due to different reactions
            const label particleShape, //1 is cylinder otherwise sphere
            const scalarField& epsilon, // [0] bed voidage, [1] ash layer voidage  
            const scalar deq //equivalent diameter to a sphere for Sh number 
        ) const = 0;

        //- Add to devolatilisation mass
        void addToCharOxidizationMass(const scalar dMass);


        // I-O

            //- Write injection info to stream
            virtual void info(Ostream& os);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeCharOxidizationModel(CloudType)                                    \
                                                                               \
    typedef Foam::CloudType::reactingMultiphaseIBMCloudType                    \
        reactingMultiphaseIBMCloudType;                                        \
    defineNamedTemplateTypeNameAndDebug                                        \
    (                                                                          \
        Foam::CharOxidizationModel<reactingMultiphaseIBMCloudType>,            \
        0                                                                      \
    );                                                                         \
    namespace Foam                                                             \
    {                                                                          \
        defineTemplateRunTimeSelectionTable                                    \
        (                                                                      \
            CharOxidizationModel<reactingMultiphaseIBMCloudType>,              \
            dictionary                                                         \
        );                                                                     \
    }


#define makeCharOxidizationModelType(SS, CloudType)                            \
                                                                               \
    typedef Foam::CloudType::reactingMultiphaseIBMCloudType                    \
        reactingMultiphaseIBMCloudType;                                        \
    defineNamedTemplateTypeNameAndDebug                                        \
        (Foam::SS<reactingMultiphaseIBMCloudType>, 0);                         \
                                                                               \
    Foam::CharOxidizationModel<reactingMultiphaseIBMCloudType>::               \
        adddictionaryConstructorToTable                                        \
        <Foam::SS<reactingMultiphaseIBMCloudType>>                             \
        add##SS##CloudType##reactingMultiphaseIBMCloudType##ConstructorToTable_;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "CharOxidizationModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
