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
    Foam::PyrolysisModel

Description
    Templated pyrolysis model class

SourceFiles
    PyrolysisModel.C
    PyrolysisModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef PyrolysisModel_H
#define PyrolysisModel_H

#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "CloudSubModelBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class PyrolysisModel Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class PyrolysisModel
:
    public CloudSubModelBase<CloudType>
{
protected:

    // Protected data

        //- Mass of lagrangian phase converted
        scalar dMass_;


public:

    //- Runtime type information
    TypeName("pyrolysisModel");

    //- Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        PyrolysisModel,
        dictionary,
        (
            const dictionary& dict,
            CloudType& owner
        ),
        (dict, owner)
    );


    // Constructors

        //- Construct null from owner
        PyrolysisModel(CloudType& owner);

        //- Construct from dictionary
        PyrolysisModel
        (
            const dictionary& dict,
            CloudType& owner,
            const word& type
        );

        //- Construct copy
        PyrolysisModel(const PyrolysisModel<CloudType>& dm);

        //- Construct and return a clone
        virtual autoPtr<PyrolysisModel<CloudType>> clone() const = 0;


    //- Destructor
    virtual ~PyrolysisModel();


    //- Selector
    static autoPtr<PyrolysisModel<CloudType>> New
    (
        const dictionary& dict,
        CloudType& owner
    );


    // Member Functions

        //- Update model
        virtual void calculate
        (
            const scalar dt,
            const scalar age,
            const scalar mass0,
            const scalar mass,
            const scalar T,
            const scalarField& YGasEff,
            const scalarField& YLiquidEff,
            const scalarField& YSolidEff,
            label& canCombust,
            scalarField& dMassDV,
            scalarField& dMassSOLID
        ) const = 0;

        //- Add to pyrolysis mass
        void addToPyrolysisMass(const scalar dMass);


        // I-O

            //- Write injection info to stream
            virtual void info(Ostream& os);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makePyrolysisModel(CloudType)                                          \
                                                                               \
    typedef Foam::CloudType::reactingMultiphaseIBMCloudType                    \
        reactingMultiphaseIBMCloudType;                                        \
    defineNamedTemplateTypeNameAndDebug                                        \
    (                                                                          \
        Foam::PyrolysisModel<reactingMultiphaseIBMCloudType>,                  \
        0                                                                      \
    );                                                                         \
    namespace Foam                                                             \
    {                                                                          \
        defineTemplateRunTimeSelectionTable                                    \
        (                                                                      \
            PyrolysisModel<reactingMultiphaseIBMCloudType>,                    \
            dictionary                                                         \
        );                                                                     \
    }


#define makePyrolysisModelType(SS, CloudType)                                  \
                                                                               \
    typedef Foam::CloudType::reactingMultiphaseIBMCloudType                    \
        reactingMultiphaseIBMCloudType;                                        \
    defineNamedTemplateTypeNameAndDebug                                        \
        (Foam::SS<reactingMultiphaseIBMCloudType>, 0);                         \
                                                                               \
    Foam::PyrolysisModel<reactingMultiphaseIBMCloudType>::                     \
        adddictionaryConstructorToTable                                        \
        <Foam::SS<reactingMultiphaseIBMCloudType>>                             \
        add##SS##CloudType##reactingMultiphaseIBMCloudType##ConstructorToTable_;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "PyrolysisModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
