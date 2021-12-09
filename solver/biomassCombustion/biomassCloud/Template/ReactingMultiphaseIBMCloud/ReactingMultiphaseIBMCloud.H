/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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
    Foam::ReactingMultiphaseIBMCloud

Description
    Templated base class for multiphase reacting cloud

    - Adds to reacting cloud
      - multiphase composition
      - devolatilisatsion
      - surface reactions

SourceFiles
    ReactingMultiphaseIBMCloudI.H
    ReactingMultiphaseIBMCloud.C

\*---------------------------------------------------------------------------*/

#ifndef ReactingMultiphaseIBMCloud_H
#define ReactingMultiphaseIBMCloud_H

#include "reactingMultiphaseIBMCloud.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes

template<class CloudType>
class PyrolysisModel;

template<class CloudType>
class CharOxidizationModel;

/*---------------------------------------------------------------------------*\
                   Class ReactingMultiphaseIBMCloud Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class ReactingMultiphaseIBMCloud
:
    public CloudType,
    public reactingMultiphaseIBMCloud
{
public:

    // Public Typedefs

        //- Type of cloud this cloud was instantiated for
        typedef CloudType cloudType;

        //- Type of parcel the cloud was instantiated for
        typedef typename CloudType::particleType parcelType;

        //- Convenience typedef for this cloud type
        typedef ReactingMultiphaseIBMCloud<CloudType> reactingMultiphaseIBMCloudType;


private:

    // Private Data

        //- Cloud copy pointer
        autoPtr<ReactingMultiphaseIBMCloud<CloudType>> cloudCopyPtr_;


protected:

    // Protected data

        //- Parcel constant properties
        typename parcelType::constantProperties constProps_;


        // References to the cloud sub-models

            //- Pyrolysis model
            autoPtr
            <
                PyrolysisModel<ReactingMultiphaseIBMCloud<CloudType>>
            >
            pyrolysisModel_;

            //- Char OxidizationModel model
            autoPtr
            <
                CharOxidizationModel<ReactingMultiphaseIBMCloud<CloudType>>
            >
            charOxidizationModel_;


        // Check

//             //- Total mass transferred to continuous phase via pyrolysis
//             scalar dMassPyrolysis_;
// 
//             //- Total mass transferred to continuous phase via surface
//             //  reactions
//             scalar dMassCharOxidization_;


    // Protected Member Functions

        // Initialisation

            //- Set cloud sub-models
            void setModels();


        // Cloud evolution functions

            //- Reset state of cloud
            void cloudReset(ReactingMultiphaseIBMCloud<CloudType>& c);


public:

    // Constructors

        //- Construct given carrier gas fields
        ReactingMultiphaseIBMCloud
        (
            const word& cloudName,
            const volScalarField& rho,
            const volVectorField& U,
            const dimensionedVector& g,
            const SLGThermo& thermo,
            bool readFields = true
        );


        //- Copy constructor with new name
        ReactingMultiphaseIBMCloud
        (
            ReactingMultiphaseIBMCloud<CloudType>& c,
            const word& name
        );

        //- Copy constructor with new name - creates bare cloud
        ReactingMultiphaseIBMCloud
        (
            const fvMesh& mesh,
            const word& name,
            const ReactingMultiphaseIBMCloud<CloudType>& c
        );

        //- Disallow default bitwise copy construction
        ReactingMultiphaseIBMCloud(const ReactingMultiphaseIBMCloud&) = delete;

        //- Construct and return clone based on (this) with new name
        virtual autoPtr<Cloud<parcelType>> clone(const word& name)
        {
            return autoPtr<Cloud<parcelType>>
            (
                new ReactingMultiphaseIBMCloud(*this, name)
            );
        }

        //- Construct and return bare clone based on (this) with new name
        virtual autoPtr<Cloud<parcelType>> cloneBare(const word& name) const
        {
            return autoPtr<Cloud<parcelType>>
            (
                new ReactingMultiphaseIBMCloud(this->mesh(), name, *this)
            );
        }


    //- Destructor
    virtual ~ReactingMultiphaseIBMCloud();


    // Member Functions

        // Access

            //- Return a reference to the cloud copy
            inline const ReactingMultiphaseIBMCloud& cloudCopy() const;

            //- Return the constant properties
            inline const typename parcelType::constantProperties&
                constProps() const;

            //- Return access to the constant properties
            inline typename parcelType::constantProperties& constProps();


            // Sub-models

                //- Return const access to pyrolysis model
                inline const PyrolysisModel
                <
                    ReactingMultiphaseIBMCloud<CloudType>
                >&
                pyrolysis() const;

                //- Return reference to pyrolysis model
                inline PyrolysisModel
                <
                    ReactingMultiphaseIBMCloud<CloudType>
                >&
                pyrolysis();

                //- Return const access to reacting surface reaction model
                inline const CharOxidizationModel
                <
                    ReactingMultiphaseIBMCloud<CloudType>
                >&
                charOxidization() const;

                //- Return reference to reacting surface reaction model
                inline CharOxidizationModel
                <
                    ReactingMultiphaseIBMCloud<CloudType>
                >&
                charOxidization();


        // Cloud evolution functions

            //- Set parcel thermo properties
            void setParcelThermoProperties
            (
                parcelType& parcel,
                const scalar lagrangianDt
            );

            //- Check parcel properties
            void checkParcelProperties
            (
                parcelType& parcel,
                const scalar lagrangianDt,
                const bool fullyDescribed
            );

            //- Store the current cloud state
            void storeState();

            //- Reset the current cloud to the previously stored state
            void restoreState();

            //- Reset the cloud source terms
            void resetSourceTerms();

            //- Evolve the cloud
            void evolve();


        // Mapping

            //- Remap the cells of particles corresponding to the
            //  mesh topology change with a default tracking data object
            virtual void autoMap(const mapPolyMesh&);


        // I-O

            //- Print cloud information
            void info();

            //- Write the field data for the cloud
            virtual void writeFields() const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const ReactingMultiphaseIBMCloud&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "ReactingMultiphaseIBMCloudI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "ReactingMultiphaseIBMCloud.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
