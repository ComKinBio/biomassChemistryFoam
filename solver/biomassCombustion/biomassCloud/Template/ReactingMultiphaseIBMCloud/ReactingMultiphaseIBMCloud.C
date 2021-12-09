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

\*---------------------------------------------------------------------------*/

#include "ReactingMultiphaseIBMCloud.H"

#include "PyrolysisModel.H"
#include "CharOxidizationModel.H"

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

template<class CloudType>
void Foam::ReactingMultiphaseIBMCloud<CloudType>::setModels()
{
    pyrolysisModel_.reset
    (
        PyrolysisModel<ReactingMultiphaseIBMCloud<CloudType>>::New
        (
            this->subModelProperties(),
            *this
        ).ptr()
    );

    charOxidizationModel_.reset
    (
        CharOxidizationModel<ReactingMultiphaseIBMCloud<CloudType>>::New
        (
            this->subModelProperties(),
            *this
        ).ptr()
    );
}


template<class CloudType>
void Foam::ReactingMultiphaseIBMCloud<CloudType>::cloudReset
(
    ReactingMultiphaseIBMCloud<CloudType>& c
)
{
    CloudType::cloudReset(c);

    pyrolysisModel_.reset(c.pyrolysisModel_.ptr());
    charOxidizationModel_.reset(c.charOxidizationModel_.ptr());
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::ReactingMultiphaseIBMCloud<CloudType>::ReactingMultiphaseIBMCloud
(
    const word& cloudName,
    const volScalarField& rho,
    const volVectorField& U,
    const dimensionedVector& g,
    const SLGThermo& thermo,
    bool readFields
)
:
    CloudType(cloudName, rho, U, g, thermo, false),
    reactingMultiphaseIBMCloud(),
    cloudCopyPtr_(nullptr),
    constProps_(this->particleProperties()),
    pyrolysisModel_(nullptr),
    charOxidizationModel_(nullptr)
{
    if (this->solution().active())
    {
        setModels();

        if (readFields)
        {
            parcelType::readFields(*this, this->composition());
            this->deleteLostParticles();
        }
    }

    if (this->solution().resetSourcesOnStartup())
    {
        resetSourceTerms();
    }
}


template<class CloudType>
Foam::ReactingMultiphaseIBMCloud<CloudType>::ReactingMultiphaseIBMCloud
(
    ReactingMultiphaseIBMCloud<CloudType>& c,
    const word& name
)
:
    CloudType(c, name),
    reactingMultiphaseIBMCloud(),
    cloudCopyPtr_(nullptr),
    constProps_(c.constProps_),
    pyrolysisModel_(c.pyrolysisModel_->clone()),
    charOxidizationModel_(c.charOxidizationModel_->clone())
{}


template<class CloudType>
Foam::ReactingMultiphaseIBMCloud<CloudType>::ReactingMultiphaseIBMCloud
(
    const fvMesh& mesh,
    const word& name,
    const ReactingMultiphaseIBMCloud<CloudType>& c
)
:
    CloudType(mesh, name, c),
    reactingMultiphaseIBMCloud(),
    cloudCopyPtr_(nullptr),
    constProps_(),
    pyrolysisModel_(nullptr),
    charOxidizationModel_(nullptr)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::ReactingMultiphaseIBMCloud<CloudType>::~ReactingMultiphaseIBMCloud()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::ReactingMultiphaseIBMCloud<CloudType>::setParcelThermoProperties
(
    parcelType& parcel,
    const scalar lagrangianDt
)
{
    CloudType::setParcelThermoProperties(parcel, lagrangianDt);
    
    //initialize layer properties from constant properties. Particle contains all layers!
    parcel.Tb0() = constProps_.Tb00();
    parcel.Tb1() = constProps_.Tb10();
    parcel.Tb2() = constProps_.Tb20();
    parcel.Tb3() = constProps_.Tb30();   
    parcel.rb0() = constProps_.rb00();
    parcel.rb1() = constProps_.rb10();
    parcel.rb2() = constProps_.rb20();
    parcel.rb3() = constProps_.rb30();   
    parcel.Tp0() = constProps_.Tp00();
    parcel.Tp1() = constProps_.Tp10();
    parcel.Tp2() = constProps_.Tp20();
    parcel.Tp3() = constProps_.Tp30();
    if(constProps_.parcelShape() == 1)
    {
        parcel.mp0() = constProps_.rho0()*constant::mathematical::pi*(2.0*pow3(parcel.rb0())-sqr(parcel.rb0())*constProps_.xi0()),
        parcel.mp1() = constProps_.rho0()*constant::mathematical::pi*(2.0*pow3(parcel.rb1()) - 2.0*pow3(parcel.rb0()) - (sqr(parcel.rb1()) - sqr(parcel.rb0()))*constProps_.xi0());
        parcel.mp2() = constProps_.rho0()*constant::mathematical::pi*(2.0*pow3(parcel.rb2()) - 2.0*pow3(parcel.rb1()) - (sqr(parcel.rb2()) - sqr(parcel.rb1()))*constProps_.xi0());
        parcel.mp3() = constProps_.rho0()*constant::mathematical::pi*(2.0*pow3(parcel.rb3()) - 2.0*pow3(parcel.rb2()) - (sqr(parcel.rb3()) - sqr(parcel.rb2()))*constProps_.xi0());
    }
    else
    {
        parcel.mp0() = constProps_.rho0()*(4.0/3.0)*constant::mathematical::pi*pow3(parcel.rb0());
        parcel.mp1() = constProps_.rho0()*(4.0/3.0)*constant::mathematical::pi*(pow3(parcel.rb1())-pow3(parcel.rb0()));
        parcel.mp2() = constProps_.rho0()*(4.0/3.0)*constant::mathematical::pi*(pow3(parcel.rb2())-pow3(parcel.rb1()));
        parcel.mp3() = constProps_.rho0()*(4.0/3.0)*constant::mathematical::pi*(pow3(parcel.rb3())-pow3(parcel.rb2()));
    }
    parcel.flagBoiling() = 1;
    parcel.flagDevo() = 1;

}


template<class CloudType>
void Foam::ReactingMultiphaseIBMCloud<CloudType>::checkParcelProperties
(
    parcelType& parcel,
    const scalar lagrangianDt,
    const bool fullyDescribed
)
{
    CloudType::checkParcelProperties(parcel, lagrangianDt, fullyDescribed);
}


template<class CloudType>
void Foam::ReactingMultiphaseIBMCloud<CloudType>::storeState()
{
    cloudCopyPtr_.reset
    (
        static_cast<ReactingMultiphaseIBMCloud<CloudType>*>
        (
            clone(this->name() + "Copy").ptr()
        )
    );
}


template<class CloudType>
void Foam::ReactingMultiphaseIBMCloud<CloudType>::restoreState()
{
    cloudReset(cloudCopyPtr_());
    cloudCopyPtr_.clear();
}


template<class CloudType>
void Foam::ReactingMultiphaseIBMCloud<CloudType>::resetSourceTerms()
{
    CloudType::resetSourceTerms();
}


template<class CloudType>
void Foam::ReactingMultiphaseIBMCloud<CloudType>::evolve()
{
    if (this->solution().canEvolve())
    {
        typename parcelType::trackingData td(*this);

        this->solve(*this, td);
    }
}


template<class CloudType>
void Foam::ReactingMultiphaseIBMCloud<CloudType>::autoMap
(
    const mapPolyMesh& mapper
)
{
    Cloud<parcelType>::autoMap(mapper);

    this->updateMesh();
}


template<class CloudType>
void Foam::ReactingMultiphaseIBMCloud<CloudType>::info()
{
    CloudType::info();

    this->pyrolysis().info(Info);
    this->charOxidization().info(Info);
}


template<class CloudType>
void Foam::ReactingMultiphaseIBMCloud<CloudType>::writeFields() const
{
    if (this->compositionModel_.valid())
    {
        CloudType::particleType::writeFields(*this, this->composition());
    }
}


// ************************************************************************* //
