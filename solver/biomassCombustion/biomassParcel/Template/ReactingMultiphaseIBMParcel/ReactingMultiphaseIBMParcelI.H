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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ParcelType>
inline Foam::ReactingMultiphaseIBMParcel<ParcelType>::constantProperties::
constantProperties()
:
    ParcelType::constantProperties(),
    parcelShape_(this->dict_, 1),
    xi0_(this->dict_, 1),
    shrinkageFactorAlpha_(this->dict_, 0.10),
    shrinkageFactorBeta_(this->dict_, 0.39),
    shrinkageFactorGamma_(this->dict_, 0.95),
    Tb00_(this->dict_, 0.0),
    Tb10_(this->dict_, 0.0),
    Tb20_(this->dict_, 0.0),
    Tb30_(this->dict_, 0.0),
    rb00_(this->dict_, 0.0),
    rb10_(this->dict_, 0.0),
    rb20_(this->dict_, 0.0),    
    rb30_(this->dict_, 0.0),
    Tp00_(this->dict_, 0.0),
    Tp10_(this->dict_, 0.0),     
    Tp20_(this->dict_, 0.0),  
    Tp30_(this->dict_, 0.0),
    deltaTime_(this->dict_, 0.0),
    hCoeff_(this->dict_, 0.0)
{}


template<class ParcelType>
inline Foam::ReactingMultiphaseIBMParcel<ParcelType>::constantProperties::
constantProperties
(
    const constantProperties& cp
)
:
    ParcelType::constantProperties(cp),
    parcelShape_(cp.parcelShape_),
    xi0_(cp.xi0_),
    shrinkageFactorAlpha_(cp.shrinkageFactorAlpha_),
    shrinkageFactorBeta_(cp.shrinkageFactorBeta_),
    shrinkageFactorGamma_(cp.shrinkageFactorGamma_),
    Tb00_(cp.Tb00_),
    Tb10_(cp.Tb10_),
    Tb20_(cp.Tb20_),
    Tb30_(cp.Tb30_),
    rb00_(cp.rb00_),
    rb10_(cp.rb10_),
    rb20_(cp.rb20_),
    rb30_(cp.rb30_),
    Tp00_(cp.Tp00_),
    Tp10_(cp.Tp10_),
    Tp20_(cp.Tp20_),
    Tp30_(cp.Tp30_),
    deltaTime_(cp.deltaTime_),
    hCoeff_(cp.hCoeff_)  
{}


template<class ParcelType>
inline Foam::ReactingMultiphaseIBMParcel<ParcelType>::constantProperties::
constantProperties
(
    const dictionary& parentDict
)
:
    ParcelType::constantProperties(parentDict),
    parcelShape_(this->dict_, word("parcelShape")),
    xi0_(this->dict_, word("xi0")),
    shrinkageFactorAlpha_(this->dict_, word("Srinkalpha")),
    shrinkageFactorBeta_(this->dict_, word("Srinkbeta")),
    shrinkageFactorGamma_(this->dict_, word("Srinkgamma")),
    Tb00_(this->dict_, "Tb00"),
    Tb10_(this->dict_, "Tb10"),
    Tb20_(this->dict_, "Tb20"),
    Tb30_(this->dict_, "Tb30"),
    rb00_(this->dict_, "rb00"),
    rb10_(this->dict_, "rb10"),
    rb20_(this->dict_, "rb20"),    
    rb30_(this->dict_, "rb30"),
    Tp00_(this->dict_, "Tp00"),
    Tp10_(this->dict_, "Tp10"),     
    Tp20_(this->dict_, "Tp20"),  
    Tp30_(this->dict_, "Tp30"),
    deltaTime_(this->dict_, "deltaTime"),
    hCoeff_(this->dict_, "hCoeff")
{}


template<class ParcelType>
inline Foam::ReactingMultiphaseIBMParcel<ParcelType>::ReactingMultiphaseIBMParcel
(
    const polyMesh& mesh,
    const barycentric& coordinates,
    const label celli,
    const label tetFacei,
    const label tetPti
)
:
    ParcelType(mesh, coordinates, celli, tetFacei, tetPti),
    Tb0_(0.0),
    Tb1_(0.0),
    Tb2_(0.0),
    Tb3_(0.0),
    rb0_(0.0),
    rb1_(0.0),
    rb2_(0.0),    
    rb3_(0.0),
    Tp0_(0.0),
    Tp1_(0.0),     
    Tp2_(0.0),  
    Tp3_(0.0),
    mp0_(0.0),
    mp1_(0.0),
    mp2_(0.0),
    mp3_(0.0),
    rDry_(0.0),
    rDevo_(0.0),
    rChar_(0.0),
    rComb_(0.0),
    ash_inchar_t_(0.0),
    QDry_(0.0),
    QComb_(0.0),
    flagBoiling_(0),
    flagDevo_(0),
    cumTime_(0.0)
{}


template<class ParcelType>
inline Foam::ReactingMultiphaseIBMParcel<ParcelType>::ReactingMultiphaseIBMParcel
(
    const polyMesh& mesh,
    const vector& position,
    const label celli
)
:
    ParcelType(mesh, position, celli),
    Tb0_(0.0),
    Tb1_(0.0),
    Tb2_(0.0),
    Tb3_(0.0),
    rb0_(0.0),
    rb1_(0.0),
    rb2_(0.0),    
    rb3_(0.0),
    Tp0_(0.0),
    Tp1_(0.0),     
    Tp2_(0.0),  
    Tp3_(0.0),
    mp0_(0.0),
    mp1_(0.0),
    mp2_(0.0),
    mp3_(0.0),
    rDry_(0.0),
    rDevo_(0.0),
    rChar_(0.0),
    rComb_(0.0),
    ash_inchar_t_(0.0),
    QDry_(0.0),
    QComb_(0.0),
    flagBoiling_(0),
    flagDevo_(0),
    cumTime_(0.0)
{}


template<class ParcelType>
inline Foam::ReactingMultiphaseIBMParcel<ParcelType>::ReactingMultiphaseIBMParcel
(
    const polyMesh& mesh,
    const barycentric& coordinates,
    const label celli,
    const label tetFacei,
    const label tetPti,
    const label typeId,
    const scalar nParticle0,
    const scalar d0,
    const scalar dTarget0,
    const vector& U0,
    const vector& f0,
    const vector& angularMomentum0,
    const vector& torque0,
    const scalarField& Y0,
    const scalarField& YGas0,
    const scalarField& YLiquid0,
    const scalarField& YSolid0,
    const constantProperties& constProps
)
:
    ParcelType
    (
        mesh,
        coordinates,
        celli,
        tetFacei,
        tetPti,
        typeId,
        nParticle0,
        d0,
        dTarget0,
        U0,
        f0,
        angularMomentum0,
        torque0,
        Y0,
        constProps
    ),
    Tb0_(constProps.Tb00()),
    Tb1_(constProps.Tb10()),
    Tb2_(constProps.Tb20()),
    Tb3_(constProps.Tb30()),
    rb0_(constProps.rb00()),
    rb1_(constProps.rb10()),
    rb2_(constProps.rb20()),    
    rb3_(constProps.rb30()),
    Tp0_(constProps.Tp00()),
    Tp1_(constProps.Tp10()),     
    Tp2_(constProps.Tp20()),  
    Tp3_(constProps.Tp30()),
    mp0_(constProps.rho0()*(4.0/3.0)*constant::mathematical::pi*pow3(rb0_)),
    mp1_(constProps.rho0()*(4.0/3.0)*constant::mathematical::pi*(pow3(rb1_)-pow3(rb0_))),
    mp2_(constProps.rho0()*(4.0/3.0)*constant::mathematical::pi*(pow3(rb2_)-pow3(rb1_))),
    mp3_(constProps.rho0()*(4.0/3.0)*constant::mathematical::pi*(pow3(rb3_)-pow3(rb2_))),
    rDry_(0.0),
    rDevo_(0.0),
    rChar_(0.0),
    rComb_(0.0),
    ash_inchar_t_(0.0),
    QDry_(0.0),
    QComb_(0.0),
    flagBoiling_(1),
    flagDevo_(1),
    cumTime_(0.0)
    
{
    if(constProps.parcelShape() == 1)
    {
        mp0_ = constProps.rho0()*constant::mathematical::pi*(2*pow3(rb0_)-sqr(rb0_)*constProps.xi0()),
        mp1_ = constProps.rho0()*constant::mathematical::pi*(2.0*pow3(rb1_) - 2.0*pow3(rb0_) - (sqr(rb1_) - sqr(rb0_))*constProps.xi0());
        mp2_ = constProps.rho0()*constant::mathematical::pi*(2.0*pow3(rb2_) - 2.0*pow3(rb1_) - (sqr(rb2_) - sqr(rb1_))*constProps.xi0());
        mp3_ = constProps.rho0()*constant::mathematical::pi*(2.0*pow3(rb3_) - 2.0*pow3(rb2_) - (sqr(rb3_) - sqr(rb2_))*constProps.xi0());
    }
}


// * * * * * * * * * constantProperties Member Functions * * * * * * * * * * //


template<class ParcelType>
inline Foam::scalar
Foam::ReactingMultiphaseIBMParcel<ParcelType>::constantProperties::Tb00() const
{
    return Tb00_.value();
}

template<class ParcelType>
inline Foam::scalar
Foam::ReactingMultiphaseIBMParcel<ParcelType>::constantProperties::Tb10() const
{
    return Tb10_.value();
}

template<class ParcelType>
inline Foam::scalar
Foam::ReactingMultiphaseIBMParcel<ParcelType>::constantProperties::Tb20() const
{
    return Tb20_.value();
}

template<class ParcelType>
inline Foam::scalar
Foam::ReactingMultiphaseIBMParcel<ParcelType>::constantProperties::Tb30() const
{
    return Tb30_.value();
}

template<class ParcelType>
inline Foam::scalar
Foam::ReactingMultiphaseIBMParcel<ParcelType>::constantProperties::rb00() const
{
    return rb00_.value();
}

template<class ParcelType>
inline Foam::scalar
Foam::ReactingMultiphaseIBMParcel<ParcelType>::constantProperties::rb10() const
{
    return rb10_.value();
}

template<class ParcelType>
inline Foam::scalar
Foam::ReactingMultiphaseIBMParcel<ParcelType>::constantProperties::rb20() const
{
    return rb20_.value();
}

template<class ParcelType>
inline Foam::scalar
Foam::ReactingMultiphaseIBMParcel<ParcelType>::constantProperties::rb30() const
{
    return rb30_.value();
}

template<class ParcelType>
inline Foam::scalar
Foam::ReactingMultiphaseIBMParcel<ParcelType>::constantProperties::Tp00() const
{
    return Tp00_.value();
}

template<class ParcelType>
inline Foam::scalar
Foam::ReactingMultiphaseIBMParcel<ParcelType>::constantProperties::Tp10() const
{
    return Tp10_.value();
}

template<class ParcelType>
inline Foam::scalar
Foam::ReactingMultiphaseIBMParcel<ParcelType>::constantProperties::Tp20() const
{
    return Tp20_.value();
}

template<class ParcelType>
inline Foam::scalar
Foam::ReactingMultiphaseIBMParcel<ParcelType>::constantProperties::Tp30() const
{
    return Tp30_.value();
}

template<class ParcelType>
inline Foam::scalar
Foam::ReactingMultiphaseIBMParcel<ParcelType>::constantProperties::deltaTime() const
{
    return deltaTime_.value();
}

template<class ParcelType>
inline Foam::scalar
Foam::ReactingMultiphaseIBMParcel<ParcelType>::constantProperties::hCoeff() const
{
    return hCoeff_.value();
}

template<class ParcelType>
inline Foam::label
Foam::ReactingMultiphaseIBMParcel<ParcelType>::constantProperties::parcelShape() const
{
    return parcelShape_.value();
}

template<class ParcelType>
inline Foam::scalar
Foam::ReactingMultiphaseIBMParcel<ParcelType>::constantProperties::xi0() const
{
    return xi0_.value();
}

template<class ParcelType>
inline Foam::scalar
Foam::ReactingMultiphaseIBMParcel<ParcelType>::constantProperties::shrinkageFactorAlpha() const
{
    return shrinkageFactorAlpha_.value();
}

template<class ParcelType>
inline Foam::scalar
Foam::ReactingMultiphaseIBMParcel<ParcelType>::constantProperties::shrinkageFactorBeta() const
{
    return shrinkageFactorBeta_.value();
}

template<class ParcelType>
inline Foam::scalar
Foam::ReactingMultiphaseIBMParcel<ParcelType>::constantProperties::shrinkageFactorGamma() const
{
    return shrinkageFactorGamma_.value();
}

// * * * * * * * * * * ThermoParcel Member Functions * * * * * * * * * * * * //
template<class ParcelType>
inline Foam::scalar
Foam::ReactingMultiphaseIBMParcel<ParcelType>::Tb0() const
{
    return Tb0_;
}

template<class ParcelType>
inline Foam::scalar
Foam::ReactingMultiphaseIBMParcel<ParcelType>::Tb1() const
{
    return Tb1_;
}

template<class ParcelType>
inline Foam::scalar
Foam::ReactingMultiphaseIBMParcel<ParcelType>::Tb2() const
{
    return Tb2_;
}

template<class ParcelType>
inline Foam::scalar
Foam::ReactingMultiphaseIBMParcel<ParcelType>::Tb3() const
{
    return Tb3_;
}

template<class ParcelType>
inline Foam::scalar
Foam::ReactingMultiphaseIBMParcel<ParcelType>::rb0() const
{
    return rb0_;
}

template<class ParcelType>
inline Foam::scalar
Foam::ReactingMultiphaseIBMParcel<ParcelType>::rb1() const
{
    return rb1_;
}

template<class ParcelType>
inline Foam::scalar
Foam::ReactingMultiphaseIBMParcel<ParcelType>::rb2() const
{
    return rb2_;
}

template<class ParcelType>
inline Foam::scalar
Foam::ReactingMultiphaseIBMParcel<ParcelType>::rb3() const
{
    return rb3_;
}

template<class ParcelType>
inline Foam::scalar
Foam::ReactingMultiphaseIBMParcel<ParcelType>::Tp0() const
{
    return Tp0_;
}

template<class ParcelType>
inline Foam::scalar
Foam::ReactingMultiphaseIBMParcel<ParcelType>::Tp1() const
{
    return Tp1_;
}

template<class ParcelType>
inline Foam::scalar
Foam::ReactingMultiphaseIBMParcel<ParcelType>::Tp2() const
{
    return Tp2_;
}

template<class ParcelType>
inline Foam::scalar
Foam::ReactingMultiphaseIBMParcel<ParcelType>::Tp3() const
{
    return Tp3_;
}

template<class ParcelType>
inline Foam::scalar
Foam::ReactingMultiphaseIBMParcel<ParcelType>::mp0() const
{
    return mp0_;
}

template<class ParcelType>
inline Foam::scalar
Foam::ReactingMultiphaseIBMParcel<ParcelType>::mp1() const
{
    return mp1_;
}

template<class ParcelType>
inline Foam::scalar
Foam::ReactingMultiphaseIBMParcel<ParcelType>::mp2() const
{
    return mp2_;
}

template<class ParcelType>
inline Foam::scalar
Foam::ReactingMultiphaseIBMParcel<ParcelType>::mp3() const
{
    return mp3_;
}

template<class ParcelType>
inline Foam::scalar
Foam::ReactingMultiphaseIBMParcel<ParcelType>::rDry() const
{
    return rDry_;
}

template<class ParcelType>
inline Foam::scalar
Foam::ReactingMultiphaseIBMParcel<ParcelType>::rDevo() const
{
    return rDevo_;
}

template<class ParcelType>
inline Foam::scalar
Foam::ReactingMultiphaseIBMParcel<ParcelType>::rChar() const
{
    return rChar_;
}

template<class ParcelType>
inline Foam::scalar
Foam::ReactingMultiphaseIBMParcel<ParcelType>::rComb() const
{
    return rComb_;
}

template<class ParcelType>
inline Foam::scalar
Foam::ReactingMultiphaseIBMParcel<ParcelType>::ash_inchar_t() const
{
    return ash_inchar_t_;
}

template<class ParcelType>
inline Foam::scalar
Foam::ReactingMultiphaseIBMParcel<ParcelType>::QDry() const
{
    return QDry_;
}

template<class ParcelType>
inline Foam::scalar
Foam::ReactingMultiphaseIBMParcel<ParcelType>::QComb() const
{
    return QComb_;
}

template<class ParcelType>
inline Foam::label
Foam::ReactingMultiphaseIBMParcel<ParcelType>::flagBoiling() const
{
    return flagBoiling_;
}

template<class ParcelType>
inline Foam::label
Foam::ReactingMultiphaseIBMParcel<ParcelType>::flagDevo() const
{
    return flagDevo_;
}

template<class ParcelType>
inline Foam::scalar
Foam::ReactingMultiphaseIBMParcel<ParcelType>::cumTime() const
{
    return cumTime_;
}

template<class ParcelType>
inline Foam::scalar&
Foam::ReactingMultiphaseIBMParcel<ParcelType>::Tb0()
{
    return Tb0_;
}

template<class ParcelType>
inline Foam::scalar&
Foam::ReactingMultiphaseIBMParcel<ParcelType>::Tb1()
{
    return Tb1_;
}

template<class ParcelType>
inline Foam::scalar&
Foam::ReactingMultiphaseIBMParcel<ParcelType>::Tb2()
{
    return Tb2_;
}

template<class ParcelType>
inline Foam::scalar&
Foam::ReactingMultiphaseIBMParcel<ParcelType>::Tb3()
{
    return Tb3_;
}

template<class ParcelType>
inline Foam::scalar&
Foam::ReactingMultiphaseIBMParcel<ParcelType>::rb0()
{
    return rb0_;
}

template<class ParcelType>
inline Foam::scalar&
Foam::ReactingMultiphaseIBMParcel<ParcelType>::rb1()
{
    return rb1_;
}

template<class ParcelType>
inline Foam::scalar&
Foam::ReactingMultiphaseIBMParcel<ParcelType>::rb2()
{
    return rb2_;
}

template<class ParcelType>
inline Foam::scalar&
Foam::ReactingMultiphaseIBMParcel<ParcelType>::rb3()
{
    return rb3_;
}

template<class ParcelType>
inline Foam::scalar&
Foam::ReactingMultiphaseIBMParcel<ParcelType>::Tp0()
{
    return Tp0_;
}

template<class ParcelType>
inline Foam::scalar&
Foam::ReactingMultiphaseIBMParcel<ParcelType>::Tp1()
{
    return Tp1_;
}

template<class ParcelType>
inline Foam::scalar&
Foam::ReactingMultiphaseIBMParcel<ParcelType>::Tp2()
{
    return Tp2_;
}

template<class ParcelType>
inline Foam::scalar&
Foam::ReactingMultiphaseIBMParcel<ParcelType>::Tp3()
{
    return Tp3_;
}

template<class ParcelType>
inline Foam::scalar&
Foam::ReactingMultiphaseIBMParcel<ParcelType>::mp0()
{
    return mp0_;
}

template<class ParcelType>
inline Foam::scalar&
Foam::ReactingMultiphaseIBMParcel<ParcelType>::mp1()
{
    return mp1_;
}

template<class ParcelType>
inline Foam::scalar&
Foam::ReactingMultiphaseIBMParcel<ParcelType>::mp2()
{
    return mp2_;
}

template<class ParcelType>
inline Foam::scalar&
Foam::ReactingMultiphaseIBMParcel<ParcelType>::mp3()
{
    return mp3_;
}

template<class ParcelType>
inline Foam::scalar&
Foam::ReactingMultiphaseIBMParcel<ParcelType>::rDry()
{
    return rDry_;
}

template<class ParcelType>
inline Foam::scalar&
Foam::ReactingMultiphaseIBMParcel<ParcelType>::rDevo()
{
    return rDevo_;
}

template<class ParcelType>
inline Foam::scalar&
Foam::ReactingMultiphaseIBMParcel<ParcelType>::rChar()
{
    return rChar_;
}

template<class ParcelType>
inline Foam::scalar&
Foam::ReactingMultiphaseIBMParcel<ParcelType>::rComb()
{
    return rComb_;
}

template<class ParcelType>
inline Foam::scalar&
Foam::ReactingMultiphaseIBMParcel<ParcelType>::ash_inchar_t()
{
    return ash_inchar_t_;
}

template<class ParcelType>
inline Foam::scalar&
Foam::ReactingMultiphaseIBMParcel<ParcelType>::QDry()
{
    return QDry_;
}

template<class ParcelType>
inline Foam::scalar&
Foam::ReactingMultiphaseIBMParcel<ParcelType>::QComb()
{
    return QComb_;
}

template<class ParcelType>
inline Foam::label&
Foam::ReactingMultiphaseIBMParcel<ParcelType>::flagBoiling()
{
    return flagBoiling_;
}

template<class ParcelType>
inline Foam::label&
Foam::ReactingMultiphaseIBMParcel<ParcelType>::flagDevo()
{
    return flagDevo_;
}

template<class ParcelType>
inline Foam::scalar&
Foam::ReactingMultiphaseIBMParcel<ParcelType>::cumTime()
{
    return cumTime_;
}

// ************************************************************************* //
