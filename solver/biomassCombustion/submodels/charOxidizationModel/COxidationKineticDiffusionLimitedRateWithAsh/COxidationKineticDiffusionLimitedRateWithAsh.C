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

\*---------------------------------------------------------------------------*/

#include "COxidationKineticDiffusionLimitedRateWithAsh.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::COxidationKineticDiffusionLimitedRateWithAsh<CloudType>::
COxidationKineticDiffusionLimitedRateWithAsh
(
    const dictionary& dict,
    CloudType& owner
)
:
    CharOxidizationModel<CloudType>(dict, owner, typeName),
    COmega_(readScalar(this->coeffDict().lookup("COmiga"))),
    ep3_(readScalar(this->coeffDict().lookup("ep3"))),
    C1_(readScalar(this->coeffDict().lookup("C1"))),
    C2_(readScalar(this->coeffDict().lookup("C2"))),
    E_(readScalar(this->coeffDict().lookup("E"))),
    CsLocalId_(-1),
    O2GlobalId_(owner.composition().carrierId("O2")),
    CO2GlobalId_(owner.composition().carrierId("CO2")),
    COGlobalId_(owner.composition().carrierId("CO")),
    WC_(0.0),
    WO2_(0.0),
    HcCO2_(0.0),
    HcCO_(0.0)
{
    // Determine Cs ids
    label idSolid = owner.composition().idSolid();
    CsLocalId_ = owner.composition().localId(idSolid, "C");

    // Set local copies of thermo properties
    WO2_ = owner.thermo().carrier().Wi(O2GlobalId_);
    const scalar WCO2 = owner.thermo().carrier().Wi(CO2GlobalId_);
    WC_ = WCO2 - WO2_;

    HcCO2_ = owner.thermo().carrier().Hc(CO2GlobalId_);
    HcCO_ = owner.thermo().carrier().Hc(COGlobalId_);

    const scalar YCloc = owner.composition().Y0(idSolid)[CsLocalId_];
    const scalar YSolidTot = owner.composition().YMixture0()[idSolid];
    Info<< "    C(s): particle mass fraction = " << YCloc*YSolidTot << endl;
}


template<class CloudType>
Foam::COxidationKineticDiffusionLimitedRateWithAsh<CloudType>::
COxidationKineticDiffusionLimitedRateWithAsh
(
    const COxidationKineticDiffusionLimitedRateWithAsh<CloudType>& srm
)
:
    CharOxidizationModel<CloudType>(srm),
    //Sb_(srm.Sb_),
    COmega_(srm.COmega_),
    ep3_(srm.ep3_),
    C1_(srm.C1_),
    C2_(srm.C2_),
    E_(srm.E_),
    CsLocalId_(srm.CsLocalId_),
    O2GlobalId_(srm.O2GlobalId_),
    CO2GlobalId_(srm.CO2GlobalId_),
    COGlobalId_(srm.COGlobalId_),
    WC_(srm.WC_),
    WO2_(srm.WO2_),
    HcCO2_(srm.HcCO2_),
    HcCO_(srm.HcCO_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::COxidationKineticDiffusionLimitedRateWithAsh<CloudType>::
~COxidationKineticDiffusionLimitedRateWithAsh()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
Foam::scalar Foam::COxidationKineticDiffusionLimitedRateWithAsh<CloudType>::calculate
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
) const
{
   // Fraction of remaining combustible material
    const label idSolid = CloudType::parcelType::SLD;
    const scalar fComb = YMixture[idSolid]*YSolid[CsLocalId_];

    // Surface combustion active combustible fraction is consumed
    if (fComb < SMALL)
    {
        return 0.0;
    }

    const SLGThermo& thermo = this->owner().thermo();

    // Local mass fraction of O2 in the carrier phase
    const scalar YO2 = thermo.carrier().Y(O2GlobalId_)[celli];

    // Local mole concentration of O2 in the carrier phase,            
    const scalar C_O2 = YO2*rhoc/WO2_;

    // Diffusion rate coefficient
    const scalar D0 = C1_*pow(0.5*(T + Tc), 2.0);
    const scalar D0ash = C1_*pow(T, 2.0);

    // Schmidt number
    const scalar Sc = muc/(rhoc*D0);

    // Sherwood number
    const scalar Sh = 2.0+1.1*pow(Sc, 0.3333)*pow(Rec, 0.6);

    // Mass convection coefficient hma
    const scalar hma = D0*Sh/(2.0*deq);

    // Mass diffusion coefficient through the ash layer hmi
    const scalar hmi = D0ash*pow(ep3_, 2.0)/(d-di);

    // General mass transfer rate for combustion
    const scalar beta_d = hmi*hma/(hmi+hma);

    // Kinetic rate beta_r
    const scalar beta_r = C2_*T*exp(-E_/T);
        
    // Correlation for the formation of CO and CO2 omega_c
    const scalar omega_c = 2.0*(1.0+4.3*exp(-COmega_/T))/(2.0+4.3*exp(-COmega_/T));
    
    // Particle surface area for the char layer
    scalar Ap;
    if (particleShape == 1)
    {
        Ap = 2.0*constant::mathematical::pi*(3.0*Foam::sqr(di)-di*Xi);
    }
    else
    {
        Ap = constant::mathematical::pi*sqr(di*2.0);
    }

    // Change in C mass [kg]
    scalar dmC = WC_*omega_c*C_O2*Ap*(beta_r*beta_d/(beta_r+beta_d))*dt;

    // Limit mass transfer by availability of C
    dmC = min(mass*fComb, dmC);

    // Molar consumption
    const scalar dOmega = dmC/omega_c/WC_;

    // Change in O2 mass [kg]
    const scalar dmO2 = dOmega*WO2_;

    // Mass of newly created CO2 [kg]
    const scalar dmCO2 = dOmega*(2.0-omega_c)*(WC_ + WO2_);

    // Mass of newly created CO [kg]
    const scalar dmCO = dOmega*(2.0*(omega_c-1.0))*(WC_ + WO2_/2.0);

    // Update local particle C mass
    dMassSolid[CsLocalId_] = dOmega*WC_*omega_c;
 
    // Update carrier O2 and CO2 mass
    dMassSRCarrier[O2GlobalId_] -= dmO2;
    dMassSRCarrier[CO2GlobalId_] += dmCO2;
    dMassSRCarrier[COGlobalId_] += dmCO;


    const scalar LHV_index_CO2 = 0.0;               /* per definition */
    const scalar LHV_index_CO = 10.25e6;            /* J/kg D&D p. 26 */
    const scalar LHV_index_C = 32.79e6;         /* J/kg Alberto's code */
    // const scalar HsC = thermo.solids().properties()[CsLocalId_].Hs(T);

    // Correlation for heatRatio CO and CO2
    //const scalar COCO2Ratio = 4.3*exp(-COmega_/T);
    //const scalar heatRatio = (COCO2Ratio+0.3)/(COCO2Ratio+1.0);

    // carrier sensible enthalpy exchange handled via change in mass
    ; 
    // Heat of reaction [J]
    return dmCO2*LHV_index_CO2 + dmCO*LHV_index_CO - dmC*LHV_index_C;
}


// ************************************************************************* //
