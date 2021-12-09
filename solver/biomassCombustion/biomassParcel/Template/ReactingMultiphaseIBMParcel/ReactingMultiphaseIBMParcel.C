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

#include "ReactingMultiphaseIBMParcel.H"
#include "mathematicalConstants.H"

#define CP_GAS 1.18698e3

using namespace Foam::constant::mathematical;

// * * * * * * * * * * * *  Helper Member Functions * * * * * * * * * * * * //

template<class ParcelType>
Foam::scalar Foam::ReactingMultiphaseIBMParcel<ParcelType>::radiusForCylinder
(
    const scalar Xi, 
    const scalar V0
)
{
    const scalar a = 2.0*constant::mathematical::pi;
    scalar b = -Xi*constant::mathematical::pi;
    //double c = 0;
    scalar d = -V0;
    scalar root1, root2, root3;
    scalar root;
    
    if (d == 0 || d>0)
    {
        Info<<"ERROR: The volume of cylinder should be positive!! Now the value: "<<V0<<"\n"<<endl;
        return 0;
    } //End if d == 0

    b /= a;
    d /= a;
    
    scalar disc, q, r, dum1, s, t, term1, r13;
    q = (-(b*b))/9.0;
    r = -(27.0*d) + b*(-2.0*(b*b));
    r /= 54.0;
    disc = q*q*q + r*r;
    term1 = (b/3.0);
    
    if (disc > 0) { // one root real, two are complex
        s = r + Foam::sqrt(disc);
        s = ((s < 0) ? -Foam::pow(-s, (1.0/3.0)) : Foam::pow(s, (1.0/3.0)));
        t = r - Foam::sqrt(disc);
        t = ((t < 0) ? -Foam::pow(-t, (1.0/3.0)) : Foam::pow(t, (1.0/3.0)));
        root = -term1 + s + t;
        return root;
    } 
    // End if (disc > 0)
    
    // The remaining options are all real
    if (disc == 0){ // All roots real, at least two are equal.
        r13 = ((r < 0) ? -Foam::pow(-r,(1.0/3.0)) : Foam::pow(r,(1.0/3.0)));
        root1 = -term1 + 2.0*r13;
        root2 = -(r13 + term1);
        root = ((root1 > root2) ? root1 : root2);
        return root;
    } // End if (disc == 0)
    
    // Only option left is that all roots are real and unequal (to get here, q < 0)
    q = -q;
    dum1 = q*q*q;
    dum1 = Foam::acos(r/Foam::sqrt(dum1));
    r13 = 2.0*Foam::sqrt(q);
    root1 = -term1 + r13*Foam::cos(dum1/3.0);
    root2 = -term1 + r13*Foam::cos((dum1 + 2.0*constant::mathematical::pi)/3.0);
    root3 = -term1 + r13*Foam::cos((dum1 + 4.0*constant::mathematical::pi)/3.0);
    if (root1 > root2){
        root = ((root1 > root3) ? root1 : root3);
    }
    else{
        root = ((root2 > root3) ? root2 : root3);
    }
    return root;
}  //End of cubicSolve


template<class ParcelType>
Foam::scalar Foam::ReactingMultiphaseIBMParcel<ParcelType>::Area_Sph
(
    const scalar radius
)
{
    return 4.0*constant::mathematical::pi*Foam::sqr(radius);
}


template<class ParcelType>
Foam::scalar Foam::ReactingMultiphaseIBMParcel<ParcelType>::Area_cylinderL
(
    const scalar radius,
    const scalar Xi
)
{
    return 2.0*constant::mathematical::pi*(3.0*Foam::sqr(radius)-radius*Xi);
}


template<class ParcelType>
Foam::scalar Foam::ReactingMultiphaseIBMParcel<ParcelType>::Volume_cylinderL
(
    const scalar radius,
    const scalar Xi
)
{
    return constant::mathematical::pi*(2.0*Foam::pow3(radius)-sqr(radius)*Xi);
}


template<class ParcelType>
Foam::scalar Foam::ReactingMultiphaseIBMParcel<ParcelType>::Vol_Rin 
(
    const scalar rin,
    const scalar rout
)
{
    
    return (4./3.)*constant::mathematical::pi*(Foam::pow3(rout) - Foam::pow3(rin));
}


template<class ParcelType>
Foam::scalar Foam::ReactingMultiphaseIBMParcel<ParcelType>::Vol_Rin_cylinderL
(
    const scalar rin,
    const scalar rout,
    const scalar Xi
)
{
    return constant::mathematical::pi*(2.0*Foam::pow3(rout) - 2.0*Foam::pow3(rin) - (Foam::sqr(rout) - Foam::sqr(rin))*Xi);
}


template<class ParcelType>
Foam::scalar Foam::ReactingMultiphaseIBMParcel<ParcelType>::R_Par
(
    const scalar rin,
    const scalar rout
)
{
    scalar temp ,r;
    
    temp = (Foam::pow3(rin)+Foam::pow3(rout))/2.0;
    
    r = Foam::cbrt(temp);
    
    return r;
}


template<class ParcelType>
Foam::scalar Foam::ReactingMultiphaseIBMParcel<ParcelType>::R_Par_cylinderL
(
    const scalar rin,
    const scalar rout,
    const scalar Xi
)
{
    scalar temp ,r;
    
    temp = (2.0*constant::mathematical::pi*(pow3(rout) + pow3(rin))-constant::mathematical::pi*(pow(rout,2.0) + pow(rin,2.0))*Xi)/2.0;
    
    r = radiusForCylinder(Xi, temp);
    
    return r;
}


template<class ParcelType>
Foam::scalar Foam::ReactingMultiphaseIBMParcel<ParcelType>::d_dr
(
    const scalar kp,
    const scalar Ab,
    const scalar ri,
    const scalar rj
)
{
    return  kp*Ab*rj/(ri*(ri-rj));  
}


template<class ParcelType>
Foam::scalar Foam::ReactingMultiphaseIBMParcel<ParcelType>::d_dr_cylinderL
(
    const scalar kp,
    const scalar Ab,
    const scalar ri,
    const scalar rj,
    const scalar Xi
)
{
   return  kp*Ab*Xi/((ri*Xi-3.0*Foam::sqr(ri))*Foam::log((3.0-Xi/rj)/(3.0-Xi/ri)));  
}


template<class ParcelType>
Foam::scalar Foam::ReactingMultiphaseIBMParcel<ParcelType>::Fb1 
(
    const scalar Tb1
)
{
        
    return  Foam::pow(10.,(8.07131-1730.63/(Tb1-39.724)))/760.0;
}


template<class ParcelType>
Foam::scalar Foam::ReactingMultiphaseIBMParcel<ParcelType>::eq4 
(
    const scalar h_coe,
    const scalar emissi,
    const scalar Ste_Bol,
    const scalar kp3,
    const scalar Ab3,
    const scalar rb3,
    const scalar rp3,
    const scalar Tg,
    const scalar Tp3,
    const scalar G,
    const scalar Source
)
{
    //the mothed used here is General formula for roots, see wikipedia, well it is quite unstable in practice.
    scalar aa, dd, ee, delta0, delta, delta1, Q, S, q, Tb3,sqrt1,sqrt2,sqrt3;
    aa = emissi*Ste_Bol*Ab3;
    dd = h_coe*Ab3+Ab3*kp3*rp3/(rb3*(rb3-rp3));
    ee = -Ab3*kp3*rp3/(rb3*(rb3-rp3))*Tp3-h_coe*Ab3*Tg-aa*(G/(4.*Ste_Bol))-Source;
//  Info<<"Tg emission:   "<<G/(4.*Ste_Bol)<<endl;
    delta = 256.0*Foam::pow3(aa)*Foam::pow3(ee)-27.0*aa*aa*Foam::pow4(dd);
    delta1 = 27.0*aa*dd*dd;
    delta0 = 12.0*aa*ee;
    sqrt1 = -27.0*delta;
    if (sqrt1  < 0)
    {
            FatalErrorInFunction
            <<"sqrt1 < 0" << exit(FatalError);
    }
                    
    Q = Foam::cbrt(0.5*(delta1+Foam::sqrt(sqrt1)));
    sqrt2 = (Q+delta0/Q)/(3.0*aa);
    if (sqrt2 < 0)
    {
            FatalErrorInFunction
            <<"sqrt2 < 0" << exit(FatalError);
    }
    S = 0.5*Foam::sqrt(sqrt2);
    q = dd/aa;
    sqrt3 = (-4*S*S+q/S);
    if (sqrt3 < 0)
    {
            FatalErrorInFunction
            <<"sqrt3 < 0" << exit(FatalError);
    }               
    Tb3 = -S+0.5*Foam::sqrt(sqrt3);
    
    return Tb3;
}


template<class ParcelType>
Foam::scalar Foam::ReactingMultiphaseIBMParcel<ParcelType>::eq7_2
(
    const scalar kp1,
    const scalar kp2,
    const scalar Ab1,
    const scalar rb1,
    const scalar rp1,
    const scalar rp2,
    const scalar Tp1,
    const scalar Tp2,
    const scalar Qb1    
)
{
    scalar temp1 = d_dr(kp2,Ab1,rb1,rp2), temp2 = d_dr(kp1,Ab1,rb1,rp1);
        
    return (-Qb1+temp1*Tp2-temp2*Tp1)/(temp1-temp2);    
}


template<class ParcelType>
Foam::scalar Foam::ReactingMultiphaseIBMParcel<ParcelType>::eq7_3
(
    const scalar kp1,
    const scalar kp2,
    const scalar Ab1,
    const scalar rb1,
    const scalar rp1,
    const scalar rp2,
    const scalar Tp1,
    const scalar Tp2,
    const scalar Fb    
)
{
    scalar temp1 = d_dr(kp2,Ab1,rb1,rp2), temp2 = d_dr(kp1,Ab1,rb1,rp1);
        
    return (temp1*(1.-Fb)*Tp2-temp2*Tp1)/(temp1*(1.-Fb)-temp2);   
}


template<class ParcelType>
Foam::scalar Foam::ReactingMultiphaseIBMParcel<ParcelType>::eq4_cylinderL
(
    const scalar h_coe,
    const scalar emissi,
    const scalar Ste_Bol,
    const scalar kp3,
    const scalar Ab3,
    const scalar rb3,
    const scalar rp3,
    const scalar Tg,
    const scalar Tp3,
    const scalar G,
    const scalar Source,
    const scalar Xi
)
{
    scalar aa, dd, ee, delta0, delta, delta1, Q, S, q, Tb3;
    aa = emissi*Ste_Bol*Ab3;
    dd = h_coe*Ab3+d_dr_cylinderL(kp3, Ab3, rb3, rp3, Xi);
    ee = -d_dr_cylinderL(kp3, Ab3, rb3, rp3, Xi)*Tp3-h_coe*Ab3*Tg-aa*(G/(4.*Ste_Bol))-Source;
    delta = 256.0*Foam::pow3(aa)*Foam::pow3(ee)-27.0*aa*aa*Foam::pow4(dd);
    delta1 = 27.0*aa*dd*dd;
    delta0 = 12.0*aa*ee;
    Q = Foam::cbrt(0.5*(delta1+Foam::sqrt(-27.0*delta)));
    S = 0.5*Foam::sqrt((Q+delta0/Q)/(3.0*aa));
    q = dd/aa;
    
    Tb3 = -S+0.5*Foam::sqrt(-4*S*S+q/S);
    
    return Tb3;
}


template<class ParcelType>
Foam::scalar Foam::ReactingMultiphaseIBMParcel<ParcelType>::eq7_2_cylinderL
(
    const scalar kp1,
    const scalar kp2,
    const scalar Ab1,
    const scalar rb1,
    const scalar rp1,
    const scalar rp2,
    const scalar Tp1,
    const scalar Tp2,
    const scalar Qb1,
    const scalar Xi
)
{
    scalar temp1 = d_dr_cylinderL(kp2,Ab1,rb1,rp2,Xi), temp2 = d_dr_cylinderL(kp1,Ab1,rb1,rp1,Xi);
    
// Info<<"source: "<<Qb1<<endl;
// Info<<"(temp1-temp2): "<<(temp1-temp2)<<endl;
// Info<<"source contribution: "<<Qb1/(temp1-temp2)<<endl;

    return (-Qb1+temp1*Tp2-temp2*Tp1)/(temp1-temp2);    
}


template<class ParcelType>
Foam::scalar Foam::ReactingMultiphaseIBMParcel<ParcelType>::eq7_3_cylinderL
(
    const scalar kp1,
    const scalar kp2,
    const scalar Ab1,
    const scalar rb1,
    const scalar rp1,
    const scalar rp2,
    const scalar Tp1,
    const scalar Tp2,
    const scalar Fb,
    const scalar Xi
)
{
    scalar temp1 = d_dr_cylinderL(kp2,Ab1,rb1,rp2,Xi), temp2 = d_dr_cylinderL(kp1,Ab1,rb1,rp1,Xi);
        
    return (temp1*(1.-Fb)*Tp2-temp2*Tp1)/(temp1*(1.-Fb)-temp2);   
}

template<class ParcelType>
Foam::scalar Foam::ReactingMultiphaseIBMParcel<ParcelType>::cp_p 
(
    const label layer,
    const scalar Tp,
    const scalar moist_WB
)
{ 
    scalar heat_capacity=0.;
    
    scalar Tp_2 = Tp;
//     if (Tp_2>=400)
//     {
//         Tp_2=400;
//     }
//     
    if (layer == 0)
    {
        /* moist wood */
        scalar c_p_wood = 4.206*Tp_2 - 37.7; /* eq. (24.i) Thunman et al., Energy & Fuels 2001 */
        scalar A = 1.0e3*((0.02355*Tp_2 - 1.320*moist_WB/(1.0 - moist_WB) - 6.191)*moist_WB/(1.0 - moist_WB));
        heat_capacity = c_p_wood*(1.0 - moist_WB) + 4185.5*moist_WB + A;    
    }
    else if (layer == 1)
    {
        /* dry wood */
        heat_capacity = 4.206*Tp_2 - 37.7; /* eq. (24.i) Thunman et al., Energy & Fuels 2001 */
    }
    else if (layer == 2)
    {
        /* char */
        heat_capacity = -334.0 + 4410.0e-3*Tp_2 - 3160.0e-6*Foam::pow(Tp_2,2.0) + 1010.0e-9*Foam::pow(Tp_2,3.0) - 119.0e-12*Foam::pow(Tp_2,4.0); /* Table 5 Thunman et al., Energy & Fuels 2001 */
    }
    else if (layer == 3)
    {
        heat_capacity = 754.0 + 0.586*(Tp_2-273.15);     /* Sven */
    }
    
    return heat_capacity;
}

template<class ParcelType>
Foam::scalar Foam::ReactingMultiphaseIBMParcel<ParcelType>::cp_p_modified
(
    const label layer,
    const scalar Tp,
    const scalar moist_WB
)
{ 
    scalar heat_capacity=0.;
    
    scalar Tp_2 = Tp;
    
    if (Tp_2 < 273)
    {
        Tp_2=273;
    }
    
    if (layer == 0)
    {
        /* moist wood */
        scalar c_p_wood = 4.206*Tp_2 - 37.7; /* eq. (24.i) Thunman et al., Energy & Fuels 2001 */
        scalar A = 1.0e3*((0.02355*Tp_2 - 1.320*moist_WB/(1.0 - moist_WB) - 6.191)*moist_WB/(1.0 - moist_WB));
        heat_capacity = c_p_wood*(1.0 - moist_WB) + 4185.5*moist_WB + A;    
    }
    else if (layer == 1)
    {
        /* dry wood */
        heat_capacity = 4.206*Tp_2 - 37.7; /* eq. (24.i) Thunman et al., Energy & Fuels 2001 */
    }
    else if (layer == 2)
    {
        /* char */
        heat_capacity = 420+2.09*Tp_2+0.000685*Foam::pow(Tp_2,2.0); /* Ramin mehrabian., Fuels process techology 2012 */
    }
    else if (layer == 3)
    {
        heat_capacity = 754.0 + 0.586*(Tp_2-273.15);     /* Sven */
    }
    
    return heat_capacity;
}

template<class ParcelType>
Foam::scalar Foam::ReactingMultiphaseIBMParcel<ParcelType>::rho_p 
(
    const label layer,
    const scalar Tp,
    const scalar moist_WB
)
{
    scalar density=0;
    
    if (layer == 0)
    {
        density = 750.0;
    }
    else if (layer == 1)
    {
        density = (1.0 - moist_WB)*750.0; 
    }
    else if (layer == 2)
    {
        density = 150.0;
    }
    else if (layer == 3)
    {
        scalar dgas = 101325.0*(2.0*14.00672*1.0e-3)/(8.3145*Tp);
        density = 2000.0*(1.0-0.65) + 0.65*dgas;          /* Sven */
    }
    
    return density;
}
    
template<class ParcelType>
Foam::scalar Foam::ReactingMultiphaseIBMParcel<ParcelType>::kp_p 
(
    const label layer,
    const scalar Tp
)
{
    scalar heat_conductivity=0;
    
    if (layer == 0)
    {
        heat_conductivity = 0.2;                 
    }
    else if (layer == 1)
    {
        heat_conductivity = 0.11;
    }
    else if (layer == 2)
    {
        heat_conductivity = 0.052;
    }
    else if (layer == 3)
    {
        scalar kgas = (-0.000000000003493)*Foam::pow3(Tp) + 0.000000003319264*Foam::sqr(Tp) + 0.000060059499759*Tp + 0.008533051948052;
        heat_conductivity = 1.03*(1.0-0.65) + 0.65*kgas;  /* Sven */
    }
    
    return heat_conductivity;
} 

template<class ParcelType>
Foam::scalar Foam::ReactingMultiphaseIBMParcel<ParcelType>::cp_water_vapor 
(
    const scalar T
)
{
    scalar A, B, C, D, E, Temp;
    A = 30.09200;
    B = 6.832514;
    C = 6.793435;
    D = -2.534480;
    E = 0.082139;
    Temp = T/1000.;
    return (A+B*Temp+C*Foam::sqr(Temp)+D*Foam::pow3(Temp)+E/Foam::sqr(Temp))/18.*1000.;    
}

template<class ParcelType>
Foam::scalar Foam::ReactingMultiphaseIBMParcel<ParcelType>::deltaHvap 
(
    const scalar T
)
{
    scalar A, B, C, dH;
    A = 51798.3515625;
    B = -10.6430501938;
    C = -0.0515099391341;
    if (T < 273.15)
    {
        dH = 45.054e3;
    }
    else if ( T > 433.15 )
    {
        dH = 37.518e3;
    }
    else
    {
        dH = A + B*T + C*Foam::sqr(T);
    }
    
    return dH/18.0e-3; /* convert J/mol to J/kg */
}


// * * * * * * * * * * *  Protected Member Functions * * * * * * * * * * * * //

template<class ParcelType>
template<class TrackCloudType>
void Foam::ReactingMultiphaseIBMParcel<ParcelType>::setCellValues
(
    TrackCloudType& cloud,
    trackingData& td
)
{
    ParcelType::setCellValues(cloud, td);
}


template<class ParcelType>
template<class TrackCloudType>
void Foam::ReactingMultiphaseIBMParcel<ParcelType>::cellValueSourceCorrection
(
    TrackCloudType& cloud,
    trackingData& td,
    const scalar dt
)
{
    // Re-use correction from reacting parcel
    ParcelType::cellValueSourceCorrection(cloud, td, dt);
}

#include "IBMParcelMainCalcFunction.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //
template<class ParcelType>
Foam::scalar Foam::ReactingMultiphaseIBMParcel<ParcelType>::updateMassFractionsProtected
(
    const scalar mass0,
    const scalarField& dMassGas,
    const scalarField& dMassLiquid,
    const scalarField& dMassSolid
)
{
    scalarField& YMix = this->Y_;

    scalar massGas =
        this->updateMassFraction(mass0*YMix[this->GAS], dMassGas, this->YGas_);
    scalar massLiquid =
        this->updateMassFraction(mass0*YMix[this->LIQ], dMassLiquid, this->YLiquid_);
    scalar massSolid =
        this->updateMassFraction(mass0*YMix[this->SLD], dMassSolid, this->YSolid_);

    scalar massNew = max(massGas + massLiquid + massSolid, rootVSmall);

    YMix[this->GAS] = massGas/massNew;
    YMix[this->LIQ] = massLiquid/massNew;
    YMix[this->SLD] = 1.0 - YMix[this->GAS] - YMix[this->LIQ];

    return massNew;
}


template<class ParcelType>
template<class TrackCloudType>
void Foam::ReactingMultiphaseIBMParcel<ParcelType>::calcSurfaceReactions
(
    TrackCloudType& cloud,
    trackingData& td,
    const scalar dt,
    const scalar d,
    const scalar di,
    const scalar T,
    const scalar mass,
    const label canCombust,
    const scalar N,
    const scalarField& YMix,
    const scalarField& YGas,
    const scalarField& YLiquid,
    const scalarField& YSolid,
    scalarField& dMassSRGas,
    scalarField& dMassSRLiquid,
    scalarField& dMassSRSolid,
    scalarField& dMassSRCarrier,
    scalar& Sh,
    scalar& dhsTrans,
    scalar& QComb,
    const scalar Re,
    const scalar Tc,
    const scalar rhoc,
    const scalar muc,
    const scalar Xi,
    scalarField& dHeatSRCarrier,
    const label particleShape,
    const scalarField& epsilons,
    const scalar deq
) const
{
    // Check that model is active
    if (!cloud.charOxidization().active())
    {
        return;
    }

    // Initialise demand-driven constants
    (void)cloud.constProps().hRetentionCoeff();
    (void)cloud.constProps().TMax();

    // Check that model is active
    if (canCombust != 1)
    {
        return;
    }

    // Update surface reactions
    // const scalar hReaction = cloud.charOxidization().calculate
    cloud.charOxidization().calculate
    (
        dt,
        this->cell(),
        d,
        T,
        td.Tc(),
        td.pc(),
        td.rhoc(),
        mass,
        YGas,
        YLiquid,
        YSolid,
        YMix,
        N,
        dMassSRGas,
        dMassSRLiquid,
        dMassSRSolid,
        dMassSRCarrier,
        di,
        td.muc(),
        Re,
        Xi,
        dHeatSRCarrier,
        particleShape,
        epsilons,
        deq
    );
    

    cloud.charOxidization().addToCharOxidizationMass
    (
        this->nParticle_
       *(sum(dMassSRGas) + sum(dMassSRLiquid) + sum(dMassSRSolid))
    );

    QComb = -(dHeatSRCarrier[0]+dHeatSRCarrier[1]+dHeatSRCarrier[2]+dHeatSRCarrier[3])/dt;

//     const scalar xsi = min(T/cloud.constProps().TMax(), 1.0);
//     const scalar coeff =
//         (1.0 - xsi*xsi)*cloud.constProps().hRetentionCoeff();
// 
//     Sh += coeff*hReaction/dt;
// 
//     dhsTrans += (1.0 - coeff)*hReaction;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ParcelType>
Foam::ReactingMultiphaseIBMParcel<ParcelType>::ReactingMultiphaseIBMParcel
(
    const ReactingMultiphaseIBMParcel<ParcelType>& p
)
:
    ParcelType(p)
{}


template<class ParcelType>
Foam::ReactingMultiphaseIBMParcel<ParcelType>::ReactingMultiphaseIBMParcel
(
    const ReactingMultiphaseIBMParcel<ParcelType>& p,
    const polyMesh& mesh
)
:
    ParcelType(p, mesh)
{}


// * * * * * * * * * * * * * * IOStream operators  * * * * * * * * * * * * * //

#include "ReactingMultiphaseIBMParcelIO.C"

// ************************************************************************* //
