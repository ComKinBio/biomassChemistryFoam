/*---------------------------------------------------------------------------*\
 * =========                 |
 * \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
 * \\    /   O peration     |
 *    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
 *    \\/     M anipulation  |
 * -------------------------------------------------------------------------------
 * License
 *    This file is part of OpenFOAM.
 * 
 *    OpenFOAM is free software: you can redistribute it and/or modify it
 *    under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 * 
 *    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
 *    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 *    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
 *    for more details.
 * 
 *    You should have received a copy of the GNU General Public License
 *    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.
 * 
 * \*---------------------------------------------------------------------------*/

#include "DiBlasiPyrolysis.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::DiBlasiPyrolysis<CloudType>::
DiBlasiPyrolysis
(
    const dictionary& dict,
    CloudType& owner
)
:
    PyrolysisModel<CloudType>(dict, owner, typeName),
    volgas1_(this->coeffDict().lookup("volgas1")),
    volgas2_(this->coeffDict().lookup("volgas2")),
    volatileToGas1Map_(volgas1_.size()),
    volatileToGas2Map_(volgas2_.size()),
    residualCoeff_(readScalar(this->coeffDict().lookup("residualCoeff"))),
    secondGasVersion_(readScalar(this->coeffDict().lookup("secondGasVersion"))),
    woodProximate_(this->coeffDict().lookup("woodProximate")),
    tarProximate_(this->coeffDict().lookup("tarProximate")),
    devolKinetic1_(this->coeffDict().lookup("devolKinetic1")),
    devolKinetic2_(this->coeffDict().lookup("devolKinetic2")),
    devolKinetic3_(this->coeffDict().lookup("devolKinetic3")),
    devolKinetic4_(this->coeffDict().lookup("devolKinetic4")),
    devolKinetic5_(this->coeffDict().lookup("devolKinetic5"))
{
    if (volgas1_.empty() || volgas2_.empty())
    {
        WarningIn
        (
            "Foam::DiBlasiPyrolysis<CloudType>::"
            "DiBlasiPyrolysis"
            "("
            "const dictionary& dict, "
            "CloudType& owner"
            ")"
        )   << "Pyrolysis model selected, but parameters are not well defined"
        << nl << endl;
    }
    else
    {
        
        const label idGas = owner.composition().idGas();	
        
        Info << "Volatile info:" << endl;
        Info << " volatile gas 1:" << endl;
        forAll(volgas1_, i)
        {
            const word& specieName = volgas1_[i].name();
            const label id = owner.composition().localId(idGas, specieName);
            volatileToGas1Map_[i] = id;
            Info<< "  " << specieName << "mass fraction = "<< volgas1_[i].W() 
            <<"  map gas ID: "<< id <<endl; 	  
        }
        Info << " volatile gas 2:" << endl;
        if (secondGasVersion_ == 0)
        {
            forAll(volgas2_, i)
            {
                const word& specieName = volgas2_[i].name();;
                const label id = owner.composition().localId(idGas, specieName);
                volatileToGas2Map_[i] = id;
                Info<< "  " << specieName << "  mass fraction = "<< volgas2_[i].W() 
                <<"  map gas ID: "<< id <<endl; 	  
            }  
        }
        else if (secondGasVersion_ == 1)
        {
            forAll(volgas2_, i)
            {
                const word& specieName = volgas2_[i].name();;
                const label id = owner.composition().localId(idGas, specieName);
                volatileToGas2Map_[i] = id;
                Info<< "  " << specieName << "  map gas ID: "<< id <<"  A:"
                << volgas2_[i].A() << "  E:"<< volgas2_[i].E()<<endl; 	  
            }  
        }
        else if (secondGasVersion_ == 2)
        {
            forAll(volgas2_, i)
            {
                const word& specieName = volgas2_[i].name();;
                const label id = owner.composition().localId(idGas, specieName);
                volatileToGas2Map_[i] = id;
                Info<< "  no second cracking in solid phase"<<endl; 	  
            }  
        }
        else if (secondGasVersion_ == 3)
        {
            forAll(volgas2_, i)
            {
                const word& specieName = volgas2_[i].name();;
                const label id = owner.composition().localId(idGas, specieName);
                volatileToGas2Map_[i] = id;
                Info<< "  second cracking in solid phase"<<endl; 	  
            }  
        }
        else if (secondGasVersion_ == 4)
        {
            forAll(volgas2_, i)
            {
                const word& specieName = volgas2_[i].name();;
                const label id = owner.composition().localId(idGas, specieName);
                volatileToGas2Map_[i] = id;
                Info<< "  second cracking in solid phase"<<endl; 	  
            }  
        }
        else
        {
            FatalErrorIn
            (
                "secondGasVersion wrong"
            )   << "secondGasVersion wrong" << exit(FatalError);
        }
        
        Info << " wood proximate analysis: C " << woodProximate_[0]*100 << "%"
        << " H " << woodProximate_[1]*100 << "%" << " O " 
        << woodProximate_[2]*100<< "%" <<endl;
        Info << " tar proximate analysis: C " << tarProximate_[0]*100 << "%"
        << " H " << tarProximate_[1]*100 << "%" << " O " 
        << tarProximate_[2]*100<< "%" <<endl;
        Info << " Kinetic data:" << endl;
        Info << "  wood -> gas1 A:" << devolKinetic1_[0] << " E:" << devolKinetic1_[1] <<endl;
        Info << "  wood -> tar A:" << devolKinetic2_[0] << " E:" << devolKinetic2_[1] <<endl;
        Info << "  wood -> char A:" << devolKinetic3_[0] << " E:" << devolKinetic3_[1] <<endl;
        if (secondGasVersion_ == 0)
        {
            Info << "  tar -> gas2 A:" << devolKinetic4_[0] << " E:" << devolKinetic4_[1] <<endl; 
        }
        else if (secondGasVersion_ == 1)
        {
            Info << "  tar -> gas2: depending on the individual gas kinetic" <<endl;
        }
        
        Info << "  tar -> char A:" << devolKinetic5_[0] << " E:" << devolKinetic5_[1] <<endl;
        
        
    }
}


// template<class CloudType>
// Foam::DiBlasiPyrolysis<CloudType>::
// DiBlasiPyrolysis
// (
//     const DiBlasiPyrolysis<CloudType>& dm
// )
// :
//     PyrolysisModel<CloudType>(dm),
//     volgas1_(dm.volgas1_),
//     volgas2_(dm.volgas2_),
//     volatileToGas1Map_(dm.volatileToGas1Map_),
//     volatileToGas2Map_(dm.volatileToGas2Map_),    
//     residualCoeff_(dm.residualCoeff_),
//     secondGasVersion_(dm.secondGasVersion_),
//     woodProximate_(dm.woodProximate_),
//     tarProximate_(dm.tarProximate_),  
//     devolKinetic1_(dm.devolKinetic1_),
//     devolKinetic2_(dm.devolKinetic2_),
//     devolKinetic3_(dm.devolKinetic3_),
//     devolKinetic4_(dm.devolKinetic4_),
//     devolKinetic5_(dm.devolKinetic5_)
// {}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::DiBlasiPyrolysis<CloudType>::
~DiBlasiPyrolysis()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::DiBlasiPyrolysis<CloudType>::calculate
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
    //add
    scalarField& dMassSOLID
) const
{
    bool done = true;
    
    // Model coefficients
    const scalar A1 = devolKinetic1_[0];
    const scalar E1 = devolKinetic1_[1];
    const scalar A2 = devolKinetic2_[0];
    const scalar E2 = devolKinetic2_[1];
    const scalar A3 = devolKinetic3_[0];
    const scalar E3 = devolKinetic3_[1];
    const scalar A4 = devolKinetic4_[0];
    const scalar E4 = devolKinetic4_[1];
    const scalar A5 = devolKinetic5_[0];
    const scalar E5 = devolKinetic5_[1];
    
    // Kinetic rate
    const scalar kappa1 = A1*exp(-E1/(RR*T));
    const scalar kappa2 = A2*exp(-E2/(RR*T));
    const scalar kappa3 = A3*exp(-E3/(RR*T));
    const scalar kappa4 = A4*exp(-E4/(RR*T));
    const scalar kappa5 = A5*exp(-E5/(RR*T));
    
    
    //- Get id for wood char tar
    const label idGas = this->owner().composition().idGas();
    const label idSolid = this->owner().composition().idSolid();
    //     const label idLiquid = this->owner().composition().idLiquid();
//     const scalar YSolidTot = this->owner().composition().YMixture0()[idSolid];
//     const scalarField& YSolid = this->owner().composition().Y0(idSolid);
    const label id_wood = this->owner().composition().localId(idSolid, "activeDryWood");
    const label id_char = this->owner().composition().localId(idSolid, "C");
    const label id_tar = this->owner().composition().localId(idGas, "tar");
    
//     const label idDryWood = this->owner().composition().localId(idSolid, "wood");
    
//     const label idAsh = this->owner().composition().localId(idSolid, "ash");
    //     const label id_H2O = this->owner().composition().localId(idLiquid, "H2O");
    
    // initial fraction
//     const scalar Ywood00 = this->owner().composition().YMixture0()[idSolid]*this->owner().composition().Y0(idSolid)[idDryWood];
//     const scalar Yash00 = this->owner().composition().YMixture0()[idSolid]*this->owner().composition().Y0(idSolid)[idAsh];
//     
//     const scalar ashDryWood = (Yash00 + Ywood00)/Ywood00;
    
    //- Inititial wood mass
//     const scalar massWood0 = mass0*YSolidTot*YSolid[this->owner().composition().localId(idSolid, "wood")];//may not work
    
    //- Current wood mass
    const scalar massWood = mass;
//     const scalar massWood = mass/ashDryWood;
    
    //- Combustion allowed once all wood evolved
    //done = done && (massWood <= residualCoeff_*massWood0);

    const scalar dmasswood_temp = dt*(kappa1+kappa2+kappa3)*massWood;   
    scalar tarprodmass = 0.0;
    if (dmasswood_temp < massWood /*&& massWood - dmasswood_temp >= massWood0*0.00005*/)
    {
        dMassSOLID[id_wood] = dmasswood_temp;
        tarprodmass = kappa2*massWood*dt;
        dMassSOLID[id_char] = kappa3*massWood*dt;
        forAll(volgas1_, i)
        {
            const label id_gas1 = volatileToGas1Map_[i];
            dMassDV[id_gas1] = volgas1_[i].W()*kappa1*massWood*dt;
        }
    }
    else
    {
        dMassSOLID[id_wood] = massWood;//-massWood0*0.00005;
        tarprodmass = kappa2*(dMassSOLID[id_wood]/(kappa1+kappa2+kappa3));    
        dMassSOLID[id_char] = kappa3*(dMassSOLID[id_wood]/(kappa1+kappa2+kappa3));
        forAll(volgas1_, i)
        {
            const label id_gas1 = volatileToGas1Map_[i];
            dMassDV[id_gas1] = volgas1_[i].W()*kappa1*(dMassSOLID[id_wood]/(kappa1+kappa2+kappa3));
        }    
        
//         canCombust = 0;//need to be changed for combustion!!!!!!!!!!!!!!!!!!!!!!!
    }

    if (secondGasVersion_ == 0)
    {
        // tar
        dMassDV[id_tar] = max((1.0-dt*(kappa4+kappa5))*tarprodmass, 0.0);
        
        if (dMassDV[id_tar] == 0.0)
        {
            // gas2
            forAll(volgas2_, i)
            {
                const label id_gas2 = volatileToGas2Map_[i];
                dMassDV[id_gas2] = volgas2_[i].W()*tarprodmass*kappa4/(kappa4+kappa5)+dMassDV[id_gas2];   
            }
            
            //char
            dMassSOLID[id_char] = -(tarprodmass*kappa5/(kappa4+kappa5)+dMassSOLID[id_char]);
        }
        
        else
        {
            forAll(volgas2_, i)
            {
                const label id_gas2 = volatileToGas2Map_[i];
                dMassDV[id_gas2] = volgas2_[i].W()*tarprodmass*kappa4*dt+dMassDV[id_gas2];   
            }
            
            //char
            dMassSOLID[id_char] = -(tarprodmass*kappa5*dt+dMassSOLID[id_char]);  
        }
        
    }
    
    else if (secondGasVersion_ == 1)
    {
        List<scalar> kappa4_2(volgas2_.size());
        scalar sumkappa4_2 = 0.0;
        
        // gas2
        forAll(volgas2_, i)
        {
            const scalar A4_2 = volgas2_[i].A();
            const scalar E4_2 = volgas2_[i].E();
            
            // Kinetic rate
            kappa4_2[i] = A4_2*exp(-E4_2/(RR*T));
            
            // sum Kinetic rate
            sumkappa4_2 = kappa4_2[i]+sumkappa4_2;
        }
        
        // tar
        dMassDV[id_tar] = max((1.0-dt*(sumkappa4_2+kappa5))*tarprodmass, 0.0); 
        
        if (dMassDV[id_tar] == 0.0)
        {
            // gas2
            forAll(volgas2_, i)
            {
                const label id_gas2 = volatileToGas2Map_[i];
                dMassDV[id_gas2] = tarprodmass*kappa4_2[i]/(sumkappa4_2+kappa5)+dMassDV[id_gas2];   
            }
            
            //char
            dMassSOLID[id_char] = -(tarprodmass*kappa5/(sumkappa4_2+kappa5)+dMassSOLID[id_char]);
        }
        
        else
        {
            forAll(volgas2_, i)
            {
                const label id_gas2 = volatileToGas2Map_[i];
                dMassDV[id_gas2] = tarprodmass*kappa4_2[i]*dt+dMassDV[id_gas2];   
            }
            
            //char
            dMassSOLID[id_char] = -(tarprodmass*kappa5*dt+dMassSOLID[id_char]);  
        }
        
    }
    
    else if (secondGasVersion_ == 2)
    {
 //Info<<"dMassDV[id_tar] in sub model before"<<dMassDV[id_tar]<<endl;
        dMassDV[id_tar] = dMassDV[id_tar] + tarprodmass ;
        dMassSOLID[id_char] = - dMassSOLID[id_char];
// Info<<"dMassDV[id_tar] in sub model"<<dMassDV[id_tar]<<endl; 
    }
    else if (secondGasVersion_ == 3)
    {
        // tar
        dMassDV[id_tar] = dMassDV[id_tar] + max((1.0-dt*(kappa4+kappa5))*tarprodmass, 0.0);
        
        if ((1.0-dt*(kappa4+kappa5)) < 0.0)
        {
            // gas2
            forAll(volgas1_, i)
            {
                const label id_gas1 = volatileToGas1Map_[i];
                dMassDV[id_gas1] = dMassDV[id_gas1] + volgas1_[i].W()*tarprodmass*kappa4/(kappa4+kappa5);
            }
            
            //char
            dMassSOLID[id_char] = -(tarprodmass*kappa5/(kappa4+kappa5)+dMassSOLID[id_char]);
        }
        
        else
        {
            forAll(volgas1_, i)
            {
                const label id_gas1 = volatileToGas1Map_[i];
                dMassDV[id_gas1] = volgas1_[i].W()*tarprodmass*kappa4*dt+dMassDV[id_gas1];   
            }
            
            //char
            dMassSOLID[id_char] = -(tarprodmass*kappa5*dt+dMassSOLID[id_char]);  
        }
        
    }
    
    else if (secondGasVersion_ == 4)
    {
        // tar
        dMassDV[id_tar] = max((1.0-dt*(kappa4+kappa5))*tarprodmass, 0.0);
            

        if (dMassDV[id_tar] == 0.0)
        {
            // gas2
            forAll(volgas1_, i)
            {
                const label id_gas1 = volatileToGas1Map_[i];
                dMassDV[id_gas1] = dMassDV[id_gas1] + volgas1_[i].W()*tarprodmass*kappa4/(kappa4+kappa5);
   
            }
            
            //char
            dMassSOLID[id_char] = -(tarprodmass*kappa5/(kappa4+kappa5)+dMassSOLID[id_char]);
        }
        
        else
        {
            forAll(volgas1_, i)
            {
                const label id_gas1 = volatileToGas1Map_[i];
                dMassDV[id_gas1] = volgas1_[i].W()*tarprodmass*kappa4*dt+dMassDV[id_gas1]; 
            }
            
            //char
            dMassSOLID[id_char] = -(tarprodmass*kappa5*dt+dMassSOLID[id_char]);  
        }
        
    }

    if (done && canCombust != -1)
    {
        canCombust = 1;
    }


}


// ************************************************************************* //
