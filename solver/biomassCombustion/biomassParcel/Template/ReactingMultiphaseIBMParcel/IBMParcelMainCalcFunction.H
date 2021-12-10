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

template<class ParcelType>
template<class TrackCloudType>
void Foam::ReactingMultiphaseIBMParcel<ParcelType>::calc
(
    TrackCloudType& cloud,
    trackingData& td,
    const scalar dt
)
{
    typedef typename TrackCloudType::reactingCloudType reactingCloudType;
    const CompositionModel<reactingCloudType>& composition =
        cloud.composition();
    
    //- 1. Initialize particle IBM
    //-    before particle sub-timestep loop
    // ================================================

        
        // Define computation settings
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            
            // solver related parameters
            // t_step time step [s], t_end conversion time [s]
            scalar t_step = cloud.constProps().deltaTime();
            
            // tolerance tolerance level for tempearture [-], t_current current time [s]
            scalar t_current = 0, tolerance = 1.0e-07;
            
            // tolerance_Tb tolerance_Tp tolerance level during the calculation [-]
            scalarField tolerance_Tp(4, GREAT), tolerance_Tb(4, GREAT);
            
            // maxIters maximum numbers of interation for Tp [-] 
            label maxIters = 100;

            
        // Define local properties at beginning of timestep
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        
            const scalar np0 = this->nParticle_;
            const vector& U0 = this->U_;
            const scalar T0 = this->T_;
            // mass at beginning [kg]
            scalar massit = mp0_+mp1_+mp2_+mp3_;
            
            // Particle shape
            const label particleShape = cloud.constProps().parcelShape();
            scalar Xi;
            if (particleShape == 1)
            {
                Xi = cloud.constProps().xi0();
            }
            else
            {
                Xi = 1;
            }
            
            // Shrinkage factors
            const scalar drySrin = cloud.constProps().shrinkageFactorAlpha();
            const scalar devoSrin = cloud.constProps().shrinkageFactorBeta();
            const scalar charSrin = cloud.constProps().shrinkageFactorGamma();
            
            // Pressure
            const scalar pc = td.pc();

            const scalarField& YMix = this->Y_;
            const label idL = composition.idLiquid();
            const label idS = composition.idSolid();
            const label idDryWood= composition.localId(idS,"wood");
            const label idActiveDryWood= composition.localId(idS,"activeDryWood");
            const label idChar= composition.localId(idS,"C");
            const label idAsh= composition.localId(idS,"ash");
            const label idwater= composition.localId(idL,"H2O");
            
            // initial fraction
            const scalar Ywood00 = composition.YMixture0()[idS]*composition.Y0(idS)[idDryWood];
            const scalar Yash00 = composition.YMixture0()[idS]*composition.Y0(idS)[idAsh];
            const scalar Ywater00 = composition.YMixture0()[idL];
            const scalar YashDB00 = Yash00/(1.-Ywater00);
            const scalar YdryDB00 = 1.-Yash00/(1.-Ywater00);
            const scalar ratioWoodMoist = Ywood00/Ywater00;
            
            //force update fraction
            this->Y_[this->GAS] = 0.0;
            this->Y_[this->LIQ] = mp0_*Ywater00/massit;
            this->Y_[this->SLD] = 1.0 - this->Y_[this->GAS] - this->Y_[this->LIQ];
            
            // ash content in char layer
            scalar ash_inchar = ash_inchar_t_;
            scalar ash_inchar_old = ash_inchar_t_; 
            
            // limit of the mass left
            const scalar p_WoodLeft = 0.003333;
            
            // thermal properties
            // constant
            // Ste_Bol      W m-2 K-4   Stefan–Boltzmann constant
            // emissi       -       emissivity
            // Lat_Heat_wat J kg-1  Latent Heat Water
            // R_gas        J K−1 mol−1 gas constant 
            const scalar emissi = cloud.constProps().epsilon0();
            const scalar Ste_Bol = physicoChemical::sigma.value(); /*Lat_Heat_wat = 2.26E06,*/ 
            const scalar T_boiling = 373.15; 
        


        // Define variables during the calc
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            
            // just for loop
            label i;
            
            // rp particle layer diameter [m]
            // Ab particle boundary surface area [m2]
            // Vp layer volume [m3]
            // mp layer mass [kg]
            // dTp/dt [K/s]
            scalarField Vp(4), mp(4), Ab(4), dTpdt(4), cp(4), rb(4), rhop(4), rp(4), 
                        alphap(4), Tp(4), Tb(4), kp(4);
            
            // records during calculation
            scalarField mp_new(4), mp_old(4), Vp_old(4), Tp_old(4), Tb_old(4),
                        Tp_lastiter(4), Tb_lastiter(4);
            
            // equivalent diameter
            scalar equivalent_d;
            
            // variable used in the thermal drying model
            scalar Fb1Value = 0;
            
            scalarField& YGas_ = this->YGas();
            scalarField& YLiquid_ = this->YLiquid();
            scalarField& YSolid_ = this->YSolid();
            
            const label nYGas = YGas_.size();
            const label nYLiquid = YLiquid_.size();
            const label nYSolid = YSolid_.size();
            // Mass transfer due to phase change
            scalarField dMassPC(nYLiquid, 0.0);
            scalarField dMassPCSoild(nYSolid, 0.0);
            
            // Mass transfer due to devolatilisation
            scalarField dMassDV(nYGas, 0.0);
            scalarField dMassDVTemp(nYGas, 0.0);
            
            // Solid mass changes due to devolatilisation
            scalarField dMassDVSoild(nYSolid, 0.0);
            scalarField dMassDVSoildTemp(nYSolid, 0.0);

            // Change in carrier phase composition due to surface reactions
            scalarField dMassSRGas(nYGas, 0.0);
            scalarField dMassSRLiquid(nYLiquid, 0.0);
            scalarField dMassSRSolid(nYSolid, 0.0);
            scalarField dMassSRCarrier(composition.carrier().species().size(), 0.0);
            scalarField dHeatSRCarrier(composition.carrier().species().size(), 0.0);//size according to surfacereaction model 

            // Change in carrier phase composition due to surface reactions in total
            scalarField dMassSRCarrierTot(composition.carrier().species().size(), 0.0);
            
            // Zero field for mass transfer
            scalarField dMassGasZero(nYGas, 0.0);
            scalarField dMassLiquidZero(nYLiquid, 0.0);
            
            // Mass transfer in gas, liquid, and solid phase
            scalarField dMassGas(nYGas, 0.0);
            scalarField dMassLiquid(nYLiquid, 0.0);
            scalarField dMassSolid(nYSolid, 0.0);
            
            // Surface concentrations of emitted species
            scalarField Cs(composition.carrier().species().size(), 0.0);

            // Molar flux of species emitted from the particle (kmol/m^2/s)

            // Sum Ni*Cpi*Wi of emission species
            scalar NCpW = 0.0;
            
            // Calc surface values
            scalar h_coe;
            scalar Ts, rhos, mus, Prs, kappas, Res;
            scalar TInifinit = td.Tc();
            
            //Radiation
            tetIndices tetIs = this->currentTetIndices();
            scalar GcInitial = Foam::pow4(273.0)*(4.*Ste_Bol);
            
            // bed voidage, should be about to get from cloud.theta() function 
            // maybe used in the Sh and Nu corelations
            // TODO make mass transfer submodel like htc
            // (current, Sh are hard coded in CharOxidization model, which need e_bed input)
            scalar e_bed;
            // minimum bed voidage, hard coded as 0.26
            // which is equal spheres the densest packing, could be smaller for unequal spheres
            const scalar e_bedMin = 0.26;
            scalarField epsilons(2, 0.0); //[0]=e_bed, [1]=e_ashLayer
            epsilons[1] = 0.65; //hard coded ashLayer voidage, TODO
            
            // Heat transfer between layers due to mass transfer, Q = m * cp * deltaT
            scalar Qp0to1, Qp1to2, Qp2to3, Qg0to1, Qg1to2; /*Qg2to3,*/ 
            scalar mp0to1, mp1to2, mp2to3, mGas0, mGas1; 
            scalar cpG00to1_temp, cpG01to2_temp, /*cpG02to3_temp, */cpG11to2_temp; /*cpG12to3_temp,*/ 
            scalar cpP0to1_temp, cpP1to2_temp, cpP2to3_temp;
            
            //Henrik layer for devolatilisation
            scalar numberSubLayer = 5.;//TODO this value should read from the properties
            label subLayer;
            scalar deltRadiusSubLayer;
            scalar deltTemperature1;
            scalar b0Temperature1;
            scalar deltTemperature2;
            scalar b0Temperature2;    
            scalar radiusSubLayerL;
            scalar radiusSubLayerM;
            scalar radiusSubLayerH;
            scalar massSubLayer;
            scalar temperatureSubLayer;
            scalar QDevo = 0.0;
        
        
        // initial the calc of the particle IBM
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        
            //- Initialise layer variables at start, only once!!
            if (this->age_ == 0.)
            {
                
                massit = mp0_+mp1_+mp2_+mp3_;     
                //force update fraction
                this->Y_[this->GAS] = 0.0;
                this->Y_[this->LIQ] = mp0_*Ywater00/massit;
                this->Y_[this->SLD] = 1.0 - this->Y_[this->GAS] - this->Y_[this->LIQ];       
                //Force correcting solid mass fraction
                this->mass0_ = massit;         
            }
                        
            Tp[0] = this->Tp0_;
            Tp[1] = this->Tp1_;
            Tp[2] = this->Tp2_;
            Tp[3] = this->Tp3_;
            
            Tb[0] = this->Tb0_;
            Tb[1] = this->Tb1_;
            Tb[2] = this->Tb2_;
            Tb[3] = this->Tb3_;

            rb[0] = this->rb0_;
            rb[1] = this->rb1_;
            rb[2] = this->rb2_;
            rb[3] = this->rb3_;

            if (particleShape == 1)
            {
                rp[0] = R_Par_cylinderL(0.,rb[0],Xi);
                rp[1] = R_Par_cylinderL(rb[0],rb[1],Xi);
                rp[2] = R_Par_cylinderL(rb[1],rb[2],Xi);
                rp[3] = R_Par_cylinderL(rb[2],rb[3],Xi); 
                
                Ab[0] = Area_cylinderL(rb[0],Xi);
                Ab[1] = Area_cylinderL(rb[1],Xi);
                Ab[2] = Area_cylinderL(rb[2],Xi);
                Ab[3] = Area_cylinderL(rb[3],Xi);

                Vp[0] = Volume_cylinderL(rb[0],Xi);
                Vp[1] = Vol_Rin_cylinderL(rb[0],rb[1],Xi);
                Vp[2] = Vol_Rin_cylinderL(rb[1],rb[2],Xi);
                Vp[3] = Vol_Rin_cylinderL(rb[2],rb[3],Xi);
            }
            else
            {
                rp[0] = R_Par(0.,rb[0]);
                rp[1] = R_Par(rb[0],rb[1]);
                rp[2] = R_Par(rb[1],rb[2]);
                rp[3] = R_Par(rb[2],rb[3]); 
                
                Ab[0] = Area_Sph(rb[0]);
                Ab[1] = Area_Sph(rb[1]);
                Ab[2] = Area_Sph(rb[2]);
                Ab[3] = Area_Sph(rb[3]);

                Vp[0] = (4.0/3.0)*constant::mathematical::pi*Foam::pow3(rb[0]);
                Vp[1] = (4.0/3.0)*constant::mathematical::pi*(Foam::pow3(rb[1])-Foam::pow3(rb[0]));
                Vp[2] = (4.0/3.0)*constant::mathematical::pi*(Foam::pow3(rb[2])-Foam::pow3(rb[1]));
                Vp[3] = (4.0/3.0)*constant::mathematical::pi*(Foam::pow3(rb[3])-Foam::pow3(rb[2]));
            }
        
            mp[0] = this->mp0_;
            mp[1] = this->mp1_;
            mp[2] = this->mp2_;
            mp[3] = this->mp3_;
            
            //Force correcting solid mass fraction
            YSolid_[idDryWood] = (mp[0])*Ywood00/(YMix[this->SLD]*massit);
            YSolid_[idActiveDryWood] = mp[1]*composition.Y0(idS)[idDryWood]/(YMix[this->SLD]*massit);
            YSolid_[idAsh] = this->mass0_*Yash00/(YMix[this->SLD]*massit);
            YSolid_[idChar] =  1.0 - YSolid_[idDryWood] - YSolid_[idActiveDryWood] - YSolid_[idAsh];
            
            
        // Sources
        //~~~~~~~~

            // Explicit momentum source for particle
            vector Su = Zero;

            // Linearised momentum source coefficient
            scalar Spu = 0.0;

            // Momentum transfer from the particle to the carrier phase
            vector dUTrans = Zero;

            // Explicit enthalpy source for particle
            scalar Sh = 0.0;

            // Linearised enthalpy source coefficient
            scalar Sph = 0.0;

            // Sensible enthalpy transfer from the particle to the carrier phase
            scalar dhsTrans = 0.0;


    
    //- 2. Start particle IBM sub-timestep loop
    // ================================================   
            
        // sub-timestep loop
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        do
        {
            // 2.1 save data last time step
            // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            for ( i=0;i<=3;i++)
            {
                Tp_old[i] = Tp[i];
                Tb_old[i] = Tb[i];
                Tp_lastiter[i] = Tp[i];
                Tb_lastiter[i] = Tb[i];
                mp_old[i] = mp[i];
                mp_new[i] =mp[i];
                Vp_old[i] = Vp[i];            
            }
            ash_inchar_old = ash_inchar;
            //iteration number 
            label iter_0 = 0;
            scalar Vp1_temp;
            scalar Vp2_temp;
                
            // Calc surface values
            this->calcSurfaceValues(cloud, td, T0, Ts, rhos, mus, Prs, kappas);
            // Reynolds number
            Res = this->Re(rhos, U0, td.Uc(), this->d_, mus);
                
            // update e_bed
            e_bed = max(1 - (Vp[0]+Vp[1]+Vp[2]+Vp[3])/cloud.pMesh().cellVolumes()[this->cell()], e_bedMin);
            epsilons[0] = e_bed;
        
            // update radiation field
            if (cloud.radiation())
            {
                GcInitial = td.GInterp().interpolate(this->coordinates(), tetIs);
                //TODO there maybe a bug in the code 
                //- that the G field is not loaded for the first time step
                //- here to avoid GcInitial=0 in the first timestep calc
                if (GcInitial < Foam::pow4(200.0)*(4.*Ste_Bol))
                {
                    GcInitial = Foam::pow4(273.0)*(4.*Ste_Bol);
                }
            }
    
            // Sensible enthalpy transfer from the particle to the carrier phase during a particle time step   
            scalar dhsTrans_dt = 0.0;    
            
            // Molar flux of species emitted from the particle (kmol/m^2/s)
            scalar Ne = 0.0;
            
            Ne = (rDevo_ - rChar_)/(Ab[3]*t_step*55.0) + rDry_/(Ab[3]*t_step*18.0);
            
        
            // 2.2 Iteration loop, finish loop to get the results for one particle IBM sub-timestep 
            // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            do
            {
                // Molar flux of species emitted from the particle (kmol/m^2/s), temp variables duing iteration
                scalar Ne_temp = 0.0;
                
                // 2.2.1 Update particle properties at t+dt 
                // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                mp_new[0] = mp_old[0]-rDry_/Ywater00;
                mp_new[1] = mp_old[1]+rDry_*ratioWoodMoist-rDevo_/YdryDB00;
                mp_new[2] = mp_old[2]+rDevo_/YdryDB00*YashDB00+rChar_-rComb_;
                mp_new[3] = mp_old[3]+(ash_inchar_old/mp_old[2])*rComb_;  
                
                //*1.5 is hard coded to increase p_WoodLeft
                if (mp_new[0] > this->mass0_*p_WoodLeft*1.5)
                {
                    flagBoiling_ = 1;
                }

                if (mp_new[1] > this->mass0_*p_WoodLeft*1.5)
                {
                    flagDevo_ = 1;
                }
        
                //update volume, radius, area at T+dt
                Vp[0] = (Vp_old[0]/mp_old[0])*mp_new[0];
                Vp1_temp = (Vp_old[1]/mp_old[1])*(mp_old[1]-rDevo_/YdryDB00);
                Vp[1] = Vp1_temp+(1. - drySrin)*(Vp_old[0] - Vp[0]);
                Vp2_temp = (Vp_old[2]/mp_old[2])*(mp_old[2]-rComb_);
                Vp[2] = Vp2_temp + (1. - devoSrin)*(Vp_old[1] - Vp1_temp);
                Vp[3] = Vp_old[3] + (1. - charSrin)*(Vp_old[2] - Vp2_temp);
        
                if (particleShape == 1)
                {
                    rb[0] = radiusForCylinder(Xi,Vp[0]);
                    rp[0] = R_Par_cylinderL(0.0,rb[0],Xi);
                    Ab[0] = Area_cylinderL(rb[0],Xi);
                    rb[1] = radiusForCylinder(Xi,Vp[0]+Vp[1]);;
                    rp[1] = R_Par_cylinderL(rb[0],rb[1],Xi);
                    Ab[1] = Area_cylinderL(rb[1],Xi);
                    rb[2] = radiusForCylinder(Xi,Vp[0]+Vp[1]+Vp[2]);;
                    rp[2] = R_Par_cylinderL(rb[1],rb[2],Xi);
                    Ab[2] = Area_cylinderL(rb[2],Xi);
                    rb[3] = radiusForCylinder(Xi,Vp[0]+Vp[1]+Vp[2]+Vp[3]);;
                    rp[3] = R_Par_cylinderL(rb[2],rb[3],Xi);
                    Ab[3] = Area_cylinderL(rb[3],Xi);    
                }
                else
                {
                    rb[0] = Foam::cbrt(0.75*Vp[0]/constant::mathematical::pi);
                    rp[0] = R_Par(0.0,rb[0]);
                    Ab[0] = Area_Sph(rb[0]);
                    rb[1] = Foam::cbrt(0.75*(Vp[0]+Vp[1])/constant::mathematical::pi);
                    rp[1] = R_Par(rb[0],rb[1]);
                    Ab[1] = Area_Sph(rb[1]);
                    rb[2] = Foam::cbrt(0.75*(Vp[0]+Vp[1]+Vp[2])/constant::mathematical::pi);
                    rp[2] = R_Par(rb[1],rb[2]);
                    Ab[2] = Area_Sph(rb[2]);
                    rb[3] = Foam::cbrt(0.75*(Vp[0]+Vp[1]+Vp[2]+Vp[3])/constant::mathematical::pi);
                    rp[3] = R_Par(rb[2],rb[3]);
                    Ab[3] = Area_Sph(rb[3]);    
                }
                
                //update properties
                rhop[0] = mp_new[0]/Vp[0];//rho_p(0,Tp[0],Ywater00);
                rhop[1] = mp_new[1]/Vp[1];//rho_p(1,Tp[1],Ywater00);
                rhop[2] = mp_new[2]/Vp[2];//rho_p(2,Tp[2],Ywater00);
                rhop[3] = mp_new[3]/Vp[3];//rho_p(3,Tp[3],Ywater00);

                kp[0] = kp_p(0,Tp[0]);
                kp[1] = kp_p(1,Tp[1]);
                kp[2] = kp_p(2,Tp[2]);
                kp[3] = kp_p(3,Tp[3]);
                
                cp[0] = cp_p_modified(0,Tp[0],Ywater00);
                cp[1] = cp_p_modified(1,Tp[1],Ywater00);
                cp[2] = cp_p_modified(2,Tp[2],Ywater00);
                cp[3] = cp_p_modified(3,Tp[3],Ywater00);

                alphap[0] = kp[0]/(rhop[0]*cp[0]);
                alphap[1] = kp[1]/(rhop[1]*cp[1]);
                alphap[2] = kp[2]/(rhop[2]*cp[2]);
                alphap[3] = kp[3]/(rhop[3]*cp[3]);  
                
                // during iteration, tempearture could be negtive value
                // to increase the solver stadbility
                if (alphap[0] < 0)
                {
                    alphap[0] = 1.0e-7;
                }
                if (alphap[1] < 0)
                {
                    alphap[1] = 7.0e-8;
                }
                if (alphap[2] < 0)
                {
                    alphap[2] = 3.0e-7;
                }
                
                // storage Tb at last iteration  
                for ( i=0;i<=3; i++)
                {
                    Tb_lastiter[i] = Tb[i];         
                }            
                
                //Heat transfer coefficient
                h_coe = cloud.heatTransfer().htc(2.0*rb[3], Res, Prs, kappas, NCpW);
                //TODO correlations need to be added as htc submodel
                //h_coe = kappas/(this->d_)*(2.0+1.1*cbrt(Prs)*pow(Res,0.6));
                //h_coe = kappas/(this->d_)*(1.77+0.29*pow(e_bed,-0.81)*sqrt(Prs)*pow(Res,0.73));
                
                
                // 2.2.2 calculate particle heat transfer and reactions
                // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
                
                    // 2.2.2.1 calculate particle surface heat transfer
                    // QComb_ should be add to Tb[2], however, seems more stable to add here
                    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 

                    if (particleShape == 1)
                    {
                        Tb[3] = eq4_cylinderL(h_coe, emissi, Ste_Bol, kp[3], Ab[3], rb[3], rp[3], TInifinit, Tp[3], GcInitial, QComb_, Xi);
                    }
                    else
                    {
                        Tb[3] = eq4(h_coe, emissi, Ste_Bol, kp[3], Ab[3], rb[3], rp[3], TInifinit, Tp[3], GcInitial, QComb_);
                    }
                    
                    // 2.2.2.2 calculate char front heat balance
                    // first check whether it is needed to calculate heat transfer for char layer
                    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
                    if (mp_old[2] <= this->mass0_*p_WoodLeft)
                    {
                        Tb[2] = Tp[3];
                        rComb_ = 0.;
                        QComb_ = 0.;       
                    }
                    else
                    {
                        if (particleShape == 1)
                        {
                            Tb[2] = eq7_2_cylinderL(kp[2], kp[3], Ab[2], rb[2], rp[2], rp[3], Tp[2], Tp[3], 0, Xi);
                        }
                        else
                        {
                            Tb[2] = eq7_2(kp[2], kp[3], Ab[2], rb[2], rp[2], rp[3], Tp[2], Tp[3], 0);
                        }
                    
                        // calculate char conversion
                        // dMassSR should be the results of one converged iteration
                        forAll(dMassSRSolid, i)
                        {
                            dMassSRSolid[i] = 0.0;
                        }
                        forAll(dMassSRCarrier, i)
                        {
                            dMassSRCarrier[i] = 0.0;
                        }            
                        forAll(dHeatSRCarrier, i)
                        {
                            dHeatSRCarrier[i] = 0.0;
                        } 
                        
                        // 2.2.2.3 calculate char conversion
                        // Calc mass and enthalpy transfer due to surface reactions
                        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
                        calcSurfaceReactions
                        (
                            cloud,
                            td,
                            t_step,
                            rb[3],
                            rb[2],
                            Tb[2],
                            massit,
                            this->canCombust(),
                            Ne,
                            YMix,
                            YGas_,
                            YLiquid_,
                            YSolid_,
                            dMassSRGas,
                            dMassSRLiquid,
                            dMassSRSolid,
                            dMassSRCarrier,
                            Sh,
                            dhsTrans_dt,
                            QComb_,
                            Res,
                            TInifinit,
                            rhos,
                            mus,
                            Xi,
                            dHeatSRCarrier,
                            particleShape,
                            epsilons,
                            this->d_
                        );
                        
                        rComb_ = mag(dMassSRSolid[idChar]);              
                    } 
                    
                    // 2.2.2.4 calculate heat balance dry wood front (Tb1)
                    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
                    if (flagDevo_ == 1 )
                    {
                        this->canCombust() = 1;
                        
                        if (particleShape == 1)
                        {
                            Tb[1] = eq7_2_cylinderL(kp[1], kp[2], Ab[1], rb[1], rp[1], rp[2], Tp[1], Tp[2], QDevo,Xi);
                        }
                        else
                        {
                            Tb[1] = eq7_2(kp[1], kp[2], Ab[1], rb[1], rp[1], rp[2], Tp[1], Tp[2], QDevo);
                        }
                
                        // 2.2.2.5 calculate devolatilisation (with subLayer)
                        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
                        forAll(dMassDV, i)
                        {
                            dMassDV[i] = 0.0;
                        }
                        forAll(dMassDVSoild, i)
                        {
                            dMassDVSoild[i] = 0.0;
                        }     
            
                        if (Tb[1] >= cloud.constProps().TDevol() && Tp[1] >= cloud.constProps().TDevol())
                        {
                            deltRadiusSubLayer = (rb[1]-rb[0])/numberSubLayer;
                            
                            deltTemperature1 = (Tb[0]-Tp[1])/(rb[0]-rp[1]);
                            b0Temperature1 = (rb[0]*Tp[1]-rp[1]*Tb[0])/(rb[0]-rp[1]);
                            
                            deltTemperature2 = (Tp[1]-Tb[1])/(rp[1]-rb[1]);
                            b0Temperature2 = (rp[1]*Tb[1]-rb[1]*Tp[1])/(rp[1]-rb[1]);
                            
                            scalar massLayerCum = 0.0;
                            
                            for (subLayer = 1; subLayer < int(numberSubLayer); subLayer++)
                            {
                                radiusSubLayerL = rb[0] + (subLayer-1)*deltRadiusSubLayer;
                                radiusSubLayerM = radiusSubLayerL + 0.5*deltRadiusSubLayer;
                                radiusSubLayerH = radiusSubLayerL + deltRadiusSubLayer;
                                forAll(dMassDVTemp, i)
                                {
                                    dMassDVTemp[i] = 0.0;
                                }
                                forAll(dMassDVSoildTemp, i)
                                {
                                    dMassDVSoildTemp[i] = 0.0;
                                }
                                
                                if (particleShape == 1)
                                {
                                    massSubLayer = Vol_Rin_cylinderL(radiusSubLayerL, radiusSubLayerH,Xi)/Vp[1]*mp_new[1];
                                }
                                else
                                {
                                    massSubLayer = Vol_Rin(radiusSubLayerL, radiusSubLayerH)/Vp[1]*mp_new[1];
                                }
                            
                                massLayerCum = massLayerCum + massSubLayer;
                                
                                if (radiusSubLayerM<=rp[1])
                                {
                                    temperatureSubLayer = radiusSubLayerM*deltTemperature1 + b0Temperature1;
                                }
                                else
                                {
                                    temperatureSubLayer = radiusSubLayerM*deltTemperature2 + b0Temperature2;
                                }
                                
                                // hard coded avoid extrem value
                                if (temperatureSubLayer > 1000.0)
                                {
                                    temperatureSubLayer = 1000.0;
                                }

                                cloud.pyrolysis().calculate
                                (
                                    t_step,
                                    this->age_,
                                    this->mass0_,
                                    massSubLayer*YdryDB00,
                                    temperatureSubLayer,
                                    YMix[this->GAS]*YGas_,
                                    YMix[this->LIQ]*YLiquid_,
                                    YMix[this->SLD]*YSolid_, 
                                    this->canCombust(),
                                    dMassDVTemp,
                                    dMassDVSoildTemp
                                );
                                dMassDV = dMassDV + dMassDVTemp;
                                dMassDVSoild = dMassDVSoild + dMassDVSoildTemp;
                            }
                        
                            // calculate the left layer mass in the last step                    
                            forAll(dMassDVTemp, i)
                            {
                                dMassDVTemp[i] = 0.0;
                            }
                            forAll(dMassDVSoildTemp, i)
                            {
                                dMassDVSoildTemp[i] = 0.0;
                            }
                
                            massSubLayer = mp_new[1] - massLayerCum;
                            temperatureSubLayer =  Tb[1];
                            cloud.pyrolysis().calculate
                            (
                                t_step,
                                this->age_,
                                this->mass0_,
                                massSubLayer*YdryDB00,
                                temperatureSubLayer,
                                YMix[this->GAS]*YGas_,
                                YMix[this->LIQ]*YLiquid_,
                                YMix[this->SLD]*YSolid_, 
                                this->canCombust(),
                                dMassDVTemp,
                                dMassDVSoildTemp
                            );               
                            dMassDV = dMassDV + dMassDVTemp;
                            dMassDVSoild = dMassDVSoild + dMassDVSoildTemp;               
                            rDevo_ = mag(dMassDVSoild[idActiveDryWood]);
                            rChar_ = mag(dMassDVSoild[idChar]);
                        
                            QDevo = -rDevo_ * cloud.constProps().LDevol() / t_step;

                    
                            if (cloud.heatTransfer().BirdCorrection())
                            {

                                forAll(dMassDV, i)
                                {
                                    const label id = composition.localToCarrierId(this->GAS, i);
                                    const scalar W = composition.carrier().Wi(id);
                                    const scalar Ni = dMassDV[i]/(Ab[3]*t_step*W);

                                    Ne_temp += Ni;
                                }
                            }
                        }
                        else
                        {
                            
                            dMassDV = 0.;
                            dMassDVSoild = 0.;
                            rDevo_ = 0.;
                            rChar_ = 0.;
                            QDevo  = 0.;
                        }
                        
                        // no dry wood left
                        if (mp_old[1]+rDry_*ratioWoodMoist-rDevo_/YdryDB00<=this->mass0_*p_WoodLeft)
                        {
                            flagDevo_ = 0;
                        }                 
                    }
                    else //flagDevo_ = 0. no devo heat and reaction calc
                    {
                        Tb[1] = Tp[2];
                        rDevo_ = 0.;
                        rChar_ = 0.;
                        QDevo  = 0.;
                    }

                    
                    // 2.2.2.6 calculate heat balance drying front (Tb0) and drying rate
                    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                    if (flagBoiling_ == 1)
                    {
                        //updat Fb, QDry_ (heat used to evaporate water )
                        Fb1Value = Fb1(Tb[0]);
                        if (Fb1Value>1.)
                        {
                            Fb1Value = 1.;
                            Tb[0] = T_boiling;
                            if (particleShape == 1)
                            {
                                QDry_ = d_dr_cylinderL(kp[1],Ab[0],rb[0],rp[1],Xi)*(Tb[0]-Tp[1])-d_dr_cylinderL(kp[0],Ab[0],rb[0],rp[0],Xi)*(Tb[0]-Tp[0]);
                            }
                            else
                            {
                                QDry_ = d_dr(kp[1],Ab[0],rb[0],rp[1])*(Tb[0]-Tp[1])-d_dr(kp[0],Ab[0],rb[0],rp[0])*(Tb[0]-Tp[0]);
                            }

                            if (QDry_ < 0)
                            {
                                // FatalErrorInFunction
                                // <<"QDry_ < 0" << exit(FatalError);
                                QDry_ = 0;
                            }
                        }
                        else
                        {
                            if (particleShape == 1)
                            {
                                Tb[0] = eq7_3_cylinderL(kp[0], kp[1], Ab[0], rb[0], rp[0], rp[1], Tp[0], Tp[1], Fb1Value,Xi);
                            }
                            else
                            {
                                Tb[0] = eq7_3(kp[0], kp[1], Ab[0], rb[0], rp[0], rp[1], Tp[0], Tp[1], Fb1Value);
                            }
                            
                            Fb1Value = Fb1(Tb[0]);
                                            
                            if (Tb[0] >= T_boiling)
                            {
                                Fb1Value = 1.;
                                Tb[0] = T_boiling;
                            }
                            
                            if (particleShape == 1)
                            {
                                QDry_ = d_dr_cylinderL(kp[1],Ab[0],rb[0],rp[1],Xi)*(Tb[0]-Tp[1])*Fb1Value;
                            }
                            else
                            {
                                QDry_ = d_dr(kp[1],Ab[0],rb[0],rp[1])*(Tb[0]-Tp[1])*Fb1Value;
                            }
                            if (QDry_ < 0)
                            {
                                QDry_ = 0;
                            }
                        }
                
                        //Calculate drying rate
                        rDry_ = QDry_*t_step/deltaHvap(Tb[0]);

                        if (mp_old[0]-rDry_/Ywater00<=this->mass0_*p_WoodLeft)
                        {
                            rDry_ = (mp_old[0]-this->mass0_*p_WoodLeft)*Ywater00;
                            flagBoiling_ = 0;
                        } 
                        
                        if (cloud.heatTransfer().BirdCorrection())
                        {
                            Ne_temp += rDry_/(Ab[3]*t_step*18.0);
                        }
                        
                    }
                    else
                    {
                        Tb[0] = Tp[1];
                        QDry_ = 0.;
                        Fb1Value = 0.;
                        rDry_ = 0.;
                    }
                    
                    
                    Ne = Ne_temp;
                    
                    //2.2.2.7 calculate heat balance for layers
                    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   
                    // mp0to1:mass transfer from wet wood to dry wood, mp1to2:mass transfer from dry wood to char
                    mp0to1 = rDry_*ratioWoodMoist;                
                    mp1to2 = rDevo_/YdryDB00*YashDB00+rChar_;
                    mp2to3 = (ash_inchar_old/mp_old[2])*rComb_;
                    // mGas0:gas released from wet wood, mGas1:gas released from dry wood
                    mGas0 = rDry_;
                    mGas1 = rDevo_/YdryDB00-mp1to2;
                    
                    // cpG00to1_temp, heat capacity of gas from wet wood, wet wood to dry wood
                    // cpG01to2_temp, heat capacity of gas from wet wood, dry wood to char
                    cpG00to1_temp = cp_water_vapor(0.5*(Tb[0]+Tb[1]));
                    cpG01to2_temp = cp_water_vapor(0.5*(Tb[1]+Tb[2]));
                    // cpG02to3_temp = cp_water_vapor(0.5*(Tb[2]+Tb[3]));

                    cpG11to2_temp = CP_GAS;
                    // cpG12to3_temp = CP_GAS;
                    
                    // cpP0to1_temp, heat capacity of dry wood, wet wood to dry wood
                    // cpP1to2_temp, heat capacity of char, dry wood to char
                    cpP0to1_temp = cp_p_modified(1,0.5*(Tb[0]+Tp[1]),Ywater00);//T=(Tb[0]+Tp[1])/2
                    cpP1to2_temp = cp_p_modified(2,0.5*(Tb[1]+Tb[2]),Ywater00);
                    cpP2to3_temp = cp_p_modified(2,0.5*(Tb[2]+Tb[3]),Ywater00);
                    
                    // Qp0to1 Heat required when dry wood transfers from wet wood to dry wood
                    // Qp1to2 Heat required when char transfers from dry wood to char
                    Qp0to1 = mp0to1*cpP0to1_temp*(Tp[1]-Tb[0]);
                    Qp1to2 = mp1to2*cpP1to2_temp*(Tp[2]-Tb[1]);
                    Qp2to3 = mp2to3*cpP2to3_temp*(Tp[3]-Tb[2]);
                    //              
                    // Heat required when gas transfers from wet wood to drywood
                    // Heat required when gas transfers from dry wood to char
                    // Heat required when gas transfers from char to ash
                    Qg0to1 = mGas0*cpG00to1_temp*(Tb[1]-Tb[0]);
                    Qg1to2 = mGas0*cpG01to2_temp*(Tb[2]-Tb[1]) + mGas1*cpG11to2_temp*(Tb[2]-Tb[1]);
                    // Qg2to3 = mGas0*cpG02to3_temp*(Tb[3]-Tb[2]) + mGas1*cpG12to3_temp*(Tb[3]-Tb[2]);    
                    
                    //Updat dTpdt at t+dt
                    if (particleShape == 1)
                    {
                        
                        dTpdt[3] = (d_dr_cylinderL(alphap[3],Ab[3],rb[3],rp[3],Xi)*(Tb[3]-Tp[3])-d_dr_cylinderL(alphap[3],Ab[2],rb[2],rp[3],Xi)*(Tb[2]-Tp[3]))/(Vp[3])-(Qp2to3)/t_step/(rhop[3]*cp[3]*Vp[3]);
                        //eq 9
                        dTpdt[2] = (d_dr_cylinderL(alphap[2],Ab[2],rb[2],rp[2],Xi)*(Tb[2]-Tp[2])-d_dr_cylinderL(alphap[2],Ab[1],rb[1],rp[2],Xi)*(Tb[1]-Tp[2]))/(Vp[2])-(Qg1to2+Qp1to2)/t_step/(rhop[2]*cp[2]*Vp[2]);
                        //eq 10
                        dTpdt[1] = (d_dr_cylinderL(alphap[1],Ab[1],rb[1],rp[1],Xi)*(Tb[1]-Tp[1])-d_dr_cylinderL(alphap[1],Ab[0],rb[0],rp[1],Xi)*(Tb[0]-Tp[1]))/(Vp[1])-(Qg0to1+Qp0to1)/t_step/(rhop[1]*cp[1]*Vp[1]);
                        //eq 11-1
                        dTpdt[0] = d_dr_cylinderL((kp[0]/(cp[0]*rhop[0])),Ab[0],rb[0],rp[0],Xi)*(Tb[0]-Tp[0])/(Vp[0]);  
                    }
                    else
                    {
                    
                        dTpdt[3] = (d_dr(alphap[3],Ab[3],rb[3],rp[3])*(Tb[3]-Tp[3])-d_dr(alphap[3],Ab[2],rb[2],rp[3])*(Tb[2]-Tp[3]))/(Vp[3])-(Qp2to3)/t_step/(rhop[3]*cp[3]*Vp[3]);
                        //eq 9
                        dTpdt[2] = (d_dr(alphap[2],Ab[2],rb[2],rp[2])*(Tb[2]-Tp[2])-d_dr(alphap[2],Ab[1],rb[1],rp[2])*(Tb[1]-Tp[2]))/(Vp[2])-(Qg1to2+Qp1to2)/t_step/(rhop[2]*cp[2]*Vp[2]);
                        //eq 10
                        dTpdt[1] = (d_dr(alphap[1],Ab[1],rb[1],rp[1])*(Tb[1]-Tp[1])-d_dr(alphap[1],Ab[0],rb[0],rp[1])*(Tb[0]-Tp[1]))/(Vp[1])-(Qg0to1+Qp0to1)/t_step/(rhop[1]*cp[1]*Vp[1]);
                        //eq 11-1
                        dTpdt[0] = d_dr((kp[0]/(cp[0]*rhop[0])),Ab[0],rb[0],rp[0])*(Tb[0]-Tp[0])/(Vp[0]); 
                    }
                    
                //2.2.3 particle IBM iteration cetia update
                //update Tp and tolerance of Tp and Tb 
                // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                for ( i=0;i<=3; i++)
                {
                    Tp_lastiter[i]=Tp[i];
                    Tp[i] = Tp_old[i] + dTpdt[i]*t_step;       
                    tolerance_Tp[i] = mag((Tp_lastiter[i] - Tp[i])/Tp[i]);
                    tolerance_Tb[i] = mag((Tb_lastiter[i] - Tb[i])/Tb[i]);  
                    
                }   
                if (flagBoiling_ == 0)
                {
                    Tb[0]=Tp[1];
                    Tp[0]=Tp[1];
                    tolerance_Tp[0] = VSMALL;
                    tolerance_Tb[0] = VSMALL;
                }
                if (flagDevo_ == 0)
                {
                    Tb[1]=Tp[2];
                    Tp[1]=Tp[2];
                    tolerance_Tp[1] = VSMALL;
                    tolerance_Tb[1] = VSMALL;
                }
                if (mp_old[2] <= this->mass0_*p_WoodLeft)
                {
                    Tb[2]=Tp[3];
                    Tp[2]=Tp[3];
                    tolerance_Tb[2] = VSMALL;
                    tolerance_Tp[2] = VSMALL;
                }
                
                if (flagBoiling_ == 1  && (mp_old[1]+rDry_*ratioWoodMoist-rDevo_/YdryDB00 > this->mass0_*p_WoodLeft))
                {
                    flagDevo_ = 1;
                }
            
                iter_0 = iter_0+1;
            
            //2.2 particle IBM iteration finished for one sub-timestep 
            }while (((tolerance_Tb[0] > tolerance || tolerance_Tb[1] > tolerance || tolerance_Tb[2] > tolerance || tolerance_Tb[3] > tolerance  || tolerance_Tp[0] > tolerance || tolerance_Tp[1] > tolerance || tolerance_Tp[2] > tolerance  || tolerance_Tp[3] > tolerance)  && iter_0 < maxIters) || iter_0<2);
                    
            // 2.3 update particle properties from last sub-timestep 
            // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            
            //- Layer mass, radius, volume
            mp[0] = mp_old[0]-rDry_/Ywater00; 
            mp[1] = mp_old[1]+rDry_*ratioWoodMoist-rDevo_/YdryDB00;
            mp[2] = mp_old[2]+rDevo_/YdryDB00*YashDB00+rChar_-rComb_;
            ash_inchar = ash_inchar_old + rDevo_/YdryDB00*YashDB00 - (ash_inchar_old/mp_old[2])*rComb_;
            mp[3] = mp_old[3]+(ash_inchar_old/mp_old[2])*rComb_; 
            
            Vp[0] = (Vp_old[0]/mp_old[0])*mp[0];
            Vp1_temp = (Vp_old[1]/mp_old[1])*(mp_old[1]-rDevo_/YdryDB00);
            Vp[1] = Vp1_temp+(1. - drySrin)*(Vp_old[0] - Vp[0]);
            Vp2_temp = (Vp_old[2]/mp_old[2])*(mp_old[2]-rComb_);
            Vp[2] = Vp2_temp + (1. - devoSrin)*(Vp_old[1] - Vp1_temp);
            Vp[3] = Vp_old[3] + (1. - charSrin)*(Vp_old[2] - Vp2_temp);
        
            if (particleShape == 1)
            {
                rb[0] = radiusForCylinder(Xi,Vp[0]);
                rp[0] = R_Par_cylinderL(0.0,rb[0],Xi);
                Ab[0] = Area_cylinderL(rb[0],Xi);
                rb[1] = radiusForCylinder(Xi,Vp[0]+Vp[1]);
                rp[1] = R_Par_cylinderL(rb[0],rb[1],Xi);
                Ab[1] = Area_cylinderL(rb[1],Xi);
                rb[2] = radiusForCylinder(Xi,Vp[0]+Vp[1]+Vp[2]);
                rp[2] = R_Par_cylinderL(rb[1],rb[2],Xi);
                Ab[2] = Area_cylinderL(rb[2],Xi);
                rb[3] = radiusForCylinder(Xi,Vp[0]+Vp[1]+Vp[2]+Vp[3]);
                rp[3] = R_Par_cylinderL(rb[2],rb[3],Xi);
                Ab[3] = Area_cylinderL(rb[3],Xi);
            }
            else
            {
                rb[0] = Foam::cbrt(0.75*(Vp[0])/constant::mathematical::pi);
                rp[0] = R_Par(0.0,rb[0]);
                Ab[0] = Area_Sph(rb[0]);
                rb[1] = Foam::cbrt(0.75*(Vp[0]+Vp[1])/constant::mathematical::pi);
                rp[1] = R_Par(rb[0],rb[1]);
                Ab[1] = Area_Sph(rb[1]);
                rb[2] = Foam::cbrt(0.75*(Vp[0]+Vp[1]+Vp[2])/constant::mathematical::pi);
                rp[2] = R_Par(rb[1],rb[2]);
                Ab[2] = Area_Sph(rb[2]);
                rb[3] = Foam::cbrt(0.75*(Vp[0]+Vp[1]+Vp[2]+Vp[3])/constant::mathematical::pi);
                rp[3] = R_Par(rb[2],rb[3]);
                Ab[3] = Area_Sph(rb[3]); 
            
            }

            //- Mass and fraction update
            //- cumulate source terms from sub-timestep
            
            //drying
            dMassPC[idwater] = rDry_;
            dMassPCSoild[idDryWood] = rDry_*ratioWoodMoist;
            dMassPCSoild[idActiveDryWood] = -rDry_*ratioWoodMoist;         
            dMassLiquid = dMassLiquid + dMassPC;      
            massit = updateMassFractionsProtected(massit, dMassGasZero, dMassPC, dMassPCSoild);
            
            //Devolatilisation
            dMassGas = dMassGas + dMassDV;
            massit = updateMassFractionsProtected(massit, dMassGasZero, dMassLiquidZero, dMassDVSoild);

            //Combustion
            dMassGas = dMassGas + dMassSRGas;
            dMassLiquid = dMassLiquid + dMassSRLiquid;
            dMassSolid = dMassSolid + dMassSRSolid;
            massit = updateMassFractionsProtected(massit, dMassGasZero, dMassLiquidZero, dMassSolid);
            dMassSRCarrierTot = dMassSRCarrierTot + dMassSRCarrier;
            
            //- Add to cumulative phase change mass
            cloud.phaseChange().addToPhaseChangeMass(np0*rDry_);//need to be corrected
            cloud.pyrolysis().addToPyrolysisMass(np0*rDevo_);

            Sph += t_step*h_coe*Ab[3];
            
            dhsTrans_dt =  t_step*h_coe*Ab[3]*((Tb[3]+Tb_old[3])/2.0 - TInifinit);

            dhsTrans = dhsTrans + dhsTrans_dt;
            
            this->T_ = /*Tb[2]*/Tb[3];
        
            //- Update time
            t_current = t_current + t_step;
            cumTime_ = cumTime_ + t_step;
                
            if (t_current + t_step >= dt)
            {
                t_step = dt - t_current;
            }
        
        // 2. sub-timestep loop of the particle IBM finisded. IBM finished
        }while (t_current<dt);

    
    //- 3. Update parcel properties from IBM calculation for one fluid timestep
    // ================================================                  
                    
        // equivalent diameter for cylinder and sphericity
        if (particleShape == 1)
        {
            equivalent_d = Foam::cbrt(6*(Vp[0]+Vp[1]+Vp[2]+Vp[3])/constant::mathematical::pi);
            this->d_ = equivalent_d;
        }
        else
        {
            this->d_ = 2*rb[3];
        }
        
        //- Update parcel class value
        this->Tp0_ = Tp[0];
        this->Tp1_ = Tp[1];
        this->Tp2_ = Tp[2];
        this->Tp3_ = Tp[3];
        this->Tb0_ = Tb[0];
        this->Tb1_ = Tb[1];
        this->Tb2_ = Tb[2];
        this->Tb3_ = Tb[3];
        this->rb0_ = rb[0];
        this->rb1_ = rb[1];
        this->rb2_ = rb[2];
        this->rb3_ = rb[3];
        this->mp0_ = mp[0];
        this->mp1_ = mp[1];
        this->mp2_ = mp[2];
        this->mp3_ = mp[3];
        this->T_ = Tb[3];
        this->rho_ = (mp[0]+mp[1]+mp[2]+mp[3])/(Vp[0]+Vp[1]+Vp[2]+Vp[3]);
        this->ash_inchar_t_ = ash_inchar;

        //- Correct surface values due to emitted species
        this->correctSurfaceValues(cloud, td, Ts, Cs, rhos, mus, Prs, kappas);   

        //- Calculate new particle velocity 
        // Reynolds number
        Res = this->Re(rhos, U0, td.Uc(), this->d_, mus);
        
        this->U_ = this->calcVelocity(cloud, td, dt, Res, mus, mp[0]+mp[1]+mp[2]+mp[3], Su, dUTrans, Spu);
                    
    // TODO Remove the particle when mass falls below minimum threshold
    
        
    //- 4.Accumulate carrier phase source terms
    // ================================================      
    
    if (cloud.solution().coupled())
    {
        // Transfer mass lost to carrier mass, momentum and enthalpy sources
        forAll(YGas_, i)
        {
            scalar dm = np0*dMassGas[i]; //in fact this is from devo in this code
            label gid = composition.localToCarrierId(this->GAS, i);
            scalar hs = composition.carrier().Hs(gid, pc, Tb2_); //Tp may cause some error
            cloud.rhoTrans(gid)[this->cell()] += dm;
            cloud.UTrans()[this->cell()] += dm*U0;
            cloud.hsTrans()[this->cell()] += dm*hs;
        }
        forAll(dMassSRCarrierTot, i) //Combustion
        {
            scalar dm = np0*dMassSRCarrierTot[i];
            scalar hs = composition.carrier().Hs(i, pc, Tb2_); //Tp may cause some error
            cloud.rhoTrans(i)[this->cell()] += dm;
            cloud.UTrans()[this->cell()] += dm*U0;
            cloud.hsTrans()[this->cell()] += dm*hs;
        }
        forAll(YLiquid_, i)
        {
            scalar dm = np0*dMassLiquid[i]; //in fact this is from drying in this code
            label gid = composition.localToCarrierId(this->LIQ, i);
            scalar hs = composition.carrier().Hs(gid, pc, T0);
            cloud.rhoTrans(gid)[this->cell()] += dm;
            cloud.UTrans()[this->cell()] += dm*U0;
            cloud.hsTrans()[this->cell()] += dm*hs;
        }

        // Update momentum transfero
        cloud.UTrans()[this->cell()] += np0*dUTrans;
        cloud.UCoeff()[this->cell()] += np0*Spu;

        // Update sensible enthalpy transfer
        cloud.hsTrans()[this->cell()] += np0*dhsTrans;
        cloud.hsCoeff()[this->cell()] += np0*Sph;

        // Update radiation fields
        if (cloud.radiation())
        {
            scalar ap;
            if (particleShape == 1)
            {
                ap = constant::mathematical::pi*rb[3]*(3*rb[3]-Xi)/2;;
            }
            else
            {
                ap = this->areaP();
            }
            const scalar T4 = pow4(Tb3_);
            cloud.radAreaP()[this->cell()] += dt*np0*ap;
            cloud.radT4()[this->cell()] += dt*np0*T4;
            cloud.radAreaPT4()[this->cell()] += dt*np0*ap*T4;
        }
    }
}

// ************************************************************************* //
