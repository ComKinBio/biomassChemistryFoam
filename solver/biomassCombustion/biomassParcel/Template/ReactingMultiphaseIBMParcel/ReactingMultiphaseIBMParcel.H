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
    Foam::ReactingMultiphaseIBMParcel

Description
    Multiphase variant of the reacting parcel class with one/two-way coupling
    with the continuous phase.

SourceFiles
    ReactingMultiphaseIBMParcelI.H
    ReactingMultiphaseIBMParcel.C
    ReactingMultiphaseIBMParcelIO.C

\*---------------------------------------------------------------------------*/

#ifndef ReactingMultiphaseIBMParcel_H
#define ReactingMultiphaseIBMParcel_H

#include "particle.H"
#include "SLGThermo.H"
#include "demandDrivenEntry.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class ParcelType>
class ReactingMultiphaseIBMParcel;

template<class ParcelType>
Ostream& operator<<
(
    Ostream&,
    const ReactingMultiphaseIBMParcel<ParcelType>&
);

/*---------------------------------------------------------------------------*\
                 Class ReactingMultiphaseIBMParcel Declaration
\*---------------------------------------------------------------------------*/

template<class ParcelType>
class ReactingMultiphaseIBMParcel
:
    public ParcelType
{
    // Private Data

        //- Size in bytes of the fields
        static const std::size_t sizeofFields_;


public:

    //- Class to hold reacting multiphase particle constant properties
    class constantProperties
    :
        public ParcelType::constantProperties
    {
        // Private Data
        
        //- parcel shape, 1 = cylinder, else = sphere
        demandDrivenEntry<label> parcelShape_;
        demandDrivenEntry<scalar> xi0_;
        
            
        //- shrinkageMode factor
        demandDrivenEntry<scalar> shrinkageFactorAlpha_;
        demandDrivenEntry<scalar> shrinkageFactorBeta_;
        demandDrivenEntry<scalar> shrinkageFactorGamma_;
        

        //- Layer initial properties
        
            //- Front temperature
            demandDrivenEntry<scalar> Tb00_;
            demandDrivenEntry<scalar> Tb10_;
            demandDrivenEntry<scalar> Tb20_;
            demandDrivenEntry<scalar> Tb30_;
            
            //- Front radius
            demandDrivenEntry<scalar> rb00_;
            demandDrivenEntry<scalar> rb10_;
            demandDrivenEntry<scalar> rb20_;    
            demandDrivenEntry<scalar> rb30_;
            
            //- Layer temperature
            demandDrivenEntry<scalar> Tp00_;
            demandDrivenEntry<scalar> Tp10_;     
            demandDrivenEntry<scalar> Tp20_;  
            demandDrivenEntry<scalar> Tp30_;
            
            //- Particle model sub-timestep
            demandDrivenEntry<scalar> deltaTime_;
            
            //- Heat coefficient
            //- (ratio between coupled to particle and carrier phase)
            demandDrivenEntry<scalar> hCoeff_;  


    public:

        // Constructors

            //- Null constructor
            constantProperties();

            //- Copy constructor
            constantProperties(const constantProperties& cp);

            //- Construct from dictionary
            constantProperties(const dictionary& parentDict);


        // Access
            
            //- Return const access to parcelShape
            inline label parcelShape() const;
            inline scalar xi0() const;
            
            //- Return const access to shrinkageMode factor
            inline scalar shrinkageFactorAlpha() const;
            inline scalar shrinkageFactorBeta() const;
            inline scalar shrinkageFactorGamma() const;
            
            
            //- Return const access to layer inital properties
               
            //- Return const access to front temperature
                inline scalar Tb00() const;
                inline scalar Tb10() const;
                inline scalar Tb20() const;
                inline scalar Tb30() const;
                
                //- Return const access to front radius
                inline scalar rb00() const;
                inline scalar rb10() const;
                inline scalar rb20() const;
                inline scalar rb30() const;
                
                //- Return const access to layer temperature
                inline scalar Tp00() const;
                inline scalar Tp10() const;
                inline scalar Tp20() const;
                inline scalar Tp30() const;

                //- Return const access to particle model sub-timestep
                inline scalar deltaTime() const;

                //- Return const access toheat transfer coefficient
                inline scalar hCoeff() const; 
    };


    //- Use base tracking data
    typedef typename ParcelType::trackingData trackingData;
    
    //- Update the mass fractions (Y, YGas, YLiquid, YSolid)
    scalar updateMassFractions
    (
        const scalar mass0,
        const scalarField& dMassGas,
        const scalarField& dMassLiquid,
        const scalarField& dMassSolid
    );

protected:

    // Protected data

        // Parcel properties

            // Thermally thick particle properties
            //- Tb_ tempearture [K], rb_ radius [m]
            //- Front or inter-face
            scalar Tb0_, Tb1_, Tb2_, Tb3_, rb0_, rb1_, rb2_, rb3_;
            
            // Tp_ tempearture [K], rp_ mass-weighted radius [m], layer mass [k]
            //- Layer
            scalar Tp0_, Tp1_, Tp2_, Tp3_, mp0_, mp1_, mp2_, mp3_;
            
            //- Mass loss/gain due to drying and devolatilization and char combustion (kg)
            scalar rDry_, rDevo_, rChar_, rComb_, ash_inchar_t_;
            
            //- Heat flux due to drying and combustion [W]
            scalar QDry_, QComb_;
            
            //- Flag for boiling and devolatilisation
            //- 1 = true
            //-( canCombust_ flage )
            //-( 0 = can devolatilise, cannot combust but can change )
            //-( 1 = can devolatilise, can combust )
            //-(-1 = cannot devolatilise or combust, and cannot change )
            label flagBoiling_, flagDevo_;
            
            //- Accumulated time [s] (used for particle model sub-loop)
            scalar cumTime_;

    // Protected Member Functions
            
        //- Update the mass fractions (Y, YGas, YLiquid, YSolid)
        scalar updateMassFractionsProtected
        (
            const scalar mass0,
            const scalarField& dMassGas,
            const scalarField& dMassLiquid,
            const scalarField& dMassSolid
        );

        //- Calculate Devolatilisation
        // coded into calc
        
        //- Calculate surface reactions
        template<class TrackCloudType>
        void calcSurfaceReactions
        (
            TrackCloudType& cloud,
            trackingData& td,
            const scalar dt,           // timestep
            const scalar d,            // radius
            const scalar di,            // inner layer radius
            const scalar T,            // temperature
            const scalar mass,         // mass
            const label canCombust,     // 'can combust' flag
            const scalar N,            // flux of species emitted from particle
            const scalarField& YMix,   // mixture mass fractions
            const scalarField& YGas,   // gas-phase mass fractions
            const scalarField& YLiquid,// liquid-phase mass fractions
            const scalarField& YSolid, // solid-phase mass fractions
            scalarField& dMassSRGas,   // gas-phase mass transfer - local
            scalarField& dMassSRLiquid,// liquid-phase mass transfer - local
            scalarField& dMassSRSolid, // solid-phase mass transfer - local
            scalarField& dMassSRCarrier, // carrier phase mass transfer
            scalar& Sh,                // explicit particle enthalpy source
            scalar& dhsTrans,           // sensible enthalpy transfer to carrier
            scalar& QComb,               // record combustion heat
            const scalar Re,            // Re number
            const scalar Tc,            // continue phase temperature
            const scalar rhoc,          // continue phase density
            const scalar muc,              // continue phase viscosity
            const scalar Xi,            //cylinder shape const
            scalarField& dHeatSRCarrier,   //SR heat
            const label particleShape,
            const scalarField& epsilons,     //bed and ashlayer void
            const scalar deq          //cylinder eq d
        ) const;


public:

    // Static Data Members

        //- Runtime type information
        TypeName("ReactingMultiphaseIBMParcel");

        //- String representation of properties
        AddToPropertyList
        (
            ParcelType,
            " Tb0"
          + " Tb1"
          + " Tb2"
          + " Tb3"
          + " rb0"
          + " rb1"
          + " rb2"
          + " rb3"
          + " Tp0"
          + " Tp1"
          + " Tp2"
          + " Tp3"  
          + " mp0"
          + " mp1"
          + " mp2"
          + " mp3" 
          + " rDry"
          + " rDevo"
          + " rChar"
          + " rComb"
          + " ash_inchar_t"
          + " QDry"  
          + " QComb"
          + " flagBoiling"
          + " flagDevo"
          + " cumTime" 
        );


    // Constructors

        //- Construct from mesh, position and topology
        //  Other properties initialised as null
        inline ReactingMultiphaseIBMParcel
        (
            const polyMesh& mesh,
            const barycentric& coordinates,
            const label celli,
            const label tetFacei,
            const label tetPti
        );

        //- Construct from a position and a cell, searching for the rest of the
        //  required topology. Other properties are initialised as null.
        inline ReactingMultiphaseIBMParcel
        (
            const polyMesh& mesh,
            const vector& position,
            const label celli
        );

        //- Construct from components
        inline ReactingMultiphaseIBMParcel
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
        );

        //- Construct from Istream
        ReactingMultiphaseIBMParcel
        (
            const polyMesh& mesh,
            Istream& is,
            bool readFields = true
        );

        //- Construct as a copy
        ReactingMultiphaseIBMParcel(const ReactingMultiphaseIBMParcel& p);

        //- Construct as a copy
        ReactingMultiphaseIBMParcel
        (
            const ReactingMultiphaseIBMParcel& p,
            const polyMesh& mesh
        );

        //- Construct and return a (basic particle) clone
        virtual autoPtr<particle> clone() const
        {
            return autoPtr<particle>(new ReactingMultiphaseIBMParcel(*this));
        }

        //- Construct and return a (basic particle) clone
        virtual autoPtr<particle> clone(const polyMesh& mesh) const
        {
            return autoPtr<particle>(new ReactingMultiphaseIBMParcel(*this, mesh));
        }

        //- Factory class to read-construct particles used for
        //  parallel transfer
        class iNew
        {
            const polyMesh& mesh_;

        public:

            iNew(const polyMesh& mesh)
            :
                mesh_(mesh)
            {}

            autoPtr<ReactingMultiphaseIBMParcel<ParcelType>> operator()
            (
                Istream& is
            ) const
            {
                return autoPtr<ReactingMultiphaseIBMParcel<ParcelType>>
                (
                    new ReactingMultiphaseIBMParcel<ParcelType>(mesh_, is, true)
                );
            }
        };


    // Member Functions

        // Access

            //- Return const access to Parcel properties
        
                //- Return const access to front tempearture
                inline scalar Tb0() const;
                inline scalar Tb1() const;
                inline scalar Tb2() const;
                inline scalar Tb3() const;
                
                //- Return const access to front radius
                inline scalar rb0() const;
                inline scalar rb1() const;
                inline scalar rb2() const;
                inline scalar rb3() const;
                
                //- Return const access to layer tempearture
                inline scalar Tp0() const;
                inline scalar Tp1() const;
                inline scalar Tp2() const;
                inline scalar Tp3() const;
                
                //- Return const access to layer mass
                inline scalar mp0() const;
                inline scalar mp1() const;
                inline scalar mp2() const;
                inline scalar mp3() const;
                
                //- Return const access to conversion rate
                inline scalar rDry() const;
                inline scalar rDevo() const;
                inline scalar rChar() const;
                inline scalar rComb() const;
                inline scalar ash_inchar_t() const;
                
                //- Return const access to conversion heat
                inline scalar QDry() const; 
                inline scalar QComb() const; 
                
                //- Return const access to conversion flag
                inline label flagBoiling() const;
                inline label flagDevo() const;
                
                //- Return const access to cummulated time
                inline scalar cumTime() const;

        // Edit

            //- Return access to Parcel properties
                
                //- Return access to front tempearture
                inline scalar& Tb0();
                inline scalar& Tb1();
                inline scalar& Tb2();
                inline scalar& Tb3();
                
                //- Return access to front radius
                inline scalar& rb0();
                inline scalar& rb1();
                inline scalar& rb2();
                inline scalar& rb3();
                
                //- Return access to layer tempearture
                inline scalar& Tp0();
                inline scalar& Tp1();
                inline scalar& Tp2();
                inline scalar& Tp3();
                
                //- Return access to layer mass
                inline scalar& mp0();
                inline scalar& mp1();
                inline scalar& mp2();
                inline scalar& mp3();
                
                //- Return access to conversion rate
                inline scalar& rDry();
                inline scalar& rDevo();
                inline scalar& rChar();
                inline scalar& rComb();
                inline scalar& ash_inchar_t();
                
                //- Return access to conversion heat
                inline scalar& QDry();
                inline scalar& QComb();
                
                //- Return access to conversion flag
                inline label& flagBoiling();
                inline label& flagDevo();
                
                //- Return access to cummulated time
                inline scalar& cumTime();

        
        // Helper functions
        
            //- get radius of circle surface from the volume of cylinderL
            scalar radiusForCylinder
            (
                const scalar Xi, 
                const scalar V0
            );
        
            //- Particle surface area sphere
            scalar Area_Sph
            (
                const scalar radius
            );
            
            
            //- Particle surface area cylinderL
            scalar Area_cylinderL
            (
                const scalar radius,
                const scalar Xi
            );
            
            //- Particle volume cylinderL
            scalar Volume_cylinderL
            (
                const scalar radius,
                const scalar Xi
            );
            
            //- Shell volume sphere
            scalar Vol_Rin 
            (
                const scalar rin,
                const scalar rout
            );
            
            //- Shell volume cylinderL
            scalar Vol_Rin_cylinderL
            (
                const scalar rin,
                const scalar rout,
                const scalar Xi
            );
            
            //- Radius of the shell sphere
            scalar R_Par
            (
                const scalar rin,
                const scalar rout
            );
            
            //- Radius of the shell cylinderL
            scalar R_Par_cylinderL
            (
                const scalar rin,
                const scalar rout,
                const scalar Xi
            );
            
            //- Compute A*kp*d/dr sphere
            scalar d_dr
            (
                const scalar kp,
                const scalar Ab,
                const scalar ri,
                const scalar rj
            );
            
            //- Compute A*kp*d/dr cylinderL
            scalar d_dr_cylinderL
            (
                const scalar kp,
                const scalar Ab,
                const scalar ri,
                const scalar rj,
                const scalar Xi
            );
            
            //- Function Fb1 based on Antoine expression 
            scalar Fb1 
            (
                const scalar Tb1
            );            
            
            //- eq4 Henrik C&F
            scalar eq4
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
            );
            
            //- Heat balance over the boundary of the layer
            //- Derived from eq7 Henrik C&F
            scalar eq7_2
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
            );
            
            scalar eq7_3
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
            );
            
            //- eq4 Henrik C&F for cylinderL
            scalar eq4_cylinderL
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
            );
            
            //- Heat balance over the boundary of the layer
            //- Derived from eq7 Henrik C&F for cylinderL
            scalar eq7_2_cylinderL
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
            );
            
            scalar eq7_3_cylinderL
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
            );
            
            //- Heat capacity of the soild
            scalar cp_p
            (
                const label layer,
                const scalar Tp,
                const scalar moist_WB
            );
            
            scalar cp_p_modified
            (
                const label layer,
                const scalar Tp,
                const scalar moist_WB
            );
            
            //- Density of the soild
            scalar rho_p
            (
                const label layer,
                const scalar Tp,
                const scalar moist_WB
            );
            
            //- Conductivity of the soild
            scalar kp_p
            (
                const label layer,
                const scalar Tp
            );
            
            //- Heat capacity of the water vapor
            scalar cp_water_vapor
            (
                const scalar T
            );
            
            //- Evaporation Heat
            scalar deltaHvap
            (
                const scalar T
            );


        // Main calculation loop

            //- Set cell values
            template<class TrackCloudType>
            void setCellValues(TrackCloudType& cloud, trackingData& td);

            //- Correct cell values using latest transfer information
            template<class TrackCloudType>
            void cellValueSourceCorrection
            (
                TrackCloudType& cloud,
                trackingData& td,
                const scalar dt
            );

            //- Update parcel properties over the time interval
            template<class TrackCloudType>
            void calc
            (
                TrackCloudType& cloud,
                trackingData& td,
                const scalar dt
            );


        // I-O

            //- Read
            template<class CloudType, class CompositionType>
            static void readFields
            (
                CloudType& c,
                const CompositionType& compModel
            );

            //- Read - no composition
            template<class CloudType>
            static void readFields(CloudType& c);

            //- Write
            template<class CloudType, class CompositionType>
            static void writeFields
            (
                const CloudType& c,
                const CompositionType& compModel
            );

            //- Read - composition supplied
            template<class CloudType>
            static void writeFields(const CloudType& c);


    // Ostream Operator

        friend Ostream& operator<< <ParcelType>
        (
            Ostream&,
            const ReactingMultiphaseIBMParcel<ParcelType>&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "ReactingMultiphaseIBMParcelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "ReactingMultiphaseIBMParcel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
