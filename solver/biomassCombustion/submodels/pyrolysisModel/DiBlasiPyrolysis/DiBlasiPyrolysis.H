/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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
    Foam::DiBlasiPyrolysis

Description
    Di Blasi devolatisation model.
    - acts on a per-specie basis
    - Rate given by Arrhenius eqn

        kappa = Ax.exp(- Ex/R.T)
		gas1	kappa1
        wood = 	tar	kappa2
		char	kappa3
		
	tar = gas2	kappa4
	      char	kappa5

    Where:
        kappa = rate constant
        A1    = activation energy (user input)
        E     = pre-exponential factor (user input)
        R     = universal gas constant
        T     = temperature

    Usage:

        DiBlasiPyrolysisCoeffs
        {
            woodProximate	(0.45 0.06 0.49) // (C H O)
            tarProximate	(0.53 0.065 0.401) // (C  H  O)
            devolKinetic1	() // (A1 E1)
            devolKinetic2	() // (A2 E2)
            devolKinetic3	() // (A3 E3)
            devolKinetic4	() // (A4 E4)
            devolKinetic5	() // (A5 E5)
            volgas1
            (
                (CH4     0.5    12/16  4/16  0.0)   // (name weight_fraction  C_mass_fraction H_mass_fraction O_mass_fraction)
                (CO2     0.5    12/44  0.0   32/44)   // (name weight_fraction C_mass_fraction H_mass_fraction O_mass_fraction)
            );
            volgas2
            (
                (CH4     0.5    12/16  4/16  0.0 A E)   // (name weight_fraction  C_mass_fraction H_mass_fraction O_mass_fraction )
                (CO2     0.5    12/44  0.0   32/44 A E)   // (name weight_fraction C_mass_fraction H_mass_fraction O_mass_fraction)
            );
            
            volatileResidualCoeff 1e-6;
        }

\*---------------------------------------------------------------------------*/

#ifndef DiBlasiPyrolysis_H
#define DiBlasiPyrolysis_H

#include "PyrolysisModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
              Class DiBlasiPyrolysis Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class DiBlasiPyrolysis
:
    public PyrolysisModel<CloudType>
{
    // Helper class to store specie-local volatile data
    class volgas1
    {
        // Private data

            //- Specie name
            word name_;

            //- Weight fraction
            scalar W_;
	    
	    //- Carbon mass fraction
	    scalar C_;
	    
	    //- Hydrogen mass fraction
	    scalar H_;
	    
	    //- Oxygen mass fraction
	    scalar O_;



    public:

        // Constructors

            //- Null constructor
            volgas1()
            :
                name_(word::null),
                W_(0.0),
                C_(0.0),
                H_(0.0),
                O_(0.0)
            {}

            //- Construct from Istream
            volgas1(Istream& is)
            :
                name_(is),
                W_(readScalar(is)),
                C_(readScalar(is)),
                H_(readScalar(is)),
                O_(readScalar(is))
            {}

//             //- Construct as copy
//             volgas1(const volgas1& vd)
//             :
//                 name_(vd.name_),
//                 W_(vd.W_),
//                 C_(vd.C_),
//                 H_(vd.H_),
//                 O_(vd.O_)
//             {}


        //- Destructor
        ~volgas1()
        {}


        // Public Member Functions

            // Access

                //- Return const access to the name
                const word& name() const
                {
                    return name_;
                }

                //- Return const access to the weight fraction
                scalar W() const
                {
                    return W_;
                }
                
                //- Return const access to the carbon mass fraction
                scalar C() const
                {
                    return C_;
                }
                
                //- Return const access to the hydrogen mass fraction
                scalar H() const
                {
                    return H_;
                }
                
                //- Return const access to the oxygen mass fraction
                scalar O() const
                {
                    return O_;
                }




        // IOstream Operators

            //- Read from Istream
            friend Istream& operator>>(Istream& is, volgas1& vd)
            {
                is.readBeginList("volgas1");
                is  >>  vd.name_  >> vd.W_ >> vd.C_ >> vd.H_ >> vd.O_;
                is.readEndList("volgas1");

                return is;
            }

            //- Write to Ostream
            friend Ostream& operator<<(Ostream& os, const volgas1& vd)
            {
                os  << token::BEGIN_LIST
                    << vd.name_ << token::SPACE
                    << vd.W_ << token::SPACE
                    << vd.C_ << token::SPACE
                    << vd.H_ << token::SPACE
                    << vd.O_
                    << token::END_LIST;

                return os;
            }
    };
    
    class volgas2
    {
        // Private data

            //- Specie name
            word name_;

            //- Weight fraction
            scalar W_;
	    
	    //- Carbon mass fraction
	    scalar C_;
	    
	    //- Hydrogen mass fraction
	    scalar H_;
	    
	    //- Oxygen mass fraction
	    scalar O_;
	    
	    //- A s-1
	    scalar A_;
	    
	    //- E J/kmol
	    scalar E_;



    public:

        // Constructors

            //- Null constructor
            volgas2()
            :
                name_(word::null),
                W_(0.0),
                C_(0.0),
                H_(0.0),
                O_(0.0),
                A_(0.0),
		E_(0.0)
            {}

            //- Construct from Istream
            volgas2(Istream& is)
            :
                name_(is),
                W_(readScalar(is)),
                C_(readScalar(is)),
                H_(readScalar(is)),
                O_(readScalar(is)),
                A_(readScalar(is)),
                E_(readScalar(is))
            {}

//             //- Construct as copy
//             volgas2(const volgas2& vd)
//             :
//                 name_(vd.name_),
//                 W_(vd.W_),
//                 C_(vd.C_),
//                 H_(vd.H_),
//                 O_(vd.O_),
//                 A_(vd.A_),
//                 E_(vd.E_)
//             {}


        //- Destructor
        ~volgas2()
        {}


        // Public Member Functions

            // Access

                //- Return const access to the name
                const word& name() const
                {
                    return name_;
                }

                //- Return const access to the weight fraction
                scalar W() const
                {
                    return W_;
                }
                
                //- Return const access to the carbon mass fraction
                scalar C() const
                {
                    return C_;
                }
                
                //- Return const access to the hydrogen mass fraction
                scalar H() const
                {
                    return H_;
                }
                
                //- Return const access to the oxygen mass fraction
                scalar O() const
                {
                    return O_;
                }

                //- Return const access to the A
                scalar A() const
                {
                    return A_;
                }

                //- Return const access to the E
                scalar E() const
                {
                    return E_;
                }                

        // IOstream Operators

            //- Read from Istream
            friend Istream& operator>>(Istream& is, volgas2& vd)
            {
                is.readBeginList("volgas2");
                is  >>  vd.name_  >> vd.W_ >> vd.C_ >> vd.H_ >> vd.O_ >> vd.A_ >> vd.E_;
                is.readEndList("volgas2");

                return is;
            }

            //- Write to Ostream
            friend Ostream& operator<<(Ostream& os, const volgas2& vd)
            {
                os  << token::BEGIN_LIST
                    << vd.name_ << token::SPACE
                    << vd.W_ << token::SPACE
                    << vd.C_ << token::SPACE
                    << vd.H_ << token::SPACE
                    << vd.O_ << token::SPACE
                    << vd.A_ << token::SPACE
                    << vd.E_  
                    << token::END_LIST;

                return os;
            }
    };


    // Private data

        // Model constants
	    
	    //- List of volgas1 data - (name weight fraction for volatile gas 1)
            List<volgas1> volgas1_;
	    
	    //- List of volgas2 data - (name weight fraction for volatile gas 2)
            List<volgas2> volgas2_;

            //- Mapping between local and cloud volatile gaseous species 1
            List<label> volatileToGas1Map_;
	    
            //- Mapping between local and cloud volatile gaseous species 2
            List<label> volatileToGas2Map_;    

            //- Volatile residual coefficient (0-1)
            //  When the fraction of volatiles are depleted below this
            //  threshold, combustion can occur
            const scalar residualCoeff_;
	    
            // version of model for secondary (0: same rate for gas + char 1: different rate for gas + char
            // 2: no second tar cracking )
	    const scalar secondGasVersion_;
	    
	    //- wood proximate analysis (C H O)
	    const vector woodProximate_;
	    
	    //- tar proximate analysis (C H O)
	    const vector tarProximate_;
	    
	    //- Kinetic parameters for the 1st pyrolysis reaction (A E)
	    const List<scalar> devolKinetic1_;
	    
	     //- Kinetic parameters for the 2nd pyrolysis reaction (A E)
	    const List<scalar> devolKinetic2_;
	    
	     //- Kinetic parameters for the 3rd pyrolysis reaction (A E)
	    const List<scalar> devolKinetic3_;
	    
	     //- Kinetic parameters for the 4th pyrolysis reaction (A E)
	    const List<scalar> devolKinetic4_;
	    
	     //- Kinetic parameters for the 5th pyrolysis reaction (A E)
	    const List<scalar> devolKinetic5_;
	 


public:

    //- Runtime type information
    TypeName("DiBlasiPyrolysis");


    // Constructors

        //- Construct from dictionary
        DiBlasiPyrolysis
        (
            const dictionary& dict,
            CloudType& owner
        );

//         //- Construct copy
//         DiBlasiPyrolysis
//         (
//             const DiBlasiPyrolysis<CloudType>& dm
//         );

        //- Construct and return a clone
        virtual autoPtr<PyrolysisModel<CloudType> > clone() const
        {
            return autoPtr<PyrolysisModel<CloudType> >
            (
                new DiBlasiPyrolysis<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~DiBlasiPyrolysis();


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
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "DiBlasiPyrolysis.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
