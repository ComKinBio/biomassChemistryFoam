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
#include "IOstreams.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<class ParcelType>
Foam::string Foam::ReactingMultiphaseIBMParcel<ParcelType>::propertyList_ =
    Foam::ReactingMultiphaseIBMParcel<ParcelType>::propertyList();

template<class ParcelType>
const std::size_t Foam::ReactingMultiphaseIBMParcel<ParcelType>::sizeofFields_
(
    24*sizeof(scalar) + 2*sizeof(label)
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
// 
template<class ParcelType>
Foam::ReactingMultiphaseIBMParcel<ParcelType>::ReactingMultiphaseIBMParcel
(
    const polyMesh& mesh,
    Istream& is,
    bool readFields
)
:
    ParcelType(mesh, is, readFields),
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
{
    if (readFields)
    {
        if (is.format() == IOstream::ASCII)
        {   

            is >>  Tb0_ >> Tb1_ >> Tb2_ >> Tb3_ >> rb0_ >> rb1_ >> rb2_ >> rb3_ 
               >> Tp0_ >> Tp1_ >> Tp2_   >> Tp3_ >> mp0_ >> mp1_ >> mp2_ >> mp3_ 
               >> rDry_ >> rDevo_ >> rChar_ >> rComb_ >> ash_inchar_t_ >> QDry_ >> QComb_ >> flagBoiling_ >> flagDevo_ >> cumTime_;
        }
        else
        {
            is.read(reinterpret_cast<char*>(&Tb0_), sizeofFields_);
        }
    }

    //Check state of Istream
    is.check
    (
        "ReactingMultiphaseIBMParcel<ParcelType>::ReactingMultiphaseIBMParcel"
        "("
            "const polyMesh&, "
            "Istream&, "
            "bool"
        ")"
    );
}


template<class ParcelType>
template<class CloudType>
void Foam::ReactingMultiphaseIBMParcel<ParcelType>::readFields(CloudType& c)
{
    ParcelType::readFields(c);
}


template<class ParcelType>
template<class CloudType, class CompositionType>
void Foam::ReactingMultiphaseIBMParcel<ParcelType>::readFields
(
    CloudType& c,
    const CompositionType& compModel
)
{
    bool valid = c.size();

    ParcelType::readFields(c, compModel);

    IOField<scalar> Tb0(c.fieldIOobject("Tb0", IOobject::MUST_READ), valid);
    c.checkFieldIOobject(c, Tb0);
    
    IOField<scalar> Tb1(c.fieldIOobject("Tb1", IOobject::MUST_READ), valid);
    c.checkFieldIOobject(c, Tb1);

    IOField<scalar> Tb2(c.fieldIOobject("Tb2", IOobject::MUST_READ), valid);
    c.checkFieldIOobject(c, Tb2);

    IOField<scalar> Tb3(c.fieldIOobject("Tb3", IOobject::MUST_READ), valid);
    c.checkFieldIOobject(c, Tb3);

    IOField<scalar> rb0(c.fieldIOobject("rb0", IOobject::MUST_READ), valid);
    c.checkFieldIOobject(c, rb0);

    IOField<scalar> rb1(c.fieldIOobject("rb1", IOobject::MUST_READ), valid);
    c.checkFieldIOobject(c, rb1);

    IOField<scalar> rb2(c.fieldIOobject("rb2", IOobject::MUST_READ), valid);
    c.checkFieldIOobject(c, rb2);

    IOField<scalar> rb3(c.fieldIOobject("rb3", IOobject::MUST_READ), valid);
    c.checkFieldIOobject(c, rb3);

    IOField<scalar> Tp0(c.fieldIOobject("Tp0", IOobject::MUST_READ), valid);
    c.checkFieldIOobject(c, Tp0);

    IOField<scalar> Tp1(c.fieldIOobject("Tp1", IOobject::MUST_READ), valid);
    c.checkFieldIOobject(c, Tp1);

    IOField<scalar> Tp2(c.fieldIOobject("Tp2", IOobject::MUST_READ), valid);
    c.checkFieldIOobject(c, Tp2);

    IOField<scalar> Tp3(c.fieldIOobject("Tp3", IOobject::MUST_READ), valid);
    c.checkFieldIOobject(c, Tp3);

    IOField<scalar> mp0(c.fieldIOobject("mp0", IOobject::MUST_READ), valid);
    c.checkFieldIOobject(c, mp0);

    IOField<scalar> mp1(c.fieldIOobject("mp1", IOobject::MUST_READ), valid);
    c.checkFieldIOobject(c, mp1);

    IOField<scalar> mp2(c.fieldIOobject("mp2", IOobject::MUST_READ), valid);
    c.checkFieldIOobject(c, mp2);

    IOField<scalar> mp3(c.fieldIOobject("mp3", IOobject::MUST_READ), valid);
    c.checkFieldIOobject(c, mp3);
    
    IOField<scalar> rDry(c.fieldIOobject("rDry", IOobject::MUST_READ), valid);
    c.checkFieldIOobject(c, rDry);

    IOField<scalar> rDevo(c.fieldIOobject("rDevo", IOobject::MUST_READ), valid);
    c.checkFieldIOobject(c, rDevo); 
    
    IOField<scalar> rChar(c.fieldIOobject("rChar", IOobject::MUST_READ), valid);
    c.checkFieldIOobject(c, rChar);  

    IOField<scalar> rComb(c.fieldIOobject("rComb", IOobject::MUST_READ), valid);
    c.checkFieldIOobject(c, rChar);  

    IOField<scalar> ash_inchar_t(c.fieldIOobject("ash_inchar_t", IOobject::MUST_READ), valid);
    c.checkFieldIOobject(c, ash_inchar_t);  
    
    IOField<scalar> QDry(c.fieldIOobject("QDry", IOobject::MUST_READ), valid);
    c.checkFieldIOobject(c, QDry);  

    IOField<scalar> QComb(c.fieldIOobject("QComb", IOobject::MUST_READ), valid);
    c.checkFieldIOobject(c, QComb);  
    
    IOField<label> flagBoiling(c.fieldIOobject("flagBoiling", IOobject::MUST_READ), valid);
    c.checkFieldIOobject(c, flagBoiling);      
    
    IOField<label> flagDevo(c.fieldIOobject("flagDevo", IOobject::MUST_READ), valid);
    c.checkFieldIOobject(c, flagDevo);      
    
    IOField<scalar> cumTime(c.fieldIOobject("cumTime", IOobject::MUST_READ), valid);
    c.checkFieldIOobject(c, cumTime); 
    
    label i = 0;
    forAllIter(typename Cloud<ReactingMultiphaseIBMParcel<ParcelType> >, c, iter)
    {
        ReactingMultiphaseIBMParcel<ParcelType>& p = iter();

        p.Tb0_ = Tb0[i];
        p.Tb1_ = Tb1[i];
        p.Tb2_ = Tb2[i];
        p.Tb3_ = Tb3[i];
        p.rb0_ = rb0[i];
        p.rb1_ = rb1[i];
        p.rb2_ = rb2[i];
        p.rb3_ = rb3[i];
        p.Tp0_ = Tp0[i];
        p.Tp1_ = Tp1[i];
        p.Tp2_ = Tp2[i];
        p.Tp3_ = Tp3[i];
        p.mp0_ = mp0[i];
        p.mp1_ = mp1[i];
        p.mp2_ = mp2[i];
        p.mp3_ = mp3[i];
        p.rDry_ = rDry[i];
        p.rDevo_ = rDevo[i];   
        p.rChar_ = rChar[i]; 
        p.rComb_ = rComb[i];
        p.ash_inchar_t_ = ash_inchar_t[i]; 
        p.QDry_ = QDry[i];
        p.QComb_ = QComb[i]; 
        p.flagBoiling_ = flagBoiling[i];
        p.flagDevo_ = flagDevo[i];
        p.cumTime_ = cumTime[i];
        i++;
    } 
}


template<class ParcelType>
template<class CloudType>
void Foam::ReactingMultiphaseIBMParcel<ParcelType>::writeFields(const CloudType& c)
{
    ParcelType::writeFields(c);
}


template<class ParcelType>
template<class CloudType, class CompositionType>
void Foam::ReactingMultiphaseIBMParcel<ParcelType>::writeFields
(
    const CloudType& c,
    const CompositionType& compModel
)
{
    ParcelType::writeFields(c, compModel);

    label np = c.size();

    // Write the composition fractions
    {
        IOField<scalar> Tb0(c.fieldIOobject("Tb0", IOobject::NO_READ), np);
        IOField<scalar> Tb1(c.fieldIOobject("Tb1", IOobject::NO_READ), np);
        IOField<scalar> Tb2(c.fieldIOobject("Tb2", IOobject::NO_READ), np);
        IOField<scalar> Tb3(c.fieldIOobject("Tb3", IOobject::NO_READ), np);
        IOField<scalar> rb0(c.fieldIOobject("rb0", IOobject::NO_READ), np);
        IOField<scalar> rb1(c.fieldIOobject("rb1", IOobject::NO_READ), np);
        IOField<scalar> rb2(c.fieldIOobject("rb2", IOobject::NO_READ), np);
        IOField<scalar> rb3(c.fieldIOobject("rb3", IOobject::NO_READ), np);
        IOField<scalar> Tp0(c.fieldIOobject("Tp0", IOobject::NO_READ), np);
        IOField<scalar> Tp1(c.fieldIOobject("Tp1", IOobject::NO_READ), np);
        IOField<scalar> Tp2(c.fieldIOobject("Tp2", IOobject::NO_READ), np);
        IOField<scalar> Tp3(c.fieldIOobject("Tp3", IOobject::NO_READ), np);
        IOField<scalar> mp0(c.fieldIOobject("mp0", IOobject::NO_READ), np);
        IOField<scalar> mp1(c.fieldIOobject("mp1", IOobject::NO_READ), np);
        IOField<scalar> mp2(c.fieldIOobject("mp2", IOobject::NO_READ), np);
        IOField<scalar> mp3(c.fieldIOobject("mp3", IOobject::NO_READ), np);
        IOField<scalar> rDry(c.fieldIOobject("rDry", IOobject::NO_READ), np);
        IOField<scalar> rDevo(c.fieldIOobject("rDevo", IOobject::NO_READ), np);
        IOField<scalar> rChar(c.fieldIOobject("rChar", IOobject::NO_READ), np);
        IOField<scalar> rComb(c.fieldIOobject("rComb", IOobject::NO_READ), np);
        IOField<scalar> ash_inchar_t(c.fieldIOobject("ash_inchar_t", IOobject::NO_READ), np);
        IOField<scalar> QDry(c.fieldIOobject("QDry", IOobject::NO_READ), np);
        IOField<scalar> QComb(c.fieldIOobject("QComb", IOobject::NO_READ), np);
        IOField<label> flagBoiling(c.fieldIOobject("flagBoiling", IOobject::NO_READ), np);
        IOField<label> flagDevo(c.fieldIOobject("flagDevo", IOobject::NO_READ), np);
        IOField<scalar> cumTime(c.fieldIOobject("cumTime", IOobject::NO_READ), np);
        
        label i = 0;
        
        forAllConstIter(typename Cloud<ReactingMultiphaseIBMParcel<ParcelType> >, c, iter)
        {
            const ReactingMultiphaseIBMParcel<ParcelType>& p = iter();

            Tb0[i] = p.Tb0_;
            Tb1[i] = p.Tb1_;
            Tb2[i] = p.Tb2_;
            Tb3[i] = p.Tb3_;
            rb0[i] = p.rb0_;
            rb1[i] = p.rb1_;
            rb2[i] = p.rb2_;
            rb3[i] = p.rb3_;
            Tp0[i] = p.Tp0_;
            Tp1[i] = p.Tp1_;
            Tp2[i] = p.Tp2_;
            Tp3[i] = p.Tp3_;
            mp0[i] = p.mp0_;
            mp1[i] = p.mp1_;
            mp2[i] = p.mp2_;
            mp3[i] = p.mp3_;
            rDry[i] = p.rDry_;
            rDevo[i] = p.rDevo_; 
            rChar[i] = p.rChar_;
            rComb[i] = p.rComb_;
            ash_inchar_t[i] = p.ash_inchar_t_;
            QDry[i] = p.QDry_;
            QComb[i] = p.QComb_;
            flagBoiling[i] = p.flagBoiling_;
            flagDevo[i] = p.flagDevo_;
            cumTime[i] = p.cumTime_; 
            i++;
        }

        Tb0.write(np > 0);
        Tb1.write(np > 0);
        Tb2.write(np > 0);
        Tb3.write(np > 0);
        rb0.write(np > 0);
        rb1.write(np > 0);
        rb2.write(np > 0);
        rb3.write(np > 0);
        Tp0.write(np > 0);
        Tp1.write(np > 0);
        Tp2.write(np > 0);
        Tp3.write(np > 0);
        mp0.write(np > 0);
        mp1.write(np > 0);
        mp2.write(np > 0);
        mp3.write(np > 0);
        rDry.write(np > 0);
        rDevo.write(np > 0);
        rChar.write(np > 0);
        rComb.write(np > 0);
        ash_inchar_t.write(np > 0);
        QDry.write(np > 0);
        QComb.write(np > 0);
        flagBoiling.write(np > 0);
        flagDevo.write(np > 0);
        cumTime.write(np > 0);     
    }
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class ParcelType>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const ReactingMultiphaseIBMParcel<ParcelType>& p
)
{
    if (os.format() == IOstream::ASCII)
    {
        os  << static_cast<const ParcelType&>(p)
            << token::SPACE << p.Tb0()
            << token::SPACE << p.Tb1()
            << token::SPACE << p.Tb2()
            << token::SPACE << p.Tb3()
            << token::SPACE << p.rb0()
            << token::SPACE << p.rb1()
            << token::SPACE << p.rb2()
            << token::SPACE << p.rb3()
            << token::SPACE << p.Tp0()
            << token::SPACE << p.Tp1()
            << token::SPACE << p.Tp2()
            << token::SPACE << p.Tp3()
            << token::SPACE << p.mp0()
            << token::SPACE << p.mp1()
            << token::SPACE << p.mp2()
            << token::SPACE << p.mp3()
            << token::SPACE << p.rDry()
            << token::SPACE << p.rDevo()
            << token::SPACE << p.rChar()
            << token::SPACE << p.rComb()
            << token::SPACE << p.ash_inchar_t()
            << token::SPACE << p.QDry()
            << token::SPACE << p.QComb()
            << token::SPACE << p.flagBoiling()
            << token::SPACE << p.flagDevo()
            << token::SPACE << p.cumTime();
    }
    else
    {
        os  << static_cast<const ParcelType&>(p);
        os.write
        (
            reinterpret_cast<const char*>(&p.Tb0_),
            ReactingMultiphaseIBMParcel<ParcelType>::sizeofFields_
        );  
    }

    // Check state of Ostream
    os.check
    (
        "Ostream& operator<<"
        "("
            "Ostream&, "
            "const ReactingMultiphaseIBMParcel<ParcelType>&"
        ")"
    );

    return os;
}


// ************************************************************************* //
