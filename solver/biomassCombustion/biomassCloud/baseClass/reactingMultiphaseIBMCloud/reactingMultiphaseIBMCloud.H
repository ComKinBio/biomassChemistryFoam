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
    Foam::reactingMultiphaseIBMCloud

Description
    Virtual abstract base class for templated reactingMultiphaseIBMCloud

SourceFiles
    reactingMultiphaseIBMCloud.C

\*---------------------------------------------------------------------------*/

#ifndef reactingMultiphaseIBMCloud_H
#define reactingMultiphaseIBMCloud_H

#include "typeInfo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class reactingMultiphaseIBMCloud Declaration
\*---------------------------------------------------------------------------*/

class reactingMultiphaseIBMCloud
{
public:

    //- Runtime type information
    TypeName("reactingMultiphaseIBMCloud");

    // Constructors

        //- Null constructor
        reactingMultiphaseIBMCloud();

        //- Disallow default bitwise copy construction
        reactingMultiphaseIBMCloud(const reactingMultiphaseIBMCloud&) = delete;


    //- Destructor
    virtual ~reactingMultiphaseIBMCloud();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const reactingMultiphaseIBMCloud&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
