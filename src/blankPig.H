/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015 OpenFOAM Foundation
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
    Foam::fv::blankPig

Description
    Object moving under the action of gravity and surrounding fluid.

    \heading Source usage

    Example usage:
    \verbatim
    modelType ; 

    translationModel blankPig;

    blankPigCoeffs 
    {
        mass 0.62;
        apparentMass 0.0961;
        objectName sphere;
    }
    \endverbatim


SourceFiles
    blankPig.C

\*---------------------------------------------------------------------------*/

#ifndef blankPig_H
#define blankPig_H

#include "translationalFrame.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam {
/*---------------------------------------------------------------------------*\
                          Class blankPig Declaration
\*---------------------------------------------------------------------------*/

class blankPig : public translationalFrame {
  private:
    // Private Member Functions

    //- Disallow default bitwise copy construct
    blankPig(const centreOfVolume &);

    //- Disallow default bitwise assignment
    void operator=(const blankPig &);

    //- The ID of the patch identifying the front
    const volScalarField &alpha_;

    const dimensionedVector xinit_;
    dimensionedVector xold_;

    const scalar lambdaF_;
    const scalar lambdad_;

    virtual vector calculate_acceleration(const volScalarField &,
                                  const volSymmTensorField &);

  public:
    TypeName("blankPig");

    blankPig(const fvMesh &, const IOdictionary &);

    virtual ~blankPig() {}
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
