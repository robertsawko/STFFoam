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
    Object is propelled through hydraulic action and experiences wall
interaction.

    \heading Source usage

    Example usage:
    \verbatim
    modelType ;

    translationModel blankPig;

    blankPigCoeffs
    {
        mass 1;
        pigPatch pig;
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
    blankPig(const blankPig &);

    //- Disallow default bitwise assignment
    void operator=(const blankPig &);

    //- The ID of the patch identifying the front
    const label pigPatchI_;
    const vector pigDirection_;
    const dimensionedScalar pigMass_;

    //- Under relaxation factor for velocity
    const scalar lambda0_;

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
