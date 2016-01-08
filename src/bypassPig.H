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
    Foam::fv::bypassPig

Description
    Object is propelled through hydraulic action and experiences wall
interaction.

    \heading Source usage

    Example usage:
    \verbatim
    modelType ;

    translationModel bypassPig;

    bypassPigCoeffs
    {
        mass 1;
        pigPatch pig;
    }
    \endverbatim


SourceFiles
    bypassPig.C

\*---------------------------------------------------------------------------*/

#ifndef bypassPig_H
#define bypassPig_H

#include "translationalFrame.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam {
/*---------------------------------------------------------------------------*\
                          Class bypassPig Declaration
\*---------------------------------------------------------------------------*/

class bypassPig : public translationalFrame {
  private:
    // Private Member Functions

    //- Disallow default bitwise copy construct
    bypassPig(const bypassPig &);

    //- Disallow default bitwise assignment
    void operator=(const bypassPig &);

    //- The ID of the patch identifying the front
    const label pigPatchI_;
    const vector pigDirection_;
    const scalar mass_;
    const vector friction_;

    //- Under relaxation factor for velocity
    const scalar lambda0_;

    vector applyFriction(vector F) const;
    virtual vector calculate_acceleration(const volScalarField &,
                                          const volSymmTensorField &);

  public:
    TypeName("bypassPig");

    bypassPig(const fvMesh &, const IOdictionary &);

    virtual ~bypassPig() {}
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
