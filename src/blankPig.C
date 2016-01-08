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

\*---------------------------------------------------------------------------*/

#include "blankPig.H"
#include "fvCFD.H"

#include "addToRunTimeSelectionTable.H"

namespace Foam {
defineTypeNameAndDebug(blankPig, 0);

addToRunTimeSelectionTable(translationalFrame, blankPig, dictionary);
}

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //
//

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::blankPig::blankPig(const fvMesh &mesh, const IOdictionary &dict)
    : translationalFrame(mesh, dict),
      pigPatchI_(mesh_.boundaryMesh().findPatchID(
          dict_.subDict(typeName + "Coeffs").lookup("pigPatch"))),
      pigDirection_(dict_.subDict(typeName + "Coeffs").lookup("pigDirection")),
      pigMass_(dict_.subDict(typeName + "Coeffs").lookup("pigMass")),
      lambda0_(
          readScalar(dict_.subDict(typeName + "Coeffs").lookup("lambda"))) {}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * *
// * //

vector Foam::blankPig::calculate_acceleration(const volScalarField &p,
                                              const volSymmTensorField &R) {
    
    vector pressure = gSum(mesh_.Sf().boundaryField()[pigPatchI_] *
                           p.boundaryField()[sphereI_]);

    vector viscous = gSum(mesh_.Sf().boundaryField()[pigPatchI_] &
                          R.boundaryField()[sphereI_]);

    return (pressure + viscous) / pigMass_.value();
}

// ************************************************************************* //
