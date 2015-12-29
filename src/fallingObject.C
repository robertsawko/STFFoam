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

#include "fallingObject.H"
#include "fvCFD.H"
#include "geometricOneField.H"
#include "wordReList.H"

#include "addToRunTimeSelectionTable.H"

namespace Foam {
defineTypeNameAndDebug(fallingObject, 0);

addToRunTimeSelectionTable(translationalFrame, fallingObject, dictionary);
}

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //
//

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fallingObject::fallingObject(const fvMesh &mesh, const IOdictionary &dict)
    : translationalFrame(mesh, dict),
      sphereI_(mesh_.boundaryMesh().findPatchID(
          dict_.subDict(typeName + "Coeffs").lookup("objectName"))),
      mass_(readScalar(dict_.subDict(typeName + "Coeffs").lookup("mass"))),
      apparentMass_(readScalar(
          dict_.subDict(typeName + "Coeffs").lookup("apparentMass"))) {}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

vector
Foam::fallingObject::calculate_acceleration(const volScalarField &p,
                                            const volSymmTensorField &R) {

    vector pressure = gSum(mesh_.Sf().boundaryField()[sphereI_] *
                           p.boundaryField()[sphereI_]);

    vector viscous = gSum(mesh_.Sf().boundaryField()[sphereI_] &
                          R.boundaryField()[sphereI_]);

    vector gravity(0, 0, -9.81);

    return (pressure + viscous + apparentMass_ * gravity) / mass_;
}

// ************************************************************************* //
