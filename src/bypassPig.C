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

#include "bypassPig.H"
#include "fvCFD.H"

#include "addToRunTimeSelectionTable.H"

namespace Foam {
defineTypeNameAndDebug(bypassPig, 0);

addToRunTimeSelectionTable(translationalFrame, bypassPig, dictionary);
}

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //
//

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::bypassPig::bypassPig(const fvMesh &mesh, const IOdictionary &dict)
    : translationalFrame(mesh, dict),
      pigPatchI_(mesh_.boundaryMesh().findPatchID(
          dict_.subDict(typeName + "Coeffs").lookup("pigPatch"))),
      pigDirection_(dict_.subDict(typeName + "Coeffs").lookup("pigDirection")),
      mass_(
          readScalar(dict_.subDict(typeName + "Coeffs").lookup("mass"))),
      friction_(dict_.subDict(typeName + "Coeffs").lookup("frictionForce")),
      lambda0_(
          readScalar(dict_.subDict(typeName + "Coeffs").lookup("lambda"))) {}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * *
// * //


vector bypassPig::applyFriction(vector F) const
{
    //case 1: |v| != 0 (pig moving)
    if (mag(VF_) != 0) {
        F -= sign(pigDirection_ & VF_) * friction_;
    }
    //case 2: v == 0 |Fsum| > |F_static| (beginning to move)
    else if (mag(pigDirection_ & F) > mag(friction_))
    {
        F -= sign(pigDirection_ & F) * friction_;
    }
    //case 3: v ==0 |Fsum| < |F_static|
    else
    {
        F = vector(0, 0, 0);
    }
    return F;
}
vector Foam::bypassPig::calculate_acceleration(const volScalarField &p,
                                              const volSymmTensorField &R) {
    
    vector pressure = gSum(mesh_.Sf().boundaryField()[pigPatchI_] *
                           p.boundaryField()[pigPatchI_]);

    vector viscous = gSum(mesh_.Sf().boundaryField()[pigPatchI_] &
                          R.boundaryField()[pigPatchI_]);

    vector F = applyFriction(pressure + viscous);

    return ((F / mass_) & pigDirection_) * pigDirection_;
}

// ************************************************************************* //
