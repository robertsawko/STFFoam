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

#include "centreOfVolume.H"
#include "fvCFD.H"

#include "addToRunTimeSelectionTable.H"

namespace Foam {
defineTypeNameAndDebug(centreOfVolume, 0);

addToRunTimeSelectionTable(translationalFrame, centreOfVolume, dictionary);
}

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //
//

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::centreOfVolume::centreOfVolume(const fvMesh &mesh,
                                     const IOdictionary &dict)
    : translationalFrame(mesh, dict),
      alpha_(mesh_.lookupObject<volScalarField>("alpha.air")),
      xinit_("x_0",
             fvc::domainIntegrate(alpha_ * mesh_.C()) /
                 fvc::domainIntegrate(alpha_)),
      xold_("x_o", xinit_),
      lambdaF_(
          readScalar(dict_.subDict(typeName + "Coeffs").lookup("lambdaF"))),
      lambdad_(
          readScalar(dict_.subDict(typeName + "Coeffs").lookup("lambdad"))) {}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

vector
Foam::centreOfVolume::calculate_acceleration(const volScalarField &p,
                                             const volSymmTensorField &R) {
    const Time &t = mesh_.time();

    dimensionedVector xcurrent(
        "x_c",
        fvc::domainIntegrate(alpha_ * mesh_.C()) // centre of mass
            /
            fvc::domainIntegrate(alpha_));

    // Velocity difference (Rusche approach)
    dimensionedVector DeltaV("DeltaV",
                             (lambdaF_ * (xinit_ - xcurrent)) / t.deltaT() -
                                 lambdad_ * (xcurrent - xold_) / t.deltaT());

    xold_.value() = xcurrent.value();

    return -(DeltaV / t.deltaT()).value();
}

// ************************************************************************* //
