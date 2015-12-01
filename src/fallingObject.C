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


Foam::wordList createFileNames()
{
    DynamicList<word> names;
    names.append("velocity");
    names.append("acceleration");
    return names;
}
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //
//

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fallingObject::fallingObject(const fvMesh &mesh)
    : functionObjectFile(mesh, "frame", createFileNames()),
      dict_(IOobject("STFProperties",
                     mesh.time().constant(),
                     mesh,
                     IOobject::MUST_READ_IF_MODIFIED,
                     IOobject::NO_WRITE)),
      mesh_(mesh), VF_(vector::zero), aF_(vector::zero),
      log_(dict_.lookupOrDefault<Switch>("log", false)),
      UName_(dict_.lookupOrDefault<word>("UName", "U")),
      sphereI_(mesh_.boundaryMesh().findPatchID("sphere")),
      mass_(readScalar(dict_.lookup("mass"))),
      apparentMass_(readScalar(dict_.lookup("apparentMass"))) {

    createFiles();
}

void Foam::fallingObject::registerVelocity(volVectorField &U) {
    volVectorField::GeometricBoundaryField &patches = U.boundaryField();

    forAll(patches, patchi) {
        fvPatchVectorField &currPatch = patches[patchi];
        if (isA<uniformFixedVelocityFvPatchField>(currPatch)) {
            Info << "Registering: " << currPatch.patch().name() << " patch of "
                 << UName_ << " field. " << endl;
            uniformFixedVelocityFvPatchField &p =
                refCast<uniformFixedVelocityFvPatchField>(currPatch);
            registeredPatches.push_back(&p);
        }
    }
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fallingObject::correctBoundaryVelocity(
    volVectorField &U) const {

    for (auto patch : registeredPatches) {
        patch->correct(VF_);
    }
}

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

void Foam::fallingObject::update(const volScalarField &p,
                                      const volSymmTensorField &R) {

    aF_ = calculate_acceleration(p, R);
    VF_ += (mesh_.time().deltaTValue() * aF_);

    if (log_ && Pstream::master()) {
        file(0) << mesh_.time().timeOutputValue() << setw(1) << " " << VF_.x()
                << setw(1) << " " << VF_.y() << setw(1) << " " << VF_.z()
                << setw(1) << endl;
        ;
        file(1) << mesh_.time().timeOutputValue() << setw(1) << " " << aF_.x()
                << setw(1) << " " << aF_.y() << setw(1) << " " << aF_.z()
                << setw(1) << endl;
    }
}

// ************************************************************************* //