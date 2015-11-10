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

#include "translationalFrame.H"
#include "fvCFD.H"
#include "geometricOneField.H"
#include "uniformFixedVelocityFvPatchField.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::translationalFrame::translationalFrame(const fvMesh &mesh)
    : dict_(IOobject("STFProperties",
                            mesh.time().constant(),
                            mesh,
                            IOobject::MUST_READ_IF_MODIFIED,
                            IOobject::NO_WRITE)),
      mesh_(mesh), VF_(dict_.lookupOrDefault<vector>("velocity", vector::zero)),
      // VF_("V", dimVelocity, vector(0, 0, 0)),
      UName_(dict_.lookupOrDefault<word>("UName", "U"))
// g0_("g0", dimAcceleration, vector::zero)
{

    // if (mesh_.foundObject<uniformDimensionedVectorField>("g")) {
    // g0_ = mesh.lookupObject<uniformDimensionedVectorField>("g");
    //}

    // Registering BCs
    const volVectorField &U = mesh.lookupObject<volVectorField>(UName_);
    const volVectorField::GeometricBoundaryField &patches = U.boundaryField();

    forAll(patches, patchi) {
        const fvPatchVectorField &currPatch = patches[patchi];
        if (isA<uniformFixedVelocityFvPatchField>(currPatch)) {
            Info << "Registering: " << currPatch.patch().name() << " patch of "
                 << UName_ << " field. " << endl;
            patchIDs.push_back(patchi);
        }
    }
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fv::translationalFrame::correct(volVectorField& field){

    for(auto id: patchIDs){
        Info<< id;
    }
/*
    volScalarField::GeometricBoundaryField& patches = field.boundaryField();
    forAll (patches, patchi)
    {
        fvPatchScalarField& currPatch = patches[patchi];
        if (isA<STFInletOutletFvPatchField>(currPatch))
        {
            (refCast<STFInletOutletFvPatchField> (currPatch))
        }
    }
    */
}

// ************************************************************************* //
