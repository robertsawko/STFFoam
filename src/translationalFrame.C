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
#include "addToRunTimeSelectionTable.H"
#include "uniformFixedVelocityFvPatchField.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(translationalFrame, 0);
    addToRunTimeSelectionTable
    (
        option,
        translationalFrame,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::translationalFrame::translationalFrame(const word &name,
                                             const word &modelType,
                                             const dictionary &dict,
                                             const fvMesh &mesh)
    : option(name, modelType, dict, mesh),
      VF_(coeffs_.lookupOrDefault<vector>("velocity", vector::zero)),
      // VF_("V", dimVelocity, vector(0, 0, 0)),
      UName_(coeffs_.lookupOrDefault<word>("UName", "U")),
      g0_("g0", dimAcceleration, vector::zero) {
    fieldNames_.setSize(1, UName_);
    applied_.setSize(1, false);

    if (mesh.foundObject<uniformDimensionedVectorField>("g")) {
        g0_ = mesh.lookupObject<uniformDimensionedVectorField>("g");
    }

    // Registering BCs
    const volVectorField& U = mesh.lookupObject<volVectorField>(UName_);
    const volVectorField::GeometricBoundaryField& patches = U.boundaryField();

    forAll (patches, patchi)
    {
        const fvPatchVectorField& currPatch = patches[patchi];
        if (isA<uniformFixedVelocityFvPatchField>(currPatch)) {
            Info<< "Registering: " << currPatch.patch().name() << " patch of "
                << UName_ << " field. " << endl;
            patchIDs.push_back(patchi);
        }
    }
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
/*
void Foam::fv::translationalFrame::addSup
(
    fvMatrix<vector>& eqn,
    const label fieldi
)
{
    addSup<geometricOneField>(geometricOneField(), eqn, fieldi);
}


void Foam::fv::translationalFrame::addSup
(
    const volScalarField& rho,
    fvMatrix<vector>& eqn,
    const label fieldi
)
{
    addSup<volScalarField>(rho, eqn, fieldi);
}
*/


bool Foam::fv::translationalFrame::read(const dictionary& dict)
{
    /*
    if (option::read(dict))
    {
        return motion_.read(coeffs_);
    }
    else
    {
        return false;
    }
    */
    return true;
}

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
