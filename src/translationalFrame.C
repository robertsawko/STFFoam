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

#include "fvCFD.H"
#include "geometricOneField.H"
#include "translationalFrame.H"
#include "wordReList.H"

namespace Foam {
defineTypeNameAndDebug(translationalFrame, 0);
defineRunTimeSelectionTable(translationalFrame, dictionary);
}

Foam::wordList createFileNames() {
    DynamicList<word> names;
    names.append("velocity");
    names.append("acceleration");
    return names;
}
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //
//

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::translationalFrame::translationalFrame(const fvMesh &mesh,
                                             const IOdictionary &dict)
    : functionObjects::writeFiles("frame", mesh.time(), dict, "frame"),
      dict_(dict), mesh_(mesh), VF_(vector::zero), aF_(vector::zero),
      log_(dict_.lookupOrDefault<Switch>("log", false)),
      UName_(dict_.lookupOrDefault<word>("UName", "U")) {

    resetNames(createFileNames());
    createFiles();
}

void Foam::translationalFrame::registerVelocity(volVectorField &U) {

    forAll(U.boundaryField(), patchi) {
        fvPatchVectorField &currentPatch = U.boundaryFieldRef()[patchi];
        if (isA<frameAwareBoundary>(currentPatch)) {
            Info << "Registering: " << currentPatch.patch().name()
                 << " patch of " << UName_ << " field. " << endl;
            frameAwareBoundary &p =
                dynamic_cast<frameAwareBoundary &>(currentPatch);
            registeredPatches.push_back(&p);
        }
    }
}

autoPtr<translationalFrame> translationalFrame::New(const fvMesh &mesh) {
    IOdictionary dict = IOdictionary(IOobject("STFProperties",
                                              mesh.time().constant(),
                                              mesh,
                                              IOobject::MUST_READ_IF_MODIFIED,
                                              IOobject::NO_WRITE));

    const word modelType(dict.lookup("translationModel"));

    Info << "Selecting translational frame: " << modelType << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(modelType);

    if (cstrIter == dictionaryConstructorTablePtr_->end()) {
        FatalErrorIn("translationalFrame::New"
                     "("
                     "const dictionary&"
                     ")")
            << "Unknown translation type " << modelType << nl << nl
            << "Valid translation types:" << endl
            << dictionaryConstructorTablePtr_->sortedToc() << exit(FatalError);
    }

    return autoPtr<translationalFrame>(cstrIter()(mesh, dict));
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::translationalFrame::correctBoundaryVelocity(
    volVectorField &U) const {

    for (auto patch : registeredPatches) {
        patch->correct(VF_);
    }
}

void Foam::translationalFrame::update(const volScalarField &p,
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


bool Foam::translationalFrame::execute()
{
    return true;
}

// ************************************************************************* //
