/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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

#include "STFScouringInflowVelocityFvPatchVectorField.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "surfaceFields.H"
#include "scalarIOList.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::STFScouringInflowVelocityFvPatchVectorField::
    STFScouringInflowVelocityFvPatchVectorField(
        const fvPatch &p, const DimensionedField<vector, volMesh> &iF)
    : fixedValueFvPatchField<vector>(p, iF), depositHeight_(), pipeDiameter_(),
      translationVelocity_(vector::zero) {}

Foam::STFScouringInflowVelocityFvPatchVectorField::
    STFScouringInflowVelocityFvPatchVectorField(
        const STFScouringInflowVelocityFvPatchVectorField &ptf,
        const fvPatch &p,
        const DimensionedField<vector, volMesh> &iF,
        const fvPatchFieldMapper &mapper)
    : fixedValueFvPatchField<vector>(ptf, p, iF, mapper),
      depositHeight_(ptf.depositHeight_().clone().ptr()),
      pipeDiameter_(ptf.pipeDiameter_().clone().ptr()),
      translationVelocity_(ptf.translationVelocity_) {}

Foam::STFScouringInflowVelocityFvPatchVectorField::
    STFScouringInflowVelocityFvPatchVectorField(
        const fvPatch &p,
        const DimensionedField<vector, volMesh> &iF,
        const dictionary &dict)
    : fixedValueFvPatchField<vector>(p, iF),
      translationVelocity_(vector::zero) {
    depositHeight_ = DataEntry<scalar>::New("depositHeight", dict);
    pipeDiameter_ = DataEntry<scalar>::New("pipeDiameter", dict);

    // Value field require if mass based
    if (dict.found("value")) {
        fvPatchField<vector>::operator=(vectorField("value", dict, p.size()));
    } else {
        evaluate(Pstream::blocking);
    }
}

Foam::STFScouringInflowVelocityFvPatchVectorField::
    STFScouringInflowVelocityFvPatchVectorField(
        const STFScouringInflowVelocityFvPatchVectorField &ptf)
    : fixedValueFvPatchField<vector>(ptf),
      depositHeight_(ptf.depositHeight_().clone().ptr()),
      pipeDiameter_(ptf.pipeDiameter_().clone().ptr()),
      translationVelocity_(ptf.translationVelocity_) {}

Foam::STFScouringInflowVelocityFvPatchVectorField::
    STFScouringInflowVelocityFvPatchVectorField(
        const STFScouringInflowVelocityFvPatchVectorField &ptf,
        const DimensionedField<vector, volMesh> &iF)
    : fixedValueFvPatchField<vector>(ptf, iF),
      depositHeight_(ptf.depositHeight_().clone().ptr()),
      pipeDiameter_(ptf.pipeDiameter_().clone().ptr()),
      translationVelocity_(ptf.translationVelocity_) {}

void Foam::STFScouringInflowVelocityFvPatchVectorField::updateCoeffs() {
    if (updated()) {
        return;
    }

    const scalar t = db().time().timeOutputValue();
    const scalar pi = Foam::constant::mathematical::pi;

    const scalar avgU = -(pi * mag(translationVelocity_) * depositHeight_->value(t) *
                          (pipeDiameter_->value(t) - depositHeight_->value(t))) /
                        gSum(patch().magSf());

    tmp<vectorField> n = patch().nf();

    // volumetric flow-rate or density not given
    operator==(n * avgU);

    fixedValueFvPatchField<vector>::updateCoeffs();
}

void Foam::STFScouringInflowVelocityFvPatchVectorField::write(
    Ostream &os) const {
    fvPatchField<vector>::write(os);
    if (depositHeight_.valid())
        depositHeight_->writeData(os);
    if (pipeDiameter_.valid())
        pipeDiameter_->writeData(os);
    writeEntry("value", os);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam {
makePatchTypeField(fvPatchVectorField,
                   STFScouringInflowVelocityFvPatchVectorField);
}

// ************************************************************************* //
