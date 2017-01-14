/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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
#include "mathematicalConstants.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::STFScouringInflowVelocityFvPatchVectorField::
    STFScouringInflowVelocityFvPatchVectorField(
        const fvPatch &p, const DimensionedField<vector, volMesh> &iF)
    : fixedValueFvPatchField<vector>(p, iF), pipeDiameter_(), depositHeight_() {
}

Foam::STFScouringInflowVelocityFvPatchVectorField::
    STFScouringInflowVelocityFvPatchVectorField(
        const STFScouringInflowVelocityFvPatchVectorField &ptf,
        const fvPatch &p,
        const DimensionedField<vector, volMesh> &iF,
        const fvPatchFieldMapper &mapper)
    : fixedValueFvPatchField<vector>(p, iF), // bypass mapper
      pipeDiameter_(ptf.pipeDiameter_().clone().ptr()),
      depositHeight_(ptf.depositHeight_().clone().ptr()) {
    // Evaluate since value not mapped
    // const scalar t = this->db().time().timeOutputValue();
    // fvPatchField<vector>::operator==(pipeDiameter_->value(t) -
    // translationVelocity_);
}

Foam::STFScouringInflowVelocityFvPatchVectorField::
    STFScouringInflowVelocityFvPatchVectorField(
        const fvPatch &p,
        const DimensionedField<vector, volMesh> &iF,
        const dictionary &dict)
    : fixedValueFvPatchField<vector>(p, iF),
      pipeDiameter_(Function1<scalar>::New("pipeDiameter", dict)),
      depositHeight_(Function1<scalar>::New("depositHeight", dict)) {
    const scalar t = this->db().time().timeOutputValue();
    fvPatchField<vector>::operator=(vectorField("value", dict, p.size()));
}

Foam::STFScouringInflowVelocityFvPatchVectorField::
    STFScouringInflowVelocityFvPatchVectorField(
        const STFScouringInflowVelocityFvPatchVectorField &ptf)
    : fixedValueFvPatchField<vector>(ptf),
      frameAwareBoundary(ptf),
      pipeDiameter_(ptf.pipeDiameter_.valid()
                        ? ptf.pipeDiameter_().clone().ptr()
                        : nullptr),
      depositHeight_(ptf.depositHeight_.valid()
                         ? ptf.depositHeight_().clone().ptr()
                         : nullptr) {}

Foam::STFScouringInflowVelocityFvPatchVectorField::
    STFScouringInflowVelocityFvPatchVectorField(
        const STFScouringInflowVelocityFvPatchVectorField &ptf,
        const DimensionedField<vector, volMesh> &iF)
    : fixedValueFvPatchField<vector>(ptf, iF),
      pipeDiameter_(ptf.pipeDiameter_.valid()
                        ? ptf.pipeDiameter_().clone().ptr()
                        : nullptr),
      depositHeight_(ptf.depositHeight_.valid()
                         ? ptf.depositHeight_().clone().ptr()
                         : nullptr) {
    /*
    // For safety re-evaluate
    const scalar t = this->db().time().timeOutputValue();

    if (ptf.pipeDiameter_.valid())
    {
        fvPatchField<vector>::operator==(pipeDiameter_->value(t));
    }
    */
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::STFScouringInflowVelocityFvPatchVectorField::updateCoeffs() {
    if (this->updated()) {
        return;
    }

    const scalar t = db().time().timeOutputValue();
    using Foam::constant::mathematical::pi;

    const scalar avgU =
        -(pi * mag(currentVelocity()) * depositHeight_->value(t) *
          (pipeDiameter_->value(t) - depositHeight_->value(t))) /
        gSum(patch().magSf());

    tmp<vectorField> n = patch().nf();

    // volumetric flow-rate or density not given
    operator==(n * avgU);

    fixedValueFvPatchVectorField::updateCoeffs();
}

void Foam::STFScouringInflowVelocityFvPatchVectorField::write(
    Ostream &os) const {
    fvPatchField<vector>::write(os);
    if (pipeDiameter_.valid())
        pipeDiameter_->writeData(os);
    if (depositHeight_.valid())
        depositHeight_->writeData(os);
    os.writeKeyword("translationVelocity") << currentVelocity()
                                           << token::END_STATEMENT << nl;
    writeEntry("value", os);
}

namespace Foam {
makePatchTypeField(fvPatchVectorField,
                   STFScouringInflowVelocityFvPatchVectorField);
}
// ************************************************************************* //
