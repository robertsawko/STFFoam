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

#include "uniformFixedVelocityFvPatchField.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::uniformFixedVelocityFvPatchField::uniformFixedVelocityFvPatchField(
    const fvPatch &p, const DimensionedField<vector, volMesh> &iF)
    : fixedValueFvPatchField<vector>(p, iF), uniformValue_(),
      translationVelocity_(vector::zero) {}

Foam::uniformFixedVelocityFvPatchField::uniformFixedVelocityFvPatchField(
    const uniformFixedVelocityFvPatchField &ptf,
    const fvPatch &p,
    const DimensionedField<vector, volMesh> &iF,
    const fvPatchFieldMapper &mapper)
    : fixedValueFvPatchField<vector>(p, iF), // bypass mapper
      uniformValue_(ptf.uniformValue_().clone().ptr()),
      translationVelocity_(ptf.translationVelocity_) {
    // Evaluate since value not mapped
    const scalar t = this->db().time().timeOutputValue();
    fvPatchField<vector>::operator==(uniformValue_->value(t) -
                                     translationVelocity_);
}

Foam::uniformFixedVelocityFvPatchField::uniformFixedVelocityFvPatchField(
    const fvPatch &p,
    const DimensionedField<vector, volMesh> &iF,
    const dictionary &dict)
    : fixedValueFvPatchField<vector>(p, iF),
      uniformValue_(DataEntry<vector>::New("uniformValue", dict)),
      translationVelocity_(vector::zero) {
    const scalar t = this->db().time().timeOutputValue();
    fvPatchField<vector>::operator=(vectorField("value", dict, p.size()));
}

Foam::uniformFixedVelocityFvPatchField::uniformFixedVelocityFvPatchField(
    const uniformFixedVelocityFvPatchField &ptf)
    : fixedValueFvPatchField<vector>(ptf),
      uniformValue_(ptf.uniformValue_.valid()
                        ? ptf.uniformValue_().clone().ptr()
                        : nullptr),
      translationVelocity_(ptf.translationVelocity_) {}

Foam::uniformFixedVelocityFvPatchField::uniformFixedVelocityFvPatchField(
    const uniformFixedVelocityFvPatchField &ptf,
    const DimensionedField<vector, volMesh> &iF)
    : fixedValueFvPatchField<vector>(ptf, iF),
      uniformValue_(ptf.uniformValue_.valid()
                        ? ptf.uniformValue_().clone().ptr()
                        : nullptr),
      translationVelocity_(ptf.translationVelocity_) {
    /*
    // For safety re-evaluate
    const scalar t = this->db().time().timeOutputValue();

    if (ptf.uniformValue_.valid())
    {
        fvPatchField<vector>::operator==(uniformValue_->value(t));
    }
    */
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::uniformFixedVelocityFvPatchField::updateCoeffs() {
    if (this->updated()) {
        return;
    }

    const scalar t = this->db().time().timeOutputValue();
    operator==(uniformValue_->value(t) - translationVelocity_);

    fixedValueFvPatchVectorField::updateCoeffs();
}

void Foam::uniformFixedVelocityFvPatchField::write(Ostream &os) const {
    fvPatchField<vector>::write(os);
    if(uniformValue_.valid())
        uniformValue_->writeData(os);
    os.writeKeyword("translationVelocity") << translationVelocity_
                                           << token::END_STATEMENT << nl;
    writeEntry("value", os);
}

void Foam::uniformFixedVelocityFvPatchField::correct(const vector &VF) {
    translationVelocity_ = VF;
}

namespace Foam {
makePatchTypeField(fvPatchVectorField, uniformFixedVelocityFvPatchField);
}
// ************************************************************************* //
