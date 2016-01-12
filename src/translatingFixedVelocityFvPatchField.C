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

#include "translatingFixedVelocityFvPatchField.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::uniformFixedVelocityFvPatchField::uniformFixedVelocityFvPatchField(
    const fvPatch &p, const DimensionedField<vector, volMesh> &iF)
    : fixedValueFvPatchField<vector>(p, iF), frameAwareBoundary(),
      uniformValue_() {}

Foam::uniformFixedVelocityFvPatchField::uniformFixedVelocityFvPatchField(
    const uniformFixedVelocityFvPatchField &ptf,
    const fvPatch &p,
    const DimensionedField<vector, volMesh> &iF,
    const fvPatchFieldMapper &mapper)
    : fixedValueFvPatchField<vector>(p, iF), // bypass mapper
      uniformValue_(ptf.uniformValue_().clone().ptr()) {
    // Evaluate since value not mapped
    const scalar t = this->db().time().timeOutputValue();
    fvPatchField<vector>::operator==(uniformValue_->value(t) -
                                     currentVelocity());
}

Foam::uniformFixedVelocityFvPatchField::uniformFixedVelocityFvPatchField(
    const fvPatch &p,
    const DimensionedField<vector, volMesh> &iF,
    const dictionary &dict)
    : fixedValueFvPatchField<vector>(p, iF),
      uniformValue_(DataEntry<vector>::New("uniformValue", dict)) {
    const scalar t = this->db().time().timeOutputValue();
    fvPatchField<vector>::operator=(vectorField("value", dict, p.size()));
}

Foam::uniformFixedVelocityFvPatchField::uniformFixedVelocityFvPatchField(
    const uniformFixedVelocityFvPatchField &ptf)
    : fixedValueFvPatchField<vector>(ptf), frameAwareBoundary(ptf),
      uniformValue_(ptf.uniformValue_.valid()
                        ? ptf.uniformValue_().clone().ptr()
                        : nullptr) {}

Foam::uniformFixedVelocityFvPatchField::uniformFixedVelocityFvPatchField(
    const uniformFixedVelocityFvPatchField &ptf,
    const DimensionedField<vector, volMesh> &iF)
    : fixedValueFvPatchField<vector>(ptf, iF), frameAwareBoundary(ptf),
      uniformValue_(ptf.uniformValue_.valid()
                        ? ptf.uniformValue_().clone().ptr()
                        : nullptr) {
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
    operator==(uniformValue_->value(t) - currentVelocity());

    fixedValueFvPatchVectorField::updateCoeffs();
}

void Foam::uniformFixedVelocityFvPatchField::write(Ostream &os) const {
    fvPatchField<vector>::write(os);
    if (uniformValue_.valid())
        uniformValue_->writeData(os);
    os.writeKeyword("translationVelocity") << currentVelocity()
                                           << token::END_STATEMENT << nl;
    writeEntry("value", os);
}

namespace Foam {
makePatchTypeField(fvPatchVectorField, uniformFixedVelocityFvPatchField);
}
// ************************************************************************* //
