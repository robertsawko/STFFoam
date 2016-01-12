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

Foam::translatingFixedVelocityFvPatchField::translatingFixedVelocityFvPatchField(
    const fvPatch &p, const DimensionedField<vector, volMesh> &iF)
    : fixedValueFvPatchField<vector>(p, iF), frameAwareBoundary(),
      fixedFrameValue_() {}

Foam::translatingFixedVelocityFvPatchField::translatingFixedVelocityFvPatchField(
    const translatingFixedVelocityFvPatchField &ptf,
    const fvPatch &p,
    const DimensionedField<vector, volMesh> &iF,
    const fvPatchFieldMapper &mapper)
    : fixedValueFvPatchField<vector>(p, iF), // bypass mapper
      fixedFrameValue_(ptf.fixedFrameValue_().clone().ptr()) {
    // Evaluate since value not mapped
    const scalar t = this->db().time().timeOutputValue();
    fvPatchField<vector>::operator==(fixedFrameValue_->value(t) -
                                     currentVelocity());
}

Foam::translatingFixedVelocityFvPatchField::translatingFixedVelocityFvPatchField(
    const fvPatch &p,
    const DimensionedField<vector, volMesh> &iF,
    const dictionary &dict)
    : fixedValueFvPatchField<vector>(p, iF),
      fixedFrameValue_(DataEntry<vector>::New("fixedFrameValue", dict)) {
    const scalar t = this->db().time().timeOutputValue();
    fvPatchField<vector>::operator=(vectorField("value", dict, p.size()));
}

Foam::translatingFixedVelocityFvPatchField::translatingFixedVelocityFvPatchField(
    const translatingFixedVelocityFvPatchField &ptf)
    : fixedValueFvPatchField<vector>(ptf), frameAwareBoundary(ptf),
      fixedFrameValue_(ptf.fixedFrameValue_.valid()
                        ? ptf.fixedFrameValue_().clone().ptr()
                        : nullptr) {}

Foam::translatingFixedVelocityFvPatchField::translatingFixedVelocityFvPatchField(
    const translatingFixedVelocityFvPatchField &ptf,
    const DimensionedField<vector, volMesh> &iF)
    : fixedValueFvPatchField<vector>(ptf, iF), frameAwareBoundary(ptf),
      fixedFrameValue_(ptf.fixedFrameValue_.valid()
                        ? ptf.fixedFrameValue_().clone().ptr()
                        : nullptr) {
    /*
    // For safety re-evaluate
    const scalar t = this->db().time().timeOutputValue();

    if (ptf.fixedFrameValue_.valid())
    {
        fvPatchField<vector>::operator==(fixedFrameValue_->value(t));
    }
    */
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::translatingFixedVelocityFvPatchField::updateCoeffs() {
    if (this->updated()) {
        return;
    }

    const scalar t = this->db().time().timeOutputValue();
    operator==(fixedFrameValue_->value(t) - currentVelocity());

    fixedValueFvPatchVectorField::updateCoeffs();
}

void Foam::translatingFixedVelocityFvPatchField::write(Ostream &os) const {
    fvPatchField<vector>::write(os);
    if (fixedFrameValue_.valid())
        fixedFrameValue_->writeData(os);
    os.writeKeyword("translationVelocity") << currentVelocity()
                                           << token::END_STATEMENT << nl;
    writeEntry("value", os);
}

namespace Foam {
makePatchTypeField(fvPatchVectorField, translatingFixedVelocityFvPatchField);
}
// ************************************************************************* //
