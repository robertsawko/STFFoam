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

Class
    Foam::translatingFixedVelocityFvPatchField

Group
    grpGenericBoundaryConditions

Description
    This boundary condition provides a uniform fixed value condition.

    \heading Patch usage

    \table
        Property     | Description             | Required    | Default value
        fixedFrameValue | uniform value           |         yes |
    \endtable

    Example of the boundary condition specification:
    \verbatim
    myPatch
    {
        type            uniformFixedValue;
        fixedFrameValue    constant 0.2;
    }
    \endverbatim

Note
    The fixedFrameValue entry is a DataEntry type, able to describe time
    varying functions.  The example above gives the usage for supplying a
    constant value.

SeeAlso
    Foam::DataEntry
    Foam::fixedValueFvPatchField

SourceFiles
    translatingFixedVelocityFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef translatingFixedVelocityFvPatchField_H
#define translatingFixedVelocityFvPatchField_H

#include "fixedValueFvPatchFields.H"
#include "frameAwareBoundary.H"
#include "DataEntry.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam {

/*---------------------------------------------------------------------------*\
                Class translatingFixedVelocityFvPatchField Declaration
\*---------------------------------------------------------------------------*/

class translatingFixedVelocityFvPatchField : public fixedValueFvPatchVectorField,
                                         public frameAwareBoundary {
    // Private data

    autoPtr<DataEntry<vector>> fixedFrameValue_;

  public:
    //- Runtime type information
    TypeName("translatingFixedVelocity");

    // Constructors

    //- Construct from patch and internal field
    translatingFixedVelocityFvPatchField(const fvPatch &,
                                     const DimensionedField<vector, volMesh> &);

    //- Construct from patch, internal field and dictionary
    translatingFixedVelocityFvPatchField(const fvPatch &,
                                     const DimensionedField<vector, volMesh> &,
                                     const dictionary &);

    //- Construct by mapping given translatingFixedVelocityFvPatchField
    //  onto a new patch
    translatingFixedVelocityFvPatchField(const translatingFixedVelocityFvPatchField &,
                                     const fvPatch &,
                                     const DimensionedField<vector, volMesh> &,
                                     const fvPatchFieldMapper &);

    //- Construct as copy
    translatingFixedVelocityFvPatchField(const translatingFixedVelocityFvPatchField &);
    virtual ~translatingFixedVelocityFvPatchField() {}

    //- Construct and return a clone
    virtual tmp<fvPatchVectorField> clone() const {
        return tmp<fvPatchVectorField>(
            new translatingFixedVelocityFvPatchField(*this));
    }

    //- Construct as copy setting internal field reference
    translatingFixedVelocityFvPatchField(const translatingFixedVelocityFvPatchField &,
                                     const DimensionedField<vector, volMesh> &);

    //- Construct and return a clone setting internal field reference
    virtual tmp<fvPatchVectorField>
    clone(const DimensionedField<vector, volMesh> &iF) const {
        return tmp<fvPatchVectorField>(
            new translatingFixedVelocityFvPatchField(*this, iF));
    }

    // Member functions

    //- Update the coefficients associated with the patch field
    virtual void updateCoeffs();

    //- Write
    virtual void write(Ostream &) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
