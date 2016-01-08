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
    Foam::STFScouringInflowVelocityFvPatchVectorField

Group
    grpGenericBoundaryConditions

Description
    This boundary condition provides a uniform fixed value condition.

    \heading Patch usage

    \table
        Property     | Description             | Required    | Default value
        uniformValue | uniform value           |         yes |
    \endtable

    Example of the boundary condition specification:
    \verbatim
    myPatch
    {
        type            uniformFixedValue;
        uniformValue    constant 0.2;
    }
    \endverbatim

Note
    The uniformValue entry is a DataEntry type, able to describe time
    varying functions.  The example above gives the usage for supplying a
    constant value.

SeeAlso
    Foam::DataEntry
    Foam::fixedValueFvPatchField

SourceFiles
    STFScouringInflowVelocityFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef STFScouringInflowVelocityFvPatchVectorField_H
#define STFScouringInflowVelocityFvPatchVectorField_H

#include "fixedValueFvPatchFields.H"
#include "DataEntry.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam {

/*---------------------------------------------------------------------------*\
                Class STFScouringInflowVelocityFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class STFScouringInflowVelocityFvPatchVectorField : public fixedValueFvPatchVectorField {
    // Private data

    autoPtr<DataEntry<scalar>> pipeDiameter_;
    autoPtr<DataEntry<scalar>> depositHeight_;
    vector translationVelocity_;

  public:
    //- Runtime type information
    TypeName("STFScouringInflow");

    // Constructors

    //- Construct from patch and internal field
    STFScouringInflowVelocityFvPatchVectorField(const fvPatch &,
                                     const DimensionedField<vector, volMesh> &);

    //- Construct from patch, internal field and dictionary
    STFScouringInflowVelocityFvPatchVectorField(const fvPatch &,
                                     const DimensionedField<vector, volMesh> &,
                                     const dictionary &);

    //- Construct by mapping given STFScouringInflowVelocityFvPatchVectorField
    //  onto a new patch
    STFScouringInflowVelocityFvPatchVectorField(const STFScouringInflowVelocityFvPatchVectorField &,
                                     const fvPatch &,
                                     const DimensionedField<vector, volMesh> &,
                                     const fvPatchFieldMapper &);

    //- Construct as copy
    STFScouringInflowVelocityFvPatchVectorField(const STFScouringInflowVelocityFvPatchVectorField &);
    virtual ~STFScouringInflowVelocityFvPatchVectorField() {}

    //- Construct and return a clone
    virtual tmp<fvPatchVectorField> clone() const {
        return tmp<fvPatchVectorField>(
            new STFScouringInflowVelocityFvPatchVectorField(*this));
    }

    //- Construct as copy setting internal field reference
    STFScouringInflowVelocityFvPatchVectorField(const STFScouringInflowVelocityFvPatchVectorField &,
                                     const DimensionedField<vector, volMesh> &);

    //- Construct and return a clone setting internal field reference
    virtual tmp<fvPatchVectorField>
    clone(const DimensionedField<vector, volMesh> &iF) const {
        return tmp<fvPatchVectorField>(
            new STFScouringInflowVelocityFvPatchVectorField(*this, iF));
    }

    // Member functions

    //- Update the coefficients associated with the patch field
    virtual void updateCoeffs();

    virtual void correct(const vector &VF);

    //- Write
    virtual void write(Ostream &) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
