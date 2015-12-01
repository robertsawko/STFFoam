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

Class
    Foam::fv::fallingObject

Description
    Solid-body 6-DoF acceleration source

    \heading Source usage

    Example usage:
    \verbatim
    SBM
    {
        type            fallingObject;
        active          yes;

        fallingObjectCoeffs
        {
            timeDataFileName "constant/acceleration-terms.dat";
        }
    }
    \endverbatim


SourceFiles
    fallingObject.C

\*---------------------------------------------------------------------------*/

#ifndef fallingObject_H
#define fallingObject_H

#include "functionObjectFile.H"
#include "fvMesh.H"
#include "dimensionedTypes.H"
#include "Switch.H"
#include <vector>
#include "uniformFixedVelocityFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam {
/*---------------------------------------------------------------------------*\
                          Class fallingObject Declaration
\*---------------------------------------------------------------------------*/

class fallingObject : public functionObjectFile {
  private:
    // Private Member Functions

    //- Disallow default bitwise copy construct
    fallingObject(const fallingObject &);

    //- Disallow default bitwise assignment
    void operator=(const fallingObject &);
    IOdictionary dict_;
    const fvMesh &mesh_;
    vector VF_;
    vector aF_;
    Switch log_;
    word UName_;

    //- The ID of the patch identifying the front
    label sphereI_;
    scalar mass_;
    scalar apparentMass_;

    std::vector<uniformFixedVelocityFvPatchField *> registeredPatches;

    vector calculate_acceleration(const volScalarField &,
                                  const volSymmTensorField &);

  public:
    fallingObject(const fvMesh &mesh);

    virtual ~fallingObject() {}

    void registerVelocity(volVectorField &U);
    void correctBoundaryVelocity(volVectorField &U) const;
    void update(const volScalarField &p, const volSymmTensorField &R);
    dimensionedVector aF() const {
        return dimensionedVector("aF", dimAcceleration, aF_);
    }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //