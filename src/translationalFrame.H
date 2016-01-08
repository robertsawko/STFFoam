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
    Foam::fv::translationalFrame

Description
    Non-inertial frame supporting non-uniform acceleration. The type of
    acceleration can be specified via a motion model. Configured through
    constant/STFProperties


SourceFiles
    translationalFrame.C

\*---------------------------------------------------------------------------*/

#ifndef translationalFrame_H
#define translationalFrame_H

#include "functionObjectFile.H"
#include "fvMesh.H"
#include "dimensionedTypes.H"
#include "Switch.H"
#include <vector>
#include "frameAwareBoundary.H"

#include "autoPtr.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam {
/*---------------------------------------------------------------------------*\
                          Class translationalFrame Declaration
\*---------------------------------------------------------------------------*/

class translationalFrame : public functionObjectFile {
  private:
    //- Disallow default bitwise copy construct
    translationalFrame(const translationalFrame &);

    //- Disallow default bitwise assignment
    void operator=(const translationalFrame &);

  protected:
    IOdictionary dict_;
    const fvMesh &mesh_;
    vector VF_;
    vector aF_;
    Switch log_;
    word UName_;

    std::vector<frameAwareBoundary *> registeredPatches;

    virtual vector calculate_acceleration(const volScalarField &,
                                  const volSymmTensorField &) = 0;

  public:
    TypeName("translationalFrame");

    declareRunTimeSelectionTable(autoPtr,
                                 translationalFrame,
                                 dictionary,
                                 (const fvMesh &mesh, const IOdictionary &dict),
                                 (mesh, dict));

    // Selector
    static autoPtr<translationalFrame>
    New(const fvMesh &mesh);
    
    translationalFrame(const fvMesh &, const IOdictionary &);

    virtual ~translationalFrame() {}

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
