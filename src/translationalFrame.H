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
    Solid-body 6-DoF acceleration source

    \heading Source usage

    Example usage:
    \verbatim
    SBM
    {
        type            translationalFrame;
        active          yes;

        translationalFrameCoeffs
        {
            timeDataFileName "constant/acceleration-terms.dat";
        }
    }
    \endverbatim


SourceFiles
    translationalFrame.C

\*---------------------------------------------------------------------------*/

#ifndef translationalFrame_H
#define translationalFrame_H

#include "fvMesh.H"
#include "dimensionedTypes.H"
#include <vector>
#include "uniformFixedVelocityFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                          Class translationalFrame Declaration
\*---------------------------------------------------------------------------*/

class translationalFrame {
private:
  // Private Member Functions

  //- Disallow default bitwise copy construct
  translationalFrame(const translationalFrame &);

  //- Disallow default bitwise assignment
  void operator=(const translationalFrame &);
  IOdictionary dict_;
  const fvMesh& mesh_;
  vector VF_;
  word UName_;

  //- The ID of the patch identifying the front
  label sphereI_;
  scalar mass_;
  scalar apparentMass_;

  std::vector<uniformFixedVelocityFvPatchField*> registeredPatches;

public:
  translationalFrame(const fvMesh &mesh);

  virtual ~translationalFrame() {}

  void registerVelocity(volVectorField& U); 
  void correctBoundaryVelocity(volVectorField& U) const; 
  void update(const volScalarField& p, const volSymmTensorField& R); 
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
