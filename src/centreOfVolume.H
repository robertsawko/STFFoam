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
    Foam::fv::centreOfVolumea

Description
    Object moving under the action of gravity and surrounding fluid.

    \heading Source usage

    Example usage:
    \verbatim
    modelType ; 

    translationModel centreOfVolumea;

    centreOfVolumeaCoeffs 
    {
        mass 0.62;
        apparentMass 0.0961;
        objectName sphere;
    }
    \endverbatim


SourceFiles
    centreOfVolumea.C

\*---------------------------------------------------------------------------*/

#ifndef centreOfVolumea_H
#define centreOfVolumea_H

#include "translationalFrame.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam {
/*---------------------------------------------------------------------------*\
                          Class centreOfVolumea Declaration
\*---------------------------------------------------------------------------*/

class centreOfVolumea : public translationalFrame {
  private:
    // Private Member Functions

    //- Disallow default bitwise copy construct
    centreOfVolumea(const fallingObject &);

    //- Disallow default bitwise assignment
    void operator=(const centreOfVolumea &);

    //- The ID of the patch identifying the front
    label sphereI_;
    scalar mass_;
    scalar apparentMass_;

    virtual vector calculate_acceleration(const volScalarField &,
                                  const volSymmTensorField &);

  public:
    TypeName("centreOfVolumea");

    centreOfVolumea(const fvMesh &, const IOdictionary &);

    virtual ~centreOfVolumea() {}

    virtual void update(const volScalarField &p, const volSymmTensorField &R);
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //