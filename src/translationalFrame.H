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

#include "fvOption.H"
#include "dimensionedTypes.H"
#include <vector>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                          Class translationalFrame Declaration
\*---------------------------------------------------------------------------*/

class translationalFrame : public option {

protected:
  // Protected data

  //- Run-time selectable acceleration model
  vector VF_;

  word UName_;

  dimensionedVector g0_;

  std::vector<label> patchIDs;

private:
  // Private Member Functions

  //- Disallow default bitwise copy construct
  translationalFrame(const translationalFrame &);

  //- Disallow default bitwise assignment
  void operator=(const translationalFrame &);

  //- Source term to momentum equation
  // template <class RhoFieldType>
  // void addSup(const RhoFieldType &rho, fvMatrix<vector> &eqn,
  //            const label fieldi);

public:
  //- Runtime type information
  TypeName("translationalFrame");

  // Constructors

  //- Construct from components
  translationalFrame(const word &name, const word &modelType,
                   const dictionary &dict, const fvMesh &mesh);

  //- Destructor
  virtual ~translationalFrame() {}

  // Member Functions
  //

  virtual void correct(volVectorField& field);

  //- Source term to momentum equation
  // virtual void addSup(fvMatrix<vector> &eqn, const label fieldi);

  //- Source term to compressible momentum equation
  // virtual void addSup(const volScalarField &rho, fvMatrix<vector> &eqn, const
  // label fieldi);

  //- Read dictionary
  virtual bool read(const dictionary &dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//#ifdef NoRepository
    //#include "translationalFrameTemplates.C"
//#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //