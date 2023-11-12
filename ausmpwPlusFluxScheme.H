/*---------------------------------------------------------------------------*\

    HiSA: High Speed Aerodynamic solver

    Copyright (C) 2014-2016 Johan Heyns - CSIR, South Africa
    Copyright (C) 2014-2016 Oliver Oxtoby - CSIR, South Africa
    Copyright (C) 1991-2008 OpenCFD Ltd.

-------------------------------------------------------------------------------
License
    This file is part of HiSA.

    HiSA is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    HiSA is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with HiSA.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::ausmpwPlusFluxScheme

Description
    AUSMPW+ Scheme

    References:
    \verbatim
        Kim K H, Kim C, Rho O H. (1998).
	Accurate Computations of Hypersonic Flows Using AUSMPW+ Scheme
	and Shock-Aligned Grid Technique.
	Albuquerque: 29th AIAA, Fluid Dynamics Conference.
    \endverbatim

Note
    This code is not suitable for the moving mesh.

SourceFiles
    ausmpwPlusFluxScheme.C

Authors
    Lei Du


\*---------------------------------------------------------------------------*/

#ifndef ausmpwPlusFluxScheme_H
#define ausmpwPlusFluxScheme_H

#include "fluxScheme.H"
#include "psiThermo.H"
//#include "upwind.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
              Class ausmpwPlusFluxScheme Declaration
\*---------------------------------------------------------------------------*/

class ausmpwPlusFluxScheme
:
    public fluxScheme
{

protected:
    // Protected data
        const fvMesh& mesh_;
        const psiThermo& thermo_;
        const volScalarField& rho_;
        const volVectorField& U_;
        const volVectorField& rhoU_;
        const volScalarField& rhoE_;
        const dictionary& dict_;

public:
    //- Runtime type information
    TypeName("AUSMPWPLUS");


	// Constructors

        //- Construct from components
        ausmpwPlusFluxScheme
        (
            const dictionary& dict,
            const psiThermo& thermo,
            const volScalarField& rho,
            const volVectorField& U,
            const volVectorField& rhoU,
            const volScalarField& rhoE
        );


    // Destructor

        ~ausmpwPlusFluxScheme();


    // Member Functions

        //- Calculate mass, momentum and energy flux
        virtual void calcFlux(surfaceScalarField& phi, surfaceVectorField& phiUp, surfaceScalarField& phiEp, surfaceVectorField& Up);

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
