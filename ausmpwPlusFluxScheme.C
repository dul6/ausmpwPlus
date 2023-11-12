/*---------------------------------------------------------------------------*\

    HiSA: High Speed Aerodynamic solver

    Copyright (C) 2014-2018 Oliver Oxtoby - CSIR, South Africa
    Copyright (C) 2014-2018 Johan Heyns - CSIR, South Africa
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

\*---------------------------------------------------------------------------*/

#include "ausmpwPlusFluxScheme.H"
#include "addToRunTimeSelectionTable.H"
#include "bound.H"
#include "fvcSurfaceReconstruct.H"
#include "cellFaceFunctions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(ausmpwPlusFluxScheme, 0);
addToRunTimeSelectionTable(fluxScheme, ausmpwPlusFluxScheme, dictionary);


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

ausmpwPlusFluxScheme::ausmpwPlusFluxScheme
(
    const dictionary& dict,
    const psiThermo& thermo,
    const volScalarField& rho,
    const volVectorField& U,
    const volVectorField& rhoU,
    const volScalarField& rhoE
)
:
    fluxScheme(typeName, dict),
    mesh_(U.mesh()),
    thermo_(thermo),
    rho_(rho),
    U_(U),
    rhoU_(rhoU),
    rhoE_(rhoE),
    dict_(dict)
{}


// * * * * * * * * * * * * * * * * Destructors * * * * * * * * * * * * * * * //

ausmpwPlusFluxScheme::~ausmpwPlusFluxScheme()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::ausmpwPlusFluxScheme::calcFlux(surfaceScalarField& phi, surfaceVectorField& phiUp, surfaceScalarField& phiEp, surfaceVectorField& Up)
{
    const volScalarField& p = thermo_.p();

    tmp<surfaceVectorField> U_L, U_R;
    fvc::surfaceReconstruct(U_, U_L, U_R, "reconstruct(U)");
    tmp<surfaceScalarField> phi_L = U_L()&mesh_.Sf();
    tmp<surfaceScalarField> phi_R = U_R()&mesh_.Sf();

    U_L->rename("U_L");
    U_R->rename("U_R");
    tmp<surfaceVectorField> U_L_rel(U_L.ref());
    tmp<surfaceVectorField> U_R_rel(U_R.ref());
/*
    // Flux relative to mesh movement
    if (mesh_.moving())
    {
        fvc::makeRelative(phi_L.ref(), U_);
        fvc::makeRelative(phi_R.ref(), U_);
    }
*/
    surfaceScalarField un_L = phi_L/mesh_.magSf();
    surfaceScalarField un_R = phi_R/mesh_.magSf();

    tmp< volScalarField > gamma = thermo_.gamma();
    tmp< volScalarField > H
    (
        (max(rhoE_/rho_,dimensionedScalar("0", rhoE_.dimensions()/rho_.dimensions(), SMALL)) +
         max(p/rho_,dimensionedScalar("0", p.dimensions()/rho_.dimensions(), SMALL)))
    );
    H->rename("H");

    tmp< volScalarField > Hrel(H.ref());

    tmp<surfaceScalarField> H_L, H_R;
    fvc::surfaceReconstruct(Hrel, H_L, H_R, "reconstruct(T)");

/*
    if (mesh_.moving())
    {
        Hrel = H() - 0.5*(U_&U_);
        volVectorField Urel(U_);
        Urel -= fvc::reconstruct(fvc::meshPhi(U_));
        Hrel.ref() += 0.5*(Urel&Urel);
    }
*/
    tmp< volScalarField > c = sqrt(2.0*(gamma()-1.0)/(gamma()+1.0)*Hrel());
    c->rename("c");
    gamma.clear();

    tmp<surfaceScalarField> c_L, c_R;
    fvc::surfaceReconstruct(c, c_L, c_R, "reconstruct(T)");
    c_L = sqr(c_L())/max(c_L(), un_L);
    c_R = sqr(c_R())/max(c_R(),-un_R);
    tmp< surfaceScalarField > c_face(min(c_L(),c_R()));
    c_L.clear();
    c_R.clear();

    // Critical Mach number
    tmp< surfaceScalarField > Mach_L(un_L/c_face());

    // Split Mach numbers
    tmp<surfaceScalarField> Mach_plus_L =
        calcForEachFace
        (
            [](const scalar& MLf)
            {
                if (mag(MLf) < 1.0)
                {
                    //scalar ML2p =  0.25*sqr(MLf+1);
                    //scalar ML2m = -0.25*sqr(MLf-1);
                    //return ML2p;             // beta = 0
                    //return ML2p*(1 - 2*ML2m);  // beta = 1/8
                    return 0.25*sqr(MLf+1);
                }
                else
                {
                    return max(MLf, 0);
                }
            },
            Mach_L()
        );

    // Pressure flux
    tmp<surfaceScalarField> p_plus_L =
        calcForEachFace
        (
            [](const scalar& MLf)
            {
                if (mag(MLf) < 1.0)
                {
                    //scalar ML2p =  0.25*sqr(MLf+1);
                    //scalar ML2m = -0.25*sqr(MLf-1);
                    //return ML2p*(2 - MLf - 3*MLf*ML2m);  //alpha = 3/16
                    return 0.25*sqr(MLf+1)*(2-MLf);
                }
                else
                {
                    return (MLf > 0 ? 1.0 : 0.0);
                }
            },
            Mach_L()
        );

    tmp< surfaceScalarField > Mach_R(un_R/c_face());

    // Split Mach numbers
    tmp<surfaceScalarField> Mach_minus_R =
        calcForEachFace
        (
            [](const scalar& MRf)
            {
                if (mag(MRf) < 1.0)
                {
                    //scalar MR2m = -0.25*sqr(MRf-1);
                    //scalar MR2p =  0.25*sqr(MRf+1);
                    //return MR2m;             // beta = 0
                    //return MR2m*(1 + 2*MR2p);  // beta = 1/8
                    return -0.25*sqr(MRf-1);
                }
                else
                {
                    return min(MRf, 0);
                }
            },
            Mach_R()
        );

    // Pressure flux
    tmp<surfaceScalarField> p_minus_R =
        calcForEachFace
        (
            [](const scalar& MRf)
            {
                if (mag(MRf) < 1.0)
                {
                    //scalar MR2m = -0.25*sqr(MRf-1);
                    //scalar MR2p =  0.25*sqr(MRf+1);
                    //return MR2m*(-2 - MRf + 3*MRf*MR2p);  //alpha = 3/16
                    return 0.25*sqr(MRf-1)*(2+MRf);
                }
                else
                {
                    return (MRf < 0 ? 1.0 : 0.0);
                }
            },
            Mach_R()
        );

    tmp<surfaceScalarField> p_L, p_R;
    fvc::surfaceReconstruct(p, p_L, p_R, "reconstruct(rho)");

    tmp<surfaceScalarField> omega = 1-pow(min(p_L()/p_R(), p_R()/p_L()), 3);

    tmp<surfaceScalarField> rho_L, rho_R;
    fvc::surfaceReconstruct(rho_, rho_L, rho_R, "reconstruct(rho)");

    tmp< surfaceScalarField > Mach_1_2 = Mach_plus_L() + Mach_minus_R();

    surfaceScalarField p_1_2 = p_plus_L()*p_L() + p_minus_R()*p_R();

    tmp<surfaceScalarField> p_L_F(p_L()/p_1_2);
    tmp<surfaceScalarField> f_L =
        calcForEachFace
        (
            [](const scalar& MLf, const scalar& PLf)
            {
                if (mag(MLf) < 1.0)
                {
                    //scalar ML2p =  0.25*sqr(MLf+1);
                    //scalar ML2m = -0.25*sqr(MLf-1);
                    //return ML2p;             // beta = 0
                    return PLf-1;  // beta = 1/8
                }
                else
                {
                    return 0.0;
                }
            },
            Mach_L(),
            p_L_F()
        );//
    p_L_F.clear();
    Mach_L.clear();

    tmp< surfaceScalarField > p_R_F(p_R()/p_1_2);
    tmp<surfaceScalarField> f_R =
        calcForEachFace
        (
            [](const scalar& MRf, const scalar& PRf)
            {
                if (mag(MRf) < 1.0)
                {
                    return PRf-1;  // beta = 1/8
                }
                else
                {
                    return 0.0;
                }
            },
            Mach_R(),
            p_R_F()
        );//
    p_R_F.clear();
    Mach_R.clear();

#if OPENFOAM >= 1712
    Mach_1_2->setOriented(true);
#endif

    tmp<surfaceScalarField> Mach_plus_L_B1 = Mach_plus_L()+Mach_minus_R()*((1-omega())*(1+f_R())-f_L());
    tmp<surfaceScalarField> Mach_plus_L_B2 = Mach_plus_L()*omega()*(1+f_L());

    tmp<surfaceScalarField> Mach_minus_R_B1 = Mach_minus_R()*omega()*(1+f_R());
    tmp<surfaceScalarField> Mach_minus_R_B2 = Mach_minus_R()+Mach_plus_L()*((1-omega())*(1+f_L())-f_R());

    p_L.clear();
    p_R.clear();
    p_plus_L.clear();
    p_minus_R.clear();
    f_L.clear();
    f_R.clear();
    Mach_plus_L.clear();
    Mach_minus_R.clear();

    tmp<surfaceVectorField> U_f = surfaceFieldSelect(U_L_rel, U_R_rel, Mach_1_2(), 0);
    tmp<surfaceScalarField> Mach_plus_L_B = surfaceFieldSelect(Mach_plus_L_B1, Mach_plus_L_B2, Mach_1_2(), 0);
    tmp<surfaceScalarField> Mach_minus_R_B = surfaceFieldSelect(Mach_minus_R_B1, Mach_minus_R_B2, Mach_1_2(), 0);
    surfaceScalarField L = Mach_plus_L_B()*rho_L()*c_face();
    surfaceScalarField R = Mach_minus_R_B()*rho_R()*c_face();

    //surfaceScalarField rhoa_LR  = Mach_1_2()*c_face()*surfaceFieldSelect(rho_L, rho_R, Mach_1_2(), 0);
    surfaceScalarField rhoa_LR  = (L + R);
    //surfaceVectorField rhoaU_LR = rhoa_LR*U_f();
    surfaceVectorField rhoaU_LR = (L*U_L() + R*U_R());
    // volScalarField ee("ee",rhoE_/rho_-0.5*magSqr(U_));
    // surfaceScalarField rhoah_LR = rhoa_LR*(fvc::surfaceReconstruct(ee, Mach_1_2(), "reconstruct(T)") + 0.5*magSqr(U_f())) + p_1_2*Mach_1_2()*c_face();
    // NOTE: According to Liou, enthalpy should be interpolated.
    //surfaceScalarField rhoah_LR = rhoa_LR*(fvc::surfaceReconstruct(H(), Mach_1_2(), "reconstruct(T)"));
    surfaceScalarField rhoah_LR = (L*H_L() + R*H_R());
    // Face velocity for sigmaDotU (turbulence term)
    Up = U_f()*mesh_.magSf();
    U_f.clear();
    H.clear();
    Mach_plus_L_B.clear();
    Mach_minus_R_B.clear();
//    c_face.clear();

    phi = rhoa_LR*mesh_.magSf();
    rhoaU_LR.setOriented(true);
    phiUp = rhoaU_LR*mesh_.magSf() + p_1_2*mesh_.Sf();
    phiEp = rhoah_LR*mesh_.magSf();

    if (mesh_.moving())
    {
        phiEp += p_1_2 * fvc::meshPhi(U_);
        //phiEp += fvc::meshPhi(U_)*fvc::surfaceReconstruct(p, Mach_1_2(), "reconstruct(T)");
        // Ensure consistent interpolation with pressure term above
        //phiEp += fvc::meshPhi(U_)*fvc::surfaceReconstruct(rho_, Mach_1_2(), "reconstruct(rho)")*fvc::surfaceReconstruct((p/rho_)(), Mach_1_2(), "reconstruct(T)");
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
