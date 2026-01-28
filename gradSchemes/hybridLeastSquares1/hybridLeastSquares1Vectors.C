/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2020-2021 OpenCFD Ltd.
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

#include "hybridLeastSquares1Vectors.H"
#include "volFields.H"
#include "regIOobject.H"
#include "cpuTime.H"
#include <exception>
#include <mutex>

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(hybridLeastSquares1Vectors, 0);
}


// * * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * //

Foam::hybridLeastSquares1Vectors::hybridLeastSquares1Vectors(const fvMesh& mesh)
:
    MeshObject_type(mesh),
    pVectors_
    (
        IOobject
        (
            "HybridLeastSquaresP",
            mesh.pointsInstance(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        mesh,
        dimensionedVector(dimless/dimLength, Zero)
    ),
    nVectors_
    (
        IOobject
        (
            "HybridLeastSquaresN",
            mesh.pointsInstance(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        mesh,
        dimensionedVector(dimless/dimLength, Zero)
    )
{
    calcLeastSquaresVectors();
}


// * * * * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * * //

Foam::hybridLeastSquares1Vectors::~hybridLeastSquares1Vectors()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::hybridLeastSquares1Vectors::calcLeastSquaresVectors()
{
    cpuTime calcTimer;
    Info << "Calculating hybrid least square gradient vectors" << nl;

    const fvMesh& mesh = this->mesh();

    // Set local references to mesh data
    const labelUList& owner = mesh.owner();
    const labelUList& neighbour = mesh.neighbour();

    const volVectorField& C = mesh.C();
    const surfaceScalarField& w = mesh.weights();
    const surfaceScalarField& magSf = mesh.magSf();
    const surfaceVectorField& Cf = mesh.Cf();
    const surfaceVectorField& Sf = mesh.Sf();

    // Optional blending field (0..1) to modify least-squares vectors
    const word blendingFieldName("gradSchemeBlending");
    const volScalarField* blendingPtr =
        mesh.findObject<volScalarField>(blendingFieldName);

    if (!blendingPtr)
    {
        const word inst =
            mesh.time().findInstance(mesh.dbDir(), blendingFieldName);

        IOobject fieldHeader
        (
            blendingFieldName,
            inst,
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE,
            IOobject::REGISTER
        );

        if (fieldHeader.typeHeaderOk<volScalarField>(true, true, false))
        {
            auto* vfPtr = new volScalarField(fieldHeader, mesh);
            regIOobject::store(vfPtr);
            blendingPtr = vfPtr;
        }
    }

    const scalarField* betaPtr =
        blendingPtr ? &blendingPtr->primitiveField() : nullptr;


    // Set up temporary storage for the dd tensor (before inversion)
    tensorField dd(mesh.nCells(), Zero);

    forAll(owner, facei)
    {
        const label own = owner[facei];
        const label nei = neighbour[facei];
        
        const vector d(C[nei] - C[own]);
        const scalar magSfByMagSqr = magSf[facei]/magSqr(d);
        const tensor wddLsq(magSfByMagSqr*sqr(d));
        
        const scalar t = mag((Cf[facei] - C[own]) & d)/(mag(d) + VSMALL);
        const tensor wddGauss(t*(Sf[facei]*d));
        
        const scalar bOwn = betaPtr
            ? max(min((*betaPtr)[own], scalar(1)), scalar(0))
            : scalar(0);
        const scalar bNei = betaPtr
            ? max(min((*betaPtr)[nei], scalar(1)), scalar(0))
            : scalar(0);                                                                                   

        dd[own] +=
            (1.0 - bOwn)*(1.0 - w[facei])*wddLsq
          + bOwn*wddGauss;
        dd[nei] +=
            (1.0 - bNei)*w[facei]*wddLsq
          + bNei*wddGauss;
    }


    surfaceVectorField::Boundary& pVectorsBf =
        pVectors_.boundaryFieldRef();

    forAll(pVectorsBf, patchi)
    {
        const fvsPatchScalarField& pw = w.boundaryField()[patchi];
        const fvsPatchScalarField& pMagSf = magSf.boundaryField()[patchi];
        const fvsPatchVectorField& pCf = Cf.boundaryField()[patchi];
        const fvsPatchVectorField& pSf = Sf.boundaryField()[patchi];

        const fvPatch& p = pw.patch();
        const labelUList& faceCells = p.patch().faceCells();

        // Build the d-vectors
        const vectorField pd(p.delta());

        if (pw.coupled())
        {
            forAll(pd, patchFacei)
            {
                const vector& d = pd[patchFacei];

                const scalar bFace = betaPtr
                    ? max
                    (
                        min((*betaPtr)[faceCells[patchFacei]], scalar(1)),
                        scalar(0)
                    )
                    : scalar(0);
                const scalar magSfByMagSqr =
                    pMagSf[patchFacei]/magSqr(d);
                const scalar t =
                    mag((pCf[patchFacei] - C[faceCells[patchFacei]]) & d)
                   /(mag(d) + VSMALL);
                const tensor wddGauss(t*(pSf[patchFacei]*d));

                dd[faceCells[patchFacei]] +=
                    (1.0 - bFace)
                   *((1 - pw[patchFacei])*magSfByMagSqr)*sqr(d)
                  + bFace*wddGauss;
            }
        }
        else
        {
            const fvsPatchScalarField& pNonOrthDelta =
                mesh.nonOrthDeltaCoeffs().boundaryField()[patchi];

            forAll(pd, patchFacei)
            {
                const vector d =
                    pSf[patchFacei]
                   /(
                        pMagSf[patchFacei]
                       *pNonOrthDelta[patchFacei]
                    );

                const scalar bFace = betaPtr
                    ? max
                    (
                        min((*betaPtr)[faceCells[patchFacei]], scalar(1)),
                        scalar(0)
                    )
                    : scalar(0);
                const scalar magSfByMagSqr =
                    pMagSf[patchFacei]/magSqr(d);
                const scalar t =
                    mag((pCf[patchFacei] - C[faceCells[patchFacei]]) & d)
                   /(mag(d) + VSMALL);
                const tensor wddGauss(t*(pSf[patchFacei]*d));

                dd[faceCells[patchFacei]] +=
                    (1.0 - bFace)*magSfByMagSqr*sqr(d)
                  + bFace*wddGauss;
            }
        }
    }


    // Invert the dd tensors - including failsafe checks
    const tensorField invDd(inv(dd));


    // Revisit all faces and calculate the pVectors_ and nVectors_ vectors
    forAll(owner, facei)
    {
        const label own = owner[facei];
        const label nei = neighbour[facei];

        const vector d(C[nei] - C[own]);
        const scalar magSfByMagSqrd = magSf[facei]/magSqr(d);
        const scalar t = mag((Cf[facei] - C[own]) & d)/(mag(d) + VSMALL);

        const scalar bOwn = betaPtr
            ? max(min((*betaPtr)[own], scalar(1)), scalar(0))
            : scalar(0);
        const scalar bNei = betaPtr
            ? max(min((*betaPtr)[nei], scalar(1)), scalar(0))
            : scalar(0);

        pVectors_[facei] =
            (1.0 - w[facei])*((1.0 - bOwn)*magSfByMagSqrd*(invDd[own] & d)
            +bOwn*t*(invDd[own] & Sf[facei]));
        nVectors_[facei] =
            -w[facei]*((1.0 - bNei)*magSfByMagSqrd*(invDd[nei] & d)
            +bNei*t*(invDd[nei] & Sf[facei]));
    }

    forAll(pVectorsBf, patchi)
    {
        fvsPatchVectorField& patchLsP = pVectorsBf[patchi];

        const fvsPatchScalarField& pw = w.boundaryField()[patchi];
        const fvsPatchScalarField& pMagSf = magSf.boundaryField()[patchi];
        const fvsPatchVectorField& pCf = Cf.boundaryField()[patchi];
        const fvsPatchVectorField& pSf = Sf.boundaryField()[patchi];

        const fvPatch& p = pw.patch();
        const labelUList& faceCells = p.faceCells();

        // Build the d-vectors
        const vectorField pd(p.delta());

        if (pw.coupled())
        {
            forAll(pd, patchFacei)
            {
                const vector& d = pd[patchFacei];

                const scalar bFace = betaPtr
                    ? max
                    (
                        min((*betaPtr)[faceCells[patchFacei]], scalar(1)),
                        scalar(0)
                    )
                    : scalar(0);
                const scalar magSfByMagSqrd =
                    pMagSf[patchFacei]/magSqr(d);
                const scalar t =
                    mag((pCf[patchFacei] - C[faceCells[patchFacei]]) & d)
                   /(mag(d) + VSMALL);

                patchLsP[patchFacei] =
                    (1.0 - pw[patchFacei])
                   *((1.0 - bFace)*magSfByMagSqrd
                   *(invDd[faceCells[patchFacei]] & d)
                   + bFace*t*(invDd[faceCells[patchFacei]] & pSf[patchFacei]));
            }
        }
        else
        {
            const fvsPatchScalarField& pNonOrthDelta =
                mesh.nonOrthDeltaCoeffs().boundaryField()[patchi];

            forAll(pd, patchFacei)
            {
                const vector d =
                    pSf[patchFacei]
                   /(
                        pMagSf[patchFacei]
                       *pNonOrthDelta[patchFacei]
                    );

                const scalar bFace = betaPtr
                    ? max
                    (
                        min((*betaPtr)[faceCells[patchFacei]], scalar(1)),
                        scalar(0)
                    )
                    : scalar(0);
                const scalar magSfByMagSqrd =
                    pMagSf[patchFacei]/magSqr(d);
                const scalar t =
                    mag((pCf[patchFacei] - C[faceCells[patchFacei]]) & d)
                   /(mag(d) + VSMALL);

                patchLsP[patchFacei] =
                    (1.0 - bFace)*magSfByMagSqrd
                   *(invDd[faceCells[patchFacei]] & d)
                  + bFace*t*(invDd[faceCells[patchFacei]] & pSf[patchFacei]);
            }
        }
    }

    Info << "Finished calculating hybrid least square gradient vectors in "
        << calcTimer.elapsedCpuTime() << " s" << nl;
}


bool Foam::hybridLeastSquares1Vectors::movePoints()
{
    calcLeastSquaresVectors();
    return true;
}


// ************************************************************************* //
