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

#include "hybridFaceLeastSquaresVectors.H"
#include "fvcSurfaceIntegrate.H"
#include "volFields.H"
#include "regIOobject.H"
#include "cpuTime.H"
#include <exception>
#include <mutex>

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(hybridFaceLeastSquaresVectors, 0);
}


// * * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * //

Foam::hybridFaceLeastSquaresVectors::hybridFaceLeastSquaresVectors
(
    const fvMesh& mesh
)
:
    MeshObject_type(mesh),
    pVectors_
    (
        IOobject
        (
            "HybridFaceLeastSquaresP",
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
            "HybridFaceLeastSquaresN",
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

Foam::hybridFaceLeastSquaresVectors::~hybridFaceLeastSquaresVectors()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::hybridFaceLeastSquaresVectors::calcLeastSquaresVectors()
{
    cpuTime calcTimer;
    Info << "Calculating hybrid face least square gradient vectors" << nl;

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

    const surfaceScalarField& weight = magSf;

    // Set up temporary storage for mean of edge centres
    vectorField meanPoints(mesh.nCells(), vector::zero);
    scalarField weightSum = fvc::surfaceSum(weight)().internalField();

    surfaceVectorField Ce
    (
        IOobject
        (
            "Ce",
            mesh.pointsInstance(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        mesh,
        dimLength
    );

    forAll(owner, facei)
    {
        const label own = owner[facei];
        const label nei = neighbour[facei];

        Ce[facei] = w[facei]*C[own] + (1 - w[facei])*C[nei];
        meanPoints[own] += weight[facei]*Ce[facei];
        meanPoints[nei] += weight[facei]*Ce[facei];
    }

    forAll(w.boundaryField(), patchi)
    {
        const fvsPatchScalarField& pw = w.boundaryField()[patchi];
        fvsPatchVectorField& pCe = Ce.boundaryFieldRef()[patchi];
        const fvsPatchScalarField& pWeight = weight.boundaryField()[patchi];
        const fvsPatchVectorField& pSf = Sf.boundaryField()[patchi];
        const fvsPatchScalarField& pMagSf = magSf.boundaryField()[patchi];
        const labelUList& faceCells = pw.patch().faceCells();

        tmp<vectorField> tpd = pw.patch().delta();
        const vectorField& pd = tpd();

        if (pw.coupled())
        {
            forAll(pw, patchFacei)
            {
                pCe[patchFacei] =
                    C[faceCells[patchFacei]]
                  + (1 - pw[patchFacei])*pd[patchFacei];
                meanPoints[faceCells[patchFacei]] +=
                    pWeight[patchFacei]*pCe[patchFacei];
            }
        }
        else
        {
            forAll(pw, patchFacei)
            {
                const vector& delta = pd[patchFacei];
                const vector unitArea = pSf[patchFacei]/pMagSf[patchFacei];
                const scalar nonOrthDelta =
                    1.0/max(unitArea & delta, 0.05*mag(delta));

                pCe[patchFacei] =
                    C[faceCells[patchFacei]]
                  + (pSf[patchFacei]/(pMagSf[patchFacei]*nonOrthDelta));

                meanPoints[faceCells[patchFacei]] +=
                    pWeight[patchFacei]*pCe[patchFacei];
            }
        }
    }

    meanPoints /= weightSum;

    // Set up temporary storage for the dd tensor (before inversion)
    tensorField dd(mesh.nCells(), Zero);

    scalarField norm =
        fvc::surfaceSum(weight*magSqr(mesh.delta()))->primitiveField();

    forAll(owner, facei)
    {
        const label own = owner[facei];
        const label nei = neighbour[facei];

        dd[own] +=
            weight[facei]/norm[own]
           *Ce[facei]*(Ce[facei] - meanPoints[own]);
        dd[nei] +=
            weight[facei]/norm[nei]
           *Ce[facei]*(Ce[facei] - meanPoints[nei]);
    }

    surfaceVectorField::Boundary& pVectorsBf =
        pVectors_.boundaryFieldRef();

    forAll(pVectorsBf, patchi)
    {
        const fvsPatchScalarField& pw = w.boundaryField()[patchi];
        const fvsPatchVectorField& pCe = Ce.boundaryField()[patchi];
        const fvsPatchScalarField& pWeight = weight.boundaryField()[patchi];

        const fvPatch& p = pw.patch();
        const labelUList& faceCells = p.patch().faceCells();

        forAll(pw, patchFacei)
        {
            dd[faceCells[patchFacei]] +=
                pWeight[patchFacei]/norm[faceCells[patchFacei]]
               *pCe[patchFacei]
               *(pCe[patchFacei] - meanPoints[faceCells[patchFacei]]);
        }
    }

    // Stabilise before inversion on 1D or 2D meshes
    if (mesh.nSolutionD() < 3)
    {
        vector emptyCmpt = (Vector<label>::one - mesh.solutionD())/2;
        dd += cmptMultiply(emptyCmpt*emptyCmpt, tensor::I);
    }

    // Invert the dd tensors - including failsafe checks
    const tensorField invDd(inv(dd));

    // Revisit all faces and calculate the pVectors_ and nVectors_ vectors
    forAll(owner, facei)
    {
        const label own = owner[facei];
        const label nei = neighbour[facei];

        const vector d(C[nei] - C[own]);
        const scalar t = mag((Cf[facei] - C[own]) & d)/(mag(d) + VSMALL);

        const scalar bOwn = betaPtr
            ? max(min((*betaPtr)[own], scalar(1)), scalar(0))
            : scalar(0);
        const scalar bNei = betaPtr
            ? max(min((*betaPtr)[nei], scalar(1)), scalar(0))
            : scalar(0);

        const vector lsOwn =
            weight[facei]/norm[own]*(invDd[own] & Ce[facei]);
        const vector lsNei =
            weight[facei]/norm[nei]*(invDd[nei] & Ce[facei]);
        const vector gaussOwn =
            (1.0 - w[facei])*t*(invDd[own] & Sf[facei]);
        const vector gaussNei =
            w[facei]*t*(invDd[nei] & Sf[facei]);

        pVectors_[facei] = (1.0 - bOwn)*lsOwn + bOwn*gaussOwn;
        nVectors_[facei] = -((1.0 - bNei)*lsNei + bNei*gaussNei);
    }

    forAll(pVectorsBf, patchi)
    {
        fvsPatchVectorField& patchLsP = pVectorsBf[patchi];

        const fvsPatchScalarField& pw = w.boundaryField()[patchi];
        const fvsPatchScalarField& pWeight = weight.boundaryField()[patchi];
        const fvsPatchVectorField& pCe = Ce.boundaryField()[patchi];
        const fvsPatchVectorField& pCf = Cf.boundaryField()[patchi];
        const fvsPatchVectorField& pSf = Sf.boundaryField()[patchi];

        const fvPatch& p = pw.patch();
        const labelUList& faceCells = p.faceCells();

        tmp<vectorField> tpd = p.delta();
        const vectorField& pd = tpd();

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
            const scalar t =
                mag((pCf[patchFacei] - C[faceCells[patchFacei]]) & d)
               /(mag(d) + VSMALL);

            const vector lsFace =
                pWeight[patchFacei]/norm[faceCells[patchFacei]]
               *(invDd[faceCells[patchFacei]] & pCe[patchFacei]);

            const vector gaussFace =
                t*(invDd[faceCells[patchFacei]] & pSf[patchFacei]);

            const vector gaussPart = pw.coupled()
                ? (1.0 - pw[patchFacei])*gaussFace
                : gaussFace;

            patchLsP[patchFacei] =
                (1.0 - bFace)*lsFace + bFace*gaussPart;
        }
    }

    Info << "Finished calculating hybrid face least square gradient vectors in "
        << calcTimer.elapsedCpuTime() << " s" << nl;
}


bool Foam::hybridFaceLeastSquaresVectors::movePoints()
{
    calcLeastSquaresVectors();
    return true;
}


// ************************************************************************* //
