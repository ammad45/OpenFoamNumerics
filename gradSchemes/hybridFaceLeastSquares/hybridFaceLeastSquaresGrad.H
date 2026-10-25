/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Class
        Foam::fv::hybridFaceLeastSquares

    Group
        grpFvGradSchemes

    Description
        Least-squares gradient scheme that uses precomputed blending factors
        (gradSchemeBlending, 0..1) when constructing least-squares vectors.
        If the blending field is missing, the standard least-squares vectors
        are used.

    Usage
        In fvSchemes:
            gradSchemes
            {
                default hybridFaceLeastSquares;
            }

        The blending field must be precomputed using computeGradSchemeBlending.

    SourceFiles
        hybridFaceLeastSquaresGrad.C

\*---------------------------------------------------------------------------*/

#ifndef hybridFaceLeastSquares_H
#define hybridFaceLeastSquares_H

#include "gradScheme.H"

namespace Foam
{
namespace fv
{

template<class Type>
class hybridFaceLeastSquares
:
    public fv::gradScheme<Type>
{
    // Private Member Functions

        //- No copy construct
        hybridFaceLeastSquares(const hybridFaceLeastSquares&) = delete;

        //- No copy assignment
        void operator=(const hybridFaceLeastSquares&) = delete;

public:

    //- Runtime type information
    TypeName("hybridFaceLeastSquares");

    // Constructors

        //- Construct from mesh
        hybridFaceLeastSquares(const fvMesh& mesh)
        :
            gradScheme<Type>(mesh)
        {}

        //- Construct from mesh and Istream (unused contents)
        hybridFaceLeastSquares(const fvMesh& mesh, Istream& schemeData);

    // Member Functions

        //- Calculate gradient
        virtual tmp< GeometricField< typename outerProduct<vector, Type>::type, fvPatchField, volMesh> > calcGrad
        (
            const GeometricField<Type, fvPatchField, volMesh>& vsf,
            const word& name
        ) const;
};

} // End namespace fv
} // End namespace Foam

#ifdef NoRepository
    #include "hybridFaceLeastSquaresGrad.C"
#endif

#endif

// ************************************************************************* //

