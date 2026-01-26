/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Class
        Foam::fv::hybridLeastSquares

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
                default hybridLeastSquares;
            }

        The blending field must be precomputed using computeGradSchemeBlending.

    SourceFiles
        hybridLeastSquaresGrad.C

\*---------------------------------------------------------------------------*/

#ifndef hybridLeastSquares_H
#define hybridLeastSquares_H

#include "gradScheme.H"

namespace Foam
{
namespace fv
{

template<class Type>
class hybridLeastSquares
:
    public fv::gradScheme<Type>
{
    // Private Member Functions

        //- No copy construct
        hybridLeastSquares(const hybridLeastSquares&) = delete;

        //- No copy assignment
        void operator=(const hybridLeastSquares&) = delete;

public:

    //- Runtime type information
    TypeName("hybridLeastSquares");

    // Constructors

        //- Construct from mesh
        hybridLeastSquares(const fvMesh& mesh)
        :
            gradScheme<Type>(mesh)
        {}

        //- Construct from mesh and Istream (unused contents)
        hybridLeastSquares(const fvMesh& mesh, Istream& schemeData);

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
    #include "hybridLeastSquaresGrad.C"
#endif

#endif

// ************************************************************************* //

