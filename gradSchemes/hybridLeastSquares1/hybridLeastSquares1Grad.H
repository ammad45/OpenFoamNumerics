/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Class
        Foam::fv::hybridLeastSquares1

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
                default hybridLeastSquares1;
            }

        The blending field must be precomputed using computeGradSchemeBlending.

    SourceFiles
        hybridLeastSquares1Grad.C

\*---------------------------------------------------------------------------*/

#ifndef hybridLeastSquares1_H
#define hybridLeastSquares1_H

#include "gradScheme.H"

namespace Foam
{
namespace fv
{

template<class Type>
class hybridLeastSquares1
:
    public fv::gradScheme<Type>
{
    // Private Member Functions

        //- No copy construct
        hybridLeastSquares1(const hybridLeastSquares1&) = delete;

        //- No copy assignment
        void operator=(const hybridLeastSquares1&) = delete;

public:

    //- Runtime type information
    TypeName("hybridLeastSquares1");

    // Constructors

        //- Construct from mesh
        hybridLeastSquares1(const fvMesh& mesh)
        :
            gradScheme<Type>(mesh)
        {}

        //- Construct from mesh and Istream (unused contents)
        hybridLeastSquares1(const fvMesh& mesh, Istream& schemeData);

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
    #include "hybridLeastSquares1Grad.C"
#endif

#endif

// ************************************************************************* //

