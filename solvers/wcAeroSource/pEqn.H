volScalarField rAU(1.0/UEqn.A());
volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p));
surfaceScalarField phiHbyA("phiHbyA", fvc::flux(HbyA));

if (pimple.ddtCorr())
{
    phiHbyA += MRF.zeroFilter(fvc::interpolate(rAU)*fvc::ddtCorr(U, phi, Uf));
}


MRF.makeRelative(phiHbyA);

if (p.needReference())
{
    fvc::makeRelative(phiHbyA, U);
    adjustPhi(phiHbyA, U, p);
    fvc::makeAbsolute(phiHbyA, U);
}

tmp<volScalarField> rAtU(rAU);

if (pimple.consistent())
{
    rAtU = 1.0/max(1.0/rAU - UEqn.H1(), 0.1/rAU);
    phiHbyA +=
        fvc::interpolate(rAtU() - rAU)*fvc::snGrad(p)*mesh.magSf();
    HbyA -= (rAU - rAtU())*fvc::grad(p);
}

if (pimple.nCorrPISO() <= 1)
{
    tUEqn.clear();
}

// Update the pressure BCs to ensure flux consistency
constrainPressure(p, U, phiHbyA, rAtU(), MRF);

// Non-orthogonal pressure corrector loop
while (pimple.correctNonOrthogonal())
{
    /// APAC Pressure equation
    fvScalarMatrix pEqn(
	-inv_c2 * ( fvm::ddt(p) + (U&fvc::grad(p)) )
	+ fvm::laplacian(rAtU(), p) == fvc::div(phiHbyA)
    );

    /// APAC buffer zone technique
    if (applyBufferZone && runTime.value() > buffStartTime_ && pRefHeader.registerObject() )
    {
        const volScalarField& pBuffRef(mesh.lookupObject<volScalarField>(pBufferRefName_));
        tmp<fvMatrix<scalar>> tBufferZonePEqn
        (
            fvm::Sp(blendFactor_, p) - blendFactor_*pBuffRef
        );

        fvMatrix<scalar>&  bufferPEqn = tBufferZonePEqn.ref();
        pEqn -= inv_c2/runTime.deltaT() * bufferPEqn; // 
    }

    pEqn.setReference(pRefCell, pRefValue);

    pEqn.solve(mesh.solver(p.select(pimple.finalInnerIter())));

    if (pimple.finalNonOrthogonalIter())
    {
        phi = phiHbyA - pEqn.flux();
    }
}

#include "continuityErrs.H"

// Explicitly relax pressure for momentum corrector
p.relax();

U = HbyA - rAtU*fvc::grad(p);
U.correctBoundaryConditions();
fvOptions.correct(U);

// Correct Uf if the mesh is moving
fvc::correctUf(Uf, U, phi);

// Make the fluxes relative to the mesh motion
fvc::makeRelative(phi, U);
