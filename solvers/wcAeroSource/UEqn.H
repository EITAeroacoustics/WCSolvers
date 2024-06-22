// Solve the Momentum equation

MRF.correctBoundaryVelocity(U);

tmp<fvVectorMatrix> tUEqn
(
    fvm::ddt(U) + fvm::div(phi, U)
  + MRF.DDt(U)
  + turbulence->divDevReff(U)
 ==
    fvOptions(U)
);
fvVectorMatrix& UEqn = tUEqn.ref();

/// APAC Buffer zone technique
//
if (applyBufferZone  && runTime.value() > buffStartTime_ && URefHeader.registerObject())
{
    const volVectorField& UBuffRef(mesh.lookupObject<volVectorField>(UBufferRefName_));
    tmp<fvMatrix<vector>> tBufferZoneUEqn
    (
        fvm::Sp(blendFactor_, U) - blendFactor_ * UBuffRef
    );
    fvMatrix<vector> &bufferUEqn = tBufferZoneUEqn.ref();
    UEqn += 1/runTime.deltaT()*bufferUEqn;
}



UEqn.relax();

fvOptions.constrain(UEqn);

if (pimple.momentumPredictor())
{
    solve(UEqn == -fvc::grad(p));

    fvOptions.correct(U);
}
