// Solve the Momentum equation

frameOfReference->correctBoundaryVelocity(U);

tmp<fvVectorMatrix> UEqn
(
    fvm::ddt(U) + fvm::div(phi, U)
  + frameOfReference->aF()
  + turbulence->divDevReff(U)
 ==
    fvOptions(U)
);

UEqn().relax();

fvOptions.constrain(UEqn());

if (pimple.momentumPredictor())
{
    solve(UEqn() == -fvc::grad(p));

    fvOptions.correct(U);
}
