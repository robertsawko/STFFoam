    // Solve the Momentum equation

    MRF.correctBoundaryVelocity(U);
    frameOfReference->correctBoundaryVelocity(U);

    fvVectorMatrix UEqn
    (
        fvm::ddt(rho, U) + fvm::div(rhoPhi, U)
      + rho * frameOfReference->aF()
      + MRF.DDt(rho, U)
      + fvc::div(UdmModel.tauDm())
      + turbulence->divDevRhoReff(U)
     ==
        fvOptions(rho, U)
    );

    UEqn.relax();

    fvOptions.constrain(UEqn);

    if (pimple.momentumPredictor())
    {
        solve
        (
            UEqn
          ==
            fvc::reconstruct
            (
                (
                  - ghf*fvc::snGrad(rho)
                  - fvc::snGrad(p_rgh)
                )*mesh.magSf()
            )
        );

        fvOptions.correct(U);
    }
