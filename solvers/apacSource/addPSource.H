if (addPointSource) 
{
    // add source term of p
    if (pSourceType == "monopole")
    {
        forAll(mesh.C(), celli) // loop all cell centres
        {
            vector r_ = mesh.C()[celli] - centre;
            scalar d = Foam::mag(r_);
            p[celli] -= ampl * Foam::sin(omeg * runTime.value()) * Foam::exp(-Foam::pow(d, 2) / sig2);
        }
    }
    else if (pSourceType == "dipole")
    {
        forAll(mesh.C(), celli) // loop all cell centres
        {
            vector r_ = mesh.C()[celli] - centre;
            scalar d = Foam::mag(r_);
            p[celli] -= ampl * Foam::sin(omeg * runTime.value()) * Foam::exp(-Foam::pow(d, 2) / sig2) * r_.y() / cc * 1.41421356237310; // sqrt(2)
        }
    }
    else if (pSourceType == "quadrupole")
    {
        forAll(mesh.C(), celli) // loop all cell centres
        {
            vector r_ = mesh.C()[celli] - centre;
            scalar d = Foam::mag(r_);
            p[celli] -= ampl * Foam::sin(omeg * runTime.value()) * Foam::exp(-Foam::pow(d, 2) / sig2)  * r_.y()  * r_.x() / cc / cc * 2;
        }
    }
    else if (pSourceType == "none")
    {
        /* do nothing  */
    }
    else
    {
        Info << "Unknown point source type: " << pSourceType << endl;
        return 0;
    }
}
