const dictionary &pointSourceDict = mesh.solutionDict().optionalSubDict("PSOURCE");
bool addPointSource = pointSourceDict.lookupOrDefault<bool>("addPointSource", false);

word pSourceType = pointSourceDict.lookupOrDefault<word>("sourceType", "none");

scalar ampl = pointSourceDict.lookupOrDefault<scalar>("amplitude", 0.0);
scalar freq = pointSourceDict.lookupOrDefault<scalar>("frequency", 0.0);
scalar cc = pointSourceDict.lookupOrDefault<scalar>("cc", 0.1);
point centre = pointSourceDict.lookupOrDefault<point>("centre", vector::zero);

Info << "[*] Reading point source parameters..." << endl;

if (addPointSource)
{
    Info << "addPointSource " << addPointSource << endl;
    Info << "sourceType " << pSourceType << endl;
    Info << "amplitude " << ampl << endl;
    Info << "frequency " << freq << endl;
    Info << "centre " << centre << endl;
    Info << "const c " << cc << endl;
    Info << "[*] Done. " << endl;
}
else
{
    Info << "[*] No point souce will be added." << endl;
}

const scalar pi = constant::mathematical::pi;
const scalar omeg = freq * 2 * pi;
const scalar sig2 = 2 * Foam::pow(cc, 2);
