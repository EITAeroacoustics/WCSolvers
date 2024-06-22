
const dictionary &bufferZoneDict = mesh.solutionDict().optionalSubDict("BUFFERZONE");

bool applyBufferZone = bufferZoneDict.lookupOrDefault<bool>("applyBufferZone", false);

volScalarField blendFactor_(
    IOobject(
        "blendFactor_",
        mesh.time().timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE),
    mesh,
    dimensionedScalar("blendFactor_", dimless, 0.0));



scalar buffStartTime_ = runTime.value()  +  bufferZoneDict.lookupOrDefault<scalar>("startTime", 0.0);
//scalar frequency_ = bufferZoneDict.lookupOrDefault<scalar>("frequency", 5000.0);
scalar alpha_ = bufferZoneDict.lookupOrDefault<scalar>("alpha", 1.0);
scalar beta_ = bufferZoneDict.lookupOrDefault<scalar>("beta", 3.0);

word zoneType_ = bufferZoneDict.lookupOrDefault<word>("zoneType", "Unknown zone type");
word pBufferRefName_ = bufferZoneDict.lookupOrDefault<word>("pRefField", "pMean");
word UBufferRefName_ = bufferZoneDict.lookupOrDefault<word>("URefField", "UMean");

IOobject pRefHeader
(
    pBufferRefName_,
    runTime.timeName(),
    mesh,
    IOobject::NO_READ
);

IOobject URefHeader
(
    UBufferRefName_,
    runTime.timeName(),
    mesh,
    IOobject::NO_READ
);

if (applyBufferZone)
{
    {
        Info << "[*] Reading bufferZone parameters..." << endl;
        // 1 timeStep makes sure that Reference field is registered.
        //buffStartTime_ = Foam::max(buffStartTime_, runTime.startTime().value()) + runTime.deltaT().value();
        Info << "    applyBufferZone: " << applyBufferZone << endl;
        Info << "    startTime: " << buffStartTime_ << endl;
        //Info << "    frequency: " << frequency_ << endl;
        Info << "    alpha: " << alpha_ << endl;
        Info << "    beta: " << beta_ << endl;
        if (zoneType_ == "sphere")
        {
            scalar r1_ = bufferZoneDict.lookupOrDefault<scalar>("innerRadius", 0.0);
            scalar r2_ = bufferZoneDict.lookupOrDefault<scalar>("outerRadius", 0.0);
            point centre_ = bufferZoneDict.lookupOrDefault<point>("centre", vector::zero);
            Info << "    zoneType: " << zoneType_ << endl;
            Info << "    innerRadius: " << r1_ << endl;
            Info << "    outerRadius: " << r2_ << endl;
            Info << "    centre: " << centre_ << endl;

            forAll(mesh.C(), celli)
            {
                scalar d = mag(mesh.C()[celli] - centre_);

                if (d < r1_)
                {
                    blendFactor_[celli] = 0.0;
                }
                else if ((d >= r1_) && (d <= r2_))
                {
                    blendFactor_[celli] = Foam::pow(Foam::mag(d - r1_) / (r2_ - r1_), beta_);
                }
                else if (d > r2_)
                {
                    blendFactor_[celli] = 1.0;
                }
            }
        }
        else if (zoneType_ == "cylinder")
        {
            Info << "Not implemented bufferZone type:" << zoneType_ << endl;
            return 0;
        }
        else if (zoneType_ == "box")
        {
            point innerMin;
            point innerMax;
            point outerMin;
            point outerMax;
            bufferZoneDict.lookup("innerBoxMin") >> innerMin;
            bufferZoneDict.lookup("innerBoxMax") >> innerMax;
            bufferZoneDict.lookup("outerBoxMin") >> outerMin;
            bufferZoneDict.lookup("outerBoxMax") >> outerMax;
            Info << "    zoneType: " << zoneType_ << endl;
            boundBox innerBox, outerBox;
            innerBox.add(innerMax); innerBox.add(innerMin);
            outerBox.add(outerMax); outerBox.add(outerMin);
            Info << "    innerBox: " << innerBox << endl;
            Info << "    outerBox: " << outerBox << endl;
            if (innerBox.empty() || outerBox.empty() || !outerBox.contains(innerBox))
            {
                Info << "Failed to initialize bufferZone. Please check box descriptions." << endl;
                return 0;
            }
            innerMin = innerBox.min();
            innerMax = innerBox.max();
            outerMin = outerBox.min();
            outerMax = outerBox.max();

            forAll(mesh.C(), celli)
            {
                // positive region
                if (mesh.C()[celli].x() >= innerMax.x() && mesh.C()[celli].x() <= outerMax.x())
                {
                    blendFactor_[celli] = Foam::max(blendFactor_[celli], Foam::pow((mesh.C()[celli].x() - innerMax.x()) / (outerMax.x() - innerMax.x()), beta_));
                }
                else if (mesh.C()[celli].x() > outerMax.x())
                {
                    blendFactor_[celli] = Foam::max(blendFactor_[celli], 1.0);
                }
                if (mesh.C()[celli].y() >= innerMax.y() && mesh.C()[celli].y() <= outerMax.y())
                {
                    blendFactor_[celli] = Foam::max(blendFactor_[celli], Foam::pow((mesh.C()[celli].y() - innerMax.y()) / (outerMax.y() - innerMax.y()), beta_));
                }
                else if (mesh.C()[celli].y() > outerMax.y())
                {
                    blendFactor_[celli] = Foam::max(blendFactor_[celli], 1.0);
                }
                if (mesh.C()[celli].z() >= innerMax.z() && mesh.C()[celli].z() <= outerMax.z())
                {
                    blendFactor_[celli] = Foam::max(blendFactor_[celli], Foam::pow((mesh.C()[celli].z() - innerMax.z()) / (outerMax.z() - innerMax.z()), beta_));
                }
                else if (mesh.C()[celli].z() > outerMax.z())
                {
                    blendFactor_[celli] = Foam::max(blendFactor_[celli], 1.0);
                }

                if (mesh.C()[celli].x() >= outerMin.x() && mesh.C()[celli].x() <= innerMin.x())
                {
                    blendFactor_[celli] = Foam::max(blendFactor_[celli], Foam::pow((innerMin.x() - mesh.C()[celli].x()) / (innerMin.x() - outerMin.x()), beta_));
                }
                else if (mesh.C()[celli].x() < outerMin.x())
                {
                    blendFactor_[celli] = Foam::max(blendFactor_[celli], 1.0);
                }
                if (mesh.C()[celli].y() >= outerMin.y() && mesh.C()[celli].y() <= innerMin.y())
                {
                    blendFactor_[celli] = Foam::max(blendFactor_[celli], Foam::pow((innerMin.y() - mesh.C()[celli].y()) / (innerMin.y() - outerMin.y()), beta_));
                }
                else if (mesh.C()[celli].y() < outerMin.y())
                {
                    blendFactor_[celli] = Foam::max(blendFactor_[celli], 1.0);
                }
                if (mesh.C()[celli].z() >= outerMin.z() && mesh.C()[celli].z() <= innerMin.z())
                {
                    blendFactor_[celli] = Foam::max(blendFactor_[celli], Foam::pow((innerMin.z() - mesh.C()[celli].z()) / (innerMin.z() - outerMin.z()), beta_));
                }
                else if (mesh.C()[celli].z() < outerMin.z())
                {
                    blendFactor_[celli] = Foam::max(blendFactor_[celli], 1.0);
                }
            }
        }
        else
        {
            Info << "Unknown bufferZone type:" << zoneType_ << endl;
            return 0;
        }

        forAll(mesh.C(), celli)
        {
            blendFactor_[celli] = alpha_ * blendFactor_[celli];
        }
        
        blendFactor_.correctBoundaryConditions();
    }
    Info << "[*] Done. " << endl
         << endl;
}
else
{
    Info << "[*] bufferZone is not applied. " << endl;
}
