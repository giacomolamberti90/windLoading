/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    pisoFoam

Description
    Transient solver for incompressible, turbulent flow, using the PISO
    algorithm.

    Sub-models include:
    - turbulence modelling, i.e. laminar, RAS or LES
    - run-time selectable MRF and finite volume options, e.g. explicit porosity

\*---------------------------------------------------------------------------*/

#include "turbulentTransportModel.H"
#include "fvOptions.H"

#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "turbulenceModel.H"

#include "IFstream.H"
#include "OFstream.H"
#include "graph.H"
#include "interpolateXY.H"
#include "PrimitivePatchInterpolation.H"
#include "Random.H"
#include "pisoControl.H"
#include "fvOptions.H"

//model constant
#define cw_ PI/4   

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "postProcess.H"
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControl.H"
    #include "createFields.H"
    #include "createFvOptions.H"
    #include "initContinuityErrs.H"

    turbulence->validate();

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;
	
	// To adjust time step 
        const bool adjustTimeStep =
        runTime.controlDict().lookupOrDefault("adjustTimeStep", false);

        scalar maxCo =
        runTime.controlDict().lookupOrDefault<scalar>("maxCo", 1.0);
        
        scalar maxDeltaT = 
        runTime.controlDict().lookupOrDefault<scalar>("maxDeltaT", GREAT);	

        #include "CourantNo.H"
        #include "setDeltaT.H"

	#include "Ugen.H"

        // Pressure-velocity PISO corrector
        {
            #include "UEqn.H"
            
            // --- PISO loop
            while (piso.correct())
            {
                #include "pEqn.H"
            }
        }
        
        laminarTransport.correct();
        turbulence->correct();

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
        
        double timeStep = runTime.deltaT().value();
        double currentTime = runTime.value();
            
        if ((currentTime - timeStep) - int(currentTime - timeStep) < 1e-12)
        {
            dimensionedScalar p0("p0", dimensionSet(0,2,-2,0,0,0,0), 1e6);
            pMin = p0;
        }        

	forAll(mesh.C(), celli)
	{ 
	    if (p[celli] < pMin[celli])
	    {
	        pMin[celli] = p[celli];
	    }
	}            
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
