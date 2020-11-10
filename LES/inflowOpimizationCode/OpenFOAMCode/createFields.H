  Info<< "Reading field p\n" << endl;
  volScalarField p
  (
      IOobject
      (
          "p",
          runTime.timeName(),
          mesh,
          IOobject::MUST_READ,
          IOobject::AUTO_WRITE
      ),
      mesh
  );

  Info<< "Reading field U\n" << endl;
  volVectorField U
  (
      IOobject
      (
          "U",
          runTime.timeName(),
          mesh,
          IOobject::MUST_READ,
          IOobject::AUTO_WRITE
      ),
      mesh
  );
  
  Info<< "Reading field pMin\n" << endl;
  volScalarField pMin
  (
      IOobject
      (
          "pMin",
          runTime.timeName(),
          mesh,
          IOobject::MUST_READ,
          IOobject::AUTO_WRITE
      ),
      mesh
  );

#include "createPhi.H"
#include "createMRF.H"
#include "inflowProperties.H"  

// read time step
double timeStep=runTime.deltaT().value();

int j,k,jj,kk,count;
int curTimeIndex=1;

double PI = 3.1415927;
double ran,uni;

// uniform grid resolution
double dy=LY/NY;    // grid size in the y-dir of the virtual uniform mesh 
double dz=LZ/NZ;    // grid size in the z-dir of the virtual uniform mesh 

// streamwise filter legth-scale
double NLX_u = lagT_u/timeStep;
double NLX_v = lagT_v/timeStep;
double NLX_w = lagT_w/timeStep;

// vertical filter legth-scale (yLu = 0.2*xLu  windTunnel slides)
double LY_u = 0.2*Uref*lagT_u;
int NLY_u = ceil(LY_u/dy);
int NLY_v = NLY_u;
int NLY_w = NLY_u;

int NLYW_u = NLY_u;
int NLYW_v = NLY_u;
int NLYW_w = NLY_u;

Info<<"\nVertical length-scales \n"<< "LY_u = "<<LY_u<<"\n"<< endl;
Info<<"\nVertical filter length-scales \n"<< "NLY_u = "<<NLY_u<<"\n"<< endl;

// spanwise filter legth-scale (zLu = 0.3*xLu  windTunnel slides)
double LZ_u = 0.3*Uref*lagT_u;
int NLZ_u = ceil(LZ_u/dz);
int NLZ_v = NLZ_u;
int NLZ_w = NLZ_u;

int NLZW_u = NLZ_u;
int NLZW_v = NLZ_u;
int NLZW_w = NLZ_u;

Info<<"\nSpanwise length-scales \n"<< "LZ_u = "<<LZ_u<<"\n"<< endl;
Info<<"\nSpanwise filter length-scales \n"<< "NLZ_u = "<<NLZ_u<<"\n"<< endl;

//N=2n and summation over -N ~ N 
int NLY2P1_u=1+ NLY_u*4;     //filter size in y-dir for u NLY_u 
int NLZ2P1_u=1+ NLZ_u*4;     //filter size in z-dir for u NLZ_u 
int NLY2P1_v=1+ NLY_v*4;     //filter size in y-dir for v NLY_v 
int NLZ2P1_v=1+ NLZ_v*4;     //filter size in z-dir for v NLZ_v 
int NLY2P1_w=1+ NLY_w*4;     //filter size in y-dir for w NLY_w 
int NLZ2P1_w=1+ NLZ_w*4;     //filter size in z-dir for w NLZ_w 

int NLYW2P1_u=1+ NLYW_u*4;   //filter size in y-dir near the wall NLYW_u 
int NLZW2P1_u=1+ NLZW_u*4;   //filter size in z-dir near the wall NLZW_u
int NLYW2P1_v=1+ NLYW_v*4;   //filter size in y-dir near the wall NLYW_v
int NLZW2P1_v=1+ NLZW_v*4;   //filter size in z-dir near the wall NLZW_v 
int NLYW2P1_w=1+ NLYW_w*4;   //filter size in y-dir near the wall NLYW_w 
int NLZW2P1_w=1+ NLZW_w*4;   //filter size in z-dir near the wall NLZW_w 
  
int NRY=3*NY;     //random array size in the y-dir NY 
int NRZ=3*NZ;     //random array size in the z-dir NZ 

//filter coeffs
double by_u[NLY2P1_u], bz_u[NLZ2P1_u], byz_u[NLY2P1_u][NLZ2P1_u]; 
double by_v[NLY2P1_v], bz_v[NLZ2P1_v], byz_v[NLY2P1_v][NLZ2P1_v]; 
double by_w[NLY2P1_w], bz_w[NLZ2P1_w], byz_w[NLY2P1_w][NLZ2P1_w];

//filter coeffs near the wall
double byw_u[NLYW2P1_u], bzw_u[NLZW2P1_u], byzw_u[NLYW2P1_u][NLZW2P1_u]; 
double byw_v[NLYW2P1_v], bzw_v[NLZW2P1_v], byzw_v[NLYW2P1_v][NLZW2P1_v]; 
double byw_w[NLYW2P1_w], bzw_w[NLZW2P1_w], byzw_w[NLYW2P1_w][NLZW2P1_w]; 

// random numbers  
double rndmx1[NRY][NRZ], rndmy1[NRY][NRZ], rndmz1[NRY][NRZ];

// correaltion functions  
double pSi_u[NY][NZ], pSi_v[NY][NZ], pSi_w[NY][NZ];

//generated velocity
double u1_[NY][NZ], v1_[NY][NZ], w1_[NY][NZ];
double u2_[NY][NZ], v2_[NY][NZ], w2_[NY][NZ];

//piso entries
label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell(p, mesh.solutionDict().subDict("PISO"), pRefCell, pRefValue);

singlePhaseTransportModel laminarTransport(U, phi);

autoPtr<incompressible::turbulenceModel> turbulence
(
    incompressible::turbulenceModel::New(U, phi, laminarTransport)
);

// find zone ID
label cellZoneID = mesh.cellZones().findZoneID(genZone);

const labelList& cells = mesh.cellZones()[cellZoneID];

// Field declaration and initialization
scalarField ycentre(cells.size(), 0);
scalarField zcentre(cells.size(), 0);
scalarField RSuu(cells.size(), 0);
scalarField RSvv(cells.size(), 0);
scalarField RSuv(cells.size(), 0);
scalarField RSww(cells.size(), 0);

vectorField Unew(mesh.cellZones()[cellZoneID].size(), vector(0,0,0));

forAll(cells, celli){

    label cellID = cells[celli];
    
    ycentre[celli] = mesh.C()[cellID].component(vector::Y);
    zcentre[celli] = mesh.C()[cellID].component(vector::Z);
}

// get coefficients, by and bz : these are generated only once whole time

  double sumy_u=0.0;
  double sumy_v=0.0;
  double sumy_w=0.0;

  double sumz_u=0.0;
  double sumz_v=0.0;
  double sumz_w=0.0;

  double sumyw_u=0.0;
  double sumyw_v=0.0;
  double sumyw_w=0.0;

  double sumzw_u=0.0;
  double sumzw_v=0.0;
  double sumzw_w=0.0;

        //sumy_i
        for (j=0; j<NLY2P1_u;j++){
      	sumy_u += Foam::exp(-0.5*PI*fabs(float(j-2*NLY_u)/float(NLY_u)) ) 
      		* Foam::exp(-0.5*PI*fabs(float(j-2*NLY_u)/float(NLY_u)) );
        }
        for (j=0; j<NLY2P1_v;j++){
      	sumy_v += Foam::exp(-0.5*PI*fabs(float(j-2*NLY_v)/float(NLY_v)) ) 
      		* Foam::exp(-0.5*PI*fabs(float(j-2*NLY_v)/float(NLY_v)) );
        }
        for (j=0; j<NLY2P1_w;j++){
      	sumy_w += Foam::exp(-0.5*PI*fabs(float(j-2*NLY_w)/float(NLY_w)) ) 
      		* Foam::exp(-0.5*PI*fabs(float(j-2*NLY_w)/float(NLY_w)) );
        }


	//sumz_i
        for (k=0; k<NLZ2P1_u;k++){
      	sumz_u += Foam::exp(-0.5*PI*fabs(float(k-2*NLZ_u)/float(NLZ_u))) 
      		* Foam::exp(-0.5*PI*fabs(float(k-2*NLZ_u)/float(NLZ_u))); 
        }
        for (k=0; k<NLZ2P1_v;k++){
      	sumz_v += Foam::exp(-0.5*PI*fabs(float(k-2*NLZ_v)/float(NLZ_v))) 
      		* Foam::exp(-0.5*PI*fabs(float(k-2*NLZ_v)/float(NLZ_v))); 
        }
        for (k=0; k<NLZ2P1_w;k++){
      	sumz_w += Foam::exp(-0.5*PI*fabs(float(k-2*NLZ_w)/float(NLZ_w))) 
      		* Foam::exp(-0.5*PI*fabs(float(k-2*NLZ_w)/float(NLZ_w))); 
        }

        //sumyw_i
        for (j=0; j<NLYW2P1_u;j++){
      	sumyw_u += Foam::exp(-0.5*PI*fabs(float(j-2*NLYW_u)/float(NLYW_u)) ) 
      		 * Foam::exp(-0.5*PI*fabs(float(j-2*NLYW_u)/float(NLYW_u)) );
        }
        for (j=0; j<NLYW2P1_v;j++){
      	sumyw_v += Foam::exp(-0.5*PI*fabs(float(j-2*NLYW_v)/float(NLYW_v)) ) 
      		 * Foam::exp(-0.5*PI*fabs(float(j-2*NLYW_v)/float(NLYW_v)) );
        }
        for (j=0; j<NLYW2P1_w;j++){
      	sumyw_w += Foam::exp(-0.5*PI*fabs(float(j-2*NLYW_w)/float(NLYW_w)) ) 
      		 * Foam::exp(-0.5*PI*fabs(float(j-2*NLYW_w)/float(NLYW_w)) );
        }

        //sumzw_i
        for (k=0; k<NLZW2P1_u;k++){
      	sumzw_u += Foam::exp(-0.5*PI*fabs(float(k-2*NLZW_u)/float(NLZW_u)) ) 
      		 * Foam::exp(-0.5*PI*fabs(float(k-2*NLZW_u)/float(NLZW_u)) );
        }
        for (k=0; k<NLZW2P1_v;k++){
      	sumzw_v += Foam::exp(-0.5*PI*fabs(float(k-2*NLZW_v)/float(NLZW_v)) ) 
      		 * Foam::exp(-0.5*PI*fabs(float(k-2*NLZW_v)/float(NLZW_v)) );
        }
        for (k=0; k<NLZW2P1_w;k++){
      	sumzw_w += Foam::exp(-0.5*PI*fabs(float(k-2*NLZW_w)/float(NLZW_w)) ) 
      		 * Foam::exp(-0.5*PI*fabs(float(k-2*NLZW_w)/float(NLZW_w)) );
        }


        sumy_u=Foam::sqrt(sumy_u);
        sumy_v=Foam::sqrt(sumy_v);
        sumy_w=Foam::sqrt(sumy_w);

        sumz_u=Foam::sqrt(sumz_u);
        sumz_v=Foam::sqrt(sumz_v);
        sumz_w=Foam::sqrt(sumz_w);

        sumyw_u=Foam::sqrt(sumyw_u);
        sumyw_v=Foam::sqrt(sumyw_v);
        sumyw_w=Foam::sqrt(sumyw_w);

        sumzw_u=Foam::sqrt(sumzw_u);
        sumzw_v=Foam::sqrt(sumzw_v);
        sumzw_w=Foam::sqrt(sumzw_w);

        // by_i
	  for (j=0;j<NLY2P1_u;j++){
      	by_u[j] = Foam::exp(-0.5*PI*fabs(float(j-2*NLY_u)/float(NLY_u)))/sumy_u;
	  }
	  for (j=0;j<NLY2P1_v;j++){
      	by_v[j] = Foam::exp(-0.5*PI*fabs(float(j-2*NLY_v)/float(NLY_v)))/sumy_v;
	  }
	  for (j=0;j<NLY2P1_w;j++){
      	by_w[j] = Foam::exp(-0.5*PI*fabs(float(j-2*NLY_w)/float(NLY_w)))/sumy_w;
	  }

        //bz_i
	  for (k=0; k<NLZ2P1_u;k++){
      	bz_u[k]= Foam::exp(-0.5*PI*fabs(float(k-2*NLZ_u)/float(NLZ_u)))/sumz_u;  
        }
	  for (k=0; k<NLZ2P1_v;k++){
      	bz_v[k]= Foam::exp(-0.5*PI*fabs(float(k-2*NLZ_v)/float(NLZ_v)))/sumz_v;  
        }
	  for (k=0; k<NLZ2P1_w;k++){
      	bz_w[k]= Foam::exp(-0.5*PI*fabs(float(k-2*NLZ_w)/float(NLZ_w)))/sumz_w;  
        }

        //byz_i
	  for (k=0; k<NLZ2P1_u;k++){
      	for (j=0;j<NLY2P1_u;j++){
      	   byz_u[j][k]=by_u[j]*bz_u[k];
      	}
        }
	  for (k=0; k<NLZ2P1_v;k++){
      	for (j=0;j<NLY2P1_v;j++){
      	   byz_v[j][k]=by_v[j]*bz_v[k];
      	}
        }
	  for (k=0; k<NLZ2P1_w;k++){
      	for (j=0;j<NLY2P1_w;j++){
      	   byz_w[j][k]=by_w[j]*bz_w[k];
      	}
        }

        //byw_i
	  for (j=0;j<NLYW2P1_u;j++){
      	byw_u[j] = Foam::exp(-0.5*PI*fabs(float(j-2*NLYW_u)/float(NLYW_u)))/sumyw_u;
	  }
	  for (j=0;j<NLYW2P1_v;j++){
      	byw_v[j] = Foam::exp(-0.5*PI*fabs(float(j-2*NLYW_v)/float(NLYW_v)))/sumyw_v;
	  }
	  for (j=0;j<NLYW2P1_w;j++){
      	byw_w[j] = Foam::exp(-0.5*PI*fabs(float(j-2*NLYW_w)/float(NLYW_w)))/sumyw_w;
	  }


        //bzw_i
	  for (k=0;k<NLZW2P1_u;k++){
      	bzw_u[k] = Foam::exp(-0.5*PI*fabs(float(k-2*NLZW_u)/float(NLZW_u)))/sumzw_u;
	  }
	  for (k=0;k<NLZW2P1_v;k++){
      	bzw_v[k] = Foam::exp(-0.5*PI*fabs(float(k-2*NLZW_v)/float(NLZW_v)))/sumzw_v;
	  }
	  for (k=0;k<NLZW2P1_w;k++){
      	bzw_w[k] = Foam::exp(-0.5*PI*fabs(float(k-2*NLZW_w)/float(NLZW_w)))/sumzw_w;
	  }


        //byzw_i
	  for (k=0; k<NLZW2P1_u;k++){
      	for (j=0;j<NLYW2P1_u;j++){
      	   byzw_u[j][k]=byw_u[j]*bzw_u[k];
      	}
        }
	  for (k=0; k<NLZW2P1_v;k++){
      	for (j=0;j<NLYW2P1_v;j++){
      	   byzw_v[j][k]=byw_v[j]*bzw_v[k];
      	}
        }
	  for (k=0; k<NLZW2P1_w;k++){
      	for (j=0;j<NLYW2P1_w;j++){
      	   byzw_w[j][k]=byw_w[j]*bzw_w[k];
      	}
        }


// filter coefficient END

// Read Exp Data
IFstream ifU(runTime.constant()/"UInlet");
IFstream ifuuBar(runTime.constant()/"uuBarInlet");
IFstream ifuvBar(runTime.constant()/"uvBarInlet");
IFstream ifvvBar(runTime.constant()/"vvBarInlet");
IFstream ifwwBar(runTime.constant()/"wwBarInlet");

graph graphU
(
      "graphU",
      "y",
      "U",
      ifU
);

graph graphuuBar
(
      "graphuuBar",
      "y",
      "uuBar",
      ifuuBar
);

graph graphuvBar
(
      "graphuvBar",
      "y",
      "uvBar",
      ifuvBar
);

graph graphvvBar
(
      "graphvvBar",
      "y",
      "vvBar",
      ifvvBar
);

graph graphwwBar
(
      "graphwwBar",
      "y",
      "wwBar",
      ifwwBar
);


// interpolate data
Info<< "\nInterpolate data \n" << endl;
scalarField expU = interpolateXY
(
      ycentre,
      graphU.x(),
      graphU.y()
);
scalarField expuuBar = interpolateXY
(
      ycentre,
      graphuuBar.x(),
      graphuuBar.y()
);
scalarField expuvBar = interpolateXY
(
      ycentre,
      graphuvBar.x(),
      graphuvBar.y()
);
scalarField expvvBar = interpolateXY
(
      ycentre,
      graphvvBar.x(),
      graphvvBar.y()
);
scalarField expwwBar = interpolateXY
(
      ycentre,
      graphwwBar.x(),
      graphwwBar.y()
);

// read experimental data 
forAll(cells, celli){

    RSuu[celli] = Foam::sqrt(SMALL+ expuuBar[celli]);
    RSuv[celli] = expuvBar[celli] / RSuu[celli];
    RSvv[celli] = Foam::sqrt(SMALL + expvvBar[celli] 
                 -Foam::sqr( RSuv[celli] ) );
    RSww[celli] = Foam::sqrt(SMALL + expwwBar[celli]);
}



