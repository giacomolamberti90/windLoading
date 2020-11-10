Info<< "Reading inflow properties\n" << endl;

    IOdictionary inflowProperties
    (
        IOobject
        (
            "inflowProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );

    word genZone  
    (                                                      	
        inflowProperties.lookupOrDefault<word>("genZone", "genPlane")
    );
    
    double  refD
    (
        inflowProperties.lookupOrDefault("refD", 1.0)
    );

    double  LY
    (
        inflowProperties.lookupOrDefault("LY", 2.0)
    );

    double  LZ
    (
        inflowProperties.lookupOrDefault("LZ", 3.5)
    );

    int NY
    (
        inflowProperties.lookupOrDefault("NY", 100)
    );

    int  NZ
    (
        inflowProperties.lookupOrDefault("NZ", 100)
    );

    double  lagT_u
    (
     	inflowProperties.lookupOrDefault("lagT_u", 0.05)
    );
    
    double  lagT_v
    (
     	inflowProperties.lookupOrDefault("lagT_v", 0.025)
    );
    
    double  lagT_w
    (
     	inflowProperties.lookupOrDefault("lagT_w", 0.025)
    );

    double  Uref
    (
     	inflowProperties.lookupOrDefault("Uref", 10.28)
    );    
    
    double  windDirection
    (
     	inflowProperties.lookupOrDefault("windDirection", 0)
    );    
