// Random Field generation
for(k=0; k < NRZ; k++){

      for(j=0;j<NRY;j++){
  
	    ran = 0.0;
          for(count=1; count<=12; ++count){
      	uni = ((double)rand() / ((double)(RAND_MAX)+(double)(1)) );
      	ran += uni;
          }
          rndmx1[j][k]=ran -6.0;


	    ran = 0.0;
          for(count=1; count<=12; ++count){
      	uni = ((double)rand() / ((double)(RAND_MAX)+(double)(1)) );
      	ran += uni;
          }
          rndmy1[j][k]=ran -6.0;


	    ran = 0.0;
          for(count=1; count<=12; ++count){
      	uni = ((double)rand() / ((double)(RAND_MAX)+(double)(1)) );
      	ran += uni;
          }
          rndmz1[j][k]=ran -6.0;
      }
} 		  
// Random Field generation End

// Spatial correlation (y, z-dir) implementation
for(kk = 0; kk < NZ; kk++){
      for( jj = 0; jj < NY; jj++){ 
  
          pSi_u[jj][kk] =0.0;
          pSi_v[jj][kk] =0.0;
          pSi_w[jj][kk] =0.0;

          if(jj*dy < 0.1*refD || jj*dy > 1.9*refD){//near wall region

      	    for(k = 0 ; k < NLZW2P1_u; k++){
      		  for (j = 0; j < NLYW2P1_u; j++) {
      		    pSi_u[jj][kk] += byzw_u[j][k]*rndmx1[j+jj+NY-2*NLYW_u][k+kk+NZ-2*NLZW_u];
      		  }
      	    }
      	    for(k = 0 ; k < NLZW2P1_v; k++){
      		  for (j = 0; j < NLYW2P1_v; j++) {
      		    pSi_v[jj][kk] += byzw_v[j][k]*rndmy1[j+jj+NY-2*NLYW_v][k+kk+NZ-2*NLZW_v];
      		  }
      	    }
      	    for(k = 0 ; k < NLZW2P1_w; k++){
      		  for (j = 0; j < NLYW2P1_w; j++) {
      		    pSi_w[jj][kk] += byzw_w[j][k]*rndmz1[j+jj+NY-2*NLYW_w][k+kk+NZ-2*NLZW_w];
      		  }
      	    }
           }

            else{
      	    for(k = 0 ; k < NLZ2P1_u; k++){
      		  for (j = 0; j < NLY2P1_u; j++) {
      		    pSi_u[jj][kk] += byz_u[j][k]*rndmx1[j+jj+NY-2*NLY_u][k+kk+NZ-2*NLZ_u];
      		  }
      	    }
      	    for(k = 0 ; k < NLZ2P1_v; k++){
      		  for (j = 0; j < NLY2P1_v; j++) {
      		    pSi_v[jj][kk] += byz_v[j][k]*rndmy1[j+jj+NY-2*NLY_v][k+kk+NZ-2*NLZ_v];
      		  }
      	    }
      	    for(k = 0 ; k < NLZ2P1_w; k++){
      		  for (j = 0; j < NLY2P1_w; j++) {
      		    pSi_w[jj][kk] += byz_w[j][k]*rndmz1[j+jj+NY-2*NLY_w][k+kk+NZ-2*NLZ_w];
      		  }
      	    }
           }
      }
}

// Spatially correlated signal is ready on 'UNIFORM' mesh


// correlation (x-dir) implementation
for(k=0; k < NZ; k++){
      for(j=0;j<NY;j++){
  
          u2_[j][k] = u1_[j][k]* Foam::exp(-cw_/NLX_u ) 
      		+ pSi_u[j][k]*Foam::sqrt(1.0-Foam::exp(-2.0*cw_/NLX_u ) );

          v2_[j][k] = v1_[j][k]* Foam::exp(-cw_/NLX_v ) 
      		+ pSi_v[j][k]*Foam::sqrt(1.0-Foam::exp(-2.0*cw_/NLX_v ) );

          w2_[j][k] = w1_[j][k]* Foam::exp(-cw_/NLX_w )  
      		+ pSi_w[j][k]*Foam::sqrt(1.0-Foam::exp(-2.0*cw_/NLX_w ) );        
      }
} 
// correlation (x-dir) implementation END

// cyclic boundary manipulation (not perfect..)
for(j=0;j<NY;j++){
  
   u2_[j][0]=u2_[j][NZ-1];
   v2_[j][0]=v2_[j][NZ-1];
   w2_[j][0]=w2_[j][NZ-1];

   u2_[j][1]=u2_[j][NZ-2];
   v2_[j][1]=v2_[j][NZ-2];
   w2_[j][1]=w2_[j][NZ-2];
} 

// Impose digital filter in cellZone
forAll(cells, celli)
{
    // Mapping from 'UNIFORM' to 'NON-UNIFORM' mesh &
    jj=(int)(0.5+ycentre[celli]/dy);
    kk=(int)(0.5+zcentre[celli]/dz+NZ/2);

    // Modify fields in the cellzone
    Unew[celli] = expU[celli]*vector(1,0,0) + 
      		  vector(RSuu[celli]*u2_[jj][kk],		
      		         RSuv[celli]*u2_[jj][kk] + RSvv[celli]*v2_[jj][kk],
      		         RSww[celli]*w2_[jj][kk]);
      		         
    // Rotate velocity field
    Unew[celli] = vector
                  ( 
                      Unew[celli][0]*Foam::cos(windDirection*PI/180) - Unew[celli][2]*Foam::sin(windDirection*PI/180),
                      Unew[celli][1],
                      Unew[celli][0]*Foam::sin(windDirection*PI/180) + Unew[celli][2]*Foam::cos(windDirection*PI/180)
                  );
}

for(k=0; k < NZ; k++){
      for(j=0;j<NY;j++){

          // Save current signal for the Next time step
          u1_[j][k] = u2_[j][k];
          v1_[j][k] = v2_[j][k];
          w1_[j][k] = w2_[j][k];   
      }
}

