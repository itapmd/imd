
/******************************************************************************
*
* IMD -- The ITAP Molecular Dynamics Program
*
* Copyright 1996-2006 Institute for Theoretical and Applied Physics,
* University of Stuttgart, D-70550 Stuttgart
*
******************************************************************************
/
/******************************************************************************
*
* imd_laser_profiles.c -- various intensity profiles for tem_xy laser modes
*
******************************************************************************/

/****************************************************************
* $Revision$
* $Date$
*****************************************************************/


/* Funkctions for the various laser modes */



double laser_intensity_profile_laguerre_00(double x, double y, double z)
{

  double rho;
  
  if( x >= laser_offset ){
    
    rho = ( y - laser_sigma_w_y ) * ( y - laser_sigma_w_y ) + ( z - laser_sigma_w_z ) * ( z - laser_sigma_w_z );
    rho *= laser_sigma_w0;
    
    return exp(-rho);
  }
  
  else
  {
    return 0.0;
  }
}

double laser_intensity_profile_laguerre_01(double x, double y, double z)
{

  double phi;
  double rho;
  
  if( x >= laser_offset ){
    
    rho = ( y - laser_sigma_w_y ) * ( y - laser_sigma_w_y ) + ( z - laser_sigma_w_z ) * ( z - laser_sigma_w_z );
    rho *= laser_sigma_w0;
    
    phi = atan2( ( y - laser_sigma_w_y ) , ( z - laser_sigma_w_z ) );

    return 2.0 * rho * cos(phi) * cos(phi) * exp(- rho) ;
    
  
    
  }
  
  else
  {
    return 0.0;
  }
}

double laser_intensity_profile_laguerre_02(double x, double y, double z)
{

  double phi;
  double rho;

  if( x >= laser_offset ){
    
    rho = ( y - laser_sigma_w_y ) * ( y - laser_sigma_w_y ) + ( z - laser_sigma_w_z ) * ( z - laser_sigma_w_z );
    rho *= 2.0 * laser_sigma_w0;
    phi = atan2( ( y - laser_sigma_w_y ) , ( z - laser_sigma_w_z ) );
  
    return rho * rho * cos(2.0 * phi) * cos(2.0 * phi) * exp(- rho) / ( 2.0 * M_PI );
    
  
    
  }
  
  else
  {
    return 0.0;
  }
}

double laser_intensity_profile_laguerre_03(double x, double y, double z)
{

  double phi;
  double rho;
  
  if( x >= laser_offset ){
    
    rho = ( y - laser_sigma_w_y ) * ( y - laser_sigma_w_y ) + ( z - laser_sigma_w_z ) * ( z - laser_sigma_w_z );
    rho *= 2.0 * laser_sigma_w0;
    phi = atan2( ( y - laser_sigma_w_y ) , ( z - laser_sigma_w_z ) );

    return rho * rho * rho * cos(3.0 * phi) * cos(3.0 * phi) * exp(- rho) / ( 6.0 * M_PI ); 
    
  }
  
  else
  {
    return 0.0;
  }
}


double laser_intensity_profile_laguerre_10(double x, double y, double z)
{

  double rho;
  double laguerre;
  
  if( x >= laser_offset ){
    
    rho = ( y - laser_sigma_w_y ) * ( y - laser_sigma_w_y ) + ( z - laser_sigma_w_z ) * ( z - laser_sigma_w_z );
    rho *= 2.0 * laser_sigma_w0;
    
    laguerre = - rho +1.0;
    laguerre *= laguerre;
    
    return laguerre * exp(- rho);
    
  }
  
  else
  {
    return 0.0;
  }
}

double laser_intensity_profile_laguerre_11(double x, double y, double z)
{

  double phi;
  double rho;
  double laguerre;
  
  if( x >= laser_offset ){
    
    rho = ( y - laser_sigma_w_y ) * ( y - laser_sigma_w_y ) + ( z - laser_sigma_w_z ) * ( z - laser_sigma_w_z );
    rho *= 2.0 * laser_sigma_w0;
    
    phi = atan2( ( y - laser_sigma_w_y ) , ( z - laser_sigma_w_z ) );

    laguerre = - rho + 2.0;
    laguerre *= laguerre;
    
    return rho * laguerre *  cos(phi) * cos(phi) * exp(- rho) / ( 2.0 * M_PI );
    
  
    
  }
  
  else
  {
    return 0.0;
  }
}

double laser_intensity_profile_laguerre_12(double x, double y, double z)
{

  double phi;
  double rho;
  double laguerre;
  
  if( x >= laser_offset ){
    
    rho = ( y - laser_sigma_w_y ) * ( y - laser_sigma_w_y ) + ( z - laser_sigma_w_z ) * ( z - laser_sigma_w_z );
    rho *= 2.0 * laser_sigma_w0;
    
    phi = atan2( ( y - laser_sigma_w_y ) , ( z - laser_sigma_w_z ) );

    laguerre = -rho + 3.0;
    laguerre *= laguerre;

    return rho * rho * laguerre * cos(2.0 * phi) * cos(2.0 * phi) * exp(- rho) / ( 6.0 * M_PI );
    
  
    
  }
  
  else
  {
    return 0.0;
  }
}


double laser_intensity_profile_laguerre_13(double x, double y, double z)
{

  double phi;
  double rho;
  double laguerre;
  
  if( x >= laser_offset ){
    
    rho = ( y - laser_sigma_w_y ) * ( y - laser_sigma_w_y ) + ( z - laser_sigma_w_z ) * ( z - laser_sigma_w_z );
    rho *= 2.0 * laser_sigma_w0;
    
    phi = atan2( ( y - laser_sigma_w_y ) , ( z - laser_sigma_w_z ) );

    laguerre = - rho + 4.0;
    laguerre *= laguerre;
    return rho * rho * rho * laguerre * cos(3.0*phi) * cos(3.0*phi) * exp(- rho) / ( 24.0 * M_PI );
    
  
    
  }
  
  else
  {
    return 0.0;
  }
}


double laser_intensity_profile_laguerre_20(double x, double y, double z)
{

  double rho;
  double laguerre;
  
  if( x >= laser_offset ){
    
    rho = ( y - laser_sigma_w_y ) * ( y - laser_sigma_w_y ) + ( z - laser_sigma_w_z ) * ( z - laser_sigma_w_z );
    rho *= 2.0 * laser_sigma_w0;
    
    laguerre = 0.5 * ( rho * rho - 4.0 * rho + 2.0 );
    laguerre *= laguerre;

    return laguerre * exp(- rho);
    
  
    
  }
  
  else
  {
    return 0.0;
  }
}

double laser_intensity_profile_laguerre_21(double x, double y, double z)
{

  double phi;
  double rho;
  double laguerre;
  
  if( x >= laser_offset ){
    
    rho = ( y - laser_sigma_w_y ) * ( y - laser_sigma_w_y ) + ( z - laser_sigma_w_z ) * ( z - laser_sigma_w_z );
    rho *= 2.0 * laser_sigma_w0;
    
    phi = atan2( ( y - laser_sigma_w_y ) , ( z - laser_sigma_w_z ) );

    laguerre = 0.5 * rho * rho - 3.0 * rho + 3.0;
    laguerre *= laguerre;
    
    return rho * laguerre * cos(phi) * cos(phi) * exp(- rho) / ( 3.0 * M_PI );
    
  
    
  }
  
  else
  {
    return 0.0;
  }
}

double laser_intensity_profile_laguerre_22(double x, double y, double z)
{

  double phi;
  double rho;
  double laguerre;
  
  if( x >= laser_offset ){
    
    rho = ( y - laser_sigma_w_y ) * ( y - laser_sigma_w_y ) + ( z - laser_sigma_w_z ) * ( z - laser_sigma_w_z );
    rho *= 2.0 * laser_sigma_w0;
    
    phi = atan2( ( y - laser_sigma_w_y ) , ( z - laser_sigma_w_z ) );

    laguerre = 0.5 * rho * rho - 4.0 * rho + 6.0;
    laguerre *= laguerre;
    
    return rho * rho * laguerre * cos(2.0 * phi) * cos(2.0 * phi) * exp(- rho) / ( 12.0  * M_PI );
    
  
    
  }
  
  else
  {
    return 0.0;
  }
}


double laser_intensity_profile_laguerre_23(double x, double y, double z)
{

  double phi;
  double rho;
  double laguerre;
  
  if( x >= laser_offset ){
    
    rho = ( y - laser_sigma_w_y ) * ( y - laser_sigma_w_y ) + ( z - laser_sigma_w_z ) * ( z - laser_sigma_w_z );
    rho *= 2.0 * laser_sigma_w0;
    
    phi = atan2( ( y - laser_sigma_w_y ) , ( z - laser_sigma_w_z ) );

    laguerre = 0.5 * rho * rho - 5.0 * rho + 10.0;
    laguerre *= laguerre;
    
    return rho * rho * rho  * laguerre * cos(3.0 * phi) * cos(3.0 * phi) * exp(- rho) / ( 60.0 * M_PI );
    
  
    
  }
  
  else
  {
    return 0.0;
  }
}


double laser_intensity_profile_laguerre_30(double x, double y, double z)
{


  double rho;
  double laguerre;
  
  if( x >= laser_offset ){
    
    rho = ( y - laser_sigma_w_y ) * ( y - laser_sigma_w_y ) + ( z - laser_sigma_w_z ) * ( z - laser_sigma_w_z );
    rho *= 2.0 * laser_sigma_w0;
    
    laguerre = (- rho * rho * rho + 9.0 * rho * rho - 18.0 * rho + 6.0) / ( 6.0 * M_PI );
    laguerre *= laguerre;

    return laguerre * exp(- rho);
        
  }
  
  else
  {
    return 0.0;
  }
}

double laser_intensity_profile_laguerre_31(double x, double y, double z)
{

  double phi;
  double rho;
  double laguerre;
  
  if( x >= laser_offset ){
    
    rho = ( y - laser_sigma_w_y ) * ( y - laser_sigma_w_y ) + ( z - laser_sigma_w_z ) * ( z - laser_sigma_w_z );
    rho *= 2.0 * laser_sigma_w0;
    
    phi = atan2( ( y - laser_sigma_w_y ) , ( z - laser_sigma_w_z ) );

    laguerre = - rho * rho * rho / 6.0 + 2.0 * rho * rho - 6.0 * rho + 4.0;
    laguerre *= laguerre;

    return rho * laguerre * cos(phi) * cos(phi) * exp(- rho) / ( 4.0 * M_PI );
    
  
    
  }
  
  else
  {
    return 0.0;
  }
}

double laser_intensity_profile_laguerre_32(double x, double y, double z)
{

  double phi;
  double rho;
  double laguerre;
  
  if( x >= laser_offset ){
    
    rho = ( y - laser_sigma_w_y ) * ( y - laser_sigma_w_y ) + ( z - laser_sigma_w_z ) * ( z - laser_sigma_w_z );
    rho *= 2.0 * laser_sigma_w0;
    
    phi = atan2( ( y - laser_sigma_w_y ) , ( z - laser_sigma_w_z ) );
    
    laguerre = - rho * rho * rho / 6.0 + 2.5 * rho  * rho - 10.0 * rho + 10.0;
    laguerre *= laguerre;
    
    return rho * rho * laguerre * cos(2.0 * phi) * cos(2.0 * phi) * exp(- rho) / ( 20.0 * M_PI ); 


  }
  
  else
  {
    return 0.0;
  }
}


double laser_intensity_profile_laguerre_33(double x, double y, double z)
{

  double phi;
  double rho;
  double laguerre;
  
  if( x >= laser_offset ){
    
    rho = ( y - laser_sigma_w_y ) * ( y - laser_sigma_w_y ) + ( z - laser_sigma_w_z ) * ( z - laser_sigma_w_z );
    rho *= 2.0 * laser_sigma_w0;
    
    phi = atan2( ( y - laser_sigma_w_y ) , ( z - laser_sigma_w_z ) );

    laguerre = - rho * rho * rho / 6.0 + 3.0 * rho * rho - 15.0 * rho + 20.0;
    laguerre *= laguerre;

    return rho * rho * rho * laguerre * cos(3.0 * phi) * cos(3.0 * phi) * exp(- rho) / ( 120.0 * M_PI );
    
      
  }
  
  else
  {
    return 0.0;
  }
}



double laser_intensity_profile_hermite_00(double x, double y, double z)
{
  
  double rho;
  
  
  if( x >= laser_offset ){
    
    rho = ( y - laser_sigma_w_y ) * ( y - laser_sigma_w_y ) + ( z - laser_sigma_w_z ) * ( z - laser_sigma_w_z ) ;
    rho *= laser_sigma_w0;
    

    return exp(- rho ) * exp (- rho );
    
  }
  
  else
  {
    return 0.0;
  }
}

double laser_intensity_profile_hermite_01(double x, double y, double z)
{
  
  double rho;
  double zher;
  double hermite;
  
  
  if( x >= laser_offset ){
    

    rho = ( y - laser_sigma_w_y ) * ( y - laser_sigma_w_y ) + ( z - laser_sigma_w_z ) * ( z - laser_sigma_w_z ) ;
    rho *= laser_sigma_w0;
   
    zher =( z - laser_sigma_w_z ) * sqrt(2.0) * sqrt(laser_sigma_w0);
 
    hermite = zher;
    hermite *= hermite;
  
    return hermite * exp( -rho ) * exp (- rho ) ;
    
  }
  
  else
  {
    return 0.0;
  }
}

double laser_intensity_profile_hermite_02(double x, double y, double z)
{
  
  double rho;
  double yher, zher;
  double hermite;
  
  
  if( x >= laser_offset ){
    

    rho = ( y - laser_sigma_w_y ) * ( y - laser_sigma_w_y ) + ( z - laser_sigma_w_z ) * ( z - laser_sigma_w_z ) ;
    rho *= laser_sigma_w0;
     
    yher =( y - laser_sigma_w_y ) * sqrt(2.0) * sqrt(laser_sigma_w0);
    zher =( z - laser_sigma_w_z ) * sqrt(2.0) * sqrt(laser_sigma_w0);

    hermite = zher * zher - 1.0;
    hermite *= hermite;
 
    return hermite * exp(-rho ) * exp (- rho );
    
  }
  
  else
  {
    return 0.0;
  }
}

double laser_intensity_profile_hermite_03(double x, double y, double z)
{
  
  double rho;
  double yher, zher;
  double hermite;
  
  
  if( x >= laser_offset ){
    

    rho = ( y - laser_sigma_w_y ) * ( y - laser_sigma_w_y ) + ( z - laser_sigma_w_z ) * ( z - laser_sigma_w_z ) ;
    rho *= laser_sigma_w0;
    
    yher =( y - laser_sigma_w_y ) * sqrt(2.0) * sqrt(laser_sigma_w0);
    zher =( z - laser_sigma_w_z ) * sqrt(2.0) * sqrt(laser_sigma_w0);
 
    hermite = zher * zher * zher - 3.0 * zher;
    hermite *= hermite;

    return hermite * exp(- rho ) * exp (- rho );
    
  }
  
  else
  {
    return 0.0;
  }
}




double laser_intensity_profile_hermite_10(double x, double y, double z)
{
  double rho;
  double yher;
  double hermite;
  
  
  if( x >= laser_offset ){
    
    rho = ( y - laser_sigma_w_y ) * ( y - laser_sigma_w_y ) + ( z - laser_sigma_w_z ) * ( z - laser_sigma_w_z ) ;
    rho *= laser_sigma_w0;
  
    yher =( y - laser_sigma_w_y ) * sqrt(2.0) * sqrt(laser_sigma_w0);
  
    hermite = yher;
    hermite *= hermite;
  
    return hermite * exp( -rho ) * exp (- rho ) ;
    
  }
  
  else
  {
    return 0.0;
  } 
}

double laser_intensity_profile_hermite_11(double x, double y, double z)
{
  
  double rho;
  double yher, zher;
  double hermite;
  
  
  if( x >= laser_offset ){
    

    rho = ( y - laser_sigma_w_y ) * ( y - laser_sigma_w_y ) + ( z - laser_sigma_w_z ) * ( z - laser_sigma_w_z ) ;
    rho *= laser_sigma_w0;
   

    yher =( y - laser_sigma_w_y ) * sqrt(2.0) * sqrt(laser_sigma_w0);
    zher =( z - laser_sigma_w_z ) * sqrt(2.0) * sqrt(laser_sigma_w0);
 
    hermite = yher * zher;
    hermite *= hermite;
  
    return hermite * exp( -rho ) * exp (- rho ) ;
    
  }
  
  else
  {
    return 0.0;
  }
}

double laser_intensity_profile_hermite_12(double x, double y, double z)
{
  
  double rho;
  double yher, zher;
  double hermite;
  
  
  if( x >= laser_offset ){
    

    rho = ( y - laser_sigma_w_y ) * ( y - laser_sigma_w_y ) + ( z - laser_sigma_w_z ) * ( z - laser_sigma_w_z ) ;
    rho *= laser_sigma_w0;
     
    yher =( y - laser_sigma_w_y ) * sqrt(2.0) * sqrt(laser_sigma_w0);
    zher =( z - laser_sigma_w_z ) * sqrt(2.0) * sqrt(laser_sigma_w0);

    hermite = yher * ( zher * zher - 1.0);
    hermite *= hermite;
 
    return hermite * exp(-rho ) * exp (- rho );
    
  }
  
  else
  {
    return 0.0;
  }
}

double laser_intensity_profile_hermite_13(double x, double y, double z)
{
  
  double rho;
  double yher, zher;
  double hermite;
  
  
  if( x >= laser_offset ){
    

    rho = ( y - laser_sigma_w_y ) * ( y - laser_sigma_w_y ) + ( z - laser_sigma_w_z ) * ( z - laser_sigma_w_z ) ;
    rho *= laser_sigma_w0;
    
    yher =( y - laser_sigma_w_y ) * sqrt(2.0) * sqrt(laser_sigma_w0);
    zher =( z - laser_sigma_w_z ) * sqrt(2.0) * sqrt(laser_sigma_w0);
 
    hermite = yher * ( zher * zher * zher - 3.0 * zher );
    hermite *= hermite;

    return hermite * exp(- rho ) * exp (- rho );
    
  }
  
  else
  {
    return 0.0;
  }
}


double laser_intensity_profile_hermite_20(double x, double y, double z)
{
  double rho;
  double yher;
  double hermite;
  
  
  if( x >= laser_offset ){
    

    rho = ( y - laser_sigma_w_y ) * ( y - laser_sigma_w_y ) + ( z - laser_sigma_w_z ) * ( z - laser_sigma_w_z ) ;
    rho *= laser_sigma_w0;
     
    yher =( y - laser_sigma_w_y ) * sqrt(2.0) * sqrt(laser_sigma_w0);
 
    hermite = yher * yher - 1.0;
    hermite *= hermite;
 
    return hermite * exp(-rho ) * exp (- rho );
    
  }
  
  else
  {
    return 0.0;
  }
}

double laser_intensity_profile_hermite_21(double x, double y, double z)
{
  double rho;
  double yher, zher;
  double hermite;
  
  
  if( x >= laser_offset ){
    

    rho = ( y - laser_sigma_w_y ) * ( y - laser_sigma_w_y ) + ( z - laser_sigma_w_z ) * ( z - laser_sigma_w_z ) ;
    rho *= laser_sigma_w0;
     
    yher =( y - laser_sigma_w_y ) * sqrt(2.0) * sqrt(laser_sigma_w0);
    zher =( z - laser_sigma_w_z ) * sqrt(2.0) * sqrt(laser_sigma_w0);

    hermite = zher * ( yher * yher - 1.0);
    hermite *= hermite;
 
    return hermite * exp(-rho ) * exp (- rho );
    
  }
  
  else
  {
    return 0.0;
  }
}

double laser_intensity_profile_hermite_22(double x, double y, double z)
{
  
  double rho;
  double yher, zher;
  double hermite;
  
  
  if( x >= laser_offset ){
    

    rho = ( y - laser_sigma_w_y ) * ( y - laser_sigma_w_y ) + ( z - laser_sigma_w_z ) * ( z - laser_sigma_w_z ) ;
    rho *= laser_sigma_w0;
     
    yher =( y - laser_sigma_w_y ) * sqrt(2.0) * sqrt(laser_sigma_w0);
    zher =( z - laser_sigma_w_z ) * sqrt(2.0) * sqrt(laser_sigma_w0);

    hermite = ( yher * yher - 1.0 ) * ( zher * zher - 1.0);
    hermite *= hermite;
 
    return hermite * exp(-rho ) * exp (- rho );
    
  }
  
  else
  {
    return 0.0;
  }
}

double laser_intensity_profile_hermite_23(double x, double y, double z)
{
  
  double rho;
  double yher, zher;
  double hermite;
  
  
  if( x >= laser_offset ){
    

    rho = ( y - laser_sigma_w_y ) * ( y - laser_sigma_w_y ) + ( z - laser_sigma_w_z ) * ( z - laser_sigma_w_z ) ;
    rho *= laser_sigma_w0;
    
    yher =( y - laser_sigma_w_y ) * sqrt(2.0) * sqrt(laser_sigma_w0);
    zher =( z - laser_sigma_w_z ) * sqrt(2.0) * sqrt(laser_sigma_w0);
 
    hermite = ( yher * yher - 1.0) * ( zher * zher * zher - 3.0 * zher );
    hermite *= hermite;

    return hermite * exp(- rho ) * exp (- rho );
    
  }
  
  else
  {
    return 0.0;
  }
}




double laser_intensity_profile_hermite_30(double x, double y, double z)
{
  double rho;
  double yher;
  double hermite;
  
  
  if( x >= laser_offset ){
    

    rho = ( y - laser_sigma_w_y ) * ( y - laser_sigma_w_y ) + ( z - laser_sigma_w_z ) * ( z - laser_sigma_w_z ) ;
    rho *= laser_sigma_w0;
    
    yher =( y - laser_sigma_w_y ) * sqrt(2.0) * sqrt(laser_sigma_w0);
 
    hermite = yher * yher * yher - 3.0 * yher;
    hermite *= hermite;

    return hermite * exp(- rho ) * exp (- rho );
    
  }
  
  else
  {
    return 0.0;
  }
}

double laser_intensity_profile_hermite_31(double x, double y, double z)
{
  double rho;
  double yher, zher;
  double hermite;
  
  
  if( x >= laser_offset ){
    

    rho = ( y - laser_sigma_w_y ) * ( y - laser_sigma_w_y ) + ( z - laser_sigma_w_z ) * ( z - laser_sigma_w_z ) ;
    rho *= laser_sigma_w0;
    
    yher =( y - laser_sigma_w_y ) * sqrt(2.0) * sqrt(laser_sigma_w0);
    zher =( z - laser_sigma_w_z ) * sqrt(2.0) * sqrt(laser_sigma_w0);
 
    hermite = zher * ( yher * yher * yher - 3.0 * yher );
    hermite *= hermite;

    return hermite * exp(- rho ) * exp (- rho );
    
  }
  
  else
  {
    return 0.0;
  }

}

double laser_intensity_profile_hermite_32(double x, double y, double z)
{
  double rho;
  double yher, zher;
  double hermite;
  
  
  if( x >= laser_offset ){
    

    rho = ( y - laser_sigma_w_y ) * ( y - laser_sigma_w_y ) + ( z - laser_sigma_w_z ) * ( z - laser_sigma_w_z ) ;
    rho *= laser_sigma_w0;
    
    yher =( y - laser_sigma_w_y ) * sqrt(2.0) * sqrt(laser_sigma_w0);
    zher =( z - laser_sigma_w_z ) * sqrt(2.0) * sqrt(laser_sigma_w0);
 
    hermite = ( zher * zher - 1.0) * ( yher * yher * yher - 3.0 * yher );
    hermite *= hermite;

    return hermite * exp(- rho ) * exp (- rho );
    
  }
  
  else
  {
    return 0.0;
  }
}

double laser_intensity_profile_hermite_33(double x, double y, double z)
{
  
  double rho;
  double yher, zher;
  double hermite;
  
  
  if( x >= laser_offset ){
    

    rho = ( y - laser_sigma_w_y ) * ( y - laser_sigma_w_y ) + ( z - laser_sigma_w_z ) * ( z - laser_sigma_w_z ) ;
    rho *= laser_sigma_w0;
    
    yher =( y - laser_sigma_w_y ) * sqrt(2.0) * sqrt(laser_sigma_w0);
    zher =( z - laser_sigma_w_z ) * sqrt(2.0) * sqrt(laser_sigma_w0);
 
    hermite = ( yher * yher * yher - 3.0 * yher ) * ( zher * zher * zher - 3.0 * zher );
    hermite *= hermite;

    return hermite * exp(- rho ) * exp (- rho );
    
  }
  
  else
  {
    return 0.0;
  }
}



