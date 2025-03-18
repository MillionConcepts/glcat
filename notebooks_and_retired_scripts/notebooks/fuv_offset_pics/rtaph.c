#include "rtaph.h"

static char iam[]     = "RTAph";
static char version[] = "$Id: rtaph.c,v 1.28 2011/03/29 21:51:17 cmillion Exp $";


/* ----------------------------------------------------------------------
  Does a file exist?  If so return "1".
*/
/*@@*/
int rtaph_FileExist( char filename[] )
{
  FILE *fu;
  int ii;
  fu = fopen(filename,"r");
  if (fu == NULL) { ii=0; } else { ii=1; fclose(fu); }
  return(ii);
}


/* ----------------------------------------------------------------------
  Read WIG2, WLK2, and CLK2  FITS calibration files.  -tab 23sep2010 
*/
/*@@*/
void rtaph_read_calibrations( WIG2type WIG2[], WLK2type WLK2[], CLK2type CLK2[], char calpath[], int eclipse )
{
  /**/
  char caldir[200];
  char ffile[200];
  /**/
  fitsfile *fptr;
  /**/
  int status = 0;
  int kk,ii,hdutype,anynull;
  int iya,iq,iyb,ixb;
  int y_as_0, y_as_inc;
  /**/
  long nelem,felem,frow,longnull;
  /**/
  short *ya, *yb, *xb, *yy;
  /**/
  double *ycor;
  /**/
  char comment[200];
  /**/
  const int imax = 30000;
  /**/

  /* Clear arrays. */
  for (iyb=0; iyb<8; ++iyb) {
    for (iya=0; iya<35; ++iya) {
      for (ixb=0; ixb<8;  ++ixb) { for (ii=0; ii<MAXWIG2; ++ii) { WIG2[ii].ycor[iyb][iya][ixb] = 0.; } }
    }
    for (iq=0; iq<32; ++iq) { for (ii=0; ii<MAXWLK2; ++ii) { WLK2[ii].ycor[iyb][iq] = 0.; } }
    for (ii=0; ii<MAXCLK2; ++ii) { CLK2[ii].ycor[iyb] = 0.; }
  }

  /* Check eclipse number. */
  if (eclipse <= 37460) {
    printf("Pre-CSP (5/4/2010) GALEX data [%d]-- not reading new (Sept. 2010) calibrations.\n",eclipse);
    return;
  }
  printf("Post-CSP (5/4/2010) GALEX data [%d]-- reading new (Sept. 2010) calibrations.\n",eclipse);

  /* Allocate. */
  ya  = (short *)calloc(imax,sizeof(short));
  yb  = (short *)calloc(imax,sizeof(short));
  xb  = (short *)calloc(imax,sizeof(short));
  yy  = (short *)calloc(imax,sizeof(short));
  ycor= (double *)calloc(imax,sizeof(double));

  /* Check calpath.  Set string 'caldir'. */
  if (rtaph_FileExist(calpath) == 1) {
    printf("NOTE:rrc: user calpath exists '%s.'\n",calpath);
    strcpy(caldir,calpath);
  } else {
    printf("NOTE:rrc: user calpath does not exist '%s', useing 'cal'.\n",calpath);
    strcpy(caldir,"cal");
  }

  /* Read WIG2 FITS calibration file. */
  sprintf(ffile,"%s/WIG2_Sep2010.fits",caldir);
  fits_open_file( &fptr, ffile, READONLY, &status );                    fits_report_error(stderr,status);if(status!=0)exit(1);
  fits_movabs_hdu( fptr, 1, &hdutype, &status );                        fits_report_error(stderr,status);if(status!=0)exit(1);
  fits_read_key( fptr, TINT, "Y_AS_0"  , &y_as_0  , comment, &status ); fits_report_error(stderr,status);if(status!=0)exit(1);
  fits_read_key( fptr, TINT, "Y_AS_INC", &y_as_inc, comment, &status ); fits_report_error(stderr,status);if(status!=0)exit(1);
  fits_movabs_hdu( fptr, 2, &hdutype, &status );                        fits_report_error(stderr,status);if(status!=0)exit(1);
  fits_read_key( fptr, TINT, "NAXIS2", &nelem, comment, &status );      fits_report_error(stderr,status);if(status!=0)exit(1);
  printf("Reading %ld rows from '%s'.\n",nelem,ffile);
  frow = 1;
  felem= 1;
  fits_read_col( fptr, TSHORT,  1, frow, felem, nelem, &longnull,   ya, &anynull, &status );  fits_report_error(stderr,status);if(status!=0)exit(1);
  fits_read_col( fptr, TSHORT,  2, frow, felem, nelem, &longnull,   yb, &anynull, &status );  fits_report_error(stderr,status);if(status!=0)exit(1);
  fits_read_col( fptr, TSHORT,  3, frow, felem, nelem, &longnull,   xb, &anynull, &status );  fits_report_error(stderr,status);if(status!=0)exit(1);
  fits_read_col( fptr, TSHORT,  4, frow, felem, nelem, &longnull,   yy, &anynull, &status );  fits_report_error(stderr,status);if(status!=0)exit(1);
  fits_read_col( fptr, TDOUBLE, 5, frow, felem, nelem, &longnull, ycor, &anynull, &status );  fits_report_error(stderr,status);if(status!=0)exit(1);
  fits_close_file( fptr, &status ); fits_report_error(stderr,status);if(status!=0)exit(1);

  /* Load WIG2[]. */
  WIG2[0].y_as_start = y_as_0;
  WIG2[0].y_as_inc   = y_as_inc;
  for (kk=0; kk<nelem; ++kk) { WIG2[yy[kk]].ycor[yb[kk]][ya[kk]][xb[kk]] = ycor[kk]; }

  /* Read WLK2 FITS calibration file. */
  sprintf(ffile,"%s/WLK2_Sep2010.fits",caldir);
  fits_open_file( &fptr, ffile, READONLY, &status );                    fits_report_error(stderr,status);if(status!=0)exit(1);
  fits_movabs_hdu( fptr, 1, &hdutype, &status );                        fits_report_error(stderr,status);if(status!=0)exit(1);
  fits_read_key( fptr, TINT, "Y_AS_0"  , &y_as_0  , comment, &status ); fits_report_error(stderr,status);if(status!=0)exit(1);
  fits_read_key( fptr, TINT, "Y_AS_INC", &y_as_inc, comment, &status ); fits_report_error(stderr,status);if(status!=0)exit(1);
  fits_movabs_hdu( fptr, 2, &hdutype, &status );                        fits_report_error(stderr,status);if(status!=0)exit(1);
  fits_read_key( fptr, TINT, "NAXIS2", &nelem, comment, &status );      fits_report_error(stderr,status);if(status!=0)exit(1);
  printf("Reading %ld rows from '%s'.\n",nelem,ffile);
  frow = 1;
  felem= 1;
  /* Here, the 'ya[]' array is used for the Q value. */
  fits_read_col( fptr, TSHORT,  1, frow, felem, nelem, &longnull,   ya, &anynull, &status );  fits_report_error(stderr,status);if(status!=0)exit(1);
  fits_read_col( fptr, TSHORT,  2, frow, felem, nelem, &longnull,   yb, &anynull, &status );  fits_report_error(stderr,status);if(status!=0)exit(1);
  fits_read_col( fptr, TSHORT,  3, frow, felem, nelem, &longnull,   yy, &anynull, &status );  fits_report_error(stderr,status);if(status!=0)exit(1);
  fits_read_col( fptr, TDOUBLE, 4, frow, felem, nelem, &longnull, ycor, &anynull, &status );  fits_report_error(stderr,status);if(status!=0)exit(1);
  fits_close_file( fptr, &status ); fits_report_error(stderr,status);if(status!=0)exit(1);

  /* Load WLK2[]. */
  WLK2[0].y_as_start = y_as_0;
  WLK2[0].y_as_inc   = y_as_inc;
  for (kk=0; kk<nelem; ++kk) { WLK2[yy[kk]].ycor[yb[kk]][ya[kk]] = ycor[kk]; }    /* Here, the 'ya[]' array is used for the Q value. */

  /* Read CLK2 FITS calibration file. */
  sprintf(ffile,"%s/CLK2_Sep2010.fits",caldir);
  fits_open_file( &fptr, ffile, READONLY, &status );                    fits_report_error(stderr,status);if(status!=0)exit(1);
  fits_movabs_hdu( fptr, 1, &hdutype, &status );                        fits_report_error(stderr,status);if(status!=0)exit(1);
  fits_read_key( fptr, TINT, "Y_AS_0"  , &y_as_0  , comment, &status ); fits_report_error(stderr,status);if(status!=0)exit(1);
  fits_read_key( fptr, TINT, "Y_AS_INC", &y_as_inc, comment, &status ); fits_report_error(stderr,status);if(status!=0)exit(1);
  fits_movabs_hdu( fptr, 2, &hdutype, &status );                        fits_report_error(stderr,status);if(status!=0)exit(1);
  fits_read_key( fptr, TINT, "NAXIS2", &nelem, comment, &status );      fits_report_error(stderr,status);if(status!=0)exit(1);
  printf("Reading %ld rows from '%s'.\n",nelem,ffile);
  frow = 1;
  felem= 1;
  fits_read_col( fptr, TSHORT,  1, frow, felem, nelem, &longnull,   yb, &anynull, &status );  fits_report_error(stderr,status);if(status!=0)exit(1);
  fits_read_col( fptr, TSHORT,  2, frow, felem, nelem, &longnull,   yy, &anynull, &status );  fits_report_error(stderr,status);if(status!=0)exit(1);
  fits_read_col( fptr, TDOUBLE, 3, frow, felem, nelem, &longnull, ycor, &anynull, &status );  fits_report_error(stderr,status);if(status!=0)exit(1);
  fits_close_file( fptr, &status ); fits_report_error(stderr,status);if(status!=0)exit(1);

  /* Load CLK2[]. */
  CLK2[0].y_as_start = y_as_0;
  CLK2[0].y_as_inc   = y_as_inc;
  for (kk=0; kk<nelem; ++kk) { CLK2[yy[kk]].ycor[yb[kk]] = ycor[kk]; }

  /* Free. */
  free(ya);   ya  =NULL;
  free(yb);   yb  =NULL;
  free(xb);   xb  =NULL;
  free(yy);   yy  =NULL;
  free(ycor); ycor=NULL;
  return;
}


/* ----------------------------------------------------------------------
  Determine dx total correction from FODC (FUV Offset and Distortion Correction) parameters 
  and data images.   -tab September 2009
*/
double rtaph_FODC_offset_plus_shift( double xp_as, double yp_as, double sct, FODCtype FODC, int band, char xory )
{
  /**/
  const double ArcSecPerPixel = 1.5;   /* 1.5 arcsec per pixel was used during distortion map creation */
  double offset,shift,ss;
  /**/
  int col,row,depth,pixno;
  /**/

  /* Apply FUV offset. */
  if (band == 2) {
    if(FODC.fdttdc <= 20. || FODC.fdttdc >= 40.) {    /* Added this check. -tab 03apr2010 */
      fprintf(stderr,"*** ERROR: fdttdc = %lf is out of range",FODC.fdttdc); exit(1);
    }
    if (xory == 'x') { offset = FODC.fodx_coef_0 - (FODC.fodx_coef_1 * (FODC.fdttdc - 29.)); }
                else { offset = FODC.fody_coef_0 - (FODC.fody_coef_1 * (FODC.fdttdc - 29.)); }
  } else {
    offset = 0.;
  }

  /* Apply distortion correction. */
  ss    = FODC.stim_coef0 + ( sct * FODC.stim_coef1 );    /* Stim separation. */
  col   = ( xp_as - FODC.cube_x0) / FODC.cube_dx;
  row   = ( yp_as - FODC.cube_y0) / FODC.cube_dy;
  depth = ( ss    - FODC.cube_d0) / FODC.cube_dd;

  /* Account for deviant stim separation values. */
  if (depth <        0     ) { depth = 0; } else {
    if (depth >= FODC.cube_nd) depth = FODC.cube_nd - 1;
  }

  /* Shift is zero for points outside column,row bounds. */
  if ((col > -1)&&(col < FODC.cube_nc)&&(row > -1)&&(row < FODC.cube_nr)) {
    pixno = col + ( row * FODC.cube_nc ) + ( depth * FODC.cube_nc * FODC.cube_nr );
    if (xory == 'x') { shift = FODC.dcdx[pixno]; }
                else { shift = FODC.dcdy[pixno]; }
  } else { shift = 0.; }

  /* Convert from pixels to arcseconds. (bug fix)... -tab 16nov2009 */
  shift = shift * ArcSecPerPixel;

  return((offset + shift));
}


/* ----------------------------------------------------------------------
  Wrap the YA value for YA=0 and YA=1 only. -tab 07jul2010
*/
int rtaph_yap( unsigned char ya, unsigned char yb, short yamc )
{
  int yap;
  yap = ya;
  if ((yb > 1)&&(yb < 5)) {
    if ((ya == 0)&&(yamc > -50)) { yap = ya + 32; } else {
    if ((ya == 1)&&(yamc > -10)) { yap = ya + 32; }      }
  }
  return(yap);
}


/* ----------------------------------------------------------------------
  Determine Y correction factor in microns given YA and YB and FODC data.  -tab 25jun2010
    Used like : y = y - rtaph_yac( FODC[0], ya, yb );  
*/
double rtaph_yac( FODCtype FODC, unsigned char ya, unsigned char yb, short yamc )
{
  /**/
  double yac;
  int iyap, iyb;
  /**/
  /* No adjustment for pre-CSP event (May 4, 2010) eclipses. */
  if (FODC.eclipse <= 37460) return(0.);
  iyap= rtaph_yap( ya, yb, yamc );       /* Wrap YA value? */
  iyb = yb;
  yac = FODC.yac[iyap][iyb];             /* Lookup correction. */
  return(yac);
}


/* ----------------------------------------------------------------------
  Determine second Y correction factor in microns given parameters. -tab 24sep2010
    Used like : y = y + rtaph_yac2( );
*/
double rtaph_yac2( WIG2type WIG2[], WLK2type WLK2[], CLK2type CLK2[], 
                   unsigned char q, unsigned char xb, unsigned char yb, unsigned char ya,
                   double y, double aspum )
{
  /**/
  double y_as, yac, yac_as;
  int ii;
  /**/
  yac=0.;
  y_as = y * aspum;
  if ((y_as > -2000.)&&(y_as < 2000.)) {
    yac_as = 0.;
    ii= ((int)y_as - WIG2[0].y_as_start) / WIG2[0].y_as_inc; yac_as=          WIG2[ii].ycor[yb][ya][xb];
    ii= ((int)y_as - WLK2[0].y_as_start) / WLK2[0].y_as_inc; yac_as= yac_as + WLK2[ii].ycor[yb][q];
    ii= ((int)y_as - CLK2[0].y_as_start) / CLK2[0].y_as_inc; yac_as= yac_as + CLK2[ii].ycor[yb];
    yac= yac_as / aspum;
  }
  return(yac);
}




/* ----------------------------------------------------------------------
  Find function and valeus to correct streaking problem... -tab 07jul2010, 14jul2010
  NUV only routine.  Run before rtaph_setup_FODC() ...
  This routine is run only for data after eclipse 37460 about when the 
  Y direction streaking started to occur.
*/
void rtaph_setup_FODC_init( PhRecRaw6Byte * ph6,  int nph6,  int band,  Cal * cal, FODCtype FODC[], int subvis, char calpath[] )
{
  /**/
  int ii, kk, i, stimok;
  int iya,iyb,ixb;
  /**/
  char line[100];
  char wrd[100];
  char svs[10];
  /**/
  short  xamc, yamc;
  short yap;
  /**/
  unsigned char q, xa, xb, yb;
  unsigned char ya;
  /**/
  double avex,avey,num,rmsx,rmsy,avermsx,avermsy,numrmsx,numrmsy;
  double yac_coef0, yac_coef1, yac_ybs[8];
  double yac, stimsep, my[4], my_as[4], mx_as[4], ny[4];
  double ry_as[4], rx_as[4];
  double x, y, x_as, y_as, y_asc, tt;
  double xraw, yraw, xraw0, yraw0;
  double y1, y2, Y1, Y2;
  double x1, x2, X1, X2;
  double rr1, rr0, slope_scale;
  double coef0_yb[8],coef1_yb[8];
  double ybmean[8],ybrms[8],ybnum[8];
  /**/
  char stimfile[80];
  /**/
  FILE *infu;
  FILE *outfu;
  /**/
  int max_STPH;
  int nSTPH[4];
  STPHtype *STPH[4];
  /**/
  /* ORIGINAL NUV average stim positions (pre-eclipse 37423 data)  -tab 14sep2009 */
  /* Data will be scaled so the stims match these positions. */
  const double onuv_avx1 = -2722.27;  const double onuv_avy1 =  2453.89;
  const double onuv_avx2 =  2468.84;  const double onuv_avy2 =  2453.78;
  const double onuv_avx3 = -2721.87;  const double onuv_avy3 = -2565.81;
  const double onuv_avx4 =  2469.85;  const double onuv_avy4 = -2567.83;
  /**/
  double nuv_avx1, nuv_avx2, nuv_avx3, nuv_avx4;
  double nuv_avy1, nuv_avy2, nuv_avy3, nuv_avy4;
  /**/
  int n8;
  double xx, *x8, *y8, *w8;  /*--*/
  double c11, c12, c13, c22, c23;
  /**/

  if (subvis > -1) { sprintf(svs,"_sv%2.2d",subvis); } else { strcpy(svs,""); }

  /* FYI: Using data from 38150 to 38268, median values and (RMS) ...  -tab 14jul2010
     yy1 = 2549.96 (0.82)   xx1 = -2722.53 (0.36)
     yy2 = 2550.10 (0.75)   xx2 =  2470.29 (0.84)
     yy3 =-2538.57 (1.30)   xx3 = -2721.98 (0.31)
     yy4 =-2538.62 (1.36)   xx4 =  2471.09 (0.86) 
  */

  /* NUV only. */
  if (band != 1) return;

  /* REFERENCE NUV average stim positions  -tab 14jul2010 */
  /* These used to find the current stim positions. */
  if (FODC[0].eclipse >= 38150) {   /* Stims moved after the NUV_YCLK change. */
    nuv_avx1 = -2722.53;  nuv_avy1 =  2549.96;
    nuv_avx2 =  2470.29;  nuv_avy2 =  2550.10;
    nuv_avx3 = -2721.98;  nuv_avy3 = -2538.57;
    nuv_avx4 =  2471.09;  nuv_avy4 = -2538.62;
  } else {
    nuv_avx1 = onuv_avx1; nuv_avy1 = onuv_avy1;
    nuv_avx2 = onuv_avx2; nuv_avy2 = onuv_avy2;
    nuv_avx3 = onuv_avx3; nuv_avy3 = onuv_avy3;
    nuv_avx4 = onuv_avx4; nuv_avy4 = onuv_avy4;
  }

  /* Version note. */
  printf("NOTE: Running rtaph_setup_FODC_init() [version 22nov2010] \n");
  printf("NOTE: Nominal NUV average stim positions (arcsec) [%d]:\n",FODC[0].eclipse);
  printf(" x: %10.2f %10.2f %10.2f %10.2f   y: %10.2f %10.2f %10.2f %10.2f \n",
     nuv_avx1, nuv_avx2, nuv_avx3, nuv_avx4, nuv_avy1, nuv_avy2, nuv_avy3, nuv_avy4 );

  /* Set. */
  for (ii=0; ii<4; ++ii) { nSTPH[ii]=0; }

  /* Allocate big STPH array. */
  max_STPH = (nph6 / 25) + 1000;
  if (max_STPH < 900000) max_STPH = 900000;
  for (ii=0; ii<4; ++ii) {
    nSTPH[ii]=0;
    STPH[ii] = (STPHtype *)calloc(max_STPH,sizeof(STPHtype));
  }
  x8 = (double *)calloc(max_STPH,sizeof(double));
  y8 = (double *)calloc(max_STPH,sizeof(double));
  w8 = (double *)calloc(max_STPH,sizeof(double));

  /* Load stim photons. */
  for(i=0; i<nph6; i++, ph6++) {

  /* ========= Unpack raw6 into timing and charge data */
    q    = ((ph6->phb4 & 0x03) << 3) + (( ph6->phb5 & 0xe0) >> 5);
    xb   = ph6->phb1 >> 5;
    xamc = (short) ((ph6->phb1 & 0x1f) << 7) +
           (short) ((ph6->phb2 & 0xfe) >> 1) -
           (short) ((ph6->phb1 & 0x10) << 8);
    yb   = ((ph6->phb2 & 0x01) << 2) + ((ph6->phb3 & 0xc0) >> 6);
    yamc = (short) ((ph6->phb3 & 0x3f) << 6) +
           (short) ((ph6->phb4 & 0xfc) >> 2) -
           (short) ((ph6->phb3 & 0x20) << 7);
    xa   = ((ph6->phb5 & 0x10) >> 4) + ((ph6->phb5 & 0x03)  << 3) +
           ((ph6->phb5 & 0x0c) >> 1);
    tt   = ph6->t;
  /* ==========  */

    xraw0 = xb* FODC[0].NUV_XCLK  + xamc ;
    yraw0 = yb* FODC[0].NUV_YCLK  + yamc ;
    ya = ((int)(((yraw0/(2* FODC[0].NUV_YCLK ) - xraw0/(2* FODC[0].NUV_XCLK )   )+10)*32)+xa) % 32;
    xraw = xraw0+((int)((xa+7)%32)-16)*FODC[0].NUV_XSLP;
    yraw = yraw0+((int)((ya+7)%32)-16)*FODC[0].NUV_YSLP;
    x = (xraw - FODC[0].NUV_XCEN)*(FODC[0].NUV_XSCL);   /* +or- 32725 to detsize */
    y = (yraw - FODC[0].NUV_YCEN)*(FODC[0].NUV_YSCL);

  /* Convert detector position from microns to arcsecs  [NUV] */
    x_as = x*cal->aspum;
    y_as = y*cal->aspum;

  /* Stims. */
    kk=-1;
    if ((x_as > (nuv_avx1 - 90.001))&&(x_as < (nuv_avx1 + 90.001))&&(y_as > (nuv_avy1 - 90.001))&&(y_as < (nuv_avy1 + 90.001))) { kk=0; }
    if ((x_as > (nuv_avx2 - 90.001))&&(x_as < (nuv_avx2 + 90.001))&&(y_as > (nuv_avy2 - 90.001))&&(y_as < (nuv_avy2 + 90.001))) { kk=1; }
    if ((x_as > (nuv_avx3 - 90.001))&&(x_as < (nuv_avx3 + 90.001))&&(y_as > (nuv_avy3 - 90.001))&&(y_as < (nuv_avy3 + 90.001))) { kk=2; }
    if ((x_as > (nuv_avx4 - 90.001))&&(x_as < (nuv_avx4 + 90.001))&&(y_as > (nuv_avy4 - 90.001))&&(y_as < (nuv_avy4 + 90.001))) { kk=3; }
    if (kk > -1) {
      STPH[kk][nSTPH[kk]].t    = tt;
      STPH[kk][nSTPH[kk]].q    = q;
      STPH[kk][nSTPH[kk]].xb   = xb;
      STPH[kk][nSTPH[kk]].xamc = xamc;
      STPH[kk][nSTPH[kk]].yb   = yb;
      STPH[kk][nSTPH[kk]].yamc = yamc;
      STPH[kk][nSTPH[kk]].xa   = xa;
      STPH[kk][nSTPH[kk]].ya   = ya;
      ++nSTPH[kk];
    } 
  }

  /* Compute initial mean positions, load STPH[][].y values and my[] values. */
  /* Now using all stim photons, rather than restricted YA range. -tab 27sep2010 */
  stimok=1;
  for (kk=0; kk<4; ++kk) {
    my[kk]=0.; my_as[kk]=0.; mx_as[kk]=0.; ny[kk]=0.;
    for (ii=0; ii<nSTPH[kk]; ++ii) {
      xraw0 = STPH[kk][ii].xb* FODC[0].NUV_XCLK + STPH[kk][ii].xamc ;
      yraw0 = STPH[kk][ii].yb* FODC[0].NUV_YCLK + STPH[kk][ii].yamc ;
      ya = ((int)(((yraw0/(2*FODC[0].NUV_YCLK) - xraw0/(2*FODC[0].NUV_XCLK)   )+10)*32)+STPH[kk][ii].xa) % 32;
      xraw = xraw0+((int)((STPH[kk][ii].xa+7)%32)-16)*FODC[0].NUV_XSLP;
      yraw = yraw0+((int)((STPH[kk][ii].ya+7)%32)-16)*FODC[0].NUV_YSLP;
      x = (xraw - FODC[0].NUV_XCEN)*(FODC[0].NUV_XSCL);   /* +or- 32725 to detsize */
      y = (yraw - FODC[0].NUV_YCEN)*(FODC[0].NUV_YSCL);
      STPH[kk][ii].y = y;
      my[kk]    = my[kk]    + y;
      my_as[kk] = my_as[kk] + (y*cal->aspum);
      mx_as[kk] = mx_as[kk] + (x*cal->aspum);
      ny[kk] = ny[kk] + 1.;
    }
    if (ny[kk] < 20.) { 
      fprintf(stderr,"===WARNING: not enough values for stim number %d [%f].\n",kk,ny[kk]); 
      stimok=0;
    } else {
      my[kk]    = my[kk]    / ny[kk];
      my_as[kk] = my_as[kk] / ny[kk];
      mx_as[kk] = mx_as[kk] / ny[kk];
    }
  }

  /* Save results (or read back in). */
  sprintf(wrd,"stim_positions_%d.dat",FODC[0].eclipse);
  if (stimok == 1) {
    printf("NOTE: Writing stim positions file '%s'.\n",wrd);
    outfu = fopen(wrd,"w");
    fprintf(outfu,"%20.12f %20.12f\n",mx_as[0],my_as[0]);
    fprintf(outfu,"%20.12f %20.12f\n",mx_as[1],my_as[1]);
    fprintf(outfu,"%20.12f %20.12f\n",mx_as[2],my_as[2]);
    fprintf(outfu,"%20.12f %20.12f\n",mx_as[3],my_as[3]);
    fclose(outfu);
  } else {
    if (rtaph_FileExist(wrd) == 1) {
      printf("NOTE: Reading stim positions file '%s'.\n",wrd);
      infu  = fopen(wrd,"r");
      fscanf(infu,"%lf %lf",&mx_as[0],&my_as[0]);
      fscanf(infu,"%lf %lf",&mx_as[1],&my_as[1]);
      fscanf(infu,"%lf %lf",&mx_as[2],&my_as[2]);
      fscanf(infu,"%lf %lf",&mx_as[3],&my_as[3]);
      fclose(infu);
      my[0] = my_as[0] / cal->aspum;
      my[1] = my_as[1] / cal->aspum;
      my[2] = my_as[2] / cal->aspum;
      my[3] = my_as[3] / cal->aspum;
    } else {
      fprintf(stderr,"***ERROR: Cannot find '%s' in order to read previous STIM positions.\n",wrd);
      exit(1);
    }
  }

  /* Compute stim separation and report. */
  stimsep = ( (mx_as[1] - mx_as[0])  + (mx_as[3] - mx_as[2]) + (my_as[0] - my_as[2]) + (my_as[1] - my_as[3]) ) / 4.;
  printf("\n");
  printf("Init: Number of stim photons: %0.1f %0.1f %0.1f %0.1f .\n",ny[0],ny[1],ny[2],ny[3]);
  printf("Init: Mean x values at stim positions (arcsec): %14.6f %14.6f %14.6f %14.6f \n", mx_as[0], mx_as[1], mx_as[2], mx_as[3]);
  printf("Init: Mean y values at stim positions (arcsec): %14.6f %14.6f %14.6f %14.6f \n", my_as[0], my_as[1], my_as[2], my_as[3]);
  printf("Init: Mean y values at stim positions(microns): %14.6f %14.6f %14.6f %14.6f [aspum=%14.6f]\n",my[0],my[1],my[2],my[3],cal->aspum);

  /* Compute RMS from mean positions. */
  for (kk=0; kk<4; ++kk) {
    ry_as[kk]=0.; rx_as[kk]=0.; ny[kk]=0.;
    for (ii=0; ii<nSTPH[kk]; ++ii) {
      xraw0 = STPH[kk][ii].xb* FODC[0].NUV_XCLK + STPH[kk][ii].xamc ;
      yraw0 = STPH[kk][ii].yb* FODC[0].NUV_YCLK + STPH[kk][ii].yamc ;
      ya = ((int)(((yraw0/(2*FODC[0].NUV_YCLK) - xraw0/(2*FODC[0].NUV_XCLK)   )+10)*32)+STPH[kk][ii].xa) % 32;
      xraw = xraw0+((int)((STPH[kk][ii].xa+7)%32)-16)*FODC[0].NUV_XSLP;
      yraw = yraw0+((int)((STPH[kk][ii].ya+7)%32)-16)*FODC[0].NUV_YSLP;
      x = (xraw - FODC[0].NUV_XCEN)*(FODC[0].NUV_XSCL);   /* +or- 32725 to detsize */
      y = (yraw - FODC[0].NUV_YCEN)*(FODC[0].NUV_YSCL);
      ry_as[kk] = ry_as[kk] + ( ((y*cal->aspum) - my_as[kk]) * ((y*cal->aspum) - my_as[kk]) );
      rx_as[kk] = rx_as[kk] + ( ((x*cal->aspum) - mx_as[kk]) * ((x*cal->aspum) - mx_as[kk]) );
      ny[kk] = ny[kk] + 1.;
    }
    if (ny[kk] > 0.) {
      ry_as[kk] = sqrt(( ry_as[kk] / ny[kk] ));
      rx_as[kk] = sqrt(( rx_as[kk] / ny[kk] ));
    } else {
      ry_as[kk] = 0.;
      rx_as[kk] = 0.;
    }
  }
  printf("Init: RMS  x values at stim positions (arcsec): %14.6f %14.6f %14.6f %14.6f \n", rx_as[0], rx_as[1], rx_as[2], rx_as[3]);
  printf("Init: RMS  y values at stim positions (arcsec): %14.6f %14.6f %14.6f %14.6f \n", ry_as[0], ry_as[1], ry_as[2], ry_as[3]);
  printf("Init: (arcsec): Stim Sep = %10.3f     Average: X RMS = %8.3f    Y RMS = %8.3f \n",
         stimsep, (rx_as[0] +  rx_as[1] +  rx_as[2] +  rx_as[3]) / 4., (ry_as[0] +  ry_as[1] +  ry_as[2] +  ry_as[3]) / 4. );

  /* Save raw stim separation value for use in distortion maps. -tab 18oct2010 */
  FODC[0].raw_stimsep = stimsep;
  printf("Raw stim separation is %10.3f\n",FODC[0].raw_stimsep);

  /* Report Mean positions and RMS values for each YA value. */
  printf("\n");
  printf("|st|ya|yb|xb| avex    | avey    | rmsx    | rmsy    | num   |\n");
  avermsx=0.; avermsy=0.;
  numrmsx=0.; numrmsy=0.;
  for (iya=0; iya<32; ++iya) {
  for (iyb=0; iyb<8;  ++iyb) {
  for (ixb=0; ixb<8;  ++ixb) {
    for (kk=0; kk<4; ++kk) {
      avex=0.; avey=0.; num=0.;
      for (ii=0; ii<nSTPH[kk]; ++ii) {
        xraw0 = STPH[kk][ii].xb* FODC[0].NUV_XCLK + STPH[kk][ii].xamc ;
        yraw0 = STPH[kk][ii].yb* FODC[0].NUV_YCLK + STPH[kk][ii].yamc ;
        ya = ((int)(((yraw0/(2*FODC[0].NUV_YCLK) - xraw0/(2*FODC[0].NUV_XCLK)   )+10)*32)+STPH[kk][ii].xa) % 32;
        if ((ya == iya)&&(STPH[kk][ii].yb == iyb)&&(STPH[kk][ii].xb == ixb)) {
          xraw = xraw0+((int)((STPH[kk][ii].xa+7)%32)-16)*FODC[0].NUV_XSLP;
          yraw = yraw0+((int)((STPH[kk][ii].ya+7)%32)-16)*FODC[0].NUV_YSLP;
          x = (xraw - FODC[0].NUV_XCEN)*(FODC[0].NUV_XSCL);   /* +or- 32725 to detsize */
          y = (yraw - FODC[0].NUV_YCEN)*(FODC[0].NUV_YSCL);
          avex = avex + x;
          avey = avey + y;
          num  = num  + 1.;
        }
      }
      if (num > 100.) {
        avex = avex / num;
        avey = avey / num;
        rmsx=0.; rmsy=0.;
        for (ii=0; ii<nSTPH[kk]; ++ii) {
          xraw0 = STPH[kk][ii].xb* FODC[0].NUV_XCLK + STPH[kk][ii].xamc ;
          yraw0 = STPH[kk][ii].yb* FODC[0].NUV_YCLK + STPH[kk][ii].yamc ;
          ya = ((int)(((yraw0/(2*FODC[0].NUV_YCLK) - xraw0/(2*FODC[0].NUV_XCLK)   )+10)*32)+STPH[kk][ii].xa) % 32;
          if ((ya == iya)&&(STPH[kk][ii].yb == iyb)&&(STPH[kk][ii].xb == ixb)) {
            xraw = xraw0+((int)((STPH[kk][ii].xa+7)%32)-16)*FODC[0].NUV_XSLP;
            yraw = yraw0+((int)((STPH[kk][ii].ya+7)%32)-16)*FODC[0].NUV_YSLP;
            x = (xraw - FODC[0].NUV_XCEN)*(FODC[0].NUV_XSCL);   /* +or- 32725 to detsize */
            y = (yraw - FODC[0].NUV_YCEN)*(FODC[0].NUV_YSCL);
            rmsx = rmsx + ( (x - avex) * (x - avex) );
            rmsy = rmsy + ( (y - avey) * (y - avey) );
          }
        }
        rmsx = sqrt(( rmsx / num ));   avermsx = avermsx + rmsx;  numrmsx = numrmsx + 1.;
        rmsy = sqrt(( rmsy / num ));   avermsy = avermsy + rmsy;  numrmsy = numrmsy + 1.;
        printf(" %2d %2d %2d %2d %9.3f %9.3f %9.4f %9.4f %7.1f \n",
                   kk,iya,iyb,ixb, avex * cal->aspum, avey * cal->aspum, rmsx * cal->aspum, rmsy * cal->aspum, num );
      }
    }
  }
  }
  }
  if (numrmsx > 0.) avermsx = avermsx / numrmsx;
  if (numrmsy > 0.) avermsy = avermsy / numrmsy;
  printf(" Average RMS X = %9.3f    Average RMS Y = %9.3f \n", avermsx * cal->aspum, avermsy * cal->aspum );
  printf("\n");


  /* Write RAW stim tables. */
  kk=1; sprintf(stimfile,"rawstims%s-%2.2d.tbl",svs,kk);
  while (rtaph_FileExist(stimfile) == 1) { ++kk; sprintf(stimfile,"rawstims%s-%2.2d.tbl",svs,kk); }
  printf("NOTE: Writing raw stims table: %s .\n",stimfile);
  outfu = fopen(stimfile,"w");
  fprintf(outfu,"| yb | yap | rely       |st| x_as    | y_as    |  ya  | xa   | xamc | yamc | q    | yac      | xb   |\n");
  /*              1234 12345 123456789012 12 123456789 123456789 123456 123456 123456 123456 123456 1234567890 123456 */
  for (kk=0; kk<4; ++kk) {
    for (ii=0; ii<nSTPH[kk]; ++ii) {
      xraw0 = STPH[kk][ii].xb* FODC[0].NUV_XCLK + STPH[kk][ii].xamc ;
      yraw0 = STPH[kk][ii].yb* FODC[0].NUV_YCLK + STPH[kk][ii].yamc ;
      ya = ((int)(((yraw0/(2*FODC[0].NUV_YCLK) - xraw0/(2*FODC[0].NUV_XCLK)   )+10)*32)+STPH[kk][ii].xa) % 32;
      xraw = xraw0+((int)((STPH[kk][ii].xa+7)%32)-16)*FODC[0].NUV_XSLP;
      yraw = yraw0+((int)((STPH[kk][ii].ya+7)%32)-16)*FODC[0].NUV_YSLP;
      x = (xraw - FODC[0].NUV_XCEN)*(FODC[0].NUV_XSCL);   /* +or- 32725 to detsize */
      y = (yraw - FODC[0].NUV_YCEN)*(FODC[0].NUV_YSCL);
      x_as = x*cal->aspum;
      y_as = y*cal->aspum;
      yap  = rtaph_yap( ya, STPH[kk][ii].yb, STPH[kk][ii].yamc );     /* possible wrap around... -tab 07jul2010 */
      fprintf(outfu," %4d %5d %12.3f %2d %9.3f %9.3f %6d %6d %6d %6d %6d %10.2f %6d \n",
           STPH[kk][ii].yb, yap, (y - my[kk]),
           kk, x_as, y_as, ya, STPH[kk][ii].xa, STPH[kk][ii].xamc, STPH[kk][ii].yamc,
           STPH[kk][ii].q, 0., STPH[kk][ii].xb );
    }
  }
  fclose(outfu);

  /* Compute Y scale and shift factors:  yprime_as = (m * y_as) + B    -tab 01jul2010 */
  y1 = (my_as[0] + my_as[1] ) / 2.;   y2 = (my_as[2] + my_as[3]) / 2;
  Y1 = (onuv_avy1+ onuv_avy2) / 2.;   Y2 = (onuv_avy3+ onuv_avy4) / 2.;
  FODC[0].My = (Y1 - Y2) / (y1 - y2);
  FODC[0].By = ( Y1 - (FODC[0].My * y1) ) / cal->aspum;
  printf("Init: FODC: Y scale and shift (microns): My=%18.12f   By=%18.12f \n",FODC[0].My,FODC[0].By);

  /* Compute X scale and shift factors:  xprime_as = (m * x_as) + B    -tab 01jul2010 */
  x1 = (mx_as[0] + mx_as[2] ) / 2.;   x2 = (mx_as[1] + mx_as[3]) / 2;
  X1 = (onuv_avx1+ onuv_avx3) / 2.;   X2 = (onuv_avx2+ onuv_avx4) / 2.;
  FODC[0].Mx = (X1 - X2) / (x1 - x2);
  FODC[0].Bx = ( X1 - (FODC[0].Mx * x1) ) / cal->aspum;
  printf("Init: FODC: X scale and shift (microns): Mx=%18.12f   Bx=%18.12f \n",FODC[0].Mx,FODC[0].Bx);

  /* SCALED: Compute mean positions, load STPH[][].y values and my[] values. */
  if (stimok == 1) {
    for (kk=0; kk<4; ++kk) {
      my[kk]=0.; my_as[kk]=0.; mx_as[kk]=0.; ny[kk]=0.;
      for (ii=0; ii<nSTPH[kk]; ++ii) {
        xraw0 = STPH[kk][ii].xb* FODC[0].NUV_XCLK + STPH[kk][ii].xamc ;
        yraw0 = STPH[kk][ii].yb* FODC[0].NUV_YCLK + STPH[kk][ii].yamc ;
        ya = ((int)(((yraw0/(2*FODC[0].NUV_YCLK) - xraw0/(2*FODC[0].NUV_XCLK)   )+10)*32)+STPH[kk][ii].xa) % 32;
        xraw = xraw0+((int)((STPH[kk][ii].xa+7)%32)-16)*FODC[0].NUV_XSLP;
        yraw = yraw0+((int)((STPH[kk][ii].ya+7)%32)-16)*FODC[0].NUV_YSLP;
        x = (xraw - FODC[0].NUV_XCEN)*(FODC[0].NUV_XSCL);   /* +or- 32725 to detsize */
        y = (yraw - FODC[0].NUV_YCEN)*(FODC[0].NUV_YSCL);
    /* Rescale x,y value. */
        x = (FODC[0].Mx * x) + FODC[0].Bx;
        y = (FODC[0].My * y) + FODC[0].By;
        STPH[kk][ii].y = y;
        my[kk]    = my[kk]    + y;
        my_as[kk] = my_as[kk] + (y*cal->aspum);
        mx_as[kk] = mx_as[kk] + (x*cal->aspum);
        ny[kk] = ny[kk] + 1.;
      }
      if (ny[kk] < 2) { fprintf(stderr,"***ERROR: no values for stim number %d .\n",kk); exit(1); }
      my[kk]    = my[kk]    / ny[kk];
      my_as[kk] = my_as[kk] / ny[kk];
      mx_as[kk] = mx_as[kk] / ny[kk];
    }
    stimsep = ( (mx_as[1] - mx_as[0])  + (mx_as[3] - mx_as[2]) + (my_as[0] - my_as[2]) + (my_as[1] - my_as[3]) ) / 4.;
    printf("\n");
    printf("Scal: Number of stim photons: %0.1f %0.1f %0.1f %0.1f .\n",ny[0],ny[1],ny[2],ny[3]);
    printf("Scal: Mean x values at stim positions (arcsec): %14.6f %14.6f %14.6f %14.6f \n", mx_as[0], mx_as[1], mx_as[2], mx_as[3]);
    printf("Scal: Mean y values at stim positions (arcsec): %14.6f %14.6f %14.6f %14.6f \n", my_as[0], my_as[1], my_as[2], my_as[3]);
    printf("Scal: Mean y values at stim positions(microns): %14.6f %14.6f %14.6f %14.6f [aspum=%14.6f]\n",my[0],my[1],my[2],my[3],cal->aspum);

    /* SCALED: Compute RMS from mean positions. */
    for (kk=0; kk<4; ++kk) {
      ry_as[kk]=0.; rx_as[kk]=0.; ny[kk]=0.;
      for (ii=0; ii<nSTPH[kk]; ++ii) {
        xraw0 = STPH[kk][ii].xb* FODC[0].NUV_XCLK + STPH[kk][ii].xamc ;
        yraw0 = STPH[kk][ii].yb* FODC[0].NUV_YCLK + STPH[kk][ii].yamc ;
        ya = ((int)(((yraw0/(2*FODC[0].NUV_YCLK) - xraw0/(2*FODC[0].NUV_XCLK)   )+10)*32)+STPH[kk][ii].xa) % 32;
        xraw = xraw0+((int)((STPH[kk][ii].xa+7)%32)-16)*FODC[0].NUV_XSLP;
        yraw = yraw0+((int)((STPH[kk][ii].ya+7)%32)-16)*FODC[0].NUV_YSLP;
        x = (xraw - FODC[0].NUV_XCEN)*(FODC[0].NUV_XSCL);   /* +or- 32725 to detsize */
        y = (yraw - FODC[0].NUV_YCEN)*(FODC[0].NUV_YSCL);
    /* Rescale x,y value. */
        x = (FODC[0].Mx * x) + FODC[0].Bx;
        y = (FODC[0].My * y) + FODC[0].By;
        ry_as[kk] = ry_as[kk] + ( ((y*cal->aspum) - my_as[kk]) * ((y*cal->aspum) - my_as[kk]) );
        rx_as[kk] = rx_as[kk] + ( ((x*cal->aspum) - mx_as[kk]) * ((x*cal->aspum) - mx_as[kk]) );
        ny[kk] = ny[kk] + 1.;
      }
      ry_as[kk] = sqrt(( ry_as[kk] / ny[kk] ));
      rx_as[kk] = sqrt(( rx_as[kk] / ny[kk] ));
    }
    printf("Scal: RMS  x values at stim positions (arcsec): %14.6f %14.6f %14.6f %14.6f \n", rx_as[0], rx_as[1], rx_as[2], rx_as[3]);
    printf("Scal: RMS  y values at stim positions (arcsec): %14.6f %14.6f %14.6f %14.6f \n", ry_as[0], ry_as[1], ry_as[2], ry_as[3]);
    printf("Scal: (arcsec): Stim Sep = %10.3f     Average: X RMS = %8.3f    Y RMS = %8.3f \n",
           stimsep, (rx_as[0] +  rx_as[1] +  rx_as[2] +  rx_as[3]) / 4., (ry_as[0] +  ry_as[1] +  ry_as[2] +  ry_as[3]) / 4. );
    printf("\n");

    /* SCALED: Load fit arrays. */
    n8=0; /*--*/
    for (kk=0; kk<4; ++kk) {
      for (ii=0; ii<nSTPH[kk]; ++ii) {
        xraw0 = STPH[kk][ii].xb* FODC[0].NUV_XCLK + STPH[kk][ii].xamc ;
        yraw0 = STPH[kk][ii].yb* FODC[0].NUV_YCLK + STPH[kk][ii].yamc ;
        ya = ((int)(((yraw0/(2*FODC[0].NUV_YCLK) - xraw0/(2*FODC[0].NUV_XCLK)   )+10)*32)+STPH[kk][ii].xa) % 32;
        xraw = xraw0+((int)((STPH[kk][ii].xa+7)%32)-16)*FODC[0].NUV_XSLP;
        yraw = yraw0+((int)((STPH[kk][ii].ya+7)%32)-16)*FODC[0].NUV_YSLP;
        x = (xraw - FODC[0].NUV_XCEN)*(FODC[0].NUV_XSCL);   /* +or- 32725 to detsize */
        y = (yraw - FODC[0].NUV_YCEN)*(FODC[0].NUV_YSCL);
        x = (FODC[0].Mx * x) + FODC[0].Bx;
        y = (FODC[0].My * y) + FODC[0].By;
        x_as = x*cal->aspum;
        y_as = y*cal->aspum;
        yap  = rtaph_yap( ya, STPH[kk][ii].yb, STPH[kk][ii].yamc );     /* possible wrap around... -tab 07jul2010 */
    /* Only fit or scale to YA > 2 points and YB=2. */
        if ((ya > 2)&&(STPH[kk][ii].yb == 2)) {   /*--*/
          x8[n8]=yap;  y8[n8]=(y - my[kk]);  w8[n8]=1.;  ++n8;  /*--*/
        }
      }
    }
  } else {
   n8 = 0;
  }

  /* Fit straight line to YB=2 points only.... -tab 30jun2010 */
  if (n8 > 200) {
    printf("NOTE: Found %d points for YA correction fit.\n",n8);
    c11=0.; c12=0.; c13=0.; c22=0.; c23=0.;
    for (ii=0; ii<n8; ++ii) {
      c11 = c11 + w8[ii];                  /* weight sum */
      c12 = c12 + x8[ii] * w8[ii];         /* x times weight */
      c13 = c13 + y8[ii] * w8[ii];         /* y times weight */
      c22 = c22 + x8[ii] * x8[ii] * w8[ii];
      c23 = c23 + y8[ii] * w8[ii] * x8[ii];
    }
    if (c11 < 3.) { fprintf(stderr,"***ERROR: not enough fit points[%d] for YB=2.\n",(int)c11); exit(1); }
    yac_coef1 = ( (c13 * c12) - (c23 * c11) ) / ( (c12 * c12) - (c22 * c11) );
    yac_coef0 = ( c13 - (c12 * yac_coef1) ) / c11;
    sprintf(wrd,"yac_coef_%d.dat",FODC[0].eclipse);
    printf("NOTE: Writing YAC coeffecients to '%s'.\n",wrd);
    outfu = fopen(wrd,"w");
    fprintf(outfu,"%20.12f %20.12f\n",yac_coef0,yac_coef1);
    fclose(outfu);
  } else {
    printf("NOTE: YAC stim fit: Not enough points for fit, n8=%d.\n",n8);
    sprintf(wrd,"yac_coef_%d.dat",FODC[0].eclipse);
    if (rtaph_FileExist(wrd) == 1) {
      printf("NOTE: Reading previous YAC coeffecients from '%s'.\n",wrd);
      infu = fopen(wrd,"r");
      fscanf(infu,"%lf %lf",&yac_coef0,&yac_coef1);
      fclose(infu);
    } else {
      fprintf(stderr,"***ERROR: Cannot find '%s' in order to read previous YAC coeffecients.\n",wrd);
      exit(1);
    }
  }
  printf("Scal: YA correction coef for YB=2: %20.12e %20.12e \n",yac_coef0,yac_coef1);


  /* Write fit points, fit, and other values. */
  if (stimok == 1) {
    kk=1; sprintf(stimfile,"yacfit%s-%2.2d.tbl",svs,kk);
    while (rtaph_FileExist(stimfile) == 1) { ++kk; sprintf(stimfile,"yacfit%s-%2.2d.tbl",svs,kk); }
    printf("NOTE: Writing YAC fit values table: %s .\n",stimfile);
    outfu = fopen(stimfile,"w");
    fprintf(outfu,"| yb | yap | rely     |st| x_as    | y_as    |  ya  | xa   | xamc | yamc | q    | yac      | xb   | rely_fit |\n");
    /*              1234 12345 1234567890 12 123456789 123456789 123456 123456 123456 123456 123456 1234567890 123456 1234567890 */
    for (kk=0; kk<4; ++kk) {
      for (ii=0; ii<nSTPH[kk]; ++ii) {
        xraw0 = STPH[kk][ii].xb* FODC[0].NUV_XCLK + STPH[kk][ii].xamc ;
        yraw0 = STPH[kk][ii].yb* FODC[0].NUV_YCLK + STPH[kk][ii].yamc ;
        ya = ((int)(((yraw0/(2*FODC[0].NUV_YCLK) - xraw0/(2*FODC[0].NUV_XCLK)   )+10)*32)+STPH[kk][ii].xa) % 32;
        xraw = xraw0+((int)((STPH[kk][ii].xa+7)%32)-16)*FODC[0].NUV_XSLP;
        yraw = yraw0+((int)((STPH[kk][ii].ya+7)%32)-16)*FODC[0].NUV_YSLP;
        x = (xraw - FODC[0].NUV_XCEN)*(FODC[0].NUV_XSCL);   /* +or- 32725 to detsize */
        y = (yraw - FODC[0].NUV_YCEN)*(FODC[0].NUV_YSCL);
        x = (FODC[0].Mx * x) + FODC[0].Bx;
        y = (FODC[0].My * y) + FODC[0].By;
        x_as = x*cal->aspum;
        y_as = y*cal->aspum;
        yap  = rtaph_yap( ya, STPH[kk][ii].yb, STPH[kk][ii].yamc );     /* possible wrap around... -tab 07jul2010 */
        if ((ya > 2)&&(STPH[kk][ii].yb == 2)) {   /*--*/
          rr1 = yac_coef0 + ((double)yap * yac_coef1);
          fprintf(outfu," %4d %5d %10.3f %2d %9.3f %9.3f %6d %6d %6d %6d %6d %10.2f %6d %10.3f \n",
             STPH[kk][ii].yb, yap, (y - my[kk]),
             kk, x_as, y_as, ya, STPH[kk][ii].xa, STPH[kk][ii].xamc, STPH[kk][ii].yamc,
             STPH[kk][ii].q, 0., STPH[kk][ii].xb, rr1 );
        }
      }
    }
    fclose(outfu);
  }

  /* Compute yb shift factors.  -tab 30jun2010 */
  /* Use zero for all. -tab 27sep2010 */
  for (iyb=0; iyb<8; ++iyb) { yac_ybs[iyb]=0.; }

  /* Default. */
  for (iyb=0; iyb<8; ++iyb) {
    coef0_yb[iyb] = yac_coef0;
    coef1_yb[iyb] = yac_coef1;
  }

  /* Set user slope adjustment.  Use best slope adjustments from September 2010. */
  /* YB=2... */
  slope_scale = 1.04;   printf("NOTE: Using slope scale of %6.3f for YB=2.\n",slope_scale);
  rr1 = yac_coef1 * slope_scale;
  rr0 = (yac_coef0 + (16. * yac_coef1)) - (16. * rr1);
  coef0_yb[2] = rr0;
  coef1_yb[2] = rr1;
  printf(" New: YA correction coef (YB=2): %20.12e %20.12e \n",coef0_yb[2], coef1_yb[2]);

  /* YB=3... */
  slope_scale = 1.06;   printf("NOTE: Using slope scale of %6.3f for YB=3.\n",slope_scale);
  rr1 = yac_coef1 * slope_scale;
  rr0 = (yac_coef0 + (16. * yac_coef1)) - (16. * rr1);
  coef0_yb[3] = rr0;
  coef1_yb[3] = rr1;
  printf(" New: YA correction coef (YB=3): %20.12e %20.12e \n",coef0_yb[3], coef1_yb[3]);

  /* YB=4... */
  slope_scale = 1.06;   printf("NOTE: Using slope scale of %6.3f for YB=4.\n",slope_scale);
  rr1 = yac_coef1 * slope_scale;
  rr0 = (yac_coef0 + (16. * yac_coef1)) - (16. * rr1);
  coef0_yb[4] = rr0;
  coef1_yb[4] = rr1;
  printf(" New: YA correction coef (YB=4): %20.12e %20.12e \n",coef0_yb[4], coef1_yb[4]);

  /* Fill in look up array. */
  for (iyb=0; iyb<8; ++iyb) {
    for (iya=0; iya<40; ++iya) {
      xx  = (double)iya;
      FODC[0].yac[iya][iyb] = (coef0_yb[iyb] + (xx * coef1_yb[iyb])) + yac_ybs[iyb];
    }
  }

  /* Replace look up array with measured values. */
  /* Not doing this anymore. -tab 27sep2010 */

  /* Print out correction factors. */
  for (iya=0; iya<40; ++iya) {
    strcpy(line,"");
    for (iyb=0; iyb<8; ++iyb) {
      sprintf(wrd," %14.8f",FODC[0].yac[iya][iyb]);
      strcat(line,wrd);
    }
    printf("%s\n",line);
  }

  /* FINAL: Write stim tables.  Include yac values. */
  if (stimok == 1) {
    kk=1; sprintf(stimfile,"finstims%s-%2.2d.tbl",svs,kk);
    while (rtaph_FileExist(stimfile) == 1) { ++kk; sprintf(stimfile,"finstims%s-%2.2d.tbl",svs,kk); }
    printf("NOTE: Writing final scaled stims table: %s .\n",stimfile);
    outfu = fopen(stimfile,"w");
    fprintf(outfu,"| yb | yap | rely       |st| x_as    | y_as    |  ya  | xa   | xamc | yamc | q    | yac        | xb   | y_asc   | t             |\n");
    /*              1234 12345 123456789012 12 123456789 123456789 123456 123456 123456 123456 123456 123456789012 123456 123456789 123456789012345 */
    for (kk=0; kk<4; ++kk) {
      for (iyb=0; iyb<8; ++iyb) { ybmean[iyb]=0.; ybnum[iyb]=0.; ybrms[iyb]=0.; }
      for (ii=0; ii<nSTPH[kk]; ++ii) {
        xraw0 = STPH[kk][ii].xb* FODC[0].NUV_XCLK + STPH[kk][ii].xamc ;
        yraw0 = STPH[kk][ii].yb* FODC[0].NUV_YCLK + STPH[kk][ii].yamc ;
        ya = ((int)(((yraw0/(2*FODC[0].NUV_YCLK) - xraw0/(2*FODC[0].NUV_XCLK)   )+10)*32)+STPH[kk][ii].xa) % 32;
        xraw = xraw0+((int)((STPH[kk][ii].xa+7)%32)-16)*FODC[0].NUV_XSLP;
        yraw = yraw0+((int)((STPH[kk][ii].ya+7)%32)-16)*FODC[0].NUV_YSLP;
        x = (xraw - FODC[0].NUV_XCEN)*(FODC[0].NUV_XSCL);   /* +or- 32725 to detsize */
        y = (yraw - FODC[0].NUV_YCEN)*(FODC[0].NUV_YSCL);
        x = (FODC[0].Mx * x) + FODC[0].Bx;
        y = (FODC[0].My * y) + FODC[0].By;
        x_as  = x*cal->aspum;
        y_as  = y*cal->aspum;
        yap   = rtaph_yap( ya, STPH[kk][ii].yb, STPH[kk][ii].yamc );     /* possible wrap around... -tab 07jul2010 */
        yac   = rtaph_yac( FODC[0], ya, STPH[kk][ii].yb, STPH[kk][ii].yamc );
        y_asc = (y - yac)*cal->aspum;
        for (iyb=0; iyb<8; ++iyb) {
          if ((STPH[kk][ii].yb == iyb)&&(ya > 4)) {
            ybmean[iyb] = ybmean[iyb]+ y_asc;
            ybnum[iyb]  = ybnum[iyb] + 1.;
          }
        }
        fprintf(outfu," %4d %5d %12.3f %2d %9.3f %9.3f %6d %6d %6d %6d %6d %12.7f %6d %9.3f %15.4f \n",
             STPH[kk][ii].yb, yap, (y - my[kk]),
             kk, x_as, y_as, ya, STPH[kk][ii].xa, STPH[kk][ii].xamc, STPH[kk][ii].yamc,
             STPH[kk][ii].q, yac, STPH[kk][ii].xb, y_asc, STPH[kk][ii].t );
      }
      for (iyb=0; iyb<8; ++iyb) {
        if (ybnum[iyb] > 20.) { ybmean[iyb] = ybmean[iyb] / ybnum[iyb]; }
      }
      for (ii=0; ii<nSTPH[kk]; ++ii) {
        xraw0 = STPH[kk][ii].xb* FODC[0].NUV_XCLK + STPH[kk][ii].xamc ;
        yraw0 = STPH[kk][ii].yb* FODC[0].NUV_YCLK + STPH[kk][ii].yamc ;
        ya = ((int)(((yraw0/(2*FODC[0].NUV_YCLK) - xraw0/(2*FODC[0].NUV_XCLK)   )+10)*32)+STPH[kk][ii].xa) % 32;
        yraw = yraw0+((int)((STPH[kk][ii].ya+7)%32)-16)*FODC[0].NUV_YSLP;
        y = (yraw - FODC[0].NUV_YCEN)*(FODC[0].NUV_YSCL);
        y = (FODC[0].My * y) + FODC[0].By;
        yac   = rtaph_yac( FODC[0], ya, STPH[kk][ii].yb, STPH[kk][ii].yamc );
        y_asc = (y - yac)*cal->aspum;
        for (iyb=0; iyb<8; ++iyb) {
          if ((STPH[kk][ii].yb == iyb)&&(ya > 4)) {
            if (ybnum[iyb] > 20.) { ybrms[iyb] = ybrms[iyb]+ ( (y_asc - ybmean[iyb]) * (y_asc - ybmean[iyb]) ); }
          }
        }
      }
      for (iyb=0; iyb<8; ++iyb) {
        if (ybnum[iyb] > 20.) {
          ybrms[iyb] = sqrt(( ybrms[iyb] / ybnum[iyb] ));
          printf(" Corrected Stim %d:  YB=%d  Num=%7.1f  Mean=%9.3f  RMS=%7.3f \n",kk,iyb,ybnum[iyb],ybmean[iyb],ybrms[iyb]);
        }
      }
    }
    fclose(outfu);
  }                         /* if (stimok == 1) { ... */

  /* Free. */
  for (kk=0; kk<4; ++kk) { free(STPH[kk]); STPH[kk] = NULL; }
  free(x8); x8 = NULL;
  free(y8); y8 = NULL;
  free(w8); w8 = NULL;
  return;
}


/* ----------------------------------------------------------------------
  Set up FUV Offset and Distortion Correction parameters and data images
  for rtaph().   -tab September 2009
*/
void rtaph_setup_FODC( PhRecRaw6Byte * ph6,  int nph6,  int band,  Cal * cal, FODCtype FODC[], int subvis, char calpath[] )
{
  /**/
  int ok, jj, ii, kk, i, anynull, max_cube;
  int old_nssd, nssd, nn,pp,pinc;
  int ngood, nreject;
  /**/
  short  xamc, yamc;
  /**/
  unsigned char q, xa, xb, yb;
  unsigned char ya;
  /**/
  const int max_ss = 2000;
  double *stim_sep, *stim_avt, *stim_num;
  double *ssdt, *ssds, *ssdn, *ssdf;
  /**/
  double avt,navt, sx1,sy1,nn1, sx2,sy2,nn2, sx3,sy3,nn3, sx4,sy4,nn4;
  double diff,rr1,rr2,rr3,rr4;
  double c11,c12,c13,c22,c23;
  double nnlim;                 /* new variable -tab 18nov2009 */
  /**/
  double x, y, x_as, y_as, tt;
  double xraw, yraw, xraw0, yraw0;
  /**/
  char svs[10];
  char comment[200];
  long nbuffer;
  long firstpixel = 1;
  int status = 0;
  /**/
  float nullval = 0;  /* don't check for null values in the image */
  /**/
  fitsfile *fptr;
  /**/
  char line[200];
  char caldir[200];
  char ssdfile[200];
  char tdcfile[200];
  char dcdxfile[200];
  char dcdyfile[200];
  /**/
  FILE *infu;
  FILE *outfu;
  /**/
  int max_STM;
  int nSTM;
  STMtype *STM;
  /**/
  /* ORIGINAL NUV average stim positions  -tab 14sep2009 */
  const double nuv_avx1 = -2722.27;  const double nuv_avy1 =  2453.89;
  const double nuv_avx2 =  2468.84;  const double nuv_avy2 =  2453.78;
  const double nuv_avx3 = -2721.87;  const double nuv_avy3 = -2565.81;
  const double nuv_avx4 =  2469.85;  const double nuv_avy4 = -2567.83;
  /* ORIGINAL FUV average stim positions -tab 15sep2009 */
  const double fuv_avx1 = -2541.88;  const double fuv_avy1 =  2455.28;
  const double fuv_avx2 =  2632.06;  const double fuv_avy2 =  2455.02;
  const double fuv_avx3 = -2541.53;  const double fuv_avy3 = -2550.89;
  const double fuv_avx4 =  2631.68;  const double fuv_avy4 = -2550.92;
  /**/ 
  
  /* Version note. */
  printf("NOTE: Running rtaph_setup_FODC() [version 22nov2010] \n");
  printf("NOTE: Nominal NUV average stim positions (arcsec) [%d]:\n",FODC[0].eclipse);
  printf(" x: %10.2f %10.2f %10.2f %10.2f   y: %10.2f %10.2f %10.2f %10.2f \n",
     nuv_avx1, nuv_avx2, nuv_avx3, nuv_avx4, nuv_avy1, nuv_avy2, nuv_avy3, nuv_avy4 );
  
  /* Check calpath.  Set string 'caldir'. */
  if (rtaph_FileExist(calpath) == 1) {
    printf("NOTE: user calpath exists '%s.'\n",calpath);
    strcpy(caldir,calpath);
  } else {
    printf("NOTE: user calpath does not exist '%s', useing 'cal'.\n",calpath);
    strcpy(caldir,"cal");
  }

  /* Defaults. */
  if (band == 1) { FODC[0].stim_coef0 = 5105.48; } else { FODC[0].stim_coef0 = 5089.75; }  /* new default. -tab 18nov2009 */
  FODC[0].stim_coef1 = 0.;
  FODC[0].fodx_coef_0= 0.;
  FODC[0].fody_coef_0= 0.;
  FODC[0].fodx_coef_1= 0.;
  FODC[0].fody_coef_1= 0.;

  /* Current working directory. */
  strcpy(svs,"");
  getcwd( line, 200 );
  printf("Current working directory = '%s'.\n",line);
  ii=0;
  while (ii < strlen(line)-8) {
  /* The following error check is obtuse, but it's looking to see if the   */
  /* working directory contains the '02-vsn' part that indicates that it's */
  /* a secondary tile and then also checks that the first two letters of   */
  /* the tile ID are '50' indicating that it is an AIS secondary tile and  */
  /* not a non-AIS secondary tile (which start at 50500).                  */
    if ((line[ii+0]=='0')&&(line[ii+1]=='2')&&(line[ii+2]=='-')&&(line[ii+3]=='v')&&(line[ii+4]=='s')&&(line[ii+5]=='n')&&(line[ii+7]=='5')&&(line[ii+8]=='0')&&(line[ii+9]!='5')&&(line[ii+9]!='6')&&(line[ii+9]!='7')&&(line[ii+9]!='8')&&(line[ii+9]!='9')) {
      if(subvis < 1) {
        fprintf(stderr,"*** %s: subvis = %d is out of range",iam,subvis);
        exit(1);
      }
      sprintf(svs,"_sg%2.2d",subvis);  
      ii=99999; 
    }
    ++ii;
  }
  if (strcmp(svs,"") != 0) { printf("NOTE: Using subvisit string = '%s'.\n",svs); }

  /* FUV Offset */
  if (band == 2 ) {
    sprintf(tdcfile,"%s/fuv_dy_fdttdc_coef_0.tbl",caldir);
    if (rtaph_FileExist(tdcfile) == 0) {
      fprintf(stderr,"===WARNING: rtaph_setup_FODC: Cannot find FUV offset coeffecient file '%s'.\n",tdcfile);
    } else {
      infu = fopen(tdcfile,"r");
      fgets(line,sizeof(line),infu);
      while(fscanf(infu,"%lf %lf",&rr1,&rr2) != EOF) {
        ii = rr1;
        if (FODC[0].eclipse == ii) { FODC[0].fody_coef_0 = rr2; }
      }
      fclose(infu);
      FODC[0].fody_coef_1 = 0.3597;
    }
    sprintf(tdcfile,"%s/fuv_dx_fdttdc_coef_0.tbl",caldir);
    if (rtaph_FileExist(tdcfile) == 0) {
      fprintf(stderr,"===WARNING: rtaph_setup_FODC: Cannot find FUV offset coeffecient file '%s'.\n",tdcfile);
    } else {
      infu = fopen(tdcfile,"r");
      fgets(line,sizeof(line),infu);
      while(fscanf(infu,"%lf %lf",&rr1,&rr2) != EOF) {
        ii = rr1;
        if (FODC[0].eclipse == ii) { FODC[0].fodx_coef_0 = rr2; }
      }
      fclose(infu);
      FODC[0].fodx_coef_1 = 0.;
    }
    printf("For FUV offset, using: eclipse=%d  fodx=%9.5f %9.5f   fody=%9.5f %9.5f  \n",FODC[0].eclipse,
            FODC[0].fodx_coef_0, FODC[0].fodx_coef_1, FODC[0].fody_coef_0, FODC[0].fody_coef_1 );
  }
  
  /* Set filename for data cube. Read modified data for post May 4th, 2010 data.  -tab 19oct2010 */
  if (FODC[0].eclipse > 37460) {
    if (band == 1) {
      printf("Reading modified distortion data cube with raw_stimsep = %10.3f .\n",FODC[0].raw_stimsep);
      if (FODC[0].raw_stimsep < 5136.30) {
        sprintf(dcdxfile,"%s/nuv_distortion_cube_dxa.fits",caldir); 
        sprintf(dcdyfile,"%s/nuv_distortion_cube_dya.fits",caldir); 
      } else { 
        if (FODC[0].raw_stimsep < 5137.25) {
          sprintf(dcdxfile,"%s/nuv_distortion_cube_dxb.fits",caldir); 
          sprintf(dcdyfile,"%s/nuv_distortion_cube_dyb.fits",caldir); 
        } else {
          sprintf(dcdxfile,"%s/nuv_distortion_cube_dxc.fits",caldir); 
          sprintf(dcdyfile,"%s/nuv_distortion_cube_dyc.fits",caldir); 
        }
      }
    } else { 
      sprintf(dcdxfile,"%s/fuv_distortion_cube_dx.fits",caldir); 
      sprintf(dcdyfile,"%s/fuv_distortion_cube_dy.fits",caldir);
    }
  } else {
    if (band == 1) { sprintf(dcdxfile,"%s/nuv_distortion_cube_dx.fits",caldir); sprintf(dcdyfile,"%s/nuv_distortion_cube_dy.fits",caldir); }
              else { sprintf(dcdxfile,"%s/fuv_distortion_cube_dx.fits",caldir); sprintf(dcdyfile,"%s/fuv_distortion_cube_dy.fits",caldir); }
  }

  /* #@# Do these lines do anything? -tab 19oct2010. */
  fits_open_file( &fptr, dcdxfile, READONLY, &status );                           fits_report_error(stderr,status); if (status != 0) exit(1);
  fits_read_card(fptr, "NAXIS", line, &status );                                  fits_report_error(stderr,status); if (status != 0) exit(1);
  fits_read_card(fptr, "NAXIS1", line, &status );                                 fits_report_error(stderr,status); if (status != 0) exit(1);
  fits_read_card(fptr, "NAXIS2", line, &status );                                 fits_report_error(stderr,status); if (status != 0) exit(1);
  fits_close_file( fptr, &status );                                               fits_report_error(stderr,status); if (status != 0) exit(1);
  
  /* Get parameters from FITS cube image header. */
  fits_open_file( &fptr, dcdxfile, READONLY, &status );                           fits_report_error(stderr,status); if (status != 0) exit(1);
  fits_read_key( fptr, TINT   , "NAXIS1" , &FODC[0].cube_nc , comment, &status ); fits_report_error(stderr,status); if (status != 0) exit(1);
  fits_read_key( fptr, TINT   , "NAXIS2" , &FODC[0].cube_nr , comment, &status ); fits_report_error(stderr,status); if (status != 0) exit(1);
  fits_read_key( fptr, TINT   , "NAXIS3" , &FODC[0].cube_nd , comment, &status ); fits_report_error(stderr,status); if (status != 0) exit(1);
  fits_read_key( fptr, TDOUBLE, "DC_X0"  , &FODC[0].cube_x0 , comment, &status ); fits_report_error(stderr,status); if (status != 0) exit(1);
  fits_read_key( fptr, TDOUBLE, "DC_Y0"  , &FODC[0].cube_y0 , comment, &status ); fits_report_error(stderr,status); if (status != 0) exit(1);
  fits_read_key( fptr, TDOUBLE, "DC_D0"  , &FODC[0].cube_d0 , comment, &status ); fits_report_error(stderr,status); if (status != 0) exit(1);
  fits_read_key( fptr, TDOUBLE, "DC_DX"  , &FODC[0].cube_dx , comment, &status ); fits_report_error(stderr,status); if (status != 0) exit(1);
  fits_read_key( fptr, TDOUBLE, "DC_DY"  , &FODC[0].cube_dy , comment, &status ); fits_report_error(stderr,status); if (status != 0) exit(1);
  fits_read_key( fptr, TDOUBLE, "DC_DD"  , &FODC[0].cube_dd , comment, &status ); fits_report_error(stderr,status); if (status != 0) exit(1);
  fits_close_file( fptr, &status );                                               fits_report_error(stderr,status); if (status != 0) exit(1);
  
  /* Allocate distortion cubes. */
  max_cube = ( FODC[0].cube_nc * FODC[0].cube_nr * FODC[0].cube_nd ) + 1000;
  printf("Allocate distortion cube memory [%d floats each cube].\n",max_cube);
  FODC[0].dcdx = (float *)calloc(max_cube,sizeof(float));
  FODC[0].dcdy = (float *)calloc(max_cube,sizeof(float));
  
  /* Read distortion cubes. */
  nbuffer = FODC[0].cube_nc * FODC[0].cube_nr * FODC[0].cube_nd;
  printf("Reading '%s'.\n",dcdxfile);
  fits_open_file( &fptr, dcdxfile, READONLY, &status );                                 fits_report_error(stderr,status); if (status != 0) exit(1);
  fits_read_img(fptr,TFLOAT,firstpixel,nbuffer,&nullval,FODC[0].dcdx,&anynull,&status); fits_report_error(stderr,status); if (status != 0) exit(1);
  fits_close_file( fptr, &status );                                                     fits_report_error(stderr,status); if (status != 0) exit(1);
  printf("Reading '%s'.\n",dcdyfile);
  fits_open_file( &fptr, dcdyfile, READONLY, &status );                                 fits_report_error(stderr,status); if (status != 0) exit(1);
  fits_read_img(fptr,TFLOAT,firstpixel,nbuffer,&nullval,FODC[0].dcdy,&anynull,&status); fits_report_error(stderr,status); if (status != 0) exit(1);
  fits_close_file( fptr, &status );                                                     fits_report_error(stderr,status); if (status != 0) exit(1);
  
  /* Set. */
  nSTM=0;
  ok = 1;

  /* Allocate big STM array. */
  if (band == 1) { max_STM = (nph6 / 25) + 1000; }
            else { max_STM = (nph6 /  5) + 1000; }
  if (max_STM < 900000) max_STM = 900000;                /* Added. -tab 18feb2010 */
  STM = (STMtype *)calloc(max_STM,sizeof(STMtype));
  for (ii=0; ii<max_STM; ++ii) { STM[ii].stim = 0; }

  /* Allocate small stim group arrays. */
  stim_sep = (double *)calloc(max_ss,sizeof(double));
  stim_avt = (double *)calloc(max_ss,sizeof(double));
  stim_num = (double *)calloc(max_ss,sizeof(double));
  ssdt     = (double *)calloc(max_ss,sizeof(double));
  ssds     = (double *)calloc(max_ss,sizeof(double));
  ssdn     = (double *)calloc(max_ss,sizeof(double));
  ssdf     = (double *)calloc(max_ss,sizeof(double));

  /* Load all stim photons. */
  printf("Load all stim photons (max = %d).\n",max_STM);
  for(i=0; i<nph6; i++, ph6++) {
  /* ========= Unpack raw6 into timing and charge data */
    q    = ((ph6->phb4 & 0x03) << 3) + (( ph6->phb5 & 0xe0) >> 5);
    xb   = ph6->phb1 >> 5;
    xamc = (short) ((ph6->phb1 & 0x1f) << 7) +
           (short) ((ph6->phb2 & 0xfe) >> 1) -
           (short) ((ph6->phb1 & 0x10) << 8);
    yb   = ((ph6->phb2 & 0x01) << 2) + ((ph6->phb3 & 0xc0) >> 6);
    yamc = (short) ((ph6->phb3 & 0x3f) << 6) +
           (short) ((ph6->phb4 & 0xfc) >> 2) -
           (short) ((ph6->phb3 & 0x20) << 7);
    xa   = ((ph6->phb5 & 0x10) >> 4) + ((ph6->phb5 & 0x03)  << 3) +
           ((ph6->phb5 & 0x0c) >> 1);
    tt   = ph6->t;
  /* ==========  */

    if (band == 1) {  /* ................... NUV ..................... */
      xraw0 = xb*FODC[0].NUV_XCLK + xamc ;
      yraw0 = yb*FODC[0].NUV_YCLK + yamc ;
      ya = ((int)(((yraw0/(2*FODC[0].NUV_YCLK) - xraw0/(2*FODC[0].NUV_XCLK)   )+10)*32)+xa) % 32;
      xraw = xraw0+((int)((xa+7)%32)-16)*FODC[0].NUV_XSLP;
      yraw = yraw0+((int)((ya+7)%32)-16)*FODC[0].NUV_YSLP;

  /* Get detector position in microns using a standard scale and offset. */
      x = (xraw - FODC[0].NUV_XCEN)*(FODC[0].NUV_XSCL);   /* +or- 32725 to detsize */
      y = (yraw - FODC[0].NUV_YCEN)*(FODC[0].NUV_YSCL);

  /* Rescale and adjust post May 4th, 2010 data only.  -tab 05oct2010 */
      if (FODC[0].eclipse > 37460) {

        /* Rescale x,y values so that STIMs match in new data... -tab 01jul2010 */
        x = (FODC[0].Mx * x) + FODC[0].Bx;
        y = (FODC[0].My * y) + FODC[0].By;

        /* Adjust y value (NUV only). -tab 24jun2010 */
        y = y - rtaph_yac( FODC[0], ya, yb, yamc );

        /* No need to to use rtaph_yac2() here since it doesn't correct stim photons. -tab 24sep2010 */

      }

  /* Convert detector position from microns to arcsecs  [NUV] */
      x_as = x*cal->aspum;
      y_as = y*cal->aspum;
  /* Load. */
      kk=-1;
      if ((x_as > (nuv_avx1 - 20.))&&(x_as < (nuv_avx1 + 20.))&&(y_as > (nuv_avy1 - 20.))&&(y_as < (nuv_avy1 + 20.))) { kk=0; }
      if ((x_as > (nuv_avx2 - 20.))&&(x_as < (nuv_avx2 + 20.))&&(y_as > (nuv_avy2 - 20.))&&(y_as < (nuv_avy2 + 20.))) { kk=1; }
      if ((x_as > (nuv_avx3 - 20.))&&(x_as < (nuv_avx3 + 20.))&&(y_as > (nuv_avy3 - 20.))&&(y_as < (nuv_avy3 + 20.))) { kk=2; }
      if ((x_as > (nuv_avx4 - 20.))&&(x_as < (nuv_avx4 + 20.))&&(y_as > (nuv_avy4 - 20.))&&(y_as < (nuv_avy4 + 20.))) { kk=3; }
      if (kk > -1) {
        STM[nSTM].stim = kk + 1;
        STM[nSTM].t    = tt;
        STM[nSTM].x_as = x_as;
        STM[nSTM].y_as = y_as;
        ++nSTM;
        if (nSTM > (max_STM - 10)) { fprintf(stderr,"***ERROR:rtaph_setup_FODC:1: Too many stim photons.\n"); exit(1); }
      }
    } else {         /* ................... FUV ..................... */
      xraw0 = xb*cal->xclk_fuv + xamc ;
      yraw0 = yb*cal->yclk_fuv + yamc ;
      ya = ((int)(((yraw0/(2*cal->yclk_fuv) -
                    xraw0/(2*cal->xclk_fuv)   )+10)*32)+xa) % 32;
      xraw = xraw0+((int)((xa+7)%32)-16)*cal->xslp_fuv;
      yraw = yraw0+((int)((ya+7)%32)-16)*cal->yslp_fuv;

  /* Get detector position in microns using a standard scale and offset. */
      x = (xraw - cal->xcen_fuv)*(cal->xscl_fuv);
      y = (yraw - cal->ycen_fuv)*(cal->yscl_fuv);
  /* Convert detector position from microns to arcsecs [FUV] */
      x_as = x*cal->aspum;
      y_as = y*cal->aspum;
  /* Load. */
      kk=-1;
      if ((x_as > (fuv_avx1 - 20.))&&(x_as < (fuv_avx1 + 20.))&&(y_as > (fuv_avy1 - 20.))&&(y_as < (fuv_avy1 + 20.))) { kk=0; }
      if ((x_as > (fuv_avx2 - 20.))&&(x_as < (fuv_avx2 + 20.))&&(y_as > (fuv_avy2 - 20.))&&(y_as < (fuv_avy2 + 20.))) { kk=1; }
      if ((x_as > (fuv_avx3 - 20.))&&(x_as < (fuv_avx3 + 20.))&&(y_as > (fuv_avy3 - 20.))&&(y_as < (fuv_avy3 + 20.))) { kk=2; }
      if ((x_as > (fuv_avx4 - 20.))&&(x_as < (fuv_avx4 + 20.))&&(y_as > (fuv_avy4 - 20.))&&(y_as < (fuv_avy4 + 20.))) { kk=3; }
      if (kk > -1) {
        STM[nSTM].stim = kk + 1;
        STM[nSTM].t    = tt;
        STM[nSTM].x_as = x_as;
        STM[nSTM].y_as = y_as;
        ++nSTM;
        if (nSTM > (max_STM - 10)) { fprintf(stderr,"***ERROR:rtaph_setup_FODC:2: Too many stim photons.\n"); exit(1); }
      }
    }
  }
  printf("Loaded %d stim photons.\n",nSTM);

  /* Check. */
  if (nSTM < 100) {
    fprintf(stderr,"===WARNING:rtaph_setup_FODC: not enough stim photons found.\n");
    ok=0;
  }

  /* Continue, if no problems. */
  nssd=0;
  if (ok == 1) {
  
  /* Compute new stim data by dividing into groups. */
    printf("Compute new stim data.\n");
    nn = 0;
    kk = ((double)nSTM / 1000.);
    if (kk < 1) kk = 1;
    pinc = (double)nSTM / (double)kk;
    pp = 0;
    while (pp < nSTM) {
      sx1=0.; sy1=0.; sx2=0.; sy2=0.; sx3=0.; sy3=0.; sx4=0.; sy4=0.; nn1=0.; nn2=0.; nn3=0.; nn4=0.;
      avt=0.; navt=0.;
      for (ii=0; ii<pinc; ++ii) {
        jj = pp + ii;
        if (jj < nSTM) {
          if (STM[jj].stim == 1) { sx1 = sx1 + STM[jj].x_as;   sy1 = sy1 + STM[jj].y_as;  nn1 = nn1 + 1.; }
          if (STM[jj].stim == 2) { sx2 = sx2 + STM[jj].x_as;   sy2 = sy2 + STM[jj].y_as;  nn2 = nn2 + 1.; }
          if (STM[jj].stim == 3) { sx3 = sx3 + STM[jj].x_as;   sy3 = sy3 + STM[jj].y_as;  nn3 = nn3 + 1.; }
          if (STM[jj].stim == 4) { sx4 = sx4 + STM[jj].x_as;   sy4 = sy4 + STM[jj].y_as;  nn4 = nn4 + 1.; }
          if (STM[jj].stim > 0 ) { avt = avt + STM[jj].t;      navt= navt + 1.;                           }
        }
      }
      if (nn > 1) { nnlim = 100.; } else {        /* default for many stim groups. -tab 18nov2009 */
        nnlim = (double)pinc / 10.;               /* Avoiding nn=0. -tab 18nov2009 */
        if (nnlim > 100.) nnlim = 100.;           /* just in case. -tab 18nov2009 */
      }
      if ((nn1 > nnlim)&&(nn2 > nnlim)&&(nn3 > nnlim)&&(nn4 > nnlim)&&(navt > nnlim)) {    /* line change... -tab 18nov2009 */
        sx1 = sx1 / nn1;  sy1 = sy1 / nn1;
        sx2 = sx2 / nn2;  sy2 = sy2 / nn2;
        sx3 = sx3 / nn3;  sy3 = sy3 / nn3;
        sx4 = sx4 / nn4;  sy4 = sy4 / nn4;
        stim_sep[nn] =  ( ( (sx2 - sx1) + (sx4 - sx3) + (sy1 - sy3) + (sy2 - sy4) ) / 4. );
        stim_avt[nn] =  ( avt / navt );
        stim_num[nn] =  ( nn1 + nn2 + nn3 + nn3 );
        ++nn;
        if (nn > max_ss - 10) { fprintf(stderr,"***ERROR:rtaph_setup_FODC: too many stim_sep[] entries, max=%d.\n",max_ss); exit(1); }
      }
      pp = pp + pinc;
    }
    printf("NOTE: rtaph_setup_FODC: %d stim groups set (pinc=%d).\n",nn,pinc);
  
  /* Read in past stim separation data. */
  /*  example:
  |sct           |stim_sep |stim_num | sep_fit |
    752426355.885  5105.713  331199.2  5105.834
   12345678901234 123456789 123456789 123456789
  */
    nssd = 0;   /* number of stim separation data */
    if (band == 1) { 
      sprintf(ssdfile,"SSD_nuv_%5.5d%s.tbl",FODC[0].eclipse,svs); 
    } else { 
      sprintf(ssdfile,"SSD_fuv_%5.5d%s.tbl",FODC[0].eclipse,svs); 
    }
/* Only include old points that occur before first current data point. */
    if (rtaph_FileExist(ssdfile) == 1) {
      infu = fopen(ssdfile,"r");
      fgets(line,sizeof(line),infu);
      while(fscanf(infu,"%lf %lf %lf %lf",&rr1,&rr2,&rr3,&rr4) != EOF) {
        if (rr1 < (stim_avt[0] - 1.)) {
          ssdt[nssd] = rr1;
          ssds[nssd] = rr2;
          ssdn[nssd] = rr3;
          ssdf[nssd] = rr4;
          ++nssd;
          if (nssd > max_ss - 10) { fprintf(stderr,"***ERROR:rtaph_setup_FODC:1: too many ssdt[] entries, max=%d.\n",max_ss); exit(1); }
        }
      }
      fclose(infu);
    }
    printf("Read in %d previous entries from '%s'.\n",nssd,ssdfile);

  /* Save old number of entries. */
    old_nssd = nssd;

  /* Load new entries. */
    for (ii=0; ii<nn; ++ii) {
      ssdt[nssd] = stim_avt[ii];
      ssds[nssd] = stim_sep[ii];
      ssdn[nssd] = stim_num[ii];
      ssdf[nssd] = 0.;
      ++nssd;
      if (nssd > max_ss - 10) { fprintf(stderr,"***ERROR:rtaph_setup_FODC:2: too many ssdt[] entries, max=%d.\n",max_ss); exit(1); }
    }
  
  /* FIRST: Fit straight line. */
    if (nssd == 0) { 
      FODC[0].stim_coef1 = 0.;
      if (band == 1) { FODC[0].stim_coef0 = 5105.48; } else { FODC[0].stim_coef0 = 5089.75; }
    } else {
    if (nssd == 1) { 
      FODC[0].stim_coef1 = 0.; 
      FODC[0].stim_coef0 = ssds[0]; 
    } else {
      c11=0.; c12=0.; c13=0.; c22=0.; c23=0.;
      for (ii=0; ii<nssd; ++ii) {
        if (ssdn[ii] > 0.) {
          c11 = c11 + ssdn[ii];                           /* weight sum */
          c12 = c12 + ssdt[ii] * ssdn[ii];                /* x times weight */
          c13 = c13 + ssds[ii] * ssdn[ii];                /* y times weight */
          c22 = c22 + ssdt[ii] * ssdt[ii] * ssdn[ii];
          c23 = c23 + ssds[ii] * ssdn[ii] * ssdt[ii];
        }
      }
      FODC[0].stim_coef1 = ( (c13 * c12) - (c23 * c11) ) / ( (c12 * c12) - (c22 * c11) );
      FODC[0].stim_coef0 = ( c13 - (c12 * FODC[0].stim_coef1) ) / c11;
    }}

  /* Good points. */
    ngood=0; for (ii=0; ii<nssd; ++ii) { if (ssdn[ii] > 0.) ++ngood; }

  /* Reject? */
    nreject=0;
    for (ii=0; ii<nssd; ++ii) { 
      diff = ( FODC[0].stim_coef0  + ( FODC[0].stim_coef1  *  ssdt[ii] ) ) - ssds[ii];   
      if (diff < 0.) diff = -1. * diff;
      if ((diff > 1.5)&&(nreject < (ngood - 3))) { ssdn[ii]=0.; ++nreject; }
    }
    if (nreject > 0.) printf("Rejected %d point(s) during fitting of stim values.\n",nreject);

  /* FINAL: Fit straight line. */
    if (nssd == 0) { 
      printf("No points to fit-- use typical stim separation values.\n");
      FODC[0].stim_coef1 = 0.;
      if (band == 1) { FODC[0].stim_coef0 = 5105.48; } else { FODC[0].stim_coef0 = 5089.75; }
    } else {
    if (nssd == 1) {
      printf("Only 1 point found-- use single value (flat) solution.\n");
      FODC[0].stim_coef1 = 0.; 
      FODC[0].stim_coef0 = ssds[0];
    } else {
      c11=0.; c12=0.; c13=0.; c22=0.; c23=0.;
      for (ii=0; ii<nssd; ++ii) {
        if (ssdn[ii] > 0.) {
          c11 = c11 + ssdn[ii];                           /* weight sum */
          c12 = c12 + ssdt[ii] * ssdn[ii];                /* x times weight */
          c13 = c13 + ssds[ii] * ssdn[ii];                /* y times weight */
          c22 = c22 + ssdt[ii] * ssdt[ii] * ssdn[ii];
          c23 = c23 + ssds[ii] * ssdn[ii] * ssdt[ii];
        }
      }
      FODC[0].stim_coef1 = ( (c13 * c12) - (c23 * c11) ) / ( (c12 * c12) - (c22 * c11) );
      FODC[0].stim_coef0 = ( c13 - (c12 * FODC[0].stim_coef1) ) / c11;
    }}
    printf("Stim solution coeffecients(band=%d): %14.7e %14.7e [n=%d]\n",band,FODC[0].stim_coef0,FODC[0].stim_coef1,nssd); 

  /* Load fit only for new data. */
    for (ii=old_nssd; ii<nssd; ++ii) { ssdf[ii] = FODC[0].stim_coef0  + ( FODC[0].stim_coef1  *  ssdt[ii] ); }

  /* Write out old and new data. */
    printf("Writing stim data to '%s' [old=%d tot=%d].\n",ssdfile,old_nssd,nssd);
    outfu = fopen(ssdfile,"w");
    fprintf(outfu,"|sct           |stim_sep |stim_num | sep_fit |\n");
    for (ii=0; ii<nssd; ++ii) {
      fprintf(outfu," %14.3f %9.3f %9.1f %9.3f \n",ssdt[ii],ssds[ii],ssdn[ii],ssdf[ii]);
    }
    fclose(outfu);

  /* Print stim data points. */
    printf("\n|sct           |stim_sep |stim_num | sep_fit |\n");
    for (ii=0; ii<nssd; ++ii) {
      printf(" %14.3f %9.3f %9.1f %9.3f \n",ssdt[ii],ssds[ii],ssdn[ii],ssdf[ii]);
    }
    printf("\n");
  
  }

  /* Use default for post-CSP(May 4, 2010) GALEX data.  -tab 27sep2010 */
  if (FODC[0].eclipse > 37460) {
    if (band == 1) { FODC[0].stim_coef0 = 5105.48; } else { FODC[0].stim_coef0 = 5089.75; }  /* new default. -tab 18nov2009 */
    FODC[0].stim_coef1 = 0.;
    printf("Default post-CSP(May 4, 2010) data: Stim coeffecients(band=%d): %14.7e %14.7e [n=%d]\n",band,FODC[0].stim_coef0,FODC[0].stim_coef1,nssd); 
  }

  /* End. */
  printf("Return to rtaph(), band=%d, ...\n",band);
  free(STM);       STM      = NULL;
  free(stim_sep);  stim_sep = NULL;
  free(stim_avt);  stim_avt = NULL;
  free(stim_num);  stim_num = NULL;
  free(ssdt);      ssdt     = NULL;
  free(ssds);      ssds     = NULL;
  free(ssdn);      ssdn     = NULL;
  free(ssdf);      ssdf     = NULL;
  return;
}


/* ----------------------------------------------------------------------
  Set some basic parameters. -tab 19jul2010  21sep2010
*/
void rtaph_basic_FODC( Cal * cal, FODCtype FODC[] )
{
/**/
/**/
FODC[0].NUV_XCLK = cal->xclk_nuv;   FODC[0].NUV_YCLK = cal->yclk_nuv;
FODC[0].NUV_XSLP = cal->xslp_nuv;   FODC[0].NUV_YSLP = cal->yslp_nuv;
FODC[0].NUV_XCEN = cal->xcen_nuv;   FODC[0].NUV_YCEN = cal->ycen_nuv;
FODC[0].NUV_XSCL = cal->xscl_nuv;   FODC[0].NUV_YSCL = cal->yscl_nuv;
/*
Email: Date: Wed, 23 Jun 2010 17:08:00 -0700:
  Change the NUV YCLK from 1992 to 2019.   -Patrick
Tom: I found that NUV_YCLK = 2016 gives the best results. -tab 23jul2010...
*/

if (FODC[0].eclipse >= 38150) { FODC[0].NUV_YCLK = 2016; }  /* reset occurred at about 2010-06-22T23:40Z -tab 25jun2010 */

/* Echo. */
printf("FODC: NUV_X/YCLK = %20.12f  %20.12f \n",FODC[0].NUV_XCLK,FODC[0].NUV_YCLK);
printf("FODC: NUV_X/YSLP = %20.12f  %20.12f \n",FODC[0].NUV_XSLP,FODC[0].NUV_YSLP);
printf("FODC: NUV_X/YCEN = %20.12f  %20.12f \n",FODC[0].NUV_XCEN,FODC[0].NUV_YCEN);
printf("FODC: NUV_X/YSCL = %20.12f  %20.12f \n",FODC[0].NUV_XSCL,FODC[0].NUV_YSCL);
printf("FODC: cal->aspum = %20.12f \n",cal->aspum);

return;
}



int rtaph(PhRecRaw6Byte * ph6, int nph6, int band,
          Aspect * asp, int nasp,
          Cal * cal,
          PhRecX * phots_x, int *nx, int *errcnt,
          NamedList *state,
          NamedList *opts)
{
  short  xamc, yamc;
  unsigned char val,q, xa, xb, yb;
  unsigned short flags,aspskip,asphvlow;
  int maxerrs,i,j,k,col,row,pix,nmasked,naspslice,nvalid,hasasp;
  int ptrx, ptry;
  int subvis, ix, stim, rawxy, nlcmask, eclipse;
  int verbose,quiet,masked;
  double t,dt,ra,dec,aspra,aspdec,detsize,tslice,fwhm,r,ang,droll;
  double ra0,dec0,roll0;
  double fdttdc;
  char *calpath=0;
  double phra,phdec,*nuvpltsol,*fuvpltsol;
  double tmp,lt,tra,tdec,lra,ldec,maskfill;
  double pltscl,aspum;
  double x,y,xp,yp,xpp,ypp,x_as,y_as,xp_as,yp_as,xpp_as,ypp_as;
  unsigned char ya;
  double xdig,ydig,xdig_as,ydig_as,wgx,wgy;
  double xsc,ysc, wkx, wky, dxi, deta, dx, dy, xraw, yraw, xraw0,yraw0,xi,eta;
  double fptrx,fptry,blt,blu,resp;
  double stimdxy;
  double *asptx, *aspty;
  Aspect *asp0;
  FODCtype FODC;

  WIG2type *WIG2;
  WLK2type *WLK2;
  CLK2type *CLK2;

  struct { double X,Y,Z; } vec, last;
  struct celprm cel;
  struct prjprm prj;
  struct linprm lin;
  double xi_xsc,xi_ysc,eta_xsc,eta_ysc;                              
  
  double yac;

  //setlinebuf(stdout);

  droll      = nm_getn(  opts,"droll",    0);
  // Store raw detector x,y in the extended photon list instead of
  // walk-corrected xp,yp.
  rawxy      = nm_getn(  opts,"rawxy",    0);
  // Mask is in post-NLC coordinates.
  nlcmask    = nm_getn(  opts,"nlcmask",  1);
  stimdxy    = nm_getn(  opts,"stimdxy",  500);
  hasasp     = ! nm_getn(opts,"noasp",    0);
  ra0        = nm_getn(  opts,"ra0",      0); // used if ! hasasp
  dec0       = nm_getn(  opts,"dec0",     0); // used if ! hasasp
  roll0      = nm_getn(  opts,"twist0",   0); // used if ! hasasp
  verbose    = nm_getn(  opts,"verbose",  0);
  quiet      = nm_getn(  opts,"quiet",    0);
  eclipse    = nm_getn(  opts,"eclipse",  0);
  fdttdc     = nm_getn(  opts,"fdttdc",   0);  
  calpath    = nm_gets(  opts,"calpath",  0);
  subvis     = nm_getn(  opts,"subvis",   0);

/* Take out this check, moved to 'rtaph_FODC_offset_plus_shift()' ... -tab 03apr2010
  if(fdttdc <= 20. || fdttdc >= 40.) {
    fprintf(stderr,"*** %s: fdttdc = %lf is out of range",iam,fdttdc); exit(1);
  }
*/
  if(eclipse < 100) {
    fprintf(stderr,"*** %s: eclipse = %d is out of range",iam,eclipse);
    exit(1);
  }
  printf("eclipse=%d  fdttdc=%9.5f  subvis=%d  calpath='%s'\n",eclipse,fdttdc,subvis,calpath);

  if(hasasp && (! asp || nasp <= 0)) {
    fprintf(stderr,"*** %s: No aspect records supplied.\n",iam);
    exit(1);
  }

  if(! hasasp) {
    /* Make fake aspect records */
    /* This is recreated for each call. Yuch. */
    asp = (Aspect *)malloc_safe(2000*sizeof(Aspect));
    for(t=ph6->t-10, i=0, nasp=0; i<2000; t+=1, ++i, ++nasp) {
      (asp+i)->t = t;
      (asp+i)->raz = ra0;
      (asp+i)->decz = dec0;
      (asp+i)->twistz = roll0;
      (asp+i)->twistzr = roll0/R2D;
      (asp+i)->flags = 0;
    }
  }

  asp0 = asp;

  if(cal == NULL || cal->wlk_nuv_x == NULL) { // || cal->maskdef.npixx <= 100) {
    fprintf(stderr,"*** %s: No cal information supplied.\n",iam);
    exit(1);
  }

  if(band < 1 || band > 2) {
    fprintf(stderr,"*** %s: No valid band supplied (band=%d).\n",iam,band);
    exit(1);
  }

  if(state == NULL) {
    fprintf(stderr,"*** %s: State not defined.\n",iam);
    exit(1);
  }

  if(hasasp) {
    if((asptx = nm_getp(state,"asptx",NULL)) != NULL) {
      aspty = nm_getp(state,"aspty",NULL);
    }
    else {
      asptx = (double *) malloc_safe(nasp*sizeof(double));
      aspty = (double *) malloc_safe(nasp*sizeof(double));
      asprates(asp, nasp, asptx, aspty, opts);
      nm_add(state,"asptx(p)",asptx);
      nm_add(state,"aspty(p)",aspty);
    }
  }

  //gnomonicsetup(&cel,&prj,&lin); // replaced by gnomfwd/rev_simple calls

  pltscl  = cal->pltscl; // arc-seconds/mm
  aspum   = cal->aspum;  // arc-seconds/um
  detsize = cal->detsize; // nominal detector size in degrees
  dt = 1.0; // Aspect file time step
  nuvpltsol = cal->nuvpltsol;
  fuvpltsol = cal->fuvpltsol;
  xi_xsc  = cal->xi_xsc;
  xi_ysc  = cal->xi_ysc;
  eta_xsc = cal->eta_xsc;
  eta_ysc = cal->eta_ysc;

  asphvlow = band == 1 ? ASPHVLOWN : ASPHVLOWF;

  if(cal && cal->mask) {
    maskfill = detsize/(cal->maskdef.npixx*fabs(cal->maskdef.pixszx));
  }
  nmasked = 0;
  nvalid  = 0;
  aspskip = 0;

  /* Allocate and read new calibrations. -tab 05oct2010 */
  if (eclipse > 37460) {
    WIG2 = (WIG2type *)calloc(MAXWIG2,sizeof(WIG2type));
    WLK2 = (WLK2type *)calloc(MAXWLK2,sizeof(WLK2type));
    CLK2 = (CLK2type *)calloc(MAXCLK2,sizeof(CLK2type));
    rtaph_read_calibrations( WIG2, WLK2, CLK2, calpath, eclipse );
  } else {
    WIG2 = NULL;
    WLK2 = NULL;
    CLK2 = NULL;
  }

  /* Set FODC values.  -tab 19jul2010  21sep2010 */
  printf("Set FODC values.\n");
  FODC.eclipse = eclipse;
  FODC.fdttdc  = fdttdc;
  rtaph_basic_FODC( cal, &FODC );
  FODC.Mx=1.;  FODC.Bx=0.; FODC.My=1.;  FODC.By=0.;  FODC.raw_stimsep=0.;

  if (eclipse > 37460) { rtaph_setup_FODC_init( ph6, nph6, band, cal, &FODC, subvis, calpath ); }

  rtaph_setup_FODC( ph6, nph6, band, cal, &FODC, subvis, calpath );

  for(*nx=0, i=0 ; i < nph6 ; i++, ph6++) {

    // ========= Unpack raw6 into timing and charge data

    q    = ((ph6->phb4 & 0x03) << 3) + (( ph6->phb5 & 0xe0) >> 5);
    xb   = ph6->phb1 >> 5;
    xamc = (short) ((ph6->phb1 & 0x1f) << 7) +
           (short) ((ph6->phb2 & 0xfe) >> 1) -
           (short) ((ph6->phb1 & 0x10) << 8);
    yb   = ((ph6->phb2 & 0x01) << 2) + ((ph6->phb3 & 0xc0) >> 6);
    yamc = (short) ((ph6->phb3 & 0x3f) << 6) +
           (short) ((ph6->phb4 & 0xfc) >> 2) -
           (short) ((ph6->phb3 & 0x20) << 7);
    xa   = ((ph6->phb5 & 0x10) >> 4) + ((ph6->phb5 & 0x03)  << 3) +
           ((ph6->phb5 & 0x0c) >> 1);

    //    XA[counter]= ((RAW6.phb5 & 0x10) >> 4) + (( RAW6.phb5 & 0x03)  << 3) +
    //                 ((RAW6.phb5 & 0x0c) >> 1);

    flags = 0;
    masked = 0;
    maxerrs = 0;

    // ========= Rectify each band separately

    if (band == 1) {

      // ========= Band 1

      xraw0 = xb*FODC.NUV_XCLK + xamc ;
      yraw0 = yb*FODC.NUV_YCLK + yamc ;

      ya = ((int)(((yraw0/(2*FODC.NUV_YCLK) - xraw0/(2*FODC.NUV_XCLK)   )+10)*32)+xa) % 32;

      xraw = xraw0+((int)((xa+7)%32)-16)*FODC.NUV_XSLP;
      yraw = yraw0+((int)((ya+7)%32)-16)*FODC.NUV_YSLP;

      // See if this is a stim
      //          band v      v scale     options v
      stim = whichstim(1,NULL,0,xraw,yraw,stimdxy,NULL);

      // Get detector position in microns using a standard scale and offset.
      x = (xraw - FODC.NUV_XCEN)*(FODC.NUV_XSCL); // +or- 32725 to detsize
      y = (yraw - FODC.NUV_YCEN)*(FODC.NUV_YSCL);

      /* Rescale and adjust post May 4th, 2010 data only. */
      if (eclipse > 37460) {

        /* Rescale x,y values so that STIMs match in new data... -tab 01jul2010 */
        x = (FODC.Mx * x) + FODC.Bx;
        y = (FODC.My * y) + FODC.By;

        /* Adjust y value (subtract) (NUV only). -tab 24jun2010 */
        yac = rtaph_yac( FODC, ya, yb, yamc );                 /* Y correction in microns. */
        y = y - yac;

        /* Adjust y value again (add) (NUV only). -tab 24sep2010 */
        yac = rtaph_yac2( WIG2, WLK2, CLK2, q, xb, yb, ya, y, aspum );  /* Y correction in microns. */
        y = y + yac;

      }

      // Convert detector position from microns to arcsecs
      x_as = x*cal->aspum;
      y_as = y*cal->aspum;

      // apply wiggle correction

      // Scale to the standard calibration image size, which is 450 10 arcsec
      // pixels centered in a 480x480 image.
      fptrx  = x_as/10.0 + 240;
      fptry  = y_as/10.0 + 240;
      ptrx = fptrx;
      ptry = fptry;

      if(ptrx<0 || ptrx>479 || ptry<0 || ptry>479) {

        // Photon out of range

        if(*errcnt<maxerrs && ! stim) {
          fprintf(stderr,"=== RTAph: NUV Wiggle X/Y out of range. i=%d.\n"
                  "===     x,yraw=%f,%f; x,y=%f,%f; ptrx,y=%d(%f\"),%d(%f\")\n"
                  "===     x,yb=%d,%d, x,yamc=%d,%d, xa,q=%d,%d\n"
                  "===     x,yclk=%f,%f; x,ycen=%f,%f; x,yscl=%f,%f; um=%f\n",
                  i,
                  xraw,yraw,x,y,ptrx,x_as,ptry,y_as,
                  xb,yb,xamc,yamc,xa,q,
                  cal->xclk_nuv,cal->yclk_nuv,
                  cal->xcen_nuv,cal->ycen_nuv,
                  cal->xscl_nuv,cal->yscl_nuv,
                  cal->aspum
                );
        }

        if(! stim) ++(*errcnt);

        // Carry on uncorrected, but set skip flag.

        flags |= PHSKIP;
        flags |= stim?PHSTIM:PHRANGE;

        xpp    = xp    = x;
        ypp    = yp    = y;
        xpp_as = xp_as = x_as;
        ypp_as = yp_as = y_as;

        xsc    = y_as*10.;
        ysc    = x_as*10.;

        goto BAD_PHOT_SKIP;

      } // Photon out of range

      // Photon in range; do correction

      blt = fptrx-ptrx;
      blu = fptry-ptry;

      // Correction is in tenths of arcseconds
      wgx=(1-blt)*cal->wig_nuv_x[ptrx+xa*480] +
          (blt)*  cal->wig_nuv_x[(ptrx+1)+xa*480];
      wgy=(1-blu)*cal->wig_nuv_y[ptry+ya*480] +
          (blu)*  cal->wig_nuv_y[(ptry+1)+ya*480];

      xdig   = x + wgx/(10.0*cal->aspum);
      ydig   = y + wgy/(10.0*cal->aspum);

      xdig_as = xdig*cal->aspum;
      ydig_as = ydig*cal->aspum;

      // apply walk correction

      // Scale to the standard calibration image size, which is 450 10 arcsec
      // pixels centered in a 480x480 image.
      fptrx  = xdig_as/10.0 + 240;
      fptry  = ydig_as/10.0 + 240;
      ptrx = fptrx;
      ptry = fptry;

      if(ptrx<0 || ptrx>479 || ptry<0 || ptry>479) {

        // Photon out of range

        if(*errcnt<maxerrs && ! stim) {
          fprintf(stderr,"=== RTAph: NUV Walk X/Y out of range. i=%d.\n"
                  "===     x,yraw=%f,%f; x,y=%f,%f; ptrx,y=%d(%f\"),%d(%f\")\n"
                  "===     x,yb=%d,%d, x,yamc=%d,%d, xa,q=%d,%d\n"
                  "===     x,yclk=%f,%f; x,ycen=%f,%f; x,yscl=%f,%f; um=%f\n",
                  i,
                  xraw,yraw,x,y,ptrx,x_as,ptry,y_as,
                  xb,yb,xamc,yamc,xa,q,
                  cal->xclk_nuv,cal->yclk_nuv,
                  cal->xcen_nuv,cal->ycen_nuv,
                  cal->xscl_nuv,cal->yscl_nuv,
                  cal->aspum
                );
        }

        if(! stim) ++(*errcnt);

        // Carry on uncorrected, but set skip flag.

        flags |= PHSKIP;
        flags |= stim?PHSTIM:PHRANGE;

        xpp    = xp    = x;
        ypp    = yp    = y;
        xpp_as = xp_as = x_as;
        ypp_as = yp_as = y_as;

        xsc    = y_as*10.;
        ysc    = x_as*10.;

        goto BAD_PHOT_SKIP;

      } // Photon out of range

      // Photon in range; do correction

      blt = fptrx-ptrx;
      blu = fptry-ptry;

      if ((cal->wlk_nuv_x[ptrx+ptry*480+q*480*480] == -999) ||
         (cal->wlk_nuv_x[(ptrx+1)+ptry*480+q*480*480] == -999) ||
         (cal->wlk_nuv_x[ptrx+(ptry+1)*480+q*480*480]  == -999) ||
         (cal->wlk_nuv_x[ptrx+1+(ptry+1)*480+q*480*480] == -999) ||
         (cal->wlk_nuv_y[ptrx+ptry*480+q*480*480] == -999) ||
         (cal->wlk_nuv_y[(ptrx+1)+ptry*480+q*480*480] == -999) ||
         (cal->wlk_nuv_y[ptrx+(ptry+1)*480+q*480*480] == -999) ||
         (cal->wlk_nuv_y[ptrx+1+(ptry+1)*480+q*480*480]== -999)) {

        // if walk table flag is encountered then flag photon for
        // skipping and set walk correction to 0.0

        flags |= PHSKIP;
        flags |= PHBADWALK;

        wkx = 0.0;
        wky = 0.0;

      }
      else {

        // walk correction should be O.K., with interpolates on all sides.

        // Correction is in microns
        wkx=(1-blt)*(1-blu)*cal->wlk_nuv_x[ptrx+ptry*480+q*480*480] +
            (blt)*(1-blu)*  cal->wlk_nuv_x[(ptrx+1)+ptry*480+q*480*480] +
            (1-blt)*(blu)*  cal->wlk_nuv_x[ptrx+(ptry+1)*480+q*480*480] +
            (blt)*(blu)*    cal->wlk_nuv_x[ptrx+1+(ptry+1)*480+q*480*480];
        wky=(1-blt)*(1-blu)*cal->wlk_nuv_y[ptrx+ptry*480+q*480*480] +
            (blt)*(1-blu)*  cal->wlk_nuv_y[(ptrx+1)+ptry*480+q*480*480] +
            (1-blt)*(blu)*  cal->wlk_nuv_y[ptrx+(ptry+1)*480+q*480*480] +
            (blt)*(blu)*    cal->wlk_nuv_y[ptrx+1+(ptry+1)*480+q*480*480];
      }

      xp    = xdig - wkx;
      yp    = ydig - wky;

      xp_as = xp*cal->aspum;
      yp_as = yp*cal->aspum;



      // Linearity correction

      fptrx = xp_as/10.0 + 240;
      fptry = yp_as/10.0 + 240;
      ptrx  = fptrx;
      ptry  = fptry;

      if(ptrx<0 || ptrx>479 || ptry<0 || ptry>479) {

        // Photon out of range

        if(*errcnt<maxerrs && ! stim) {
          fprintf(stderr,"=== RTAph: NUV Linearity X/Y out of range. i=%d.\n"
                  "===     x,yraw=%f,%f; x,yp=%f,%f; ptrx,y=%d(%f\"),%d(%f\")\n"
                  "===     wkx,y=%f,%f\n"
                  "===     x,yb=%d,%d, x,yamc=%d,%d, xa,q=%d,%d\n"
                  "===     x,yclk=%f,%f; x,ycen=%f,%f;x,yscl=%f,%f; um=%f\n",
                  i,
                  xraw,yraw,xp,yp,ptrx,xp_as,ptry,yp_as,
                  wkx,wky,
                  xb,yb,xamc,yamc,xa,q,
                  cal->xclk_nuv,cal->yclk_nuv,
                  cal->xcen_nuv,cal->ycen_nuv,
                  cal->xscl_nuv,cal->yscl_nuv,
                  cal->aspum
                  );
        }

        if(! stim) ++(*errcnt);

        flags |= PHSKIP;
        flags |= stim?PHSTIM:PHRANGE;

        xpp    = xp;
        ypp    = yp;
        xpp_as = xp_as;
        ypp_as = yp_as;

        xsc    = yp_as*10.;
        ysc    = xp_as*10.;

        goto BAD_PHOT_SKIP;

      } // Photon out of range

      // apply linearity correction

      blt = fptrx-ptrx;
      blu = fptry-ptry;

      // these units are in asecs
      dx=(1-blt)*(1-blu)*cal->lin_nuv_x[ptrx+ptry*480] +
         (blt)*(1-blu)*  cal->lin_nuv_x[(ptrx+1)+ptry*480] +
         (1-blt)*(blu)*  cal->lin_nuv_x[ptrx+(ptry+1)*480] +
         (blt)*(blu)*    cal->lin_nuv_x[ptrx+1+(ptry+1)*480];
      dy=(1-blt)*(1-blu)*cal->lin_nuv_y[ptrx+ptry*480] +
         (blt)*(1-blu)*  cal->lin_nuv_y[(ptrx+1)+ptry*480] +
         (1-blt)*(blu)*  cal->lin_nuv_y[ptrx+(ptry+1)*480] +
         (blt)*(blu)*    cal->lin_nuv_y[ptrx+1+(ptry+1)*480];

//      xpp_as = xp_as + dx;
//      ypp_as = yp_as + dy;

      // apply linearity correction and NUV-FUV offset and distortion x,y correction
      xpp_as = xp_as + dx + rtaph_FODC_offset_plus_shift( xp_as, yp_as, ph6->t, FODC, band, 'x' );
      ypp_as = yp_as + dy + rtaph_FODC_offset_plus_shift( xp_as, yp_as, ph6->t, FODC, band, 'y' );

      xpp    = xpp_as/cal->aspum;
      ypp    = ypp_as/cal->aspum;

      // If a plate solution has been provided, apply it
      if(nuvpltsol) {
        double tmp;
        tmp = nuvpltsol[0] + nuvpltsol[1]*xpp + nuvpltsol[2]*ypp;
        ypp = nuvpltsol[3] + nuvpltsol[4]*xpp + nuvpltsol[5]*ypp;
        xpp = tmp;
        xpp_as = xpp*cal->aspum;
        ypp_as = ypp*cal->aspum;
      }

      // Now rotate/reflect to get into S/C axes in deci-arcseconds.
      // Note that NUV/FUV transformations are different!
      xsc = ypp_as*10.;
      ysc = xpp_as*10.;
      // Here's the remapping to xi and eta; sky map oreientation when roll=0
      xi  = (xi_xsc *xsc + xi_ysc *ysc);
      eta = (eta_xsc*xsc + eta_ysc*ysc);



    } else {


      // ========= Band 2

      xraw0 = xb*cal->xclk_fuv + xamc ;
      yraw0 = yb*cal->yclk_fuv + yamc ;

      ya = ((int)(((yraw0/(2*cal->yclk_fuv) -
                    xraw0/(2*cal->xclk_fuv)   )+10)*32)+xa) % 32;

      xraw = xraw0+((int)((xa+7)%32)-16)*cal->xslp_fuv;
      yraw = yraw0+((int)((ya+7)%32)-16)*cal->yslp_fuv;

      x = (xraw - cal->xcen_fuv)*(cal->xscl_fuv);
      y = (yraw - cal->ycen_fuv)*(cal->yscl_fuv);

      // See if this is a stim
      //          band v      v scale  window v   v options
      stim = whichstim(2,NULL,0,xraw,yraw,stimdxy,NULL);

      // Convert detector position from microns to arcsecs
      x_as = x*cal->aspum;
      y_as = y*cal->aspum;

      // Scale to the standard calibration image size, which is 450 10 arcsec
      // pixels centered in a 480x480 image.
      fptrx  = x_as/10.0 + 240;
      fptry  = y_as/10.0 + 240;
      ptrx = fptrx;
      ptry = fptry;

      if(ptrx<0 || ptrx>479 || ptry<0 || ptry>479) {

        // Photon out of range

        if(*errcnt<maxerrs && ! stim) {
          fprintf(stderr,"=== RTAph: FUV Wiggle X/Y out of range. i=%d.\n"
                  "===     x,yraw=%f,%f; x,y=%f,%f; ptrx,y=%d(%f\"),%d(%f\")\n"
                  "===     x,yb=%d,%d, x,yamc=%d,%d, xa,q=%d,%d\n"
                  "===     x,yclk=%f,%f; x,ycen=%f,%f; x,yscl=%f,%f; um=%f\n",
                  i,
                  xraw,yraw,x,y,ptrx,x_as,ptry,y_as,
                  xb,yb,xamc,yamc,xa,q,
                  cal->xclk_fuv,cal->yclk_fuv,
                  cal->xcen_fuv,cal->ycen_fuv,
                  cal->xscl_fuv,cal->yscl_fuv,
                  cal->aspum
                );
        }

        if(! stim) ++(*errcnt);

        // Carry on uncorrected, but set skip flag.

        flags |= PHSKIP;
        flags |= stim?PHSTIM:PHRANGE;

        xpp    = xp    = x;
        ypp    = yp    = y;
        xpp_as = xp_as = x_as;
        ypp_as = yp_as = y_as;

        xsc    = y_as*10.;
        ysc    = x_as*10.;

        goto BAD_PHOT_SKIP;

      } // Photon out of range

      // Photon in range; do correction

      blt = fptrx-ptrx;
      blu = fptry-ptry;

      // Correction is in 10ths of arcseconds
      wgx=(1-blt)*cal->wig_fuv_x[ptrx+xa*480] +
          (blt)*  cal->wig_fuv_x[(ptrx+1)+xa*480];
      wgy=(1-blu)*cal->wig_fuv_y[ptry+ya*480] +
          (blu)*  cal->wig_fuv_y[(ptry+1)+ya*480];

      xdig   = x + wgx/(10.0*cal->aspum);
      ydig   = y + wgy/(10.0*cal->aspum);

      xdig_as = xdig*cal->aspum;
      ydig_as = ydig*cal->aspum;

      // apply walk correction

      fptrx  = xdig_as/10.0 + 240;
      fptry  = ydig_as/10.0 + 240;
      ptrx = fptrx;
      ptry = fptry;

      if(ptrx<0 || ptrx>479 || ptry<0 || ptry>479) {

         // Photon out of range

        if(*errcnt<maxerrs && ! stim) {
          fprintf(stderr,"=== RTAph: FUV Walk X/Y out of range. i=%d.\n"
                  "===     x,yraw=%f,%f; x,y=%f,%f; ptrx,y=%d(%f\"),%d(%f\")\n"
                  "===     x,yb=%d,%d, x,yamc=%d,%d, xa,q=%d,%d\n"
                  "===     x,yclk=%f,%f; x,ycen=%f,%f; x,yscl=%f,%f; um=%f\n",
                  i,
                  xraw,yraw,x,y,ptrx,x_as,ptry,y_as,
                  xb,yb,xamc,yamc,xa,q,
                  cal->xclk_fuv,cal->yclk_fuv,
                  cal->xcen_fuv,cal->ycen_fuv,
                  cal->xscl_fuv,cal->yscl_fuv,
                  cal->aspum
                  );
        }

        if(! stim) ++(*errcnt);

        flags |= PHSKIP;
        flags |= stim?PHSTIM:PHRANGE;

        xpp    = xp    = x;
        ypp    = yp    = y;
        xpp_as = xp_as = x_as;
        ypp_as = yp_as = y_as;

        xsc = -y_as*10.;
        ysc = -x_as*10.;

        goto BAD_PHOT_SKIP;
      } // Photon out of range

      blt = fptrx-ptrx;
      blu = fptry-ptry;


      if ((cal->wlk_fuv_x[ptrx+ptry*480+q*480*480] == -999) ||
         (cal->wlk_fuv_x[(ptrx+1)+ptry*480+q*480*480] == -999) ||
         (cal->wlk_fuv_x[ptrx+(ptry+1)*480+q*480*480]  == -999) ||
         (cal->wlk_fuv_x[ptrx+1+(ptry+1)*480+q*480*480] == -999) ||
         (cal->wlk_fuv_y[ptrx+ptry*480+q*480*480] == -999) ||
         (cal->wlk_fuv_y[(ptrx+1)+ptry*480+q*480*480] == -999) ||
         (cal->wlk_fuv_y[ptrx+(ptry+1)*480+q*480*480] == -999) ||
         (cal->wlk_fuv_y[ptrx+1+(ptry+1)*480+q*480*480]== -999)) {

        // if walk table flag is encountered then flag photon for
        // skipping and set walk correction to 0.0

        flags |= PHSKIP;
        flags |= PHBADWALK;

        wkx = 0.0;
        wky = 0.0;

      } else {

        // walk correction should be O.K., with interpolates on all sides.


      wkx=(1-blt)*(1-blu)*cal->wlk_fuv_x[ptrx+ptry*480+q*480*480] +
          (blt)*(1-blu)*  cal->wlk_fuv_x[(ptrx+1)+ptry*480+q*480*480] +
          (1-blt)*(blu)*  cal->wlk_fuv_x[ptrx+(ptry+1)*480+q*480*480] +
          (blt)*(blu)*    cal->wlk_fuv_x[ptrx+1+(ptry+1)*480+q*480*480];
      wky=(1-blt)*(1-blu)*cal->wlk_fuv_y[ptrx+ptry*480+q*480*480] +
          (blt)*(1-blu)*  cal->wlk_fuv_y[(ptrx+1)+ptry*480+q*480*480] +
          (1-blt)*(blu)*  cal->wlk_fuv_y[ptrx+(ptry+1)*480+q*480*480] +
            (blt)*(blu)*  cal->wlk_fuv_y[ptrx+1+(ptry+1)*480+q*480*480];

      }

      xp    = xdig - wkx;
      yp    = ydig - wky;
      xp_as = xp*cal->aspum;
      yp_as = yp*cal->aspum;


      // Linearity correction

      fptrx = xp_as/10.0 + 240;
      fptry = yp_as/10.0 + 240;
      ptrx  = fptrx;
      ptry  = fptry;

      if(ptrx<0 || ptrx>479 || ptry<0 || ptry>479) {

        // Photon out of range

        if(*errcnt<maxerrs && ! stim) {
          fprintf(stderr,"=== RTAph: FUV Linearity X/Y out of range. i=%d.\n"
                  "===     x,yraw=%f,%f; x,yp=%f,%f; ptrx,y=%d(%f\"),%d(%f\")\n"
                  "===     wkx,y=%f,%f\n"
                  "===     x,yb=%d,%d, x,yamc=%d,%d, xa,q=%d,%d\n"
                  "===     x,yclk=%f,%f; x,ycen=%f,%f; x,yscl=%f,%f; um=%f\n",
                  i,
                  xraw,yraw,x,y,ptrx,xp_as,ptry,yp_as,
                  wkx,wky,
                  xb,yb,xamc,yamc,xa,q,
                  cal->xclk_fuv,cal->yclk_fuv,
                  cal->xcen_fuv,cal->ycen_fuv,
                  cal->xscl_fuv,cal->yscl_fuv,
                  cal->aspum
                  );
        }

        if(! stim) ++(*errcnt);

        flags |= PHSKIP;
        flags |= stim?PHSTIM:PHRANGE;

        xpp    = xp;
        ypp    = yp;
        xpp_as = xp_as;
        ypp_as = yp_as;
        xsc = -yp_as*10.;
        ysc = -xp_as*10.;

        goto BAD_PHOT_SKIP;
      } // Photon out of range

      blt = fptrx-ptrx;
      blu = fptry-ptry;

      // these units are in microns
      dx=(1-blt)*(1-blu)*cal->lin_fuv_x[ptrx+ptry*480] +
         (blt)*(1-blu)*  cal->lin_fuv_x[(ptrx+1)+ptry*480] +
         (1-blt)*(blu)*  cal->lin_fuv_x[ptrx+(ptry+1)*480] +
         (blt)*(blu)*    cal->lin_fuv_x[ptrx+1+(ptry+1)*480];
      dy=(1-blt)*(1-blu)*cal->lin_fuv_y[ptrx+ptry*480] +
         (blt)*(1-blu)*  cal->lin_fuv_y[(ptrx+1)+ptry*480] +
         (1-blt)*(blu)*  cal->lin_fuv_y[ptrx+(ptry+1)*480] +
         (blt)*(blu)*    cal->lin_fuv_y[ptrx+1+(ptry+1)*480];

//      xpp_as = xp_as + dx;
//      ypp_as = yp_as + dy;

      // apply linearity correction and NUV-FUV offset and distortion x,y correction
      xpp_as = xp_as + dx + rtaph_FODC_offset_plus_shift( xp_as, yp_as, ph6->t, FODC, band, 'x' );
      ypp_as = yp_as + dy + rtaph_FODC_offset_plus_shift( xp_as, yp_as, ph6->t, FODC, band, 'y' );

      xpp    = xpp_as/cal->aspum;
      ypp    = ypp_as/cal->aspum;

      // If a plate solution has been provided, apply it
      if(fuvpltsol) {
        double tmp;
        tmp = fuvpltsol[0] + fuvpltsol[1]*xpp + fuvpltsol[2]*ypp;
        ypp = fuvpltsol[3] + fuvpltsol[4]*xpp + fuvpltsol[5]*ypp;
        xpp = tmp;
        xpp_as = xpp*cal->aspum;
        ypp_as = ypp*cal->aspum;
      }

      // Now rotate/reflect to get into S/C axes in deci-arcseconds.
      // Note that NUV/FUV transformations are different!
      xsc = -ypp_as*10.;
      ysc = -xpp_as*10.;
      // Here's the remapping to xi and eta; sky map orientation when roll=0
      xi  = (xi_xsc *xsc + xi_ysc *ysc);
      eta = (eta_xsc*xsc + eta_ysc*ysc);

    } // Band 2

  BAD_PHOT_SKIP:

    // Populate extended photon record with everything but equatorial position.

    phots_x->t      = ph6->t;
    // Make full detector width = +/- 1.0 scaled by PHSCL(x) = 21750
    phots_x->x      = ROUND((rawxy?x_as:xp_as)/(3600*detsize/2.)*PHSCL(x));
    phots_x->y      = ROUND((rawxy?y_as:yp_as)/(3600*detsize/2.)*PHSCL(x));
    phots_x->q      = q;
    phots_x->val    = xa;
    phots_x->thetax = ROUND(xsc);
    phots_x->thetay = ROUND(ysc);

    if(cal && cal->mask) {
      double xi_det,eta_det;
      // Now get a xi and eta appropriate for masking. These must come
      // from the detector coordinates, so have different axis
      // transformations depending on band.  Assume either
      // values depending on selected options
      if(nlcmask) {
        // NLC-corrected photon position
        xi_det   = xi/36000;
        eta_det  = eta/36000;
      }
      else {
        // Use pre-NLC positions in either raw or walk/wiggle-corrected
        // versions
        double detx = rawxy ? x_as:xp_as;
        double dety = rawxy ? y_as:yp_as;
        if(band == 1) {
          // NUV: X_det -> Y_s/c,  Y_det -> X_s/c
          // (In degrees)
          xi_det  = xi_xsc*(dety/3600)  + xi_ysc*(detx/3600);
          eta_det = eta_xsc*(dety/3600) + eta_ysc*(detx/3600);
        }
        else {
          // FUV: X_det -> -Y_s/c, Y_det-> -X_s/c
          // (In degrees)
          xi_det  = xi_xsc*(-dety/3600)  + xi_ysc*(-detx/3600);
          eta_det = eta_xsc*(-dety/3600) + eta_ysc*(-detx/3600);
        }
      }
      /* Get photon col/row in the mask image */
      col = ( xi_det/(detsize/2.)*maskfill + 1.)/2. * cal->maskdef.npixx;
      row = (eta_det/(detsize/2.)*maskfill + 1.)/2. * cal->maskdef.npixy;
      pix = XYTOPIX(col,row,cal->maskdef.npixx); /* Pixel number in image */
      //printf("col,row,pix=%d,%d,%d\n",col,row,pix);
      if(col < cal->maskdef.npixx && col >= 0 &&
         row < cal->maskdef.npixy && row >= 0) {
        if(*((float *)cal->mask+pix) == 0) {
          ++nmasked;
          masked = 1;
          flags |= PHMASKED;
          flags |= PHSKIP;
        }
      }
    }


    // Rectification complete.


    // Assign equatorial positions from aspect data.

    // (Note that skipped photons (such as stims) will end up with bogus
    //  sky positions, but we shouldn't care, since it's to be skipped.)

    // Photon time
    t = phots_x->t;

    // Look up the next (or 1st) aspect record closest in time but smaller
    if(*nx==0 || t > (asp+1)->t) {

      //  printf("Next record: %d asptx %f, aspty %f daspra %f daspdec %f\n",
      //         ix,asptx[ix],aspty[ix],
      //         ((asp->raz)-ra0)*3600,((asp->decz)-dec0)*3600);

      // Find next pair of Aspect records to bracket the photon time
      if(*nx == 0) ix = 0;
      while(ix < nasp && asp->t >= 0 && t >= asp->t) { ++asp; ++ix; }
      if(ix>0) {
        --asp; --ix;
      }
      if(ix < 0 || ix >= nasp || t+.1 < asp->t || t-.1 > (asp+1)->t) {
        // The desired time is missing from the aspect file.
        //if(! aspskip) fprintf(stderr,"=== RTAph: "
        //                      "Photon time %f not bracketed; skip. "
        //                      "ix=%d, nasp=%d, asptime=%f\n",
        //                      t,ix,nasp,asp->t);
        aspskip = 1;
        //exit(-1);
      }
      else {
        //double tmpra,tmpdec;

        aspskip = 0;

        // Update sky->focal plane conversion for this time slice.
        // Use dec-arcsecond pixels with pixel (0,0) corresponding to
        // the boresite centered on the image.
        //gnomonicinit(0., asp->raz, -1.0/36000.0,
        //             0., asp->decz, 1.0/36000.0,
        //             -asp->twistz + droll,
        //             &cel,&prj,&lin);

        if(verbose > 1) {
          printf("### Aspect now at index %d, t=%f/%f, raz/decz=%f/%f, "
                 "flags=%d\n",
                 ix, asp->t,t, asp->raz,asp->decz,asp->flags);
        }
      
        // Doing this once seems to be necessary after the init above. I think
        // two inits in a row must cause corruption somewhere.
        // Other than avoiding that problem, this is a NO-OP.
        //gnomonicrev(&tmpra,&tmpdec,asp->raz,asp->decz,&cel,&prj,&lin);
        
        naspslice = 0;
      }
    }


    // Perform aspect correction (i.e. dedither)

    // Check that the photon is bracketed by valid aspect entries.
    // Need good entries on both sides so we know the interpolation is OK.
    if(! aspskip                             &&  // Found an aspect record 
       ! (asp->flags&(ASPSKIP|asphvlow))     &&  // Current asp OK
       ! ((asp+1)->flags&(ASPSKIP|asphvlow)) &&  // Next asp OK
       (t - asp->t) <= dt                        // In range of an aspect record
       ) { 
      // Currently assume twist is negligeably variable over one time slice.
      // More complex code will be required if not.
      // asptx/y are the rate of motion arrays which are aligned
      // x->RA, y->Dec when twist=0.
      if(hasasp) {
        dxi  =  asptx[ix]*(t - asp->t)/dt;
        deta =  aspty[ix]*(t - asp->t)/dt;
      }
      else {
        dxi = deta = 0;
      }
      xi  += dxi;
      eta += deta;

      // Here's the conversion to equatorial coordinates.
      //gnomonicrev(&phra,&phdec,xi,eta,&cel,&prj,&lin);
      // Use the faster, simpler version.
      gnomrev_simple(xi,eta,
                     asp->raz,asp->decz, -asp->twistz + droll,
                     1.0/36000.0, 0.0, 
                     &phra,&phdec);

      // Some debug printout. Nothing in this block really matters.
      if (naspslice==0 && verbose) { 
        double tra0,tdec0,traul,tdecul;
        //printf("\nxpp/xsc=%0.1f/%0.1f, xi,dxi=%0.1f/%0.1f,  "
        //       "ypp/ysc=%0.1f/%0.1f, eta,deta=%0.1f/%0.1f, "
        //       "asptx/y=%0.1f/%0.1f,  delt/dt=%0.1f/%0.1f\n",
        //       xpp,xsc, xi,dxi,
        //       ypp,ysc, eta,deta,
        //       hasasp?asptx[ix]:0,hasasp?aspty[ix]:0,t-asp->t,dt);
        //gnomonicrev(&tra0, &tdec0,  0.,0.,&cel,&prj,&lin);
        //gnomonicrev(&traul,&tdecul,-1.,1.,&cel,&prj,&lin);
        gnomrev_simple(0.,0.,
                       asp->raz,asp->decz, -asp->twistz + droll,
                       1.0/36000.0, 0.0,
                       &tra0, &tdec0);
        gnomrev_simple(-1.,1.,
                       asp->raz,asp->decz, -asp->twistz + droll,
                       1.0/36000.0, 0.0,
                       &traul,&tdecul);
        
        //printf("phra/dec = %f/%f, raz/decz   = %f/%f, roll=%f\n"
        //       "Ra0/dec0 = %f/%f, RaUL/decUL = %f/%f\n"
        //       //"cdelt = %g/%g, crpix=%f/%f, crval=%f/%f\n"
        //       ,phra,phdec,asp->raz,asp->decz, asp->twistz,
        //       tra0,tdec0,traul,tdecul
        //       //,lin.cdelt[0],lin.cdelt[1],lin.crpix[0],lin.crpix[1],
        //       //cel.ref[0],cel.ref[1]
        //       );
        // printf("\tdra,ddec=%f,%f\n",
        //        (tra-lra)*cos(tdec/R2D)*36000,(tdec-ldec)*36000);
        lt = t; lra = tra; ldec = tdec;
      }

      // Fill in the final components of the extended photon record;
      // the equatorial position.
      sph2vec(phra/R2D,phdec/R2D,&(vec.X));
      phots_x->X = vec.X*PHSCL(X);
      phots_x->Y = vec.Y*PHSCL(Y);
      phots_x->Z = vec.Z*PHSCL(Z);

    } // Do aspect correction

    else {
      // Deal with photons with no acceptable entry in the aspect file.
      // Skip aspect correction.

      phots_x->X = -1*PHSCL(X); // (will resolve to the south eq. pole if not
      phots_x->Y = -1*PHSCL(Y); //  tested for and rejected as not normalizing
      phots_x->Z = -1*PHSCL(Z); //  to unity.)
      flags |= PHBADASP;
      flags |= PHSKIP; // Skip photons w/ no acceptable entry in the aspect file

    }

    // Set photon flags
    // Carry along raw6 flags, except format specifier bits
    phots_x->flags = (ph6->flags&~PHFMTBITS);
    // Add new flags
    phots_x->flags |= flags;
    // Add format specifier bits
    PHSETFMT(phots_x,PHFMTX);

    // COunt unskipped photons, i.e. valid photons
    if(! (phots_x->flags&PHSKIP)) ++nvalid;
    // COunt everything else
    ++phots_x;
    ++*nx;
    ++naspslice;

  } // for i

  //gnomonicfree(&lin); // replaced by gnomfwd/rev_simple calls

  if(cal && cal->mask) {
    if(! quiet) printf("Masked %d photons.\n",nmasked);
  }

  if(*errcnt > 0) {
    //fprintf(stderr,"=== RTAph: Have seen %d out of range photons.\n",
    //        *errcnt);
  }

  if(! hasasp && asp0 != NULL) { free_safe(asp0); }

  free(FODC.dcdx); FODC.dcdx = NULL;
  free(FODC.dcdy); FODC.dcdy = NULL;
  if (eclipse > 37460) {
    free(WIG2); WIG2 = NULL;
    free(WLK2); WLK2 = NULL;
    free(CLK2); CLK2 = NULL;
  }

  return nvalid;

}

// Compute boresite dither rates and store for later use.
int asprates(Aspect * asp, int nasp, 
             double * asptx, double * aspty, NamedList *opts)
{
  Aspect *aspt;
  int i; //,verbose;
  double dx,dy,droll; //,t,dt,ra,dec,aspra,aspdec;
  //double ra0,dec0;

  //struct { double X,Y,Z; } vec, last;
  //struct celprm cel;
  //struct prjprm prj;
  //struct linprm lin;

  droll   = nm_getn(opts,"droll",  0);

  //gnomonicsetup(&cel,&prj,&lin); // replaced by gnomfwd/rev_simple calls

  for (i=0; i<nasp; i++, asp++){

    if(i > 0) {
      // Sky focal plane center position for the current time slice
      // on the previous focal plane position. So asptx/y approximatly contain
      // asptx,aspty(t) = x,y(t+1) - x,y(t).

      //gnomonicfwd(asp->raz,asp->decz,&dx,&dy,&cel,&prj,&lin);
      gnomfwd_simple(asp->raz, asp->decz,
                     aspt->raz,aspt->decz, -aspt->twistz + droll,
                     1.0/36000.0, 0.,
                     &dx,&dy);

      // Record essentially the rate info; relative x/y motion vs. time
      // for each time slice
      *asptx++ = dx;
      *aspty++ = dy;
    }

    // Set transformation to this focal plane position/orientation.
    //gnomonicinit(0.,asp->raz,-1.0/36000.0,
    //             0.,asp->decz,1.0/36000.0,
    //             -asp->twistz + droll,
    //             &cel,&prj,&lin);

    aspt = asp;
  }

  //gnomonicfree(&lin); // replaced by gnomfwd/rev_simple calls

  return 0;
}


