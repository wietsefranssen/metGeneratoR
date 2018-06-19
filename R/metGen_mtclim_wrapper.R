mtclim_wrapper <- function(){
  
  # /******************************************************************************
  #   mtclim_wrapper: interface between VIC and MTCLIM.
  # ******************************************************************************/
  
  #     /* initialize the mtclim data structures */ 
  #       mtclim_init(have_dewpt, have_shortwave, elevation, slope, aspect, ehoriz, whoriz,
  #                   annual_prcp, lat, Ndays, dmy, prec,
  #                   tmax, tmin, vp, hourlyrad, tiny_radfract, &ctrl, &p,
  #                   &mtclim_data);  
  
  #     /* calculate daily air temperatures */
  #       if (calc_tair(&ctrl, &p, &mtclim_data)) {
  #         nrerror("Error in calc_tair()... exiting\n");
  #       }
  calc_tair(ctrl, p, mtclim_data)
  
}
# void mtclim_wrapper(int have_dewpt, int have_shortwave, double hour_offset,
#                     double elevation, double slope, double aspect,
#                     double ehoriz, double whoriz,
#                     double annual_prcp, double lat, 
#                     int Ndays, dmy_struct *dmy, 
#                     double *prec, double *tmax, double *tmin, double *tskc,
#                     double *vp, double *hourlyrad, double *fdir) 
# /******************************************************************************
#   mtclim_wrapper: interface between VIC and MTCLIM.
# 
# Modifications:
#   2012-Feb-16 Cleaned up commented code.					TJB
# ******************************************************************************/
# {
#   control_struct ctrl;
#   parameter_struct p;
#   data_struct mtclim_data;
#   double **tiny_radfract;
#   int i;
#   
#   /* allocate space for the tiny_radfract array */
#     tiny_radfract = (double **) calloc(366, sizeof(double*));
#     if (tiny_radfract == NULL) {
#       nrerror("Memory allocation error in mtclim_init() ...\n");
#     }
#     for (i=0; i<366; i++) {
#       tiny_radfract[i] = (double *) calloc(86400, sizeof(double));
#       if (tiny_radfract[i] == NULL) {
#         nrerror("Memory allocation error in mtclim_init() ...\n");
#       }
#     }
#     
#     /* initialize the mtclim data structures */ 
#       mtclim_init(have_dewpt, have_shortwave, elevation, slope, aspect, ehoriz, whoriz,
#                   annual_prcp, lat, Ndays, dmy, prec,
#                   tmax, tmin, vp, hourlyrad, tiny_radfract, &ctrl, &p,
#                   &mtclim_data);  
#     
#     /* calculate daily air temperatures */
#       if (calc_tair(&ctrl, &p, &mtclim_data)) {
#         nrerror("Error in calc_tair()... exiting\n");
#       }
#     
#     /* calculate daily precipitation */
#       if (calc_prcp(&ctrl, &p, &mtclim_data)) {
#         nrerror("Error in calc_prcp()... exiting\n");
#       }
#     
#     /* calculate daily snowpack using simple model (this is only for radiation correction, *not* the same as the VIC snowpack estimate) */
#       if (snowpack(&ctrl, &p, &mtclim_data)) {
#         nrerror("Error in snowpack()... exiting\n");
#       }
#     
#     /* calculate srad and humidity with iterative algorithm */
#       if (calc_srad_humidity_iterative(&ctrl, &p, &mtclim_data, tiny_radfract)) { 
#         nrerror("Error in calc_srad_humidity_iterative()... exiting\n");
#       }
#     
#     /* translate the mtclim structures back to the VIC data structures */
#       mtclim_to_vic(hour_offset, Ndays,
#                     dmy, tiny_radfract, &ctrl,&mtclim_data, tskc, vp,
#                     hourlyrad, fdir);
#     
#     /* clean up */
#       if (data_free(&ctrl, &mtclim_data)) {
#         nrerror("Error in data_free()... exiting\n");
#       }
#     for (i=0; i<366; i++) {
#       free(tiny_radfract[i]);
#     }
#     free(tiny_radfract);
# }
