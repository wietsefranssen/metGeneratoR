// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// set_vp_cr
NumericVector set_vp_cr(NumericVector tair_r, NumericVector relhum_r, int nx, int ny, int nrec);
RcppExport SEXP _metGeneratoR_set_vp_cr(SEXP tair_rSEXP, SEXP relhum_rSEXP, SEXP nxSEXP, SEXP nySEXP, SEXP nrecSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type tair_r(tair_rSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type relhum_r(relhum_rSEXP);
    Rcpp::traits::input_parameter< int >::type nx(nxSEXP);
    Rcpp::traits::input_parameter< int >::type ny(nySEXP);
    Rcpp::traits::input_parameter< int >::type nrec(nrecSEXP);
    rcpp_result_gen = Rcpp::wrap(set_vp_cr(tair_r, relhum_r, nx, ny, nrec));
    return rcpp_result_gen;
END_RCPP
}
// set_min_max_hour_cr
NumericVector set_min_max_hour_cr(NumericVector radfrac, int nx);
RcppExport SEXP _metGeneratoR_set_min_max_hour_cr(SEXP radfracSEXP, SEXP nxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type radfrac(radfracSEXP);
    Rcpp::traits::input_parameter< int >::type nx(nxSEXP);
    rcpp_result_gen = Rcpp::wrap(set_min_max_hour_cr(radfrac, nx));
    return rcpp_result_gen;
END_RCPP
}
// set_max_min_lonlat_cr
NumericVector set_max_min_lonlat_cr(NumericVector tmin_map, NumericVector tmax_map, int yday, int nrec);
RcppExport SEXP _metGeneratoR_set_max_min_lonlat_cr(SEXP tmin_mapSEXP, SEXP tmax_mapSEXP, SEXP ydaySEXP, SEXP nrecSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type tmin_map(tmin_mapSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type tmax_map(tmax_mapSEXP);
    Rcpp::traits::input_parameter< int >::type yday(ydaySEXP);
    Rcpp::traits::input_parameter< int >::type nrec(nrecSEXP);
    rcpp_result_gen = Rcpp::wrap(set_max_min_lonlat_cr(tmin_map, tmax_map, yday, nrec));
    return rcpp_result_gen;
END_RCPP
}
// rad_map_final_cr
NumericVector rad_map_final_cr(int nrec, int yday, int nx_parts, double gmt_float);
RcppExport SEXP _metGeneratoR_rad_map_final_cr(SEXP nrecSEXP, SEXP ydaySEXP, SEXP nx_partsSEXP, SEXP gmt_floatSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type nrec(nrecSEXP);
    Rcpp::traits::input_parameter< int >::type yday(ydaySEXP);
    Rcpp::traits::input_parameter< int >::type nx_parts(nx_partsSEXP);
    Rcpp::traits::input_parameter< double >::type gmt_float(gmt_floatSEXP);
    rcpp_result_gen = Rcpp::wrap(rad_map_final_cr(nrec, yday, nx_parts, gmt_float));
    return rcpp_result_gen;
END_RCPP
}
// rad_map_lats_cr
NumericVector rad_map_lats_cr(int nt, int yday);
RcppExport SEXP _metGeneratoR_rad_map_lats_cr(SEXP ntSEXP, SEXP ydaySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type nt(ntSEXP);
    Rcpp::traits::input_parameter< int >::type yday(ydaySEXP);
    rcpp_result_gen = Rcpp::wrap(rad_map_lats_cr(nt, yday));
    return rcpp_result_gen;
END_RCPP
}
// solar_geom_cr
NumericVector solar_geom_cr(float lat, int yday, int timesteps_per_day);
RcppExport SEXP _metGeneratoR_solar_geom_cr(SEXP latSEXP, SEXP ydaySEXP, SEXP timesteps_per_daySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< float >::type lat(latSEXP);
    Rcpp::traits::input_parameter< int >::type yday(ydaySEXP);
    Rcpp::traits::input_parameter< int >::type timesteps_per_day(timesteps_per_daySEXP);
    rcpp_result_gen = Rcpp::wrap(solar_geom_cr(lat, yday, timesteps_per_day));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_metGeneratoR_set_vp_cr", (DL_FUNC) &_metGeneratoR_set_vp_cr, 5},
    {"_metGeneratoR_set_min_max_hour_cr", (DL_FUNC) &_metGeneratoR_set_min_max_hour_cr, 2},
    {"_metGeneratoR_set_max_min_lonlat_cr", (DL_FUNC) &_metGeneratoR_set_max_min_lonlat_cr, 4},
    {"_metGeneratoR_rad_map_final_cr", (DL_FUNC) &_metGeneratoR_rad_map_final_cr, 4},
    {"_metGeneratoR_rad_map_lats_cr", (DL_FUNC) &_metGeneratoR_rad_map_lats_cr, 2},
    {"_metGeneratoR_solar_geom_cr", (DL_FUNC) &_metGeneratoR_solar_geom_cr, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_metGeneratoR(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
