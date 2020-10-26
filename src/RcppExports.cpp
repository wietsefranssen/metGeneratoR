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
// sh2vp
NumericVector sh2vp(NumericVector q, NumericVector p);
RcppExport SEXP _metGeneratoR_sh2vp(SEXP qSEXP, SEXP pSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type q(qSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type p(pSEXP);
    rcpp_result_gen = Rcpp::wrap(sh2vp(q, p));
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
// calc_tas_cr
NumericVector calc_tas_cr(NumericVector rad_fract_map, NumericVector tmin_map, NumericVector tmax_map, int yday, int nrec, NumericVector xybox);
RcppExport SEXP _metGeneratoR_calc_tas_cr(SEXP rad_fract_mapSEXP, SEXP tmin_mapSEXP, SEXP tmax_mapSEXP, SEXP ydaySEXP, SEXP nrecSEXP, SEXP xyboxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type rad_fract_map(rad_fract_mapSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type tmin_map(tmin_mapSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type tmax_map(tmax_mapSEXP);
    Rcpp::traits::input_parameter< int >::type yday(ydaySEXP);
    Rcpp::traits::input_parameter< int >::type nrec(nrecSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type xybox(xyboxSEXP);
    rcpp_result_gen = Rcpp::wrap(calc_tas_cr(rad_fract_map, tmin_map, tmax_map, yday, nrec, xybox));
    return rcpp_result_gen;
END_RCPP
}
// rad_map_final_cr
NumericVector rad_map_final_cr(int nrec, int yday, double gmt_float, NumericVector xybox, NumericVector lats, NumericVector lons, float gmt_offset);
RcppExport SEXP _metGeneratoR_rad_map_final_cr(SEXP nrecSEXP, SEXP ydaySEXP, SEXP gmt_floatSEXP, SEXP xyboxSEXP, SEXP latsSEXP, SEXP lonsSEXP, SEXP gmt_offsetSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type nrec(nrecSEXP);
    Rcpp::traits::input_parameter< int >::type yday(ydaySEXP);
    Rcpp::traits::input_parameter< double >::type gmt_float(gmt_floatSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type xybox(xyboxSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type lats(latsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type lons(lonsSEXP);
    Rcpp::traits::input_parameter< float >::type gmt_offset(gmt_offsetSEXP);
    rcpp_result_gen = Rcpp::wrap(rad_map_final_cr(nrec, yday, gmt_float, xybox, lats, lons, gmt_offset));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_metGeneratoR_set_vp_cr", (DL_FUNC) &_metGeneratoR_set_vp_cr, 5},
    {"_metGeneratoR_sh2vp", (DL_FUNC) &_metGeneratoR_sh2vp, 2},
    {"_metGeneratoR_set_min_max_hour_cr", (DL_FUNC) &_metGeneratoR_set_min_max_hour_cr, 2},
    {"_metGeneratoR_calc_tas_cr", (DL_FUNC) &_metGeneratoR_calc_tas_cr, 6},
    {"_metGeneratoR_rad_map_final_cr", (DL_FUNC) &_metGeneratoR_rad_map_final_cr, 7},
    {NULL, NULL, 0}
};

RcppExport void R_init_metGeneratoR(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
