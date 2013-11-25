/**
 @file
 @author Ken Knowles Kenneth.Knowles@NOAA.gov
 @date 29-Aug-2012
 
 @defgroup geosearch
 */

#ifndef GEOSEARCH_H
#define GEOSEARCH_H

#include <stdexcept>
#include <string>
#include <vector>
#include <list>

#include "GeoSearchUtil.h"
#include "GeoLatRange.h"
#include "GeoLonRange.h"
#include "GeoOrbitFinder.h"
#include "GeoOrbitSegment.h"
#include "GeoSpatialExtent.h"

/*
 The following structure encapsulates parameters for overlap detection,
 which vary by orbit quadrant. To simplify the overlap detection algorithm,
 the orbit is divided into quadrants by direction of travel and hemisphere.
 Starting north bound at the equator and in order as the satellite travels:
 Q0 is Ascending, Northern Hemisphere
 Q1 is Descending, Northern Hemisphere
 Q2 is Descending, Southern Hemisphere
 Q3 is Ascending, Southern Hemisphere
 */
typedef struct {
  double nlat;
  double slat;
  GeoSpatialExtent *extent;
  bool in_play;
  SatellDir dir;
} RegionParamsByOrbitQuadrant;

const int NUM_ORBIT_QUADRANTS_GEOSEARCH = 4;

class GeoSearch : public GeoSearchUtil {
  
public:  
  GeoSearch( double nlat, double slat, double wlon, double elon, const char* accepted_dirs_str );
  ~GeoSearch();
  void initDatatypeParams( double xl_step, const GeoOrbitFinder& orbit_ );
  static std::string eqCrossCondition( std::string asc_desc_flag, float nLat, float sLat, float wLon, float eLon);
  static std::vector<GeoOrbitSegment> orbitalPassCoverage( int orbit_number, const GeoOrbitFinder& orbit, 
                                                          double lon0, double orbit_start_minutes, 
                                                          double coverage_start_minutes, double coverage_end_minutes );
  bool findOverlap( std::vector<GeoOrbitSegment> data_set_pass );
  std::string xmlOverlapList();
  bool useGeosearch() const;
  bool isInitialized() const;
  bool isOverlapEmpty() const;
  std::vector<GeoOrbitSegment> overlapList() const;

private:
  // criteria
  GeoLatRange* search_bounds_lat;
  GeoLonRange* search_bounds_lon;
  char accepted_dirs;
  GeoSpatialExtent* search_area;
  GeoOrbitFinder* orbit;
  // results
  std::vector<GeoOrbitSegment> overlap_list;
  // control
  bool use_geosearch;
  bool is_initialized;
  struct {
    double nlat;
    double slat;
    GeoSpatialExtent *extent;
    bool in_play;
    SatellDir dir;
  } search_region_quadrant[NUM_ORBIT_QUADRANTS_GEOSEARCH];

};

inline bool GeoSearch::useGeosearch() const { return use_geosearch; }
inline bool GeoSearch::isInitialized() const { return is_initialized; }
inline bool GeoSearch::isOverlapEmpty() const { return overlap_list.empty(); }
inline std::vector<GeoOrbitSegment> GeoSearch::overlapList() const { return overlap_list; }

/**
 * get orbital parameters common to a data type
 *
 * @param DBObjectSubtype specific type of DBObject we're dealing with
 *
 * @param Ref database values structure
 *
 * @param[out] inclination_deg orbit inclination
 *
 * @param[out] swath_width_m distance from one edge to the other
 *
 * @param[out] period_minutes time for one complete orbit
 *
 * @param[out] xl_step accuracy versus speed tuning
 *
 * @return 0 == normal, 1 == all default values, which is also OK
 *
 */
template< typename DBObjectSubtype >
int getOrbitParams( const DBObjectSubtype& Ref, double& inclination_deg, double& swath_width_m, 
                    double& period_minutes, double& xl_step ) {
  
  int status = 1;
  inclination_deg = geosearch::typical_inclination_deg;
  swath_width_m = geosearch::typical_swath_width_m;
  period_minutes = geosearch::typical_period_minutes;
  xl_step = Ref.getXl_step( geosearch::default_xl_step );
  
  if ( ! Ref.angla_degIsNull() ) { 
    swath_width_m = ( 2 * Ref.getAngla_deg() / 360 ) * ( 2 * M_PI * defaultRadius_m );
    status = 0;
  }
  
  if ( ! Ref.orbital_period_msIsNull() ) { 
    period_minutes = msecToMinutes( Ref.getOrbital_period_ms() );
    status = 0;
  }
  
  return status;
}

/**
 * get coverage parameters for particular row
 *
 * @tparam DBObjectSubtype specific type of DBObject we're dealing with
 *
 * @param Ref database values structure
 *
 * @param gsExec GeoSearch executive (solely to access time conversions)
 *
 * @param[out] lon0 longitude of ascending crossing for reference orbit of pass
 *
 * @param[out] direction "A" Ascending only, "D" Descending only, "B" pass contains Both directions
 *
 * @param[out] orbit_start_minutes time of ascending crossing at lon0
 *
 * @param[out] coverage_start_minutes time of start of data
 *
 * @param[out] coverage_end_minutes time of end of data
 *
 * @return 0 == all parameters, 1 == coverage times only, -1 == failure
 *
 */
template< typename DBObjectSubtype >
int getPassParams( const DBObjectSubtype& Ref, double& lon0, char& direction, int& orbit_number, 
                  double& orbit_start_minutes, double& coverage_start_minutes, double& coverage_end_minutes ) {
  
  if ( Ref.ds_start_dateIsNull() || Ref.ds_start_timeIsNull() 
      || Ref.ds_end_dateIsNull() || Ref.ds_end_timeIsNull() ) {
    return -1;
  }
 
  // should be able to at least get coverage times
  
  coverage_start_minutes = GeoSearchUtil::t_rel_minutes( Ref.getDs_start_date(), Ref.getDs_start_time() );
  
  coverage_end_minutes = GeoSearchUtil::t_rel_minutes( Ref.getDs_end_date(), Ref.getDs_end_time() );
  
  // try for crucial geographic information
  
  if ( Ref.eq_x_lon_100th_degIsNull() || Ref.eq_x_dateIsNull() || Ref.eq_x_timeIsNull() ) {
    return 1;
  }
  
  lon0 = 0.01 * Ref.getEq_x_lon_100th_deg();
  
  orbit_start_minutes = GeoSearchUtil::t_rel_minutes( Ref.getEq_x_date(), Ref.getEq_x_time()  );
  
  // other information can use default values without much trouble
  
  direction = *( Ref.getAsc_desc_flag( "B" ) );
  
  orbit_number = Ref.getOrbit_number( 0 );
    
  return 0;
}

/**
 * output a result to any kind of stream
 *
 * @tparam inv_geo sub-type of DBObject
 * 
 * @param Ref instance variable
 *
 * @param useFile style for printing column values
 *
 * @param outFile output stream to write to
 *
 */
template< typename inv_geo >
void printFileResult( inv_geo& Ref, int useFile, fstream& outFile ) {
  
  if ( useFile == 1 ) {
    Ref._PrintFile( outFile );
  } else {
    Ref._XMLFile( outFile );
  }
  
  int first_orbit = Ref.getOrbit_number();
  int last_orbit = first_orbit;
  
  if ( ! Ref.geoSearchExec->isOverlapEmpty() ) {
    outFile << Ref.geoSearchExec->xmlOverlapList();
    last_orbit = Ref.geoSearchExec->overlapList().back().orbitNumber();
  }
  
  outFile << "<first_orbit>" << first_orbit << "</first_orbit>";
  outFile << "<last_orbit>" << last_orbit << "</last_orbit>";
  outFile << "<number_orbits>" << last_orbit - first_orbit + 1 << "</number_orbits>";
  
	outFile << endl;
}


#endif // GEOSEARCH_H
