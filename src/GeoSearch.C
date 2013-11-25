//! $Id$
/**
 @file
 @author Ken Knowles Kenneth.Knowles@NOAA.gov
 @date 30-Sep-2012
 
 @defgroup geosearch
 */

#include <string>
#include <vector>
#include <sstream>
#include <ctime>

#include "GeoSearch.h"

/**
 * set up a search executive with config obtained from user input criteria
 *
 * @param nlat northern bound
 *
 * @param slat southern bound
 *
 * @param wlon western bound (regardless of min/max longitude)
 *
 * @param elon eastern bound (regardless of min/max longitude)
 *
 * @param accepted_dirs_str "A"== Asc only, "D"== Desc only, "B"== either or Both
 *
 */
GeoSearch::GeoSearch( double nlat, double slat, double wlon, double elon, const char* accepted_dirs_str ) {
  
  use_geosearch = false;
  is_initialized = false;
  
  if ( ( nlat == 0 && slat == 0 && wlon == 0 && elon == 0 )
      || ( nlat == 90 && slat == -90 && wlon == -180 && elon == 180 ) ) return; // no geosearch for global criteria
  
  use_geosearch = true;

  search_bounds_lat = new GeoLatRange( nlat, slat );
  
  search_bounds_lon = new GeoLonRange( wlon, elon );
  
  if ( strchr( "ADB", accepted_dirs_str[0] ) ) {
    accepted_dirs = accepted_dirs_str[0];
  } else {
    accepted_dirs = 'B';
  }
  
}

/**
 * set parameters that are generic to all rows in this data type
 *
 * @param xl_step accuracy adjustment parameter extracted from first row
 *
 * @param orbit orbit finder constructed from values extracted from first row
 *
 */
void GeoSearch::initDatatypeParams( double xl_step, const GeoOrbitFinder& orbit_ ) {
  
  search_area = new GeoSpatialExtent( *search_bounds_lat, *search_bounds_lon, xl_step, xl_step );
  orbit = new GeoOrbitFinder( orbit_ );
  
  // set quadrant params
  // adjust region of interest to fit within bounds of lat_range_cover
  double nlat = search_area->lat_north();
  double slat = search_area->lat_south();
  
  if ( orbit->latCover().lat_north() - accuracy_deg <= nlat ) {
    nlat = orbit->latCover().lat_north() - accuracy_deg;
    if ( nlat - accuracy_deg <= slat ) slat = nlat - accuracy_deg;
  }
  if ( slat <= orbit->latCover().lat_south() + accuracy_deg ) {
    slat = orbit->latCover().lat_south() + accuracy_deg;
    if (nlat <= slat + accuracy_deg ) nlat = slat + accuracy_deg;
  }
  
  // set parameters of search region for each orbit quadrant 
  
  for ( int q = 0; q < NUM_ORBIT_QUADRANTS_GEOSEARCH; q++ ) {
    search_region_quadrant[q].nlat = nlat;
    search_region_quadrant[q].slat = slat;
    search_region_quadrant[q].extent = NULL;
    search_region_quadrant[q].dir = ASC;
    search_region_quadrant[q].in_play = false;
  }
  
  search_region_quadrant[1].dir = search_region_quadrant[2].dir = DES;
  
  if ( slat < 0 ) search_region_quadrant[0].slat = search_region_quadrant[1].slat = latEquator + accuracy_deg;
  if ( 0 <= nlat ) search_region_quadrant[2].nlat = search_region_quadrant[3].nlat = latEquator - accuracy_deg;
  
  if ( accepted_dirs != 'D' && 0 < nlat ) search_region_quadrant[0].in_play = true;
  if ( accepted_dirs != 'A' && 0 < nlat ) search_region_quadrant[1].in_play = true;
  if ( accepted_dirs != 'A' && slat < 0 ) search_region_quadrant[2].in_play = true;
  if ( accepted_dirs != 'D' && slat < 0 ) search_region_quadrant[3].in_play = true;
  
  // create oversampled representations of the region of interest by quadrant
  
  for ( int q = 0; q < NUM_ORBIT_QUADRANTS_GEOSEARCH; q++ ) {
    
    if ( search_region_quadrant[q].in_play ) {
      search_region_quadrant[q].extent = new GeoSpatialExtent( GeoLatRange( search_region_quadrant[q].nlat, search_region_quadrant[q].slat ), 
                                                     GeoLonRange( search_area->lon_west(), search_area->lon_east() ), 
                                                     xl_step, xl_step );
    }
  }
  
  is_initialized = true;
  
}

/**
 * Destructor
 *
 */
GeoSearch::~GeoSearch() {
  
  delete search_area;
  delete orbit;
  delete search_bounds_lat;
  delete search_bounds_lon;
  for ( int q = 0; q < NUM_ORBIT_QUADRANTS_GEOSEARCH; q++ ) delete search_region_quadrant[q].extent;
}

/**
 * Determine a list of orbit segments that cover a data set
 * 
 * @param orbit_number sequential identifier for orbit starting at lon0
 *
 * @param orbit finder for swath coverage
 *
 * @param lon0 ascending longitude in degrees for orbit pass of reference segment
 * 
 * @param orbit_start_minutes time at lon0 equator crossing
 * 
 * @param start_coverage_minutes start of actual data coverage
 *
 * @param end_coverage_minutes end of data 
 *
 * @return pointer to list of new orbit segments (zero length on failure)
 *
 */
std::vector<GeoOrbitSegment> GeoSearch::orbitalPassCoverage( int orbit_number, const GeoOrbitFinder& orbit, 
                                                            double lon0, double orbit_start_minutes, 
                                                            double coverage_start_minutes, double coverage_end_minutes ) {
  
  double total_minutes = coverage_end_minutes - coverage_start_minutes;
  std::vector<GeoOrbitSegment> list;
  list.clear();
  
  if ( total_minutes <= 0 || geosearch::coverage_max_minutes < total_minutes ) return list;
  
  double one_orbit_dlon = geosearch::EarthRotation_deg( 360, orbit.full_minutes() );
  
  double start_deg = ( ( coverage_start_minutes - orbit_start_minutes ) / orbit.full_minutes() ) * 360;
  
  double end_deg = ( ( coverage_end_minutes - orbit_start_minutes ) / orbit.full_minutes() ) * 360;
  
  // normalize start of orbit progress
  // start and end are adjusted together to preserve the interval
  
  while ( start_deg < 0 ) {
    orbit_number -= 1;
    orbit_start_minutes -= orbit.full_minutes();
    lon0 += one_orbit_dlon;
    start_deg += 360;
    end_deg += 360;
  }
  
  while ( 360 <= start_deg ) {
    orbit_number += 1;
    orbit_start_minutes += orbit.full_minutes();
    lon0 -= one_orbit_dlon;
    start_deg -= 360;
    end_deg -= 360;
  }
  
  int num_segments = static_cast<int>( ceil( end_deg / 360. ) );
  
  for ( int isegment = 0; isegment < num_segments - 1; isegment++ ) {
    
    list.push_back( GeoOrbitSegment( orbit_number, orbit, lon0, coverage_start_minutes, 
                                    orbit_start_minutes, start_deg, 360. ) );
    
    start_deg = 0;
    end_deg -= 360;
    orbit_start_minutes += orbit.full_minutes();
    orbit_number += 1;
    lon0 -= one_orbit_dlon;
  }
  
  list.push_back( GeoOrbitSegment( orbit_number, orbit, lon0, coverage_start_minutes, 
                                  orbit_start_minutes, start_deg, end_deg ) );
  
  return list;
}

/**
 * Identify sections of orbit segments that overlap the search area
 *
 * @param data_set_pass full coverage of data set
 *
 * @param whether or not intersection was found
 *        if found then list of overlapping segments is stored
 */
bool GeoSearch::findOverlap( std::vector<GeoOrbitSegment> data_set_pass ) {
  
  bool found = false;
  overlap_list.clear();
  
  if ( data_set_pass.empty() ) return false;
  
  // check that there is a chance for overlap with this orbit finder
  
  if ( orbit->latCover().lat_north() <= search_area->lat_south() ||
      search_area->lat_north() <= orbit->latCover().lat_south() ) return false;
  
  // test each orbit segment for overlap in each quadrant
  
  for ( std::vector<GeoOrbitSegment>::iterator segment = data_set_pass.begin(); segment != data_set_pass.end(); segment++ ) {
    
    if ( segment->isSameOrbit( *orbit ) ) {
      
      for ( int q = 0; q < NUM_ORBIT_QUADRANTS_GEOSEARCH; q++ ) {
        
        if ( search_region_quadrant[q].in_play ) {
          GeoOrbitSegment candidate = segment->overlap( *(search_region_quadrant[q].extent), search_region_quadrant[q].dir );
          if ( ! candidate.isEmpty() ) {
            found = true;
            overlap_list.push_back( candidate );
          }
        }
      }
    }
  }
  
  return found;
}

/**
 *   Express equator crossing constraints as a where clause condition
 *
 *  @param asc_desc_flag "A" == means find orbits that pass over the region of 
 *  interest in the ascending direction only. "D" == descending only. "B" == both.
 *
 *  @note All latitudes and longitudes are expressed in decimal degrees.
 *
 *  @param nLat latitude of northern limit of region of interest
 *
 *  @param sLat latitude of southern limit of region of interest
 *
 *  @param wLon longitude of western limit of region of interest
 *
 *  @param eLon longitude of eastern limit of region of interest
 *
 */
std::string GeoSearch::eqCrossCondition( std::string asc_desc_flag, float nLat, float sLat, float wLon, float eLon) {
  
  // there's only one typical orbit, so cache it for further use
  static GeoOrbitFinder *typicalOrbit = NULL;
  if ( NULL == typicalOrbit ) {
    typicalOrbit = new GeoOrbitFinder( geosearch::typical_inclination_deg, 
                                      geosearch::typical_period_minutes, 
                                      geosearch::typical_swath_width_m );
  }
  
  GeoSpatialExtent region_tmp( GeoLatRange( nLat, sLat ), GeoLonRange( wLon, eLon ), 
                              geosearch::latStepEqCross, geosearch::lonStepEqCross );
  
  std::list<GeoLonRange> lon_xranges;
  
  if (asc_desc_flag == "A" || asc_desc_flag == "B") {      
    lon_xranges.push_back( typicalOrbit->crossingLonRange( region_tmp, ASC ) );
  }
  
  if (asc_desc_flag == "D" || asc_desc_flag == "B") {
    lon_xranges.push_back( typicalOrbit->crossingLonRange( region_tmp, DES ) );
  }
  
  return geosearch::rangeListCondition( lon_xranges, geosearch::eqCrossLonName, geosearch::latlonIntScale );
}

/**
 *  xml representation of overlap list
 *
 *  @return xml string
 *
 */
std::string GeoSearch::xmlOverlapList() {

  std::stringstream xmlStream;
  int map_index = 0;
  int firstOrbit = 999999;
  int lastOrbit = -1;
  

  for ( std::vector<GeoOrbitSegment>::iterator segment = overlap_list.begin();
       segment != overlap_list.end(); segment++ ) { 
    
    ++map_index;

    if ( segment->orbitNumber() > lastOrbit ) { lastOrbit = segment->orbitNumber(); }
    if ( segment->orbitNumber() < firstOrbit ) { firstOrbit = segment->orbitNumber(); }
    
    xmlStream << segment->xmlOverlap( map_index );
  }
  
  xmlStream << "<number_orbits>" << lastOrbit - firstOrbit + 1 << "</number_orbits>";
  
  return xmlStream.str();
}


