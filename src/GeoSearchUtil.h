/**
 @file
 @author Ken Knowles Kenneth.Knowles@NOAA.gov
 @date 04-Aug-2012
 @version $Id: GeoSearchUtil.h 16300 2013-05-08 20:15:55Z kknowles $
 @ingroup geosearch
 */

/**
 * This file contains the base class for Geosearch
 * and defines constants and auxiliary functions
 * used throughout.
 */

#ifndef GEOSEARCHUTIL_H
#define GEOSEARCHUTIL_H


#include <string>
#include <cmath>
#include <float.h>
#include <ctime>

#include "DBObject.hxx"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// satellite direction with respect to equator crossing
enum SatellDir { DES, ASC };

// swath position facing down range
static const int GEOSEARCH_NUM_SWATH_POS = 3;
enum PosIndex { 
  STBD_EDGE, ///< starboard (right)
  PORT_EDGE, ///< port (left)
  SWATH_CENTER ///< center ( may be offset from nadir )
};

// hemisphere
static const int GEOSEARCH_NUM_HEMIS = 2;
enum HemisIndex { SOUTH_HEMIS, NORTH_HEMIS };

// units for variables
// letters (or ending in _m) are in meters
// lat lon (or ending in _deg) are in decimal degrees
// greek letters (spelled out: phi, lam, beta, sigma...) are in radians
// underscore means subscript for mathematical variables
enum unitsType { METERS, DEGREES, RADIANS };

// lat/lon quadrilateral corners
enum corner_index { NW, NE, SW, SE };
static const int GEOSEARCH_NUM_CORNERS = 4;

//
// @note constants that are true even outside of Geosearch
//
namespace {
  
  /// Usefull symbolic constants
  static const double latEquator = 0;
  static const double latNorthPole = 90;
  static const double latSouthPole = -90;
  static const double latGlobalSpan = 180;
  static const double lonWestIDL = -180;
  static const double lonEastIDL = 180;
  static const double lonGlobalSpan = 360;
  
  /// spherical approximation for Earth radius in meters
  static const double defaultRadius_m = 6371007.2;
  
  // unit conversions
  inline double toRadians( double angle_deg ) { return angle_deg * M_PI/180; }
  inline double toDegrees( double angle_rad ) { return angle_rad * 180/M_PI; }
  inline long int minutesToMsec( double t_minutes ) { return static_cast<long int>( t_minutes * 60 * 1000 ); }
  inline double msecToMinutes( long int t_msec ) { return static_cast<double>( t_msec ) / 1000 / 60; }
  inline int scaledLatLon( double latlon_deg ) { return static_cast<int>( latlon_deg * 100 ); }
  
}

//
// @note operations that may have wider use 
// but the implementation here is specific to Geosearch
//
class GeoSearchUtil {
  
public:
  
  /// Tunable constant for dismissing roundoff error
  /// typically used for trig functions
  static const double insigValue = 10*FLT_EPSILON; 
  
  /// tolerance for latitude/longitude coordinates
  static const double accuracy_deg = 0.01;
  
  /// tolerance for rectangular coordinates 
  static double const accuracy_m = 1000;
  
  /// tolerance for time intervals
  static const double accuracy_minutes = 0.25;
  
  /// out of bounds indicator
  static const double UNDEFINED = 999;
  
  /// significant value tests
  static bool isZero( double quantity );
  static bool isZero_deg( double latlon_deg );
  static bool isZero_m( double distance_m );
  static bool isZero_minutes( double t_minutes );
  
  static double t_rel_minutes( time_t t_sec, double fraction_sec );
  static double t_rel_minutes( std::string dsdt_str );
  static double t_rel_minutes( std::string date_str, int time_int );
  static std::string t_rel_dsdt( double t_minutes );
  
};

/// Simple test for small numbers close enough to be considered zero
inline bool GeoSearchUtil::isZero( const double quantity ) { return fabs( quantity ) <= insigValue; }

/// Simple test for small distances close enough to be considered zero
inline bool GeoSearchUtil::isZero_m( double distance_m ) { return fabs( distance_m ) < accuracy_m; }

/// Simple test for small lat or lon close enough to be considered zero
inline bool GeoSearchUtil::isZero_deg( double latlon_deg ) { return fabs( latlon_deg ) < accuracy_deg; }

/// Simple test for small time intervals close enough to be considered zero
inline bool GeoSearchUtil::isZero_minutes( double t_minutes ) { return fabs( t_minutes ) < accuracy_minutes; }

/**
 * Convert system time to relative time in minutes
 *
 * @param t_sec system time in integer seconds
 *
 * @param fraction_sec fraction after decimal point
 *
 * @return minutes relative time within Geosearch
 *
 */
inline double GeoSearchUtil::t_rel_minutes( time_t t_sec, double fraction_sec ) {
  
  return ( static_cast<double>( t_sec ) + fraction_sec ) / 60.;
}

/**
 * Convert data set date time string to relative time in minutes
 *
 * @param dsdt_str with form "yyyy-mm-dd hh:mm:ss.xxx"
 *
 * @return minutes relative time within Geosearch
 *
 */
inline double GeoSearchUtil::t_rel_minutes( std::string dsdt_str ) {
  
  struct tm t_struct;
  
  if ( dsdt_str.length() < 10 ) return 0;
  
  t_struct.tm_year = atoi( dsdt_str.substr( 0, 4 ).c_str() ) - 1900;
  t_struct.tm_mon = atoi( dsdt_str.substr( 5, 2 ).c_str() ) - 1;
  t_struct.tm_mday = atoi( dsdt_str.substr( 8, 2).c_str() );
  
  if ( dsdt_str.length() > 10 ) {
    t_struct.tm_hour = atoi( dsdt_str.substr( 11, 2 ).c_str() );
    t_struct.tm_min = atoi( dsdt_str.substr( 14, 2 ).c_str() );
    t_struct.tm_sec = atoi( dsdt_str.substr( 17, 2 ).c_str() );
  } else {
    t_struct.tm_hour = t_struct.tm_min = t_struct.tm_sec = 0;
  }      
  
  double fraction_sec = ( dsdt_str.length() > 19 ) ? atof( dsdt_str.substr( 19 ).c_str() ) : 0;
  
  time_t t_sec = timegm( &t_struct );
  
  return t_rel_minutes( t_sec, fraction_sec );
}

/**
 * Convert date and time to relative time in minutes
 *
 * @param date_str with form "yyyymmdd"
 *
 * @param time_int with form hhmmss
 *
 * @return minutes relative time within Geosearch
 *
 */
inline double GeoSearchUtil::t_rel_minutes( std::string date_str, int time_int ) {
  
  struct tm t_struct;
  
  if ( date_str.length() < 10 ) return 0;
  
  t_struct.tm_year = atoi( date_str.substr( 0, 4 ).c_str() ) - 1900;
  t_struct.tm_mon = atoi( date_str.substr( 5, 2 ).c_str() ) - 1;
  t_struct.tm_mday = atoi( date_str.substr( 8, 2).c_str() );
  
  t_struct.tm_hour = time_int / 10000;
  t_struct.tm_min = ( time_int / 100 ) % 100;
  t_struct.tm_sec = time_int % 100;
  
  time_t t_sec = timegm( &t_struct );
  
  return t_rel_minutes( t_sec, 0.0 );
}

/**
 * Convert relative time in minutes to data set date time (dsdt) formatted string
 *
 * @param t_minutes relative time within Geosearch
 *
 * @return "yyyy-mm-dd hh:mm:ss.xxx"
 *
 */
inline std::string GeoSearchUtil::t_rel_dsdt( double t_minutes ) {
  
  time_t t_sec = static_cast<time_t>( t_minutes * 60 );
  int remainder_msec = static_cast<int>( 1000 * ( t_minutes * 60 - t_sec ) );
  if ( remainder_msec < 0 ) { t_sec -= 1; remainder_msec += 1000; }
  
  struct tm* t_struct = gmtime( &t_sec );
  
  const int T_CHAR_SIZE = 32;
  char t_char[T_CHAR_SIZE];
  
  snprintf( t_char, T_CHAR_SIZE - 1, "%4.4d-%02d-%02d %02d:%02d:%02d.%03d",
           t_struct->tm_year+1900, t_struct->tm_mon+1, t_struct->tm_mday,
           t_struct->tm_hour, t_struct->tm_min, t_struct->tm_sec, remainder_msec );
  
  return std::string( t_char );
}


//
// @note operations and constants specific to Geosearch
//
namespace geosearch {
  
  // maximum duration of an orbital swath coverage
  const int coverage_max_minutes = 24 * 60;
  
  /**
   * combine range intervals as SQL string
   *
   * @tparam Container type of ranges
   *
   * @param ranges list of ranges that have a "condition" method
   * 
   * @param field_name argument to "condition" function
   * 
   * @param scale argument to "condition" function
   * 
   * @return SQL string that ORs together the list of conditions
   */
  template< typename Container >
  std::string rangeListCondition( Container ranges, std::string field_name, int scale ) {
    
    std::string *SQL = new std::string("");
    
    for ( typename Container::iterator range = ranges.begin(); range != ranges.end(); range++ ) {
      
      std::string next = range->condition( field_name, scale );
      
      if ( next.size() > 0 ) {
        *SQL += SQL->size() > 0 ? " OR " : "(";
        *SQL += next;
      }
    }
    
    if ( SQL->size() > 0 ) *SQL += ")";
    
    return *SQL;
  }
  
  /**
   * test for inclusion within any of a group of ranges
   *
   * @tparam Container type of ranges
   * 
   * @param ranges list of ranges that have a "contains" method
   * 
   * @param value argument to "contains" function
   * 
   * @return true if at least one range contains value
   */
  template< typename Container >
  bool rangeListContains( Container ranges, double value ) {
    
    for ( typename Container::iterator range = ranges.begin(); range != ranges.end(); ranges++ ) {
      if ( range->contains( value ) ) return true;
    }
    
    return false;
  }
  
  /**
   * Distance from origin to surface of Earth at specific location.
   * 
   * @note This is the geocentric radius, which may differ from the 
   * geodetic radius that elevations are based on depending on the figure
   * used to represent the Earth. For now, assuming a spherical Earth,
   * there is no difference btw geodetic and geocentric.
   *  
   * @param lat latitude in degrees
   * 
   * @param lon longitude in degrees
   * 
   * @return distance from origin to surface in meters at the given location
   */
  inline double surfaceRadius_m( double lat, double lon ) { return defaultRadius_m; }
  
  /**
   * find amount Earth rotates for a given progress along the orbit
   *
   * @note 24 * 60 minutes per revolution is only correct here for sun-synchronous orbits.
   * Technically, the length of a rotation should be offset to compensate for
   * the movement of the Earth around the Sun and then corrected for the
   * precession of the orbit of the satellite around the Earth. In the case of
   * sun-synchronous satellites, by design, these offsets exactly cancel each other.
   * It's safe to make this assumption because all polar orbiting environmental
   * satellites are in sun-synchronous orbits. Regardless, the offsets are small.
   *
   * @param progress angle from equator crossing
   *
   * @param period_minutes time for one full orbit
   *
   * @return angle rotated
   */
  inline double EarthRotation_rad( double progress_rad, double period_minutes) { 
    return progress_rad * period_minutes / ( 24*60 ); 
  }
  
  inline double EarthRotation_deg( double progress_deg, double period_minutes ) {
    return toDegrees( EarthRotation_rad(  toRadians( progress_deg ), period_minutes ) );
  }
  
  /**
   * Central angle
   *
   * @param dist_m distance along surface of Earth
   *
   * @return angle traversed measured at center of Earth
   */
  inline double centralAngle_deg( double dist_m ) { return toDegrees( dist_m / defaultRadius_m ); }
  
  /// field name
  const std::string eqCrossLonName = "eq_x_lon_100th_deg";
  
  /// multiplier from program's double to database's int
  const int latlonIntScale = 100;
  
  /// adjust algorithm speed versus accuracy
  const double latStepEqCross = 10;
  const double lonStepEqCross = 20;
  
  /// orbital parameters
  
  // copied from old eqtrig.c
  // const double typical_inclination_deg = 98.922;
  // const double typical_period_minutes = 102.1277;
  
  // NOAA 19
  // const double typical_inclination_deg = 98.73;
  // const double typical_period_minutes = 102.14;
  
  // NPP
  const double typical_inclination_deg = 98.7172;
  const double typical_period_minutes = 101.44;
  
  // equivalent to full coverage at the equator
  const double typical_swath_width_m = 2840000;
  
  // default accuracy
  const double default_xl_step = 0.1;
  
}

#endif // GEOSEARCHUTIL_H

