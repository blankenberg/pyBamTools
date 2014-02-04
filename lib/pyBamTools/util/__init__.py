#Dan Blankenberg
import sys

import numpy

NULL_CHAR = '\x00'

MAX_INT = sys.maxint

MAX_NEG_INT = -MAX_INT

NUMPY_DTYPES = dict( uint8=numpy.uint8, uint16=numpy.uint16, uint32=numpy.uint32, uint64=numpy.uint64 )

INTEGER_TYPES = tuple( NUMPY_DTYPES.values() + [ int, long ] )

def is_integer( value ):
    """Check if value is a valid unsigned int, int, or long"""
    return isinstance( value, INTEGER_TYPES )

def get_region_overlap( region_1, region_2 ):
    return max( 0, min( region_1[1], region_2[1] ) - max( region_1[0], region_2[0] ) )

def get_region_overlap_with_positions( region_1, region_2 ):
    start = max( region_1[0], region_2[0] )
    end = min( region_1[1], region_2[1] )
    return ( max( 0, end - start ), start, end )
