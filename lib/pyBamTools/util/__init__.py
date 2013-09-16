#Dan Blankenberg
import sys

import numpy

NULL_CHAR = '\x00'

MAX_NEG_INT = -sys.maxint

NUMPY_DTYPES = dict( uint8=numpy.uint8, uint16=numpy.uint16, uint32=numpy.uint32, uint64=numpy.uint64 )

INTEGER_TYPES = tuple( NUMPY_DTYPES.values() + [ int, long ] )

def is_integer( value ):
    """Check if value is a valid unsigned int, int, or long"""
    return isinstance( value, INTEGER_TYPES )
