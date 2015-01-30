#Dan Blankenberg
import sys

import numpy

NULL_CHAR = '\x00'

MAX_INT = sys.maxint

MAX_NEG_INT = -MAX_INT

NUMPY_DTYPES = dict( uint8=numpy.uint8, uint16=numpy.uint16, uint32=numpy.uint32, uint64=numpy.uint64 )

NUMPY_DTYPE_SIZES = [ ( numpy.uint8, 255 ),
                      ( numpy.uint16, 65535 ),
                      ( numpy.uint32, 4294967295 ),
                      ( numpy.uint64, 18446744073709551615 )
                      ]

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

def guess_numpy_dtypes_from_idxstats( bams, default=None ):
    if not isinstance( bams, list ):
        bams = [ bams ]
    rval = None
    for bam_index in bams:
        if hasattr( bam_index , '_bam_index' ):
            bam_index = bam_index._bam_index
        refs = bam_index._references
        assert isinstance( refs, list ), ValueError( "Invalid object provided to guess_numpy_dtypes_from_idxstats()" )
        if rval is None:
            rval = [ None ] * len( refs )
        for i, val in enumerate( refs ):
            assert isinstance( val, dict ), ValueError( "Invalid object provided as reference element to guess_numpy_dtypes_from_idxstats()" )
            val = val.get( 'n_mapped', None )
            if val is not None:
                cur_val = rval[ i ] or 0
                rval[i] = cur_val + val
    for i, val in enumerate( rval ):
        if val is None:
            # Default to default
            dtype = default
        else:
            for ( dtype, size ) in NUMPY_DTYPE_SIZES:
                if val < size:
                    break
        rval[i] = dtype
    return rval

def guess_array_memory_usage( bam_readers, dtype, use_strand=False ):
    """Returns an estimate for the maximum amount of memory to be consumed by numpy arrays."""
    ARRAY_COUNT = 5
    if not isinstance( bam_readers, list ):
        bam_readers = [ bam_readers ]
    if isinstance( dtype, basestring ):
        dtype = NUMPY_DTYPES.get( dtype, None )
    use_strand = use_strand + 1 #if false, factor of 1, if true, factor of 2
    dtypes = guess_numpy_dtypes_from_idxstats( bam_readers, default=None )
    if not [ dt for dt in dtypes if dt is not None ]:
        dtypes = guess_numpy_dtypes_from_idxstats( bam_readers, default=dtype )
    read_groups = []
    no_read_group = False
    max_ref_size = 0
    for bam in bam_readers:
        rgs = bam.get_read_groups()
        if rgs:
            for rg in rgs:
                if rg not in read_groups:
                    read_groups.append( rg )
        else:
            no_read_group = True
        for i, ( name, length ) in enumerate( bam.get_references() ):
            if dtypes[i] is not None:
                max_ref_size = max( max_ref_size, ( length * dtypes[i]().nbytes) )
    read_groups = len( read_groups ) + no_read_group
    return max_ref_size * ARRAY_COUNT * use_strand * read_groups
