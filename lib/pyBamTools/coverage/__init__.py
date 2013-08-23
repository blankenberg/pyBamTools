#Dan Blankenberg
import sys

import numpy

from ..util.nuc import DELETIONS_KWD
from ..util.nuc import INSERTIONS_KWD

from ..util.nuc import NUCLEOTIDES_UPPER
from ..util.nuc import UNKNOWN_NUCLEOTIDE_NAME

DEFAULT_SEQUENCE_LENGTH = 255000000

class NamedRegionOverlap( object ):
    
    def __init__( self, sequence_lengths, regions=None ):
        self._sequence_lengths = sequence_lengths
        self._regions = {}
        self._complete_regions = []
        if regions:
            for region in regions:
                self.add_region( region )
        
    def add_region( self, region ):
        if isinstance( region, basestring ):
            if region not in self._complete_regions:
                self._complete_regions.append( region )
            if region in self._regions:
                del self._regions[ region ]
        else:
            region_name, region_start, region_end = region
            if region_start >= self._sequence_lengths[ region_name ]:
                return #start is before end of sequence, can't add anything
            region_start = max( region_start, 0 ) #no neg region starts
            region_end = min( region_end, self._sequence_lengths[ region_name ] ) #region maxes out at length of sequence
            if region_name not in self._complete_regions:
                if region_name not in self._regions:
                    self._regions[ region_name ] = self._coverage_array = numpy.zeros( ( self._sequence_lengths[ region_name ] ), dtype=numpy.bool )
                self._regions[ region_name ][ region_start:region_end ] = numpy.ones( region_end - region_start, dtype=numpy.bool )
    def is_empty( self ):
        return not self._complete_regions and not self._regions
    def region_overlaps( self, region_name, region_start=None, region_end=None, empty_default=False ):
        if self.is_empty():
            return empty_default
        if isinstance( region_name, tuple ):
            region_name, region_start, region_end = region_name
        if region_name in self._complete_regions:
            return True
        if region_name in self._regions:
            if region_start is None:
                if region_end is None:
                    return True
                region_start = region_end
                region_end += 1
            elif region_end is None:
                region_end = region_start + 1
            return True in self._regions[ region_name ][ region_start:region_end ]
        return False
    def iter_covered_regions( self, region_names=None ):
        if not region_names:
            region_names = self._sequence_lengths.keys()
        elif isinstance( region_names, basestring ):
            region_names = [ region_names ]
        for region_name in region_names:
            if region_name in self._complete_regions:
                yield region_name, 0, self._sequence_lengths[ region_name ]
            elif region_name in self._regions:
                try:
                    start_pos = numpy.where( self._regions[ region_name ] == True )[0][0]
                except IndexError:
                    continue
                sequence_size = self._sequence_lengths[ region_name ]
                #max_index = sequence_size - 1
                while True:
                    try:
                        end_pos = numpy.where( self._regions[ region_name ][ start_pos: ] == False )[0][0] + start_pos
                    except IndexError:
                        yield region_name, start_pos, sequence_size
                        break
                    yield region_name, start_pos, end_pos
                    try:
                        start_pos = numpy.where( self._regions[ region_name ][ end_pos: ] == True )[0][0] + end_pos
                    except IndexError:
                        break

class SequenceCoverage( object ):
    def __init__( self, size=DEFAULT_SEQUENCE_LENGTH, dtype=numpy.uint8 ):
        self._coverage_array = numpy.zeros( ( size ), dtype=dtype )
        self._size = size
    def __str__( self ):
        return str( self._coverage_array )
    def __iter__( self ):
        return iter( self._coverage_array )
    def get( self, position, default=None ):
        if 0 < position >= self._size:
            print >>sys.stderr, 'Warning: coverage of out of bounds (size=%s) position has been requested: %s' % ( self._size, position )
            return default #fixme: 0 or None for default here?
        return self._coverage_array[ position ]
    def set( self, position, value ):
        if 0 < position >= self._size:
            print >>sys.stderr, 'Warning: setting coverage of out of bounds (size=%s) position has been requested: %s to be set to %s. This information has been ignored.' % ( self._size, position, value )
            return None
        try:
            self._coverage_array[ position ] = value
        except OverflowError:
            print >> sys.stderr, "Value of %i at position %i is too large to fit into array, keeping current value (%i)." % ( value, position, self._coverage_array[ position ] )
        return self._coverage_array[ position ]
    def sum( self, beg=None, end=None ):
        if beg == end == None:
            return self._coverage_array.sum()
        else:
            return self._coverage_array[beg:end].sum()

class NucleotideCoverage( object ):
    
    def __init__( self, size=DEFAULT_SEQUENCE_LENGTH, dtype=numpy.uint8, reference_name=None ): 
        self._nucleotide_dict = {}
        self.size = size
        for nuc in NUCLEOTIDES_UPPER + [ UNKNOWN_NUCLEOTIDE_NAME ]:
            self._nucleotide_dict[ nuc ] = SequenceCoverage( size, dtype=dtype )
        self._insert_dict = {}
        self._delete_dict = {}
        self._reference_name = reference_name
    def add_read( self, read, min_base_quality=None, die_on_error=False ):
        position = read.get_position( one_based=False )
        sequence = list( read.get_seq() )
        quality = list(read.get_qual_tuple())
        cigar = read.get_cigar()
        position_offset = 0
        while cigar:
            cigar_size, cigar_op = cigar.pop( 0 )
            if cigar_op in [ 0, 7, 8 ]: 
                #M alignment match (can be a sequence match or mismatch); = sequence match; x sequence mismatch
                #= sequence match
                #X sequence mismatch
                for i in range( cigar_size ):
                    nuc = sequence.pop( 0 ).upper()
                    cov = self._nucleotide_dict.get( nuc, None )
                    if not cov:
                        cov = self._nucleotide_dict[ UNKNOWN_NUCLEOTIDE_NAME ]
                    j = position + position_offset + i
                    base_quality = quality.pop(0)
                    if min_base_quality is None or min_base_quality >= base_quality:
                        if cov.set( j, cov.get( j, 0 ) + 1 ) is None:
                            print >> sys.stderr, "Error causing read:", read.to_sam()
                position_offset += cigar_size
            elif cigar_op == 1: #I insertion to the reference
                #remove inserted seqs and qualities
                #no change in position
                i = position + position_offset
                if i < 0:
                    #TODO: see if it is ok to wrap around...eg, -2 index as option for circ chromes
                    if die_on_error:
                        raise Exception( "Insertion representation requires a position less than one:\n%s" % ( read.to_sam() ) )
                    else:
                        sys.stderr.write( "Insertion representation requires a position less than one:\n%s\n" % ( read.to_sam() ) )
                        return
                base_quality = quality[:cigar_size]
                if min_base_quality is None or numpy.mean( base_quality ) >= min_base_quality:
                    if i not in self._insert_dict:
                        self._insert_dict[i] = []
                    self._insert_dict[i].append( ''.join( sequence[:cigar_size] ) )
                sequence = sequence[cigar_size:]
                quality = quality[cigar_size:]
            elif cigar_op == 2: #D deletion from the reference
                #update position
                i = position + position_offset
                if i < 0:
                    #TODO: see if it is ok to wrap around...eg, -2 index as option for circ chromes
                    if die_on_error:
                        raise Exception( "Deletion representation requires a position less than one:\n%s" % ( read.to_sam() ) )
                    else:
                        sys.stderr.write( "Deletion representation requires a position less than one:\n%s\n" % ( read.to_sam() ) )
                        return
                if i not in self._delete_dict:
                    self._delete_dict[i] = []
                self._delete_dict[i].append( cigar_size )
                position_offset += cigar_size
            elif cigar_op == 3: #N skipped region from the reference
                #just move along the position
                position_offset += cigar_size
            elif cigar_op == 4: #S soft clipping (clipped sequences present in SEQ)
                #Do not set coverage here
                #Remove soft clipped sequence
                #Do not move position
                sequence = sequence[cigar_size:]
                quality = quality[cigar_size:]
            elif cigar_op == 5: #H hard clipping (clipped sequences NOT present in SEQ)
                #Do not move position
                pass
            elif cigar_op == 6: #P padding (silent deletion from padded reference)
                pass
            else: #unknown cigar_op
                if die_on_error:
                    raise ValueError( 'The cigar operation "%s" of size "%s" for the read below is unknown.\n%s)' % ( cigar_op, cigar_size, read ) )
                else:
                    sys.stderr.write( 'The cigar operation "%s" of size "%s" for the read below is unknown.\n%s)\n' % ( cigar_op, cigar_size, read ) )
                    return
    
    def __str__( self ):
        return str( self._coverage_array )
    def __iter__( self ):
        return iter( self._coverage_array )
    def iteritems( self ):
        return self._nucleotide_dict.iteritems()
    def get( self, position, nucleotides=True, insertions=True, deletions=True ):
        nucs = {}
        if nucleotides:
            for nuc, cov in self.iteritems():
                cov = cov.get( position )
                if cov:
                    nucs[nuc] = cov
        if deletions:
            if position in self._delete_dict:
                nucs[ DELETIONS_KWD ] = self._delete_dict[position]
        if insertions:
            if position in self._insert_dict:
                nucs[ INSERTIONS_KWD ] = self._insert_dict[position]
        return nucs
    def set( self, position, value ):
        try:
            self._coverage_array[ position ] = value
        except OverflowError:
            print >> sys.stderr, "Value of %i at position %i is too large to fit into array, keeping current value (%i)." % ( value, position, self._coverage_array[ position ] )
    def iter_coverage( self ):
        for i in xrange( self.size ):
            nucs = self.get( i )
            if nucs:
                yield i, nucs
