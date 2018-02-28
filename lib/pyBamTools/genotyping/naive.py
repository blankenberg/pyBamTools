#Dan Blankenberg
import string
from sys import stderr
from datetime import date

import numpy

from pyBamParser.fasta import IndexedReferenceSequences

from ..util.odict import odict

from ..util import MAX_NEG_INT
from ..util import NUMPY_DTYPES
from ..util import is_integer
from ..util import get_region_overlap_with_positions
from ..util import guess_numpy_dtypes_from_idxstats

from ..util.sam import SAM_HEADER_NON_TAB_RECORDS
from ..util.sam import SAM_READ_GROUP_RECORD_CODE
from ..util.sam import SAM_NO_READ_GROUP_NAME
from ..util.sam import SAM_READ_GROUP_ID_STR
from ..util.sam import SAM_READ_GROUP_SAMPLE_STR

from ..util.vcf import VCF_NO_VALUE
from ..util.vcf import VCF_UNPHASED_SEPARATOR
from ..util.vcf import VCF_PHASED_SEPARATOR
from ..util.vcf import VCF_GT_PHASED_SEPARATOR 
from ..util.vcf import VCF_INFO_FIELD_SEPARATOR
from ..util.vcf import VCF_FORMAT_FIELD_SEPARATOR

from ..util.nuc import DELETIONS_KWD
from ..util.nuc import INSERTIONS_KWD

from ..coverage import NucleotideCoverage
from ..coverage import NamedRegionOverlap

TAB_CHAR = "\t"
COMMA_CHAR = ","

PROGRAM_NAME = "Naive Variant Caller"
PROGRAM_VERSION = "0.0.3"

class ReadGroupGenotyper( object ):
    __INDEL_OFFSET__ = 0
    __DEFAULT_NUMPY_DTYPE__ = NUMPY_DTYPES.get( 'uint64' )
    
    #TODO: make these filters defined generically
    def __init__( self, bam_readers=None, reference_sequence_filename=None, dtype=None, min_support_depth=None, 
                  min_base_quality=None, min_mapping_quality=None, restrict_regions=None, use_strand=False, 
                  allow_out_of_bounds_positions=False, safe=True ):
        self._read_groups = odict() #order important?
        self._sample_names = []
        self._sequence_lengths = odict()
        self._read_group_coverage = {}#odcit()
        self._read_group_coverage_reverse = {}
        self._no_read_group_coverage = None#odcit()
        self._no_read_group_coverage_reverse = None
        self._has_no_read_groups = False
        self._reference_sequence_filename = reference_sequence_filename
        self._reference_sequences = IndexedReferenceSequences( reference_sequence_filename, sequence_filter=string.upper )
        #add seq lengths from fasta index
        for ref_seq_name in self._reference_sequences.get_sequence_names():
            self._sequence_lengths[ ref_seq_name ] = self._reference_sequences.get_sequence_size_by_name( ref_seq_name )
        if bam_readers is None:
            bam_readers = []
        elif not isinstance( bam_readers, list ):
            bam_readers = [ bam_readers ]
        self._readers = []
        for bam_reader in bam_readers:
            self.add_reader( bam_reader )
        self._use_strand = use_strand
        if dtype and isinstance( dtype, basestring ):
            dtype = NUMPY_DTYPES.get( dtype, None )
        if dtype:
            force_dtype = True
        else:
            dtype = self.__DEFAULT_NUMPY_DTYPE__
            force_dtype = False
        self._dtypes = guess_numpy_dtypes_from_idxstats( bam_readers, dtype, force_dtype=force_dtype )
        #filters:
        self._min_support_depth = min_support_depth or 0
        self._min_base_quality = min_base_quality
        self._min_mapping_quality = min_mapping_quality
        if self._min_mapping_quality is None:
            self._min_mapping_quality = MAX_NEG_INT
        if not restrict_regions:
            restrict_regions = [ ( ref_seq_name, 0, ref_seq_length ) for ref_seq_name, ref_seq_length in self._sequence_lengths.iteritems() ]
        else:
            for i in xrange( len( restrict_regions ) ):
                if isinstance( restrict_regions[i], basestring ):
                    restrict_regions[i] = ( restrict_regions[i], 0, self._sequence_lengths[ restrict_regions[i] ] )
            restrict_regions = self._reference_sequences.sort_region_list( restrict_regions )
        self._restrict_regions_list = restrict_regions
        self._restrict_regions = NamedRegionOverlap( self._sequence_lengths, regions=self._restrict_regions_list )
        self._covered_regions = NamedRegionOverlap( self._sequence_lengths )
        self._allow_out_of_bounds_positions = allow_out_of_bounds_positions
        self._safe = safe
        
        #issue warnings about zero length sequences
        for ref_seq_name, ref_seq_size in self._sequence_lengths.iteritems():
            if ref_seq_size < 1:
                print >>stderr, "Warning: %s was initialized with suspect sequence length (%i)." % ( ref_seq_name, ref_seq_size )
        
    def add_reader( self, bam_reader ):
        self._readers.append( bam_reader )
        for seq_name, seq_len in bam_reader.get_references():
            if seq_name not in self._sequence_lengths:
                self._sequence_lengths[ seq_name ] = seq_len
            else:
                assert self._sequence_lengths[ seq_name ] == seq_len, "Sequence length mismatch for '%s': %s != %s" % ( seq_name, self._sequence_lengths[ seq_name ], seq_len ) #TODO: print out offending reader name
        #add read group information
        headers = bam_reader.get_sam_header_dict()
        if SAM_READ_GROUP_RECORD_CODE in headers:
            for rg in headers[SAM_READ_GROUP_RECORD_CODE]: #{'LB': 'reads_library', 'ID': 'read_group_id', 'PL': 'ILLUMINA', 'SM': 'read_sample'}
                rg_name = rg[SAM_READ_GROUP_ID_STR]
                if rg_name not in self._read_groups:
                    self._read_groups[ rg_name ] = rg #TODO make get_readgroups method
                    #TODO: make sample name default to read group id
                    if rg[SAM_READ_GROUP_SAMPLE_STR] not in self._sample_names:
                        self._sample_names.append( rg[SAM_READ_GROUP_SAMPLE_STR] )
                else:
                    assert self._read_groups[ rg_name ] == rg, "Read Group info mismatch: '%s' != '%s'." % ( self._read_groups[ rg_name ], rg )
        else:
            self._has_no_read_groups = True #for now all non-read groups are in one category, merges multiple inputs
    
    def _process_read( self, read, seq_name, region_start, region_end ):
        region_overlap = 0
        if read.get_reference_name() == seq_name:
            position = read.get_position( one_based=False )
            end_position = read.get_end_position( one_based=False )
            if position == end_position:
                region_overlap = -1 # this read has no length, most likely due to not having a cigar available (*)
            elif read.get_mapq() >= self._min_mapping_quality:
                region_overlap, start, end = get_region_overlap_with_positions( ( position, end_position ), ( region_start, region_end ) ) #can speed this up by not using get_region_overlap method above, instead doing it in-line
                if region_overlap:
                    if self._read_groups:
                        rg_name = read.get_read_group()
                    else:
                        rg_name = None #if no readgroups in header, assume no read groups in read
                    if rg_name:
                        if self._use_strand:
                            if read.is_seq_reverse_complement():
                                coverage = self._read_group_coverage_reverse.get( rg_name, None )
                                if coverage is None:
                                    coverage = self._read_group_coverage_reverse[ rg_name ] = NucleotideCoverage( self._sequence_lengths[ seq_name ], dtype=self._dtypes[ read.get_reference_id() ] or self.__DEFAULT_NUMPY_DTYPE__, safe=self._safe )
                            else:
                                coverage = self._read_group_coverage.get( rg_name, None )
                                if coverage is None:
                                    coverage = self._read_group_coverage[ rg_name ] = NucleotideCoverage( self._sequence_lengths[ seq_name ], dtype=self._dtypes[ read.get_reference_id() ] or self.__DEFAULT_NUMPY_DTYPE__, safe=self._safe )
                        else:
                            coverage = self._read_group_coverage.get( rg_name, None )
                            if coverage is None:
                                coverage = self._read_group_coverage[ rg_name ] = NucleotideCoverage( self._sequence_lengths[ seq_name ], dtype=self._dtypes[ read.get_reference_id() ] or self.__DEFAULT_NUMPY_DTYPE__, safe=self._safe )
                    else:
                        if self._use_strand:
                            if read.is_seq_reverse_complement():
                                coverage = self._no_read_group_coverage_reverse
                                if coverage is None:
                                    coverage = self._no_read_group_coverage_reverse = NucleotideCoverage( self._sequence_lengths[ seq_name ], dtype=self._dtypes[ read.get_reference_id() ] or self.__DEFAULT_NUMPY_DTYPE__, safe=self._safe )
                            else:
                                coverage = self._no_read_group_coverage
                                if coverage is None:
                                    coverage = self._no_read_group_coverage = NucleotideCoverage( self._sequence_lengths[ seq_name ], dtype=self._dtypes[ read.get_reference_id() ] or self.__DEFAULT_NUMPY_DTYPE__, safe=self._safe )
                        else:
                            coverage = self._no_read_group_coverage
                            if coverage is None:
                                coverage = self._no_read_group_coverage = NucleotideCoverage( self._sequence_lengths[ seq_name ], dtype=self._dtypes[ read.get_reference_id() ] or self.__DEFAULT_NUMPY_DTYPE__, safe=self._safe )
                    coverage.add_read( read, min_base_quality=self._min_base_quality ) #coverage will be determined and stored for non-overlapping bits here
                    self._covered_regions.add_region( ( seq_name, start, end ) ) #store overlap for faster iteration over coverage
                elif end_position <= region_start:
                    region_overlap = -1 # ends before region starts, since sorted, we'll keep processing BAM
            else:
                region_overlap = -1 #-1 indicates that overlap was not checked
        return region_overlap
    
    def _iter_current_coverage( self ):
        for region_name, region_start, region_end in self._covered_regions.iter_covered_regions():
            for i in xrange( region_start, region_end ):
                coverage = odict()
                coverage_reverse = odict()
                if self._read_group_coverage:
                    for rg_name in self._read_groups.keys():
                        cov = self._read_group_coverage.get( rg_name, None )
                        if cov:
                            cov_val = cov.get( i, insertions=False, deletions=False )
                            cov_val.update( cov.get( i + self.__INDEL_OFFSET__, nucleotides=False ) )
                            if cov_val:
                                coverage[ rg_name ] = cov_val
                        cov_reverse = self._read_group_coverage_reverse.get( rg_name, None )
                        if cov_reverse:
                            cov_reverse_val = cov_reverse.get( i, insertions=False, deletions=False )
                            cov_reverse_val.update( cov_reverse.get( i + self.__INDEL_OFFSET__, nucleotides=False ) )
                            if cov_reverse_val:
                                coverage_reverse[ rg_name ] = cov_reverse_val
                if self._no_read_group_coverage:
                    cov = self._no_read_group_coverage.get( i, insertions=False, deletions=False )
                    cov.update( self._no_read_group_coverage.get( i + self.__INDEL_OFFSET__, nucleotides=False ) )
                    if cov:
                        coverage[ SAM_NO_READ_GROUP_NAME ] = cov
                if self._no_read_group_coverage_reverse:
                    cov = self._no_read_group_coverage_reverse.get( i, insertions=False, deletions=False )
                    cov.update( self._no_read_group_coverage_reverse.get( i + self.__INDEL_OFFSET__, nucleotides=False ) )
                    if cov:
                        coverage_reverse[ SAM_NO_READ_GROUP_NAME ] = cov
                if coverage or coverage_reverse:
                    yield region_name, i, ( coverage, coverage_reverse )
    
    def iter_coverage( self ):
        for seq_name, start, end in self._restrict_regions.iter_covered_regions():
            for reader in self._readers:
                bam_read = reader.jump( seq_name, start, next=True )
                if bam_read:
                    overlap = self._process_read( bam_read, seq_name, start, end )
                    while overlap:
                        try:
                            bam_read = reader.next()
                        except StopIteration:
                            break
                        overlap = self._process_read( bam_read, seq_name, start, end )
            
            for v in self._iter_current_coverage():
                yield v
            #reset values
            self._read_group_coverage = {}#odict()
            self._no_read_group_coverage = None
            self._read_group_coverage_reverse = {}#odict()
            self._no_read_group_coverage_reverse = None
            #CONFIRM NOT NEED TO RESET THIS: self._restrict_regions = NamedRegionOverlap( self._sequence_lengths, regions=self._restrict_regions_list ) # reset to initial values
            self._covered_regions = NamedRegionOverlap( self._sequence_lengths )
    def iter_vcf( self, ploidy=2, variants_only=False, command_line=None, info_fields=['AC', 'AF'], format_fields=['GT', 'AC', 'AF', 'NC']  ):
        #http://www.1000genomes.org/wiki/Analysis/Variant%20Call%20Format/vcf-variant-call-format-version-41
        # move these consts
        #output header
        header = '##fileformat=%s\n##fileDate=%s\n##source=%s version %s\n##reference=file://%s\n' % ( 'VCFv4.1', date.today().strftime( '%Y%m%d' ), PROGRAM_NAME, PROGRAM_VERSION, self._reference_sequence_filename )
        if 'AC' in info_fields:
            header += '##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes, for each ALT allele, in the same order as listed">\n'
        if 'AF' in info_fields:
            header += '##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency, for each ALT allele, in the same order as listed">\n'
        #how much of header is needed?
        if 'GT' in format_fields:
            header += '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        if 'AC' in format_fields:
            header += '##FORMAT=<ID=AC,Number=.,Type=Integer,Description="Allele count in genotypes, for each ALT allele, in the same order as listed">\n'
        if 'AF' in format_fields:
            header += '##FORMAT=<ID=AF,Number=.,Type=Float,Description="Allele Frequency, for each ALT allele, in the same order as listed">\n'
        if 'NC' in format_fields:
            header += '##FORMAT=<ID=NC,Number=.,Type=String,Description="Nucleotide and indel counts">\n'
        
        if command_line:
            header += '##COMMAND=%s' % ( command_line )
        
        header += '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT'
        samples = list( self._sample_names )
        if self._has_no_read_groups:
            samples.insert( 0, SAM_NO_READ_GROUP_NAME )
        for sample in samples:
            header = "%s\t%s" % ( header, sample )
        yield header
        alt_lambda = lambda x: x[0]
        vcf_id = VCF_NO_VALUE #move id back into loop when it can be provided via e.g. external file
        #move qual and filter back into loop when they are calculated / provided
        qual = VCF_NO_VALUE 
        filter = VCF_NO_VALUE 
        _calculate_allele_coverage = self._calculate_allele_coverage
        _get_ref_allele_for_position = self._get_ref_allele_for_position
        _coverage_by_sample_from_dict = self._coverage_by_sample_from_dict
        _get_info_format_value_fields = self._get_info_format_value_fields
        for seq_name, pos, ( nucs, nucs_reverse ) in self.iter_coverage():
            #vcf_id = VCF_NO_VALUE #add id back here, when list of ids can be provided via e.g. external file
            ref = _get_ref_allele_for_position( seq_name, pos, die_on_error = not self._allow_out_of_bounds_positions )
            indeled_ref, alt_tuples = _calculate_allele_coverage( nucs.values(), nucs_reverse.values(), ref, position=pos, sequence_name=seq_name, skip_list=ref, min_support_depth=self._min_support_depth )
            #for now, lets to a special check for anything not having a string for a name:
            # FIXME: clean this up after confirming fix
            #alt = []
            #for alt_name, alt_count in alt_tuples:
            #    if isinstance( alt_name, basestring ):
            #        alt.append( alt_name )
            #    else:
            #        print >>stderr, 'An alternate allele of invalid name "%s" was encountered having count "%s" at position "%s" and has been removed.' % ( alt_name, alt_count, pos )
            alt = map( alt_lambda, alt_tuples )
            
            if variants_only:
                no_alt = True
                for an_alt in alt:
                    if an_alt != VCF_NO_VALUE:
                        no_alt = False
                        break
                if no_alt:
                    continue #no alts here, don't return region. TODO: Make this configurable
            #qual = VCF_NO_VALUE 
            #filter = VCF_NO_VALUE 
            fields = [ seq_name, str( pos + 1 ), vcf_id, indeled_ref, COMMA_CHAR.join( alt ) or VCF_NO_VALUE, qual, filter ]
            
            samples_dict = _coverage_by_sample_from_dict( samples, nucs, ref )
            if nucs_reverse:
                sample_dict_reverse = _coverage_by_sample_from_dict( samples, nucs_reverse, ref )
            else:
                sample_dict_reverse = {}#odict()#real value is odict, but dict is fine when not 'real'
            
            info, format = _get_info_format_value_fields( samples, samples_dict, sample_dict_reverse, ref, alt, indeled_ref, info_fields=info_fields, format_fields=format_fields, ploidy=ploidy, min_support_depth=self._min_support_depth )
            
            
            fields.append( VCF_INFO_FIELD_SEPARATOR.join( info ) or VCF_NO_VALUE )
            if format_fields: #format names list
                fields.append( VCF_FORMAT_FIELD_SEPARATOR.join( format_fields ) )
                if format:
                    for sample_format in format:
                        fields.append( VCF_FORMAT_FIELD_SEPARATOR.join( sample_format ) )
                else:
                    for sample in samples:
                        fields.append( VCF_FORMAT_FIELD_SEPARATOR.join([ VCF_NO_VALUE ] * len( format_fields ) ) )
            else:
                fields.append( VCF_NO_VALUE )
                for sample in samples:
                    fields.append( VCF_NO_VALUE )
            
            yield TAB_CHAR.join( fields )
    def _coverage_by_sample_from_dict( self, samples, rg_dict, reference_nucleotide ):
        rval = odict()
        for sample in samples:
            sample_dict = {}
            for rg_name, rg in rg_dict.iteritems():
                if rg_name in self._read_groups and self._read_groups[rg_name][SAM_READ_GROUP_SAMPLE_STR] == sample:
                    coverage = rg_dict[ rg_name ]
                elif rg_name not in self._read_groups and rg_name == sample:
                    coverage = rg_dict[ SAM_NO_READ_GROUP_NAME ]
                else:
                    continue
                for name, value in coverage.iteritems():                    
                    if isinstance( value, list ):
                        variants = {}
                        for v in value:
                            #what todo for deletions
                            if name == INSERTIONS_KWD:
                                v = reference_nucleotide + v
                            variants[v] = variants.get( v, 0 ) + 1
                        for v, cov in variants.iteritems():
                            sample_dict[v] = sample_dict.get( v, 0 ) + cov
                        
                    else:
                        sample_dict[name] = sample_dict.get( name, 0 ) + value
            rval[sample] = sample_dict
        return rval
    
    def _rework_indels( self, coverage, reference_nucleotide, insertion_count, deletion_count, position, sequence_name ):
        max_delete_size = 0
        if deletion_count:
            delete_remaining = deletion_count
            for name, value in coverage:
                if not delete_remaining:
                    break
                if is_integer( name ):
                    delete_remaining -= 1
                    max_delete_size = max( max_delete_size, name )
        len_old_reference_nucleotide = len( reference_nucleotide )
        if max_delete_size:
            reference_nucleotide = self._get_ref_allele_for_position( sequence_name, position, length=max_delete_size+1, die_on_error = not self._allow_out_of_bounds_positions )
        rval = []
        for name, value in coverage:
            if is_integer( name ):
                #deletion
                name = reference_nucleotide[:len_old_reference_nucleotide] + reference_nucleotide[ len_old_reference_nucleotide + name :] 
            elif len( name ) > 1:
                #insertion
                name = reference_nucleotide[:len_old_reference_nucleotide] + name[len_old_reference_nucleotide:] + reference_nucleotide[len_old_reference_nucleotide:]
            else:
                name = name + reference_nucleotide[len_old_reference_nucleotide:]
            rval.append( ( name, value ) )
        return reference_nucleotide, rval
    
    def _calculate_allele_coverage( self, coverages, coverages_reverse, reference_nucleotide, position=None, sequence_name=None, skip_list = None, min_support_depth=None ):
        if not isinstance( coverages, list ):
            coverages = [ coverages ]
        else:
            coverages = list( coverages )
        if coverages_reverse:
            if not isinstance( coverages_reverse, list ):
                coverages_reverse = [ coverages_reverse ]
            coverages.extend( coverages_reverse )
        if not isinstance( skip_list, list ):
            if skip_list is None:
                skip_list = []
            else:
                skip_list = [skip_list]
        coverage_dict = {}
        insertion_count = 0
        deletion_count = 0
        for coverage in coverages:
            for nuc, cov in coverage.iteritems():
                if isinstance( cov, list ):
                    #this is an indel
                    variants = {}
                    for v in cov:
                        #what todo for deletions
                        if nuc == INSERTIONS_KWD:
                            v = reference_nucleotide + v
                            if v not in variants:
                                insertion_count += 1
                        elif nuc == DELETIONS_KWD:
                            if v not in variants:
                                deletion_count += 1
                        else:
                            print >>stderr, "An unreachable condition has been reached in _calculate_allele_coverage at position '%s'. coverages=%s. coverages_reverse=%s" % ( position, coverages, coverages_reverse )
                        variants[v] = variants.get( v, 0 ) + 1
                    for v, cov in variants.iteritems():
                        coverage_dict[v] = coverage_dict.get( v, 0 ) + cov
                else:
                    #this is a standard nucleotide
                    coverage_dict[ nuc ] = coverage_dict.get( nuc, 0 ) + cov
        #filter by depth
        if min_support_depth:
            coverage_dict2 = {}
            for name, value in coverage_dict.iteritems():
                if value >= min_support_depth:
                    coverage_dict2[name] = value
                elif is_integer( name ):
                    deletion_count -= 1
                elif len( name ) > 1:
                    insertion_count -= 1
            #coverage_dict = dict( filter( lambda x: x[1] >= min_support_depth, coverage_dict.iteritems() ) )
            coverage_dict = coverage_dict2
            del coverage_dict2
        
        #do some filter on min coverage?
        coverage = [ ( x[1], x[0] ) for x in sorted( map( lambda x: ( x[1], x[0] ), coverage_dict.iteritems() ), reverse=True ) if x[1] not in skip_list ] #is ordering important?
        if position is not None and sequence_name and ( insertion_count or deletion_count ):
            reference_nucleotide, coverage = self._rework_indels( coverage, reference_nucleotide, insertion_count, deletion_count, position, sequence_name )
        return ( reference_nucleotide, coverage )
        
    def _get_ref_allele_for_position( self, sequence_name, position, length=1, die_on_error=True ):
        return self._reference_sequences.get_sequence_by_position( sequence_name, position, length=length, unknown_sequence_character=VCF_NO_VALUE, die_on_error=die_on_error )
    def _get_info_format_value_fields( self, samples, samples_dict, samples_dict_reverse, ref, alt, indeled_ref, info_fields=['AC', 'AF'], format_fields=['GT', 'AC', 'AF', 'NC'], ploidy=2, min_support_depth=None ):
        gt_no_value = VCF_GT_PHASED_SEPARATOR[False].join( [ VCF_NO_VALUE for i in range( ploidy ) ] )
        ac_af_no_value = COMMA_CHAR.join( [ VCF_NO_VALUE for i in range( len( alt ) ) ] )
        ref_list = [indeled_ref]
        if not isinstance( alt, list ):
            alt_list = [alt]
        else:
            alt_list = alt
        all_alt_nucs_count = [0] * len( alt_list )
        all_nucs_count = 0
        format = []
        
        #move these 3 out, for faster?
        len_old_reference_nucleotide = len( ref )
        len_reference_nucleotide = len( indeled_ref )
        max_deletion_len = len_reference_nucleotide - len_old_reference_nucleotide
        
        for sample in samples:
            sample_format = []
            if sample in samples_dict or sample in samples_dict_reverse:
                coverage_dict = samples_dict.get( sample, {} )
                coverage_dict_reverse = samples_dict_reverse.get( sample, {} )
                
                genotyping_coverage_dict = dict( coverage_dict )
                for nuc, count in coverage_dict_reverse.iteritems():
                    genotyping_coverage_dict[ nuc ] = genotyping_coverage_dict.get( nuc, 0 ) + count
                
                if min_support_depth is None:
                    min_support_depth = 0
                nucs = sorted( map( lambda x: ( x[1], x[0] ), genotyping_coverage_dict.iteritems() ), reverse=True )#list of count, nuc
                #nucs_reverse = sorted( map( lambda x: ( x[1], x[0] ), coverage_dict_reverse.iteritems() ), reverse=True )#list of count, nuc
                #nucs = filter( lambda x: not isinstance( x[0], list ), nucs ) #filter out indels for now. FIXME: add sort by len or int count to allow indels here
                #resolve indels
                indeled_nucs = []
                
                
                
                #### below here should make a combined dict for genotyping here, then do both for the raw nuc counts
                
                nucs_sum = 0
                for value, name in nucs:
                    if value >= min_support_depth:
                        if is_integer( name ):
                            #deletion
                            name = indeled_ref[:len_old_reference_nucleotide] + indeled_ref[ len_old_reference_nucleotide + name :]
                        elif len( name ) > 1:
                            #insertion
                            name = indeled_ref[:len_old_reference_nucleotide] + name[len_old_reference_nucleotide:] + indeled_ref[len_old_reference_nucleotide:]
                        else:
                            name = name + indeled_ref[len_old_reference_nucleotide:]
                        indeled_nucs.append( ( value, name) )
                        nucs_sum += value
                
                all_nucs_count += nucs_sum
            
                
    
                
                #format
                #AC and AF
                ac_list = [ 0 for i in range( len( alt_list ) ) ]
                af_list = list( ac_list )
                for value, name in indeled_nucs:
                    try:
                        i = alt_list.index( name )
                        ac_list[ i ] = int( value )
                        all_alt_nucs_count[ i ] = all_alt_nucs_count[ i ] + ac_list[ i ] #needed for info field
                        af_list[ i ] = float( value ) / nucs_sum
                    except ValueError:
                        continue #this is ref allele
                #NC
                nc_field = "" #"%s:" % fields[-1]
                if self._use_strand:
                    prefix = '+'
                else:
                    prefix = ''
                for c, count in coverage_dict.iteritems():
                    if count:#should we filter by self._min_support_dept here? or display all
                        if is_integer( c ):
                            c = 'd%s' % ( c )
                        nc_field = "%s%s%s=%d," % ( nc_field, prefix, c, count )
                prefix = '-'
                for c, count in coverage_dict_reverse.iteritems():
                    if count:#should we filter by self._min_support_dept here? or display all
                        if is_integer( c ):
                            c = 'd%s' % ( c )
                        nc_field = "%s%s%s=%d," % ( nc_field, prefix, c, count )
                
                gt_possible = ref_list + alt_list
                calls = []
                len_nucs = len( indeled_nucs )
                count_per_allele = float( nucs_sum ) / ploidy
                nucs = map( lambda x: ( x[0] / count_per_allele, x[1] ), indeled_nucs )
                remain_ploidy = ploidy
                while remain_ploidy > 0:
                    nucs = sorted( nucs, reverse=True )
                    #Need to handle case where we have e.g. 2 coverage in A and G
                    #check called alleles list when read counts match
                    if len_nucs == 0:
                        nuc = None
                    else:
                        for i, ( weight, nuc ) in enumerate( nucs ):
                            j = i+1
                            if nuc not in calls or j >= len_nucs or weight > nucs[j][0] :
                                break
                        nucs[i] = ( nucs[i][0] - 1, nucs[i][1] ) 
                    calls.append( nuc )
                    remain_ploidy -= 1
                gt_list = []
                for nuc in calls:
                    try:
                        call = gt_possible.index( nuc )
                    except ValueError:
                        call = VCF_NO_VALUE #value is not one of the ref or alt alleles
                    gt_list.append( call )
                
                sample_format.append( VCF_GT_PHASED_SEPARATOR[False].join( map( str, gt_list ) ) )
                if 'AC' in format_fields:
                    sample_format.append( COMMA_CHAR.join( map( str, ac_list ) ) )
                if 'AF' in format_fields:
                    sample_format.append( COMMA_CHAR.join( map( str, af_list ) ) )
                #NC
                if 'NC' in format_fields:
                    sample_format.append( nc_field )
            else:
                sample_format.append( gt_no_value )
                if 'AC' in format:
                    sample_format.append( ac_af_no_value ) #AC
                if 'AF' in format:
                    sample_format.append( ac_af_no_value ) #AF
                sample_format.append( VCF_NO_VALUE )
            format.append( sample_format )
        info = []
        ac = COMMA_CHAR.join( map( str, all_alt_nucs_count ) )
        if 'AC' in info_fields:
            info.append( "AC=%s" % ( ac ) )
        if 'AF' in info_fields:
            if all_nucs_count:
                af = COMMA_CHAR.join( map( lambda x: str( float( x ) / float( all_nucs_count ) ), all_alt_nucs_count ) )
            else:
                af = ac
            info.append( "AF=%s" % ( af ) )
        return info, format
        
class VCFReadGroupGenotyper( ReadGroupGenotyper ):
    __INDEL_OFFSET__ = 1

