#Dan Blankenberg
import string

import numpy

from pyBamParser.fasta import IndexedReferenceSequences

from ..util.odict import odict

from ..util import MAX_NEG_INT

from ..util import NUMPY_DTYPES

from ..util.sam import SAM_HEADER_NON_TAB_RECORDS
from ..util.sam import SAM_READ_GROUP_RECORD_CODE
from ..util.sam import SAM_NO_READ_GROUP_NAME

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

class ReadGroupGenotyper( object ):
    __INDEL_OFFSET__ = 0
    __DEFAULT_NUMPY_DTYPE__ = NUMPY_DTYPES.get( 'uint8' )
    
    #TODO: make these filters defined generically
    def __init__( self, bam_readers=None, reference_sequence_filename=None, dtype=None, min_support_depth=None, 
                  min_base_quality=None, min_mapping_quality=None, restrict_regions=None, use_strand=False ):
        self._read_groups = odict() #order important?
        self._sample_names = []
        self._sequence_lengths = odict()
        self._read_group_coverage = {}#odcit()
        self._read_group_coverage_reverse = {}
        self._no_read_group_coverage = None#odcit()
        self._no_read_group_coverage_reverse = None
        self._has_no_read_groups = False
        self._reference_name = None
        self._reference_sequence_filename = reference_sequence_filename
        self._reference_sequences = IndexedReferenceSequences( reference_sequence_filename, sequence_filter=string.upper )
        if dtype and isinstance( dtype, basestring ):
            dtype = NUMPY_DTYPES.get( dtype, self.__DEFAULT_NUMPY_DTYPE__ )
        if not dtype:
            dtype = self.__DEFAULT_NUMPY_DTYPE__
        self._dtype = dtype
        if bam_readers is None:
            bam_readers = []
        elif not isinstance( bam_readers, list ):
            bam_readers = [ bam_readers ]
        self._readers = []
        for bam_reader in bam_readers:
            self.add_reader( bam_reader ) #(reader, buffered block)
        self._use_strand = use_strand
        #filters:
        self._min_support_depth = min_support_depth or 0
        self._min_base_quality = min_base_quality
        self._min_mapping_quality = min_mapping_quality
        if self._min_mapping_quality is None:
            self._min_mapping_quality = MAX_NEG_INT
        self._restrict_regions_list = restrict_regions
        self._restrict_regions = NamedRegionOverlap( self._sequence_lengths, regions=self._restrict_regions_list )
        self._covered_regions = NamedRegionOverlap( self._sequence_lengths )
            
        
    def add_reader( self, bam_reader ):
        self._readers.append( ( bam_reader, None ) ) #(reader, buffered read)
        for seq_name, seq_len in bam_reader.get_references():
            if seq_name not in self._sequence_lengths:
                self._sequence_lengths[ seq_name ] = seq_len
                #print 'found', seq_name, seq_len
            else:
                assert self._sequence_lengths[ seq_name ] == seq_len, "Sequence length mismatch for '%s': %s != %s" % ( seq_name, self._sequence_lengths[ seq_name ], seq_len ) #TODO: print out offending reader name
        #add read group information
        headers = bam_reader.get_sam_header_dict()
        if SAM_READ_GROUP_RECORD_CODE in headers:
            #print 'rg', headers[SAM_READ_GROUP_RECORD_CODE]
            for rg in headers[SAM_READ_GROUP_RECORD_CODE]: #{'LB': 'reads_library', 'ID': 'read_group_id', 'PL': 'ILLUMINA', 'SM': 'read_sample'}
                rg_name = rg['ID']
                if rg_name not in self._read_groups:
                    self._read_groups[ rg_name ] = rg #TODO make get_readgroups method
                    #TODO: make sample name default to read group id
                    if rg['SM'] not in self._sample_names:
                        self._sample_names.append( rg['SM'] )
                else:
                    assert self._read_groups[ rg_name ] == rg, "Read Group info mismatch: '%s' != '%s'." % ( self._read_groups[ rg_name ], rg )
        else:
            self._has_no_read_groups = True #for now all non-read groups are in one category, merges multiple inputs
    def _process_read( self, read ):
        position = read.get_position( one_based=False )
        end_position = read.get_end_position( one_based=False )
        if self._restrict_regions.region_overlaps( self._reference_name, position, end_position, empty_default=True ):
            if read.get_mapq >= self._min_mapping_quality:
                rg_name = read.get_read_group()
                if rg_name:
                    if self._use_strand:
                        if read.is_seq_reverse_complement():
                            coverage = self._read_group_coverage_reverse.get( rg_name, None )
                            if coverage is None:
                                coverage = self._read_group_coverage_reverse[ rg_name ] = NucleotideCoverage( self._sequence_lengths[ self._reference_name ], dtype=self._dtype, indel_offset=self.__INDEL_OFFSET__ )
                        else:
                            coverage = self._read_group_coverage.get( rg_name, None )
                            if coverage is None:
                                coverage = self._read_group_coverage[ rg_name ] = NucleotideCoverage( self._sequence_lengths[ self._reference_name ], dtype=self._dtype, indel_offset=self.__INDEL_OFFSET__ )
                    else:
                        coverage = self._read_group_coverage.get( rg_name, None )
                        if coverage is None:
                            coverage = self._read_group_coverage[ rg_name ] = NucleotideCoverage( self._sequence_lengths[ self._reference_name ], dtype=self._dtype, indel_offset=self.__INDEL_OFFSET__ )
                else:
                    if self._use_strand:
                        if read.is_seq_reverse_complement():
                            coverage = self._no_read_group_coverage_reverse
                            if coverage is None:
                                coverage = self._no_read_group_coverage_reverse = NucleotideCoverage( self._sequence_lengths[ self._reference_name ], dtype=self._dtype, indel_offset=self.__INDEL_OFFSET__ )
                        else:
                            coverage = self._no_read_group_coverage
                            if coverage is None:
                                coverage = self._no_read_group_coverage = NucleotideCoverage( self._sequence_lengths[ self._reference_name ], dtype=self._dtype, indel_offset=self.__INDEL_OFFSET__ )
                    else:
                        coverage = self._no_read_group_coverage
                        if coverage is None:
                            coverage = self._no_read_group_coverage = NucleotideCoverage( self._sequence_lengths[ self._reference_name ], dtype=self._dtype, indel_offset=self.__INDEL_OFFSET__ )
                coverage.add_read( read, min_base_quality=self._min_base_quality )
                
                self._covered_regions.add_region( ( self._reference_name, position, end_position ) ) #store overlap for faster iteration over coverage
    
    def _buffer_readers( self, condition=False ):
        readers = []
        for bam_reader, buffered_read in self._readers:
            empty_reader = False
            while not buffered_read or ( condition and not condition( buffered_read ) ):
                try:
                    buffered_read = bam_reader.next()
                except StopIteration:
                    empty_reader = True
                    break
            if empty_reader:
                continue
            readers.append( ( bam_reader, buffered_read ) )
        self._readers = readers
        return readers
    def _iter_current_coverage( self ):
        if self._reference_name and self._restrict_regions.region_overlaps( self._reference_name, empty_default=True ):
            for region_name, region_start, region_end in self._covered_regions.iter_covered_regions():
                print 'coverage', region_name, region_start, region_end
                for i in xrange( region_start, region_end ):
                    coverage = odict()
                    coverage_reverse = odict()
                    if self._read_group_coverage:
                        for rg_name in self._read_groups.keys():
                            cov = self._read_group_coverage.get( rg_name, None )
                            if cov:
                                cov = cov.get( i )
                                if cov:
                                    coverage[ rg_name ] = cov
                            cov_reverse = self._read_group_coverage_reverse.get( rg_name, None )
                            if cov_reverse:
                                cov_reverse = cov_reverse.get( i )
                                if cov_reverse:
                                    coverage_reverse[ rg_name ] = cov_reverse
                    if self._no_read_group_coverage:
                        cov = self._no_read_group_coverage.get( i )
                        if cov:
                            coverage[ SAM_NO_READ_GROUP_NAME ] = cov
                    if self._no_read_group_coverage_reverse:
                        cov = self._no_read_group_coverage_reverse.get( i )
                        if cov:
                            coverage_reverse[ SAM_NO_READ_GROUP_NAME ] = cov
                    if coverage or coverage_reverse:
                        yield self._reference_name, i, ( coverage, coverage_reverse )
    def iter_coverage( self ):
        sequence_names = self._sequence_lengths.keys()
        read_no_ref_filter = lambda x: True if x.get_reference_id() >= 0 else False
        while self._readers:
            buffered_readers = self._buffer_readers( condition=read_no_ref_filter )
            if not buffered_readers:
                break
            reference_name_index = len( sequence_names )
            for bam_reader, read in buffered_readers:
                reference_name_index = min( reference_name_index, sequence_names.index( read.get_reference_name() ) )
            next_reference_name = sequence_names[ reference_name_index ]
            if next_reference_name != self._reference_name:
                for v in self._iter_current_coverage():
                    yield v
                #reset values
                self._read_group_coverage = {}#odict()
                self._no_read_group_coverage = None
                self._read_group_coverage_reverse = {}#odict()
                self._no_read_group_coverage_reverse = None
                #CONFIRM NOT NEED TO RESET THIS: self._restrict_regions = NamedRegionOverlap( self._sequence_lengths, regions=self._restrict_regions_list ) # reset to initial values
                self._covered_regions = NamedRegionOverlap( self._sequence_lengths )
            self._reference_name = next_reference_name
            used_reads = []
            for i, ( bam_reader, read ) in enumerate( buffered_readers ):
                if read.get_reference_name() == self._reference_name:
                    used_reads.append( i )
                    self._process_read( read )
            for i in used_reads:
                self._readers[i] = ( self._readers[i][0], None )
        #flush remaining coverage
        for v in self._iter_current_coverage():
            yield v
        #reset values
        self._read_group_coverage = {}#odict()
        self._no_read_group_coverage = None
        self._read_group_coverage_reverse = {}#odict()
        self._no_read_group_coverage_reverse = None
    def iter_vcf( self, ploidy=2, variants_only=False, command_line=None, info_fields=['AC','AF'], format_fields=['GT', 'AC', 'AF', 'NC']  ):
        #http://www.1000genomes.org/wiki/Analysis/Variant%20Call%20Format/vcf-variant-call-format-version-41
        # move these consts
        
        #output header
        header = '##fileformat=%s\n##fileDate=%s\n##source=%s\n##reference=file://%s\n' % ( 'VCFv4.1', '?0', 'Dan', self._reference_sequence_filename )
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
        for seq_name, pos, ( nucs, nucs_reverse ) in self.iter_coverage():
            id = VCF_NO_VALUE
            ref = self._get_ref_allele_for_position( seq_name, pos )
            indeled_ref, alt_tuples = self._calculate_allele_coverage( nucs.values(), nucs_reverse.values() , ref, position=pos, sequence_name=seq_name, skip_list=ref, min_support_depth=self._min_support_depth )
            alt = map( lambda x: x[0], alt_tuples )
            if variants_only:
                no_alt = True
                for an_alt in alt:
                    if an_alt != VCF_NO_VALUE:
                        no_alt = False
                        break
                if no_alt:
                    continue #no alts here, don't return region. TODO: Make this configurable
            qual = VCF_NO_VALUE 
            filter = VCF_NO_VALUE 
            fields = [ seq_name, str( pos + 1 ), id, indeled_ref, ','.join( alt ) or VCF_NO_VALUE, qual, filter ]
            
            samples_dict = self._coverage_by_sample_from_dict( samples, nucs, ref )
            if nucs_reverse:
                sample_dict_reverse = self._coverage_by_sample_from_dict( samples, nucs_reverse, ref )
            else:
                sample_dict_reverse = odict()
            
            info, format = self._get_info_format_value_fields( samples, samples_dict, sample_dict_reverse, ref, alt, indeled_ref, info_fields=info_fields, format_fields=format_fields, ploidy=ploidy, min_support_depth=self._min_support_depth )
            
            
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
            
            yield "\t".join( fields )
    def _coverage_by_sample_from_dict( self, samples, rg_dict, reference_nucleotide ):
        rval = odict()
        for sample in samples:
            sample_dict = {}
            for rg_name, rg in rg_dict.iteritems():
                if rg_name in self._read_groups and self._read_groups[rg_name]['SM'] == sample:
                    coverage = rg_dict[ rg_name ]
                elif rg_name not in self._read_groups and rg_name == sample: #rg_name == sample == SAM_NO_READ_GROUP_NAME:
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
#                                has_insertion += True
                            #elif name == DELETIONS_KWD:
                                #has_deletion += True
                            #if nuc == __DELETIONS_KWD__:
                            #    v = 'd:%v' % v
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
                if isinstance( name, int ):
                    delete_remaining -= 1
                    max_delete_size = max( max_delete_size, name )
        len_old_reference_nucleotide = len( reference_nucleotide )
        if max_delete_size:
            reference_nucleotide = self._get_ref_allele_for_position( sequence_name, position, length=max_delete_size+1 )
        rval = []
        for name, value in coverage:
            if isinstance( name, int ):
                #deletion
                name = reference_nucleotide[:len_old_reference_nucleotide] + reference_nucleotide[ len_old_reference_nucleotide + name :] 
            elif len( name ) > 1:
                #insertion
                name = reference_nucleotide[:len_old_reference_nucleotide] + name[len_old_reference_nucleotide:] + reference_nucleotide[len_old_reference_nucleotide:]
            else:
                name = name + reference_nucleotide[len_old_reference_nucleotide:]
            rval.append( ( name, value) )
        return reference_nucleotide, rval
    
    def _calculate_allele_coverage( self, coverages, coverages_reverse, reference_nucleotide, position=None, sequence_name=None, skip_list = None, min_support_depth=None ):
        if not isinstance( coverages, list ):
            coverages = [ coverages ]
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
                    variants = {}
                    for v in cov:
                        #what todo for deletions
                        if nuc == INSERTIONS_KWD:
                            v = reference_nucleotide + v
                            if v not in variants:
                                insertion_count += True
                        elif nuc == DELETIONS_KWD:
                            if v not in variants:
                                deletion_count += True
                        #if nuc == __DELETIONS_KWD__:
                        #    v = 'd:%v' % v
                        variants[v] = variants.get( v, 0 ) + 1
                    for v, cov in variants.iteritems():
                        coverage_dict[v] = coverage_dict.get( v, 0 ) + cov
                else:
                #if not isinstance( cov, list ): #remove indels
                    coverage_dict[nuc] = coverage_dict.get( nuc, 0 ) + cov
        #filter by depth
        if min_support_depth:
            coverage_dict2 = {}
            for name, value in coverage_dict.iteritems():
                if value >= min_support_depth:
                    coverage_dict2[name] = value
                elif isinstance( name, int ):
                    deletion_count -= 1
                elif len( name ) > 1:
                    insertion_count -= 1
            #coverage_dict = dict( filter( lambda x: x[1] >= min_support_depth, coverage_dict.iteritems() ) )
            coverage_dict = coverage_dict2
            del coverage_dict2
        #do some filter on min coverage?
        coverage = [ ( x[1], x[0] ) for x in sorted( map( lambda x: ( x[1], x[0] ), coverage_dict.iteritems() ), reverse=True ) if x[1] not in skip_list ] #is ordering important?
        if position and sequence_name and ( insertion_count or deletion_count ):
            reference_nucleotide, coverage = self._rework_indels( coverage, reference_nucleotide, insertion_count, deletion_count, position, sequence_name )
        return ( reference_nucleotide, coverage )
        
    def _get_ref_allele_for_position( self, sequence_name, position, length=1 ):
        return self._reference_sequences.get_sequence_by_position( sequence_name, position, length=length, unknown_sequence_character=VCF_NO_VALUE )
    def _get_info_format_value_fields( self, samples, samples_dict, samples_dict_reverse, ref, alt, indeled_ref, info_fields=['AC','AF'], format_fields=['GT', 'AC', 'AF', 'NC'], ploidy=2, min_support_depth=None ):
        gt_no_value = VCF_GT_PHASED_SEPARATOR[False].join( [ VCF_NO_VALUE for i in range( ploidy ) ] )
        ac_af_no_value = ",".join( [ VCF_NO_VALUE for i in range( len( alt ) ) ] )
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
                        if isinstance( name, int ):
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
                        if isinstance( c, int ):
                            c = '%sd%s' % ( prefix, c )
                        nc_field = "%s%s%s-%s," % ( nc_field, prefix, c, count )
                prefix = '-'
                for c, count in coverage_dict_reverse.iteritems():
                    if count:#should we filter by self._min_support_dept here? or display all
                        if isinstance( c, int ):
                            c = '%sd%s' % ( prefix, c )
                        nc_field = "%s%s%s-%s," % ( nc_field, prefix, c, count )
                
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
                
                #print 'rval', rval
                sample_format.append( VCF_GT_PHASED_SEPARATOR[False].join( map( str, gt_list ) ) )
                if 'AC' in format_fields:
                    sample_format.append( ",".join( map( str, ac_list ) ) )
                if 'AF' in format_fields:
                    sample_format.append( ",".join( map( str, af_list ) ) )
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
        ac = ','.join( map( str, all_alt_nucs_count ) )
        if 'AC' in info_fields:
            info.append( "AC=%s" % ( ac ) )
        if 'AF' in info_fields:
            if all_nucs_count:
                af = ','.join( map( lambda x: str( float( x ) / float( all_nucs_count ) ), all_alt_nucs_count ) )
            else:
                af = ac
            info.append( "AF=%s" % ( af ) )
        return info, format
        
class VCFReadGroupGenotyper( ReadGroupGenotyper ):
    __INDEL_OFFSET__ = -1

