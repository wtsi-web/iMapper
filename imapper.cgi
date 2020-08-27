#!/usr/local/bin/perl -T
# imapper.cgi
# 
# Extracts the tr -tagged genomic sequence from uploaded insert sequences
# and maps to genome.
#
# Author:        jk4
# Maintainer:    jk4
# Created:       2008-01-29


#     Copyright (C) 2008  Kong, J., Zhu, F., Stalker, J. and Adams, D.J.

#     This program is free software: you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.

#     This program is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.

#     You should have received a copy of the GNU General Public License
#     along with this program.  If not, see <https://www.gnu.org/licenses/>.




use strict;
use warnings;
no warnings 'uninitialized';
use SangerPaths qw(core ensembl68 bioperl123 team113);


use SangerWeb;
use Bio::EnsEMBL::Registry;

use Storable qw(freeze thaw);

use Ssaha_wrapper;

$ENV{PATH}= '/usr/local/bin';   # solely to stop taint from barfing

$|++;


our $VERSION = do { my @r = (q$Revision: 1.13 $ =~ /\d+/mxg); sprintf '%d.'.'%03d' x $#r, @r };
our $SEQ_LIMIT = 10_000;
our $SIG_COLOUR = "#CCFF00";
our $SEQ_COLOUR = "#FFFF00";
our $CUT_COLOUR = "#FF9900";
our $CONT_COLOUR = "#CC99FF";
my $javascript = <<JS ;
<!-- 
function toggleAdv() {
    var block = document.getElementById("advopt");
    var link = document.getElementById("advtoggleanchor");
    if (block.style.display == "none"){
        block.style.display = "block";
	link.firstChild.nodeValue="Hide advanced options";
    } else {
        block.style.display = "none";  
	link.firstChild.nodeValue="Show advanced options";
    }
}

function setSigtext(dropdown) {
    var sigtext = document.getElementById("signature_seq");
    sigtext.value = dropdown.options[dropdown.selectedIndex].value;
}

// -->
JS

my $css = <<CSS ;
table.imapper th {
    padding-top: 10px;
    padding-bottom: 10px;
    padding-right: 1em;
    border-bottom: #83a4c3 dotted 1px;
    text-align: right;
    vertical-align:middle;
}

table.imapper td {
    padding-top: 10px;
    padding-bottom: 10px;
    border-bottom: #83a4c3 dotted 1px;
}

table.imap_summary th {
    padding-top: 2px;
    padding-bottom: 2px;
    padding-right: 1em;
    border-bottom: #83a4c3 dotted 1px;
    text-align: right;
    vertical-align:middle;
}

table.imap_summary td {
    padding-top: 2px;
    padding-bottom: 2px;
    border-bottom: #83a4c3 dotted 1px;
}

#key {
    background-color: #CBDCED;
}

table.imap_output {
    border-collapse:collapse;
    border-spacing:0px;
    margin-bottom: 10px;
    width: 100%;
}

table.imap_output th {
    border:1px solid #22408F;
    font-weight:normal;
    padding:0.1em 0.5em;
    text-align:left;
    width: 200px;
}

table.imap_output td {
    border:1px solid #22408F;
    padding:0.1em 0.5em;
    text-align:left;
    vertical-align:top;
}


CSS


my $sw  = SangerWeb->new({
    'title'   => q(iMapper),
    'banner'  => q(),
    'inifile' => SangerWeb->document_root() . q(/Info/header.ini),
    'script'  => $javascript,
    'style'  => $css,
});

# Output dispatch table
my %output_formats = (	'tab'	=> {header   => \&output_tabular_header,
				    each_seq => \&output_tabular,
				    footer   => \&output_tabular_footer},
			'fasta'	=> {header   => \&output_fasta_header,
				    each_seq => \&output_fasta,
				    footer   => \&output_fasta_footer},
			'gff'	=> {header   => \&output_gff_header,
				    each_seq => \&output_gff,
				    footer   => \&output_gff_footer},
			'loc'	=> {header   => \&output_loc_header,
				    each_seq => \&output_loc,
				    footer   => \&output_loc_footer},
		      );
			
my %allowed_spp = ( 'Rattus norvegicus'  => 1,
		    'Mus musculus'  => 1,
		    'Homo sapiens'  => 1,
		    'Drosophila melanogaster'  => 1,
		    'Danio rerio'  => 1,
		    'Saccharomyces cerevisiae'  => 1,
		    'Gallus gallus' => 1,
		    'Sus scrofa' => 1,
		    'Bos taurus' => 1,
		    
		    
		  );

my $default_spp = 'Mus musculus';

my $cgi = $sw->cgi();

# Three entry points into the script:
# 1) post from form - upload or paste seq, all form parameters
#	This should read the data and freeze the params, then redirect to get
#	with key.
# 2) get with key to cached data which is either 
#	    - uploaded sequence and params
#	    (in which case, parse the data and run default output)
#	    - parsed data, with output type and output_option as GET params
# 3) initial entry - no parameters, empty form
#

my $upload_seq	=  $cgi->upload('upload_seq');
my $paste_seq	= $cgi->param('paste_seq');
my $key		= $cgi->param('key');
my $format	= $cgi->param('format');
$format = 'tab' unless exists $output_formats{$format};

if ($upload_seq or $paste_seq){ # entry point 1 - post from form
    # store the data and parameters, then exit
    my $fastaHash = read_fasta($upload_seq, $paste_seq);
    unless(scalar keys(%{$fastaHash->{fasta}})){
	print $sw->header();
	print_empty_form();
	print $sw->footer();
	exit;
    }


    # detaint and defaults
    my $species			    = $cgi->param('species');
    $species = $default_spp unless exists $allowed_spp{$species};
    $species =~ s/ /_/g;	# hash in Ssaha_wrapper has underscores

    my $condition	    = $cgi->param('seq_format');
    my $signature_seq	    = uc ($cgi->param('signature_seq'));
    my $restriction_seq	    = uc ($cgi->param('restriction_seq'));
    my $threshold_sig	    = $cgi->param('threshold_sig');
    my $threshold_cont	    = $cgi->param('threshold_cont');
    my $output_option	    = $cgi->param('output1');
    my $contamination_seq1  = uc $cgi->param('contamination_seq1');
    my $contamination_seq2  = uc $cgi->param('contamination_seq2');
    my $five_flank	    = $cgi->param('five_flank');
    my $three_flank	    = $cgi->param('three_flank');
    my $contamination_bp    = $cgi->param('contamination_bp');
    my $gap_penalty	    = $cgi->param('gap_penalty');
    my $match_score	    = $cgi->param('match_score');
    my $mismatch_score	    = $cgi->param('mismatch_score');
    my $format		    = $cgi->param('format');
    my $ssaha_score         = $cgi->param('ssaha_score'); 
    my $ssaha_fold         = $cgi->param('ssaha_fold');

    if ($ssaha_score =~ /(\S+)/){
	$fastaHash->{ssaha_score} = $1;
    }else {
	$fastaHash->{ssaha_score} = 35;
    }

    if ($ssaha_fold =~ /(\S+)/){
	$fastaHash->{ssaha_fold} = $1;
    }else {
	$fastaHash->{ssaha_fold} = 1.5;
    }

    if ($condition =~ /(\w+)/){
	$fastaHash->{condition} = $1;
    }else {
	$fastaHash->{condition} = "unknown";
    }

    if ($signature_seq =~ /([ATGCN]+)/){
	$fastaHash->{signature_seq} = $1;
    }

    if ($restriction_seq =~ /([ATGC]+)/){
	$fastaHash->{restriction_seq} = $1;
    }
    if ($threshold_sig =~ /(\d+)/){
	$fastaHash->{threshold_sig} = $1;
	$fastaHash->{threshold_sig} = 100 if $fastaHash->{threshold_sig} > 100;
    }
    if ($threshold_cont =~ /(\d+)/){
	$fastaHash->{threshold_cont} = $1;
	$fastaHash->{threshold_cont} = 100 if $fastaHash->{threshold_cont} > 100;
    } 

    if ($output_option =~ /(all|good)/){
	$fastaHash->{output_option} = $1;
    }

    if ($contamination_seq1 =~ /([ATGCN]+)/){
	$fastaHash->{contamination_seq1} = $1;
    }
    if ($contamination_seq2 =~ /([ATGCN]+)/){
	$fastaHash->{contamination_seq2} = $2;
    }

    if ($five_flank =~ /(\d+)/){
	$fastaHash->{five_flank} = $1;
	$fastaHash->{five_flank} = 1_000_000 if $fastaHash->{five_flank} >  1_000_000;
    }
    else {
	$fastaHash->{five_flank} = 0;
    }

    if ($three_flank =~ /(\d+)/){
	$fastaHash->{three_flank} = $1;
	$fastaHash->{three_flank} = 1_000_000 if $fastaHash->{three_flank} >  1_000_000;
    }
    else {
	$fastaHash->{three_flank} = 0;
    }

    if ($contamination_bp =~ /(\d+)/){
	$fastaHash->{contamination_bp} = $1;
    } 
    else {
	$fastaHash->{contamination_bp} = 0;
    } 

    if ($gap_penalty =~ /([-\d]+)/){
	$fastaHash->{gap_penalty} = int($1);
    } 
    if ($match_score =~ /([-\d]+)/){
	$fastaHash->{match_score} = int($1);
    } 
    if ($mismatch_score =~ /([-\d]+)/){
	$fastaHash->{mismatch_score} = int($1);
    } 
    
    if ($format =~ /(\w+)/){
	$fastaHash->{format} = $1;
    }

    $fastaHash->{species}	    = $species;

    my $frozed = freeze($fastaHash);
    my $dbstore = $sw->dbstore();
    my $key = $dbstore->set($frozed, undef, 24) # 24 hour storage
		or warn "Problem dbstoring data\n";

    # Redirect to a get with the key.  This makes the returned browser page
    # immediately bookmarkable and navigable (back buttons don't try to
    # rePOST , etc).
    #
    # In future, this is a good point to return a "your key is x, keep clicking
    # to retrieve your results" page instead, while spawning a process to do
    # the analysis in the background, rather than in-process as now.

    print $cgi->redirect( -location =>"imapper.cgi?key=$key",
			  -method   => 'GET',
			  -status   => 303);
}
elsif ($key){ # entry point 2 - cached data
    my $dbstore = $sw->dbstore();
    my $data = $dbstore->get($key);
    my $cached_data = thaw($data);

    if ($cached_data){
	# Is cached data the raw sequence (needs processing)
	# or the processed sequence (needs outputting in appropriate format)
	if ($cached_data->{fasta}){ # raw sequence

	    # analyse the data, calling print each_seq as we go
	    $output_formats{$format}{header}->();
	    my $fastadata = map_fasta($cached_data, $output_formats{$format}{each_seq}, $key);
	    $output_formats{$format}{footer}->($fastadata, $key);
	    
	    # overwrite the raw data with the processed data
	    my $frozed = freeze($fastadata);
	    $dbstore->set($frozed, $key, 24) or warn "Problem storing $key";
	}
	else { # processed data, just output it

	    $output_formats{$format}{header}->();
	    foreach my $seq (sort keys %$cached_data){
		next if $seq =~ /^__/;	# internal keys for metadata
		$output_formats{$format}{each_seq}->($cached_data->{$seq}, $key);
	    }
	    $output_formats{$format}{footer}->($cached_data, $key);
	}
    }
    else { # no data in cache
	print $sw->header();
	print_empty_form();
	print $sw->footer();
    }
}
else {	# entry point 3 - initial entry, empty form
    print $sw->header();
    print_empty_form();
    print $sw->footer();
}

###############################################################################
################################### END #######################################
###############################################################################


sub map_fasta {
    my $dataHash = shift;
    my $output_coderef = shift;
    my $cachekey = shift;
    my $ssaha = Ssaha_wrapper->new();

    my $condition	    = $dataHash->{condition};
    my $signature_seq	    = $dataHash->{signature_seq};
    my $restriction_seq	    = $dataHash->{restriction_seq};
    my $threshold_sig   = $dataHash->{threshold_sig};
    my $threshold_cont   = $dataHash->{threshold_cont};
    my $output_option	    = $dataHash->{output_option};
    my $contamination_seq1   = uc $dataHash->{contamination_seq1};
    my $contamination_seq2   = uc $dataHash->{contamination_seq2};
    my $contamination_bp    = $dataHash->{contamination_bp};
    my $gap_penalty	    = $dataHash->{gap_penalty};
    my $match_score	    = $dataHash->{match_score};
    my $mismatch_score	    = $dataHash->{mismatch_score};
    my $format		    = $dataHash->{format};
    my $species		    = $dataHash->{species};
    my $ssaha_score         = $dataHash->{ssaha_score};
    my $ssaha_fold          = $dataHash->{ssaha_fold};




    my ($reg, $sa);
    eval {
	$reg = 'Bio::EnsEMBL::Registry';

	$reg->load_registry_from_db (	-host => 'ensdb-archive',
					-user => 'ensro',
					-port => 5304,
				    );

	$sa = $reg->get_adaptor($species, "core", "Slice");
    };
    warn $@ if $@;

    # loop over each sequence in the FASTA, align signature, identify genomic
    # sequence and map it to the genome.

    my %map_for_seq;
  
    $map_for_seq{'__species'}	= $species;

    $map_for_seq{'__total'}	= 0;    # total sequences processed
    $map_for_seq{'__sig'}	= 0;    # seqs with signature
    $map_for_seq{'__map'}	= 0;    # seqs which map to genome
    $map_for_seq{'__contam'}	= 0;    # seqs with contamination
    $map_for_seq{'__good'}	= 0;    # 'good' seqs - sig & no contamination
    $map_for_seq{'__gene'}      = 0;    # seqs mapped on genes
    my $fastaHash = $dataHash->{fasta};
    
    foreach my $name (sort keys %{$fastaHash}){

	if (exists $map_for_seq{$name}){
	    warn "Duplicate FASTA entry for $name";
	    next;
	}
	$map_for_seq{$name}{name} = $name;  # output subs get the subhash
					    # so won't know the name key
	$map_for_seq{$name}{species} = $species;
$map_for_seq{$name}{ssaha_fold} = $ssaha_fold; 
	$map_for_seq{'__total'}++;
	my $seq_1 = $fastaHash->{$name};
	my $len = length ($seq_1);
	$map_for_seq{$name}{orig_length} = $len;
	my $seq;
	my $seq_2;

	if ($condition eq "forward"){ 
	    $seq = $seq_1; 
	}
	elsif ($condition eq "reverse"){ 
	    $seq = rev_comp ($seq_1); 
	}
	else { 
	    # if unknown, cat the forward & revcomp seqs together so both are
	    # searched
	    $seq_2 = rev_comp ($seq_1);
	    $seq = $seq_1.$seq_2; 
	}
	
	$map_for_seq{$name}{sequence} = $seq;

	my @user_seq = split //, uc($seq);

	my @sig_seq = split //, uc($signature_seq);

	my @matrix = @{build_matrix(\@user_seq, \@sig_seq)};
       
	########################################################################
	#Get alignment returns length of the local alignment, position of the
	#alignment and the percentage of matches within the alignment.
	########################################################################

	my ($alignlength,$max_x,$max_y, $endup, $enddown, $align_percent,$up_align,$down_align,$align_symbol) = get_alignment(\@user_seq, \@sig_seq, \@matrix, $match_score, $mismatch_score);

	if ($align_percent > $threshold_sig) {
	    $map_for_seq{'__sig'}++;
	    # Take first XXXbp from curated seq, search for contamination seq
	    # If present, skip further mapping on this sequence

	    $map_for_seq{$name}{sig_from} = $max_x;
	    $map_for_seq{$name}{sig_to} = $endup;
	    $map_for_seq{$name}{align}{length} = $alignlength;
	    $map_for_seq{$name}{align}{percent} = $align_percent;
	    $map_for_seq{$name}{align}{upstring} = join "",@$up_align;
	    $map_for_seq{$name}{align}{downstring} = join "",@$down_align;
	    $map_for_seq{$name}{align}{symbolstring} = join "",@$align_symbol;


	    my $ori;
	    if ($condition eq "forward"){
		$ori = 'FWD';
	    }
	    elsif ($condition eq "reverse"){
		$ori = 'REV';
	    }
	    else {
		if ($endup < $len){
		    $ori = 'FWD';
		}
		else {
		    $ori = 'REV';
		}
	    }
	    $map_for_seq{$name}{orientation} = $ori;



	########################################################################
	#find up to two contaminating sequences using LSA subrotine 
	########################################################################

	    my $pattern_seq = substr ($seq,$endup,$contamination_bp); 
	    $contamination_seq1 ||= "!!!!";
	    $contamination_seq2 ||= "!!!!";
            
                if ($contamination_seq1 ne "!!!!") {

	my @user_seq1 = split //, uc($pattern_seq);

	my @sig_seq1 = split //, uc($contamination_seq1);

	my @matrix1 = @{build_matrix(\@user_seq1, \@sig_seq1)};

	my ($alignlength1,$max_x1,$max_y1,$endup1,$enddown1,$align_percent1, $up_align1,$down_align1,$align_symbol1) 
= get_alignment(\@user_seq1, \@sig_seq1, \@matrix1, $match_score, $mismatch_score);


	    if ($align_percent1 > $threshold_cont){


		$map_for_seq{$name}{contaminated} = 1;             
   		$map_for_seq{$name}{contaminate_from} = $endup+$max_x1+1;
   		$map_for_seq{$name}{contaminate_to} = $endup+$endup1+1;                                                            $map_for_seq{'__contam'}++;
		delete $map_for_seq{$name} if $output_option eq 'good';
	        $output_coderef->($map_for_seq{$name}, $cachekey);
                next;       #jump out if contaminating seq identified
	       }
            } 


                if ($contamination_seq2 ne "!!!!" and $map_for_seq{$name}{contaminated} != 1) {

	my @user_seq2 = split //, uc($pattern_seq);

	my @sig_seq2 = split //, uc($contamination_seq2);

	my @matrix2 = @{build_matrix(\@user_seq2, \@sig_seq2)};

	my ($alignlength2,$max_x2,$max_y2,$endup2,$enddown2,$align_percent2, $up_align2,$down_align2,$align_symbol2) 
= get_alignment(\@user_seq2, \@sig_seq2, \@matrix2, $match_score, $mismatch_score);


	        if ($align_percent2 > $threshold_cont){


		$map_for_seq{$name}{contaminated} = 1;             
   		$map_for_seq{$name}{contaminate_from} = $endup+$max_x2+1;
   		$map_for_seq{$name}{contaminate_to} = $endup+$endup2+1;                                                            $map_for_seq{'__contam'}++;
		delete $map_for_seq{$name} if $output_option eq 'good';
	        $output_coderef->($map_for_seq{$name}, $cachekey);
                next;
	            }
                }


	    # Definition of 'good': has signature and no contamination
            if ($map_for_seq{$name}{contaminated} != 1 and $map_for_seq{$name}{name} ne ""){
	    $map_for_seq{$name}{good} = 1;
	    $map_for_seq{'__good'}++;
               }



	    ###################################################################
	    #Find the first cut site immediately downstream of the IR and
	    #output genomic sequence.  If the cutsite is not found, output the
	    #full sequence downstream of the IR signature sequence.
	    ###################################################################
	
	    my @s = find_cut($seq, $endup, $len, $restriction_seq);
	    my $cutlen = length($restriction_seq);
	    my $seq_to_ssaha;

	    if (scalar @s){     # we have cutsites 
		$seq_to_ssaha=substr ($seq, $endup, $s[0]-$endup+$cutlen - 1);
		$map_for_seq{$name}{cutsite}{start} = $s[0];
		$map_for_seq{$name}{cutsite}{end} = $s[0] + $cutlen - 1;
	    }
	    else {
		if ($endup < $len) {
		    $seq_to_ssaha=substr ($seq, $endup, $len-$endup);
		}
		else {
		    $seq_to_ssaha=substr ($seq, $endup, $len*2-$endup);
		}
	    }


             
	    $map_for_seq{$name}{curated_sequence} = $seq_to_ssaha;


	    my @fasta_list;
	    push @fasta_list, ">user_seq\n$seq_to_ssaha";

	    my ($hit_chr, $hit_start, $hit_end, $hit_dir, $score) = ssaha_map(\@fasta_list, $ssaha, $species, $ssaha_score, $ssaha_fold);
	    if ($hit_chr){
		$map_for_seq{$name}{hit}{chr} = $hit_chr;
		$map_for_seq{$name}{hit}{start} = $hit_start;
		$map_for_seq{$name}{hit}{end} = $hit_end;
		$map_for_seq{$name}{hit}{dir} = $hit_dir;
		$map_for_seq{$name}{hit}{score} = $score;
		unless($hit_chr eq 'ERR' or $hit_chr eq 'MULTI' or $hit_chr eq 'SCORE'){
		    $map_for_seq{'__map'}++;

		    # check for overlapping genes
		    my $five_flank = $dataHash->{five_flank};
		    my $three_flank = $dataHash->{three_flank};
		    unless ($hit_dir eq 'F'){
                ($five_flank,$three_flank) = ($three_flank,$five_flank);
		    }
		    eval {
                
                
                
                my $hit_start_gene_ov;
                my $hit_stop_gene_ov;
                
                if($hit_start < $hit_end){
                
                    #hit on F strand
                    $hit_start_gene_ov= $hit_start;
                    $hit_stop_gene_ov= $hit_end;
                
                }
                else{
                    #hit on R strand
                    $hit_start_gene_ov= $hit_end;
                    $hit_stop_gene_ov= $hit_start;
                
                
                }
                
                my $slice = $sa->fetch_by_region('chromosome',$hit_chr, $hit_start_gene_ov - $five_flank, $hit_stop_gene_ov + $three_flank);
                
                
                
                foreach my $gene (@{ $slice->get_all_Genes }) {
                    if ($map_for_seq{$name}{gene_no} != 1){
                        $map_for_seq{'__gene'}++;
                    }
                    $map_for_seq{$name}{gene_no} = 1;
                    my $genename = $gene->external_name || $gene->stable_id;
                    push @{$map_for_seq{$name}{hit}{genes}}, {'name' => $genename, 'id' => $gene->stable_id};
                }
		    };
		    warn $@ if $@;
		}
	    }
	    
	} 
	else { # no signature above threshold
	    delete $map_for_seq{$name} if $output_option eq 'good';
	}

	# Stream output if required.
	if ($output_coderef){
	    $output_coderef->($map_for_seq{$name}, $cachekey);
	}

    }# end: foreach fasta entry
 
    return \%map_for_seq;
}


#######################################
# 
# Subroutine to find out the indices of the cell 
# among matrix with maximum score
#
#######################################

sub find_index_matrix {
    my ($user_seq,$sig_seq,$matrix) = @_;
    my  $m =@$user_seq;
    my  $n =@$sig_seq;
    my $index1;
    my $index2;
    my $maxcell=$matrix->[0]->[0]; # start the first element
    for (my $i=0;$i<=$m;$i++){
	for (my $j=0;$j<=$n;$j++){
	
	    if ( $maxcell<= $matrix->[$i]->[$j]){ # if find another max value
		$maxcell =$matrix->[$i]->[$j];
		$index1 = $i;
		$index2 = $j;
	    }
	}
    }
    return ($index1,$index2);
}

###########################################
# a subroutine to get score for two nucleotides,
# here defines the score for match and mismatch
#################################33

sub get_score_for_AB {
    my ($aa1,$aa2,$match_score, $mismatch_score) =  @_;
    my $DNAscore;

    $aa1 =uc $aa1;
    $aa2 = uc $aa2;
   
    if ($aa1 eq $aa2) {
	$DNAscore=$match_score;
    }
    else {
	$DNAscore=$mismatch_score;
    }   
    return $DNAscore;
}



#################################################################
# 
# A subroutine for backtracking the matrix to find out the best local alignment
# For this purpose, first to pass into the subroutine two sequences and the matrix already build.
#
##################################################################

sub get_alignment {
    my  ($sequence1,$sequence2,$matrix, $match_score, $mismatch_score)=@_;

    ############################################################
    #
    # Make two new arrays to store the alignment for two sequence
    #
    #############################################################
      my @up;
      my @down;
      my @align_symbol;
      my $alignlength;
      my $len_sigseq=@$sequence2+0;

    #################################################################
    # 
    # Call for subroutine to find out the index of cell among marix
    # which has the maximum score. Then backtracking from this cell
    #
    ##################################################################
    my ($m,$n)=find_index_matrix($sequence1,$sequence2,$matrix);
    #################################################################
    # 
    # Mark backing tracking starting position for two strands
    #
    ##################################################################

    my $endup = $m;
    my $enddown = $n;
    my $countMatch=0;

    ############################################################
    # 
    # Backtracking the alignment from cell with maximum score. 
    # The algorithm for backtracking local alignment is exactly the same 
    # with backtracking of global alignment.
    # Stop backtracking when reached cell with value zero.
    #
    #############################################################

    while( $matrix->[$m]->[$n]!=0 ){
	$alignlength++;
	my $score = get_score_for_AB($sequence1->[$m-1],$sequence2->[$n-1], $match_score, $mismatch_score);
	if($matrix->[$m]->[$n] == $matrix->[$m-1]->[$n-1] + $score){
	   unshift @up,$sequence1->[$m-1];
	   unshift @down, $sequence2->[$n-1];
	   
	   if ($sequence1->[$m-1] eq $sequence2->[$n-1]) {
	   unshift @align_symbol, "*";
	   $countMatch++;
	   }
	   else {
	   unshift @align_symbol, " ";
	   }
	   $m=$m-1;
	   $n=$n-1;
	}
	elsif($matrix->[$m]->[$n]==$matrix->[$m-1]->[$n]-1){
	    unshift @up, $sequence1->[$m-1];
	    unshift @down, "-";
	    unshift @align_symbol, " ";
	    $m=$m-1;
	}
	else {
	    unshift @up, "-";
	    unshift @down, $sequence2->[$n-1];
	    unshift @align_symbol, " ";
	    $n=$n-1;
	}
    }

my $align_percent;

    if ($alignlength > $len_sigseq) {
    $align_percent=$countMatch/$alignlength*100;
       }

    else {
    $align_percent=$countMatch/$len_sigseq*100;
    }

    return ($alignlength, $m, $n, $endup, $enddown, $align_percent, \@up, \@down, \@align_symbol);
}

#######################################
#
# Find maximum in a list 
#
#######################################
sub find_max {
  my @data = @_;
  my $max = $data[0];
  foreach my $x (@data) {
    if ($x > $max) {
      $max = $x;
    }
  }
  return $max;
}


#####################################################
#
#Subroutine read fasta files
#
#####################################
sub read_fasta {
    my ($upload_seq, $paste_seq) = @_;
    my (%name2seq, $name, $seq);

    if ($upload_seq){
	while (<$upload_seq>){
	    $paste_seq .= $_;
	}
    }

    # guess line ending and go with that:
    $paste_seq =~ /(\r\n|\r|\n)/;
    my $line_end = $1;

    my @lines = split $line_end, $paste_seq;
    
    my $seq_count = 0;
    foreach my $line(@lines){
	if ($line =~ /^>/) {
	    $seq_count++;
	    last if $seq_count > $SEQ_LIMIT;
	   $line =~ s/>//;
	    if (defined $name) {
		$name2seq{fasta}{$name} = $seq;
	    }
	    $name = $line;
	    if ($name =~ /^__/){
		# keys starting __ are reserved, so prepend a '!'
		$name = "!$name";
	    }
	    $seq = "";
	}
	else{
	    $line =~ s/[^ATGCNatgcn]//g;  # clean up sequence
	    $seq .= uc($line);	# make sequence uppercase
	}
    }
    $name2seq{fasta}{$name} = uc($seq);
    
    return \%name2seq;
}

####################################################
#
# Subroutine: do reverse complementary sequence
#
#####################################################
sub rev_comp {
    my ($seq) = @_;
    my $revcomp=reverse($seq);
    $revcomp=~ tr/atgcATGC/tacgTACG/;
    return $revcomp;
}

####################################################
#
# Subroutine: find restriction cutting site, return an array with position
# information.
#
####################################################

sub find_cut{
    my ($seq, $endup, $len, $restriction_seq) = @_;
    my @site;
    my $s;
    if ($endup < $len){              ###Only search cutting site in the sequence which is downstream the tr
        $s = substr ($seq, $endup, $len-$endup);
    }elsif ($endup > $len){
        my $f_len = $len*2;
        $s = substr ($seq, $endup,$f_len-$endup); 
    }
    while ($s =~ /$restriction_seq/g){          ###Set pattern match with restriction enzyme cutting site using regular expression
        push @site, pos($s)-length($restriction_seq)+1+$endup;
    }
    return @site;
}



sub print_empty_form {
print qq(
<h1 style="font: italic normal 900 3em arial black">iMapper</h1>
<h2 style="font: italic normal 900 1.5em arial">Insertional Mutagenesis Mapper</h2>
<br />
<fieldset>
<legend>About</legend>

<p>iMapper is a sequence analysis tool designed specifically for large-scale
analysis of insertional mutagenesis tag sequences against vertebrate and
invertebrate genomes. It trims real genomic segments from linker-based PCR
sequence input and automatically maps insertion sites onto an assembled genome.
<a href="/resources/software/imapper/help.html">Learn more</a></p> 

<p>Submissions are limited to $SEQ_LIMIT sequences.</p>

<p><a href="/resources/software/imapper/data/pb_test.fa">Download</a> a file of 96 PiggyBac traces for test use.</p>

<p>Comments & Questions: <a href="mailto:imapper\@sanger.ac.uk">iMapper\@sanger.ac.uk</a></p>

<p>
<b>Citing iMapper</b><br />
If you use iMapper to analyze your insertional datasets and find it useful, we would be appreciate if you could consider citing the reference that describes this work:<br /><br />
Kong, J., Zhu, F., Stalker, J. and Adams, D.J. (2008) <a href="http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2639305/">iMapper: a web application for the automated analysis and mapping of insertional mutagenesis sequence data against Ensembl genomes.</a><br />
<i>Bioinformatics (Oxford, England)</i>, <b>24</b>, 2923-2925.
</p>

</fieldset>

<br /><br />

<fieldset>
<legend>Submit your sequence to iMapper</legend>
<br />
<form action="imapper.cgi" method="post" enctype="multipart/form-data">
<table align="center" class="imapper">

<tr>
    <th>Species:</th>
        <td><select name="species">
);
foreach my $spp(sort keys %allowed_spp){
    my $selected = $spp eq $default_spp ? "selected" : "";
    print qq(<option $selected value="$spp">$spp</option>);
}
print qq(
    </select></td>
</tr>

<tr>
    <th>Sequencing from:</th>
        <td><select name="seq_format">
    <option value="unknown">Unknown </option>
    <option value="forward">Forward (tag seq end)</option>
    <option value="reverse">Reverse (linker seq end)</option>
    </select></td>
</tr>

<tr>
    <th>Output option:</th>
        <td><input type="radio" value="all" checked name="output1">Output all sequences<br />
<input type="radio" name="output1" value="good">Output good sequences (containing tag sequence)</td>
</tr>

<tr>
    <th>Mutagen tag sequence:</th>
        <td><input name="signature_seq" id="signature_seq" size=40 value="tatctttctagggttaa">  <br /> <br />
        Or choose from preset: <select         name="signatures" onchange="setSigtext(this)">
	<option value="TATCTTTCTAGGGTTAA">PiggyBac 5'&3'TR</option>
	<option value="TTCCGACTTCAACTGTA">Sleeping Beauty 5'&3'TR</option>
	<option value="ACAGGTGGGGTCTTTCA">Retrovirus U3LTR</option>
	</select></td>
</tr>

<tr>
    <th>Restriction site:</th>
        <td><input name="restriction_seq" size=15 value="gatc"></td>
</tr>

<tr>
    <th>Advanced options:</th>
        <td><a id="advtoggleanchor" href="javascript:toggleAdv();void(0);">Show advanced options</a><br\>
    <div id="advopt" style="display:none">
	<table style="border-top:#83a4c3 dotted 1px;border-left:#83a4c3 dotted 1px;border-right:#83a4c3 dotted 1px;">

	    <th>Alignment criteria:</th>
	    	    <td> Alignment threshold: &gt;<input name="threshold_sig" size=3 value="80">  % <br />
            Gap penalty: <input name="gap_penalty" size=4 value="-1"><br />
	    Match score: <input name="match_score" size=4 value="1"><br />
	    Mismatch score: <input name="mismatch_score" size=4 value="-1">
	    </td>
	</tr>
	<tr>
	    <th>Contaminating<br />sequences:</th>
	    	    <td> seq1: <input name="contamination_seq1" size=25><br /><br />

	    	    seq2: <input name="contamination_seq2" size=25><br /><br />
                    Alignment percentage: &gt;<input name="threshold_cont" size=3 value="90"> % <br /><br />
	    Search in first: <input name="contamination_bp" size=4 value="100"> residues
	    </td>
	</tr>
	<tr>
	    <th>Overlapping genes:</th>
	    	    <td>A gene overlaps a hit if it falls<br \>
			within the following window:
		    <br />5' flanking bp: <input name="five_flank" size=6 value=10000><br /><br />
	    	    3' flanking bp: <input name="three_flank" size=6 value=10000><br /><br />
	    </td>
	</tr>
	<tr>
	    <th>SSAHA mapping parameters:</th>
	    	    <td>A unique mapping will be determined if:
		    <br />SSAHA mapping score > <input name="ssaha_score" size=6 value="35"><br />
                    <br />and<br />
	    	    The score for best hit > <input name="ssaha_fold" size=3 value="1.5"> fold of other scores<br /><br />
	    </td>
	</tr>
	</table>
    </div></td>
</tr>
<tr>
    <th>Please upload your sequence file:</th>
        <td><input type="file" name="upload_seq" size="40"></td>
    </td>
</tr>

<tr>
    <th>Or<br />
    paste here in FASTA format:</th>
        <td><textarea name="paste_seq" rows=12 cols=50 wrap="virtual" style="float:left"></textarea></td>
    </td>
</tr>

<tr>
    <th>&nbsp;</th>
        <td><input type="submit"> <input type="reset"></td>
    </td>
</tr>
</table>
</fieldset>
</FORM>
);


}

sub wrap_sequence{
    
    #split up into blocks
    my ($sequence)= @_;
    my $num_bases_per_line= 80;
    my $sequence_chunky= "";
  
    
    while($sequence =~/(.{$num_bases_per_line})/gc){	# keep match pos
        $sequence_chunky .= "$1\n";
    }
    
    if($sequence =~/\G(.+)/){	# match from last match pos
        $sequence_chunky .= $1;
    }
    
    return $sequence_chunky;

}


sub markup_wrap_sequence {
    
    #split up into blocks
    my $seqdata= shift;
    my $num_bases_per_line= 60;
    my $seq_chunky= "";
    my $sequence = $seqdata->{sequence};
    

    my $sig_from    = $seqdata->{sig_from};
    my $sig_to	    = $seqdata->{sig_to};
    my $cut_from    = $seqdata->{cutsite}{start};	
    my $cut_to	    = $seqdata->{cutsite}{end};
    my $len	    = $seqdata->{orig_length};
    my $cutlen	    = $cut_to - $cut_from + 1;
    my $contaminate_from = $seqdata->{contaminate_from};
    my $contaminate_to = $seqdata->{contaminate_to};
    my $contaminate_len = $contaminate_to - $contaminate_from;

    if ($cut_from eq "") {  
	if ($sig_to < $len) {
	    $cut_from = $len;
	}
	else {
	    $cut_from = 2*$len;  
	}
    }

    
if ($sig_to > $len) {
	$sequence = substr($sequence,length($sequence)/2,length($sequence)/2);
        $sig_from = $sig_from-$len;
        $sig_to = $sig_to-$len;
        $cut_from = $cut_from-$len;
        $cut_to = $cut_to-$len;

        if ($contaminate_from){
        $contaminate_from  = $contaminate_from-$len;
        $contaminate_to  = $contaminate_to-$len;
           }
	}


    while($sequence =~/(.{$num_bases_per_line})/gc){	# keep match pos
	my $subseq = $1;
	$subseq =~ s/N$/1/;
	$subseq =~ s/A$/2/;
	$subseq =~ s/T$/3/;
	$subseq =~ s/G$/4/;
	$subseq =~ s/C$/5/;
        $seq_chunky .= $subseq;
    }
    
    if($sequence =~/\G(.+)/){	# match from last match pos
        $seq_chunky .= "$1<br />";
    }

    # print sequence colour-coded in blocks:
    # seq start to seq_from (max_x) : upstream stuff
    # seq_from to seq_to    : sig sequence
    # seq_to to cutsite_start : genomic sequence
    # cutsite_start to cutsite_end: cut site
    # cutsite_end to end of seq : downstream stuff


my $marked_up_seq;


    $marked_up_seq =  substr ($seq_chunky,0,$sig_from);

    $marked_up_seq .= "<span style=\"background-color: $SIG_COLOUR\">". substr ($seq_chunky,$sig_from,$sig_to-$sig_from)."</span>";

    if ($contaminate_from eq "") {     #output color-annotated sequences 

    $marked_up_seq .=  "<span style=\"background-color: $SEQ_COLOUR\">". substr ($seq_chunky,$sig_to,$cut_from-$sig_to-1). "</span>";

    if ($cut_from != $len){
	$marked_up_seq .= "<span style=\"background-color: $CUT_COLOUR\">". substr ($seq_chunky,$cut_from-1,$cutlen). "</span>";
	$marked_up_seq .= substr ($seq_chunky,$cut_from-1+$cutlen,$len-$cut_from-$cutlen+1). "</span>"; 	    
    }
    }

    else{        # output color-annotated contamination sequences
        $marked_up_seq .=  substr ($seq_chunky,$sig_to,$contaminate_from-$sig_to-1);

	$marked_up_seq .= "<span style=\"background-color: $CONT_COLOUR\">". substr ($seq_chunky,$contaminate_from-1,$contaminate_len). "</span>";

	$marked_up_seq .= substr ($seq_chunky,$contaminate_to-1,$len-$contaminate_to+1);
    }   

    $marked_up_seq =~ s/1/N\n/g;
    $marked_up_seq =~ s/2/A\n/g;
    $marked_up_seq =~ s/3/T\n/g;
    $marked_up_seq =~ s/4/G\n/g;
    $marked_up_seq =~ s/5/C\n/g;
    return $marked_up_seq;

}


sub build_matrix {
    my ($user_seq, $sig_seq) = @_;
    my @matrix;

    my $m = scalar(@$user_seq);
    my $n = scalar(@$sig_seq);

    my $gapPenalty = -1;
    my $match = 1;
    my $mismatch = -1;

    # Initialise the scoring matrix with 0

    for (my $i=0;$i<=$m;$i++){
      $matrix[$i][0] = 0;
    }
    for (my $j=0;$j<=$n;$j++){
      $matrix[0][$j] = 0;
    }


    #######################################
    # For all cells in the  matrix
    #######################################

    for (my $i=1;$i<=$m;$i++){
      for (my $j=1;$j<=$n;$j++){

    ##################################
    # Assigning MATRIX(i,j)
    ##################################

      my $score = get_score_for_AB($user_seq->[$i-1],$sig_seq->[$j-1], $match, $mismatch);

	$matrix[$i][$j] = find_max(0,
                                   $matrix[$i-1][$j] + $gapPenalty,
				   $matrix[$i][$j-1]+ $gapPenalty,
				   $matrix[$i-1][$j-1]+ $score);
	                             
      }
    }

    return \@matrix;
}


sub ssaha_map {
    my ($fasta_list, $ssaha, $species, $ssaha_score, $ssaha_fold) = @_;

    my ($chr, $start, $end, $dir, $score);
    my %ssaha_results= %{$ssaha->do_ssaha($species, $fasta_list)};

    if(%ssaha_results){
	
	if(! defined($ssaha_results{'user_seq'})){
	    #no hits
	}
	else {
	    # hits
	    my @results_list= @{$ssaha_results{'user_seq'}};
	    
	    # best hit must be score > $ssaha_score and score >= $ssaha_fold* next hit score
	    if ($results_list[0]->{Score} >= $ssaha_score){
		my $best_hit;
		if(defined $results_list[1]){
		    if ($results_list[0]->{Score} >= $ssaha_fold*$results_list[1]->{Score}){
		    $best_hit=$results_list[0];
		    }
		    else {
			$chr = "MULTI";
		    }
		}
		else {
		    # unique and above threshold
		    $best_hit=$results_list[0];
		}
	    
		$chr    = $best_hit->{'S_Name'} if $chr eq "";
		$start  = $best_hit->{'S_Start'};
		$end    = $best_hit->{'S_End'};
		$dir    = $best_hit->{'Direction'};
		$score  = $best_hit->{'Score'};
	    }
	    else {
		$chr = "SCORE";
	    }
	}
	
    }
    else {
	# problem with ssaha server e.g. server down, network problem
	# wrong port number, etc.
	$chr = "ERR";
    }

    return ($chr, $start, $end, $dir, $score);
}

sub output_tabular_header {
    print $sw->header();
    print qq(<h1 style="font: italic normal 900 3em arial black">iMapper</h1><br />);
}

sub output_tabular {
    my ($seqdata, $cachekey) = @_;

#my $ssaha_fold = $seqdata->{ssaha_fold};
#print "$ssaha_fold,@@@@@@@@@@@@@";
#print $seqdata->{ssaha_fold};


    unless (defined $seqdata->{name}){
    return;
    }

    print qq(<table class="imap_output">
	<col id="key"/>
	<col id="value"/>
    );

#    print "<tr>";
#    print "<p><input type='checkbox' name='$seqdata->{name}' value='NO' checked></p>";
#    print "</tr>";

    print "<tr>\n";
    print "<th>Name</th>\n";
    print "<td>".$seqdata->{name}."</td>\n";
    print "</tr>\n";

    print "<tr>\n";
    print "<th>Tag Sequence</th>\n";
    unless (defined $seqdata->{sig_from}){
	print "<td>No tag sequence found</td>\n";
	print "</tr>\n";
	print "</table>\n";
	return;
    }
    print "<td>From ",$seqdata->{sig_from} + 1," to ",$seqdata->{sig_to};
    my $ori = $seqdata->{orientation};
    print " on the $ori strand </td>\n";
    print "</tr>\n";

    print "<tr>\n";
    print "<th>Alignment</th>\n";
    print qq(<td><pre style="margin-top: 0px; margin-bottom: 0px;">);
    print $seqdata->{align}{upstring}."\n";
    print $seqdata->{align}{downstring}."\n";
    print $seqdata->{align}{symbolstring}."</pre></td>\n";
    print "</tr>\n";

   if (	$seqdata->{contaminated} == 1) {  # output annotated contamination sequence
    print "<tr>\n";
    print "<th>Contaminating sequence</th>\n";
    print "<td> This sequence contains contaminating sequence</td>\n";
    print "<tr>\n";
    print "<th>Annotated sequence output</th>\n";
    print qq(<td><pre style="margin-top: 0px; margin-bottom: 1em;">).markup_wrap_sequence($seqdata)."</pre>";
    print "<span style=\"background-color: $CONT_COLOUR\">NNN</span>"," contaminating sequence";  
    print "</td>\n</tr>\n";
    return;
   }


    print "<tr>\n";
    print "<th>Restriction site</th>\n";
    if ($seqdata->{cutsite}){
	print "<td>At position ".$seqdata->{cutsite}{start}."</td>\n";
    }
    else {
	print "<td>No restriction site found</td>\n";
    }
    print "</tr>\n";
    print "<tr>\n";

    print "<th>Map to genome</th>\n";
    my $has_genomic_hit = 0;
    if ($seqdata->{hit}{chr}){
	if ($seqdata->{hit}{chr} eq 'ERR'){
	    print "<td>SSAHA mapping error</td>\n";
	}
	elsif ($seqdata->{hit}{chr} eq 'MULTI'){
	    print "<td>Sequence maps to multiple locations</td>\n";
	}
	elsif ($seqdata->{hit}{chr} eq 'SCORE'){
	    print "<td>Sequence mapping score is below threshold</td>\n";
	}
	else {
	    $has_genomic_hit = 1;
	    my $chr=$seqdata->{hit}{chr};
	    my $start = $seqdata->{hit}{start};
	    my $end = $seqdata->{hit}{end};
	    my $dir = $seqdata->{hit}{dir};
	    $dir = $dir eq 'F' ? 'fwd' : 'rev';
	    print "<td>Chr: $chr Start: $start End: $end Dir: $dir";

	    # these will give exact location
	    # my $centre = ($start + $end)/2;
	    # my $hit_len = $end - $start + 1;

	    # these will give 100bp around the insertion point
	    my $centre = $dir eq 'fwd' ? $start : $end;
	    my $hit_len = 100;

	    my $species = $seqdata->{species};

        print qq(&nbsp;&nbsp;&nbsp;[ <a href="http://www.ensembl.org/$species/Location/View?r=${chr}:${centre}-${centre}" target="_blank">Ensembl Location View</a> ]);
        
	    print qq(</td>\n);
	}
    }
    else {
	    print "<td>No genome match found</td>\n";
    }
    print "</tr>\n";

    print "<tr>\n";
    print "<th>Overlapping genes</th>\n";
    print qq(<td>);
    if ($seqdata->{hit}{genes}){
	print qq(The mapped location overlaps: );
	my @urls;
	foreach my $gene(@{$seqdata->{hit}{genes}}){
	    my $name	= $gene->{'name'};
	    my $id	= $gene->{'id'};
	    my $species = $seqdata->{species};
	    
        push @urls, qq(<a href="http://www.ensembl.org/$species/Gene/Summary?g=${id}" target="_blank">$name</a>);
        
	}
	print join ", ", @urls;
    }
    else {
	if ($has_genomic_hit){
	    print qq(The mapped location does not overlap any genes in Ensembl.);
	}
	else {
	    print qq(No mapped location to check for genes);
	}

    }
    print "</td>\n</tr>\n";
    print "<tr>\n";
    print "<th>Annotated sequence output</th>\n";
    print qq(<td><pre style="margin-top: 0px; margin-bottom: 1em;">).markup_wrap_sequence($seqdata)."</pre>";
    print "<span style=\"background-color: $SIG_COLOUR\">NNN</span>"," tag sequence ";
   print "<span style=\"background-color: $SEQ_COLOUR\">NNN</span>"," genomic sequence ";
   print "<span style=\"background-color: $CUT_COLOUR\">NNN</span>"," restriction site";

    
    print "</td>\n</tr>\n";
    print "</table>";
}

sub output_tabular_footer {
    my ($seqdata, $cachekey) = @_;
    my $total_seqs  = $seqdata->{'__total'};
    my $sig_seqs    = $seqdata->{'__sig'};
    my $mapped_seqs = $seqdata->{'__map'};
    my $good_seqs   = $seqdata->{'__good'};
    my $contam_seqs = $seqdata->{'__contam'};
    my $species	    = $seqdata->{'__species'};
    my $gene_no      = $seqdata->{'__gene'};
    print qq(<br /><fieldset>
    <legend>Summary</legend>
    <table align="left" class="imap_summary">
	<tr>
	    <th>Total sequences submitted:</th>
	    <td>$total_seqs</td>
	</tr>
	<tr>
	    <th>Sequences containing tag sequence:</th>
	    <td>$sig_seqs</td>
	</tr>
	<tr>
	    <th>Sequences containing contaminant sequence:</th>
	    <td>$contam_seqs</td>
	</tr>
	<tr>
	    <th>Sequences mapped to reference genome:</th>
	    <td>$mapped_seqs</td>
	</tr>
	<tr>
	    <th>Sequences mapped overlapping genes:</th>
	    <td>$gene_no</td>
	</tr>
    </table>
    </fieldset>
    );
    print qq(<br /><br />\n);
    # Only link to other formats if there are appropriate data
    if ($sig_seqs || $mapped_seqs){
	print qq(
	<fieldset>
	<legend>Data formats</legend>
	<br />);
    if($sig_seqs){
	print qq(<a href="imapper.cgi?key=$cachekey&format=fasta">Generate FASTA file (FASTA file of clipped sequences for sequence reads containing tag and no contamination) </a><br /><br />);
    }
    if ($mapped_seqs){
	print qq(<a href="imapper.cgi?key=$cachekey&format=gff">Generate GFF file (GFF file of features for all the mapped sequences)</a><br /><br />
	<a href="http://www.ensembl.org/$species/karyoview?url_file_1=http://www.sanger.ac.uk/cgi-bin/teams/team113/imapper.cgi\%3Fkey=$cachekey\%26format=loc">Generate chromosome side view graph (Chromosome side view graph of all the insertion sites via Ensembl KaryoView)</a><br />);
    }
    print qq(</fieldset><br /><br />);
    }
    print $sw->footer();
}

sub output_fasta_header {
    print $cgi->header( -type => 'text/plain');
}

sub output_fasta {
    my ($seqdata, $cachekey) = @_;
    if ($seqdata->{curated_sequence}){
	print ">".$seqdata->{name}."\n";
	print wrap_sequence($seqdata->{curated_sequence})."\n";
    }
}

sub output_fasta_footer {
    my ($seqdata, $cachekey) = @_;
    # no footer
}

sub output_gff_header {
    print $cgi->header( -type => 'text/plain');
}

sub output_gff {
    my ($seqdata, $cachekey) = @_;
    my $dataHash = shift;

    # <seq> <source> <feature> <start> <end> <score> <strand> <frame> [attribs]
    return unless $seqdata->{hit}{chr};
    return if (   $seqdata->{hit}{chr} eq 'ERR' 
		 || $seqdata->{hit}{chr} eq 'MULTI' 
		 || $seqdata->{hit}{chr} eq 'SCORE'
		 );

    my $chr	= $seqdata->{hit}{chr};
    my $start	= $seqdata->{hit}{start};
    my $end	= $seqdata->{hit}{end};
    my $dir	= $seqdata->{hit}{dir};
    $dir = $dir eq 'F' ? '+' : '-';
    my $score	= $seqdata->{hit}{score};
    my $name	= $seqdata->{name};
    my $ori = $seqdata->{orientation}." strand";
    my $hit_gene;
    my $gene_name;

    if ($seqdata->{hit}{genes}){
	$hit_gene = "Hit on gene";
	foreach my $gene(@{$seqdata->{hit}{genes}}){
	    $gene_name	.= $gene->{'id'}."\t".$gene->{'name'}."\t";
	    }
        } 
    else
        {
         $hit_gene = "No hit gene";
         $gene_name = ".";
        }


    print join "\t", ($chr, 'ssaha', $ori, $start, $end, $score, $dir, '.', $name, $hit_gene, $gene_name);
    print "\n";
}

sub output_gff_footer {
    my ($seqdata, $cachekey) = @_;
    # no footer
}


sub output_loc_header {
    print $cgi->header( -type => 'text/plain');
}

# loc is simple chr, start, end, name output format for karyoview upload
sub output_loc {
    my ($seqdata, $cachekey) = @_;
    return unless $seqdata->{hit}{chr};
    return if (   $seqdata->{hit}{chr} eq 'ERR' 
		 || $seqdata->{hit}{chr} eq 'MULTI' 
		 || $seqdata->{hit}{chr} eq 'SCORE' 
		 );

    my $chr	= $seqdata->{hit}{chr};
    my $start	= $seqdata->{hit}{start};
    my $end	= $seqdata->{hit}{end};
    my $name	= $seqdata->{name};

    print join "\t", ($chr, $start, $end, $name);
    print "\n";
}

sub output_loc_footer {
    my ($seqdata, $cachekey) = @_;
    # no footer
}
