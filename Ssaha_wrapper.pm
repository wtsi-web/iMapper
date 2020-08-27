package Ssaha_wrapper;
# wrapper module around the ssaha client script ssaha2Client.pl that is
# provied with ssaha2server
#
# this is used for client/server based ssaha
#
# place any port/host config info and ssaha parameters into the constructor
#
# Author:        sr7
# Created:       2008-02

#     Copyright (C) 2008  Stephen Rice

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




use IPC::Open2;

use strict;
use warnings;
no warnings 'uninitialized';

$ENV{PATH}= '/usr/local/bin';   # solely to stop taint from barfing

sub new {
    my ($class)= @_;
    
    my $self= {};
    
    #Put host and ports for the ssaha servers here
    
    $self->{'Homo_sapiens'}	= {'port' => 40011, 'host' => 'ssaha03'};
    $self->{'Mus_musculus'}	= {'port' => 40014, 'host' => 'ssaha02'};
    $self->{'Rattus_norvegicus'} = {'port' => 40016, 'host' => 'ssaha15'};
    $self->{'Drosophila_melanogaster'} = {'port' => 40008, 'host' => 'ssaha08'};
    $self->{'Danio_rerio'} = {'port' => 40007, 'host' => 'ssaha02'};
    $self->{'Saccharomyces_cerevisiae'} = {'port' => 40017, 'host' => 'ssaha02'};

    
    #Put commandline params for ssaha2 here:
    
    #Usage: ./ssaha2Client.pl -seeds <5> -score <30> -depth <50> -centre <0> -species <0> -trace <0> -identity <50.0> -port <60000> -server <localhost> -align <0> -copyright <0> -output <ssaha2>\n\n");
  
  $self->{'params'}={
  
    '-seeds' => 2,
    '-score' => 20,
    '-depth'=> 5,
    '-identity' => 30.0,
    '-align' => 0,
    '-output' =>'ssaha2',
  
  };
    
    bless($self, $class);
    
    return $self;


}

sub do_ssaha{
    
    #input the species name to search against
    # and a array-ref of sequences in fasta format
    
    #on success returns a hash of results
    
    #fasta_header => [
    #           {score => $score, Q_name => $q_name, ...}
    #           {}
    #           ...
    #       ]
            
    # on failure returns undef
    
    
    my ($self, $species, $aref)= @_;
    
    my $host= $self->{$species}->{'host'};
    my $port= $self->{$species}->{'port'};
    
    if(!defined($host)){
        print STDERR "No ssaha host defined for this species\n";
        return undef;
    
    }
    
    if(!defined($port)){
        print STDERR "No ssaha port defined for this species\n";
        return undef;
        
    
    }
    
    # convert the array of fasta files into a hash
    my $fasta_string= "";
    
    foreach my $fasta(@{$aref}){
        
        
        $fasta_string .= "$fasta\n";
    
    }
    
   
    # process each fasta file in the list
    
    # construct the command
    my $comm= "/WWW/SANGER_docs/bin-offline/SSAHA/ssaha2Client.pl -port $port -server $host ";
    
    foreach my $option(keys %{$self->{params}} ){
        $comm .= " " . $option . " " . $self->{params}->{$option};
    
    }


    #file handles
    
    my $w_fh;
    my $r_fh;
    
    my %ssaha_results;
    
    eval{
    
        open2($r_fh, $w_fh, $comm);
        print $w_fh "$fasta_string";
        close($w_fh);
        
        # get the ssaha output and process
        
        while(<$r_fh>){
        
            my $line= $_;
            chomp $line;
            
            my $q_name;
            
            
            if($line =~ /^Matches For Query .+: (\S+)/){
                $q_name= $1;
                $ssaha_results{$q_name}= undef;
            }
            
            if($line =~/^\d/){
                # data line starts with a digit (score)
                # -tags does not work with client/server version of ssaha
                
                
                my @fields= split(/\s+/, $line);
                
                my $q_name= $fields[1];
                
                my $href= {'Score' => $fields[0], 'Q_Name' => $fields[1], 'S_Name' => $fields[2], 'Q_Start' => $fields[3], 'Q_End' => $fields[4], 'S_Start' => $fields[5], 'S_End' => $fields[6], 'Direction' => $fields[7], 'Num_bases' => $fields[8], 'Identity' => $fields[9]};
                
                push @{$ssaha_results{$q_name}}, $href;
            
            }
        
        }
        
        close($r_fh);
    
    };
    
    if($@){
        #error occurred
	warn "SSAHA: $@";
        return undef;
    }
    else{
        return \%ssaha_results;
    
    }
    
}

1;
