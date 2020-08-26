# iMapper
iMapper (Insertional Mutagenesis Mapper) was a web-based tool for the efficient analysis of insertion site sequence reads against vertebrate and invertebrate Ensembl genomes.
Until recently, this tool was available on the Sanger Institute's web site. The code is provided here for anyone that wishes to reconstruct it.

The code is comprised of a Perl CGI script: imapper.cgi and a Perl module: Ssaha_wrapper.pm.

immaper.cgi contains the core code and Ssaha_wrapper.pm contains code that interfaces with ssaha_server.

