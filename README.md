# iMapper
iMapper (Insertional Mutagenesis Mapper) was a web-based tool for the efficient analysis of insertion site sequence reads against vertebrate and invertebrate Ensembl genomes.
Until recently, this tool was available on the Sanger Institute's web site. The code is provided here for anyone that wishes to reconstruct it.

The code is comprised of a Perl CGI script: imapper.cgi and a Perl module: Ssaha_wrapper.pm.

immaper.cgi contains the core code and Ssaha_wrapper.pm contains code that interfaces with ssaha_server.

If you use iMapper to analyze your insertional datasets and find it useful, we would be appreciate if you could consider citing the reference that describes this work:
Kong, J., Zhu, F., Stalker, J. and Adams, D.J. (2008) iMapper: a web application for the automated analysis and mapping of insertional mutagenesis sequence data against Ensembl genomes. Bioinformatics (Oxford, England), 24, 2923-2925.

