




                    SSEalign Version 1.1 - Usage Guide
                 ===========================================

  1. Prerequisites for SSEalgin

1.1 Linux platform

1.2 Perl package

1.3 SSpro8 software

1.4 EMBOSS toolkits


  2. Detail of the required tools
  
2.1 Linux platform

SSEalign package should be compatible with any Linux operating system. If a problem occurs during the installation or if SSEalign is not runningproperly, please read the information reported below about the possible problems.
  
2.2 Perl package

The scripts in SSEalign are written in perl. These scripts assume perl is
installed in the default location for linux operating systems: /usr/bin/perl
If perl is not installed in this folder of your system, please create a
symlink to your perl installation in /usr/bin.  
  
2.3 SSpro8 software

The SSpro8 software is a tool for prediction of protein secondary structure and can be download in http://scratch.proteomics.ics.uci.edu/.

2.4 EMBOSS toolkits

Our alignment tool EMBOSS-strecher is included in EMBOSS toolkits and this tool can be download in http://emboss.sourceforge.net/download/. 

2.5 other suggestions

Note also that the 'tmp' folders in the package must not be removed or renamed as they are used by the various tools to store intermediate files during a run. Temporary files are removed after each run and several instances of SSEalign can run simultaneously.


  3. Testing SSEalign usage

To test the usage of SSEalign, please change directory to SSEalign storage area and run SSEalign on the provided example datasets. if it runs with no errors, means that the SSEalign have been installed successfully.


    perl SSEalign_two_groups.pl example_SSE_group1.fasta example_SSE_group1.fasta example_output.txt 59.18 2 

    perl SSEalign_provided_list.pl example_total_SSE.fasta example_provided_list.txt example_output_table.txt example_output_detail.txt



  4. Citation

If you have applied the SSEalign in your studies ,please cite our paper:
Yang, Zhiyuan and Tsui, Stephen. "SSEalign: accurate function prediction of bacterial unannotated protein, based on effective training dataset" Bioinformatics (under revision). 


  
