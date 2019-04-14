




                    SSEalign Version 1.1 - Usage Guide
                 ===========================================

1. Prerequisites for SSEalgin

1.1 Linux platform

1.2 Perl package

1.3 SSpro8 software

1.4 Dynamic programming algorithm package


2. Detail of the required tools
  
2.1 Linux platform

SSEalign package should be compatible with any Linux operating system. If a problem occurs during the installation or if SSEalign is not runningproperly, please read the information reported below about the possible problems.
  
2.2 Perl package

The scripts in SSEalign are written in perl. These scripts assume perl is
installed in the default location for linux operating systems: /usr/bin/perl
If perl is not installed in this folder of your system, please create a
symlink to your perl installation in /usr/bin.  
  
2.3 SSpro8 software

The SSpro8 software is a tool for prediction of protein secondary structure and can be download in http://scratch.proteomics.ics.uci.edu/. Other tools for protein secondary structure, such JPRED4 and PsiPred is also acceptable, but their results will not be not as good as SSpro8.

2.4 Dynamic programming algorithm package

The package of Dynamic programming algorithm package is recommended to download in https://github.com/noporpoise/seq-align. Please make sure to configure this package and add the path into the .bashrc before using SSEalign. The EMBOSS tool is also acceptable for dynamic programming algorithm but its file size is much large. 

2.5 other suggestions

Note also that the 'tmp' folders in the package must not be removed or renamed as they are used by the various tools to store intermediate files during a run. Temporary files are removed after each run and several instances of SSEalign can run simultaneously.

3. Testing SSEalign usage

To test the usage of SSEalign, please change directory to SSEalign storage area and run SSEalign on the provided example datasets. if it runs with no errors, means that the SSEalign have been installed successfully.


    perl SSEalign_two_groups.pl example_SSE_group1.fasta example_SSE_group2.fasta example_output.txt 

    perl SSEalign_provided_list.pl example_total_SSE.fasta example_provided_list.txt example_output_table.txt example_output_detail.txt


4. Citations

If you applied the SSEalign in your studies ,please cite our paper:

Zhiyuan Yang, Stephen Kwok-Wing Tsui. Functional Annotation of Proteins Encoded by the Minimal Bacterial Genome Based on Secondary Structure Element Alignment [J]. Journal of Proteome Research, 2018, 17(7):2511-2520.


  
