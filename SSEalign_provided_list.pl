#!/usr/bin/perl
# Description: This script is for Secondary Structure Element Alignment for provided list of protein-protein pair
# usage example: SSEalign_provided_list.pl <total_SSE.fasta> <provided_list.txt> <output_table.txt> <output_align.txt>
# Author: YANG Zhiyuan ( yangzhiyuan@link.cuhk.edu.hk )
# Institute: The Chinese University of Hong Kong
# Date: 2017-03-23

use strict;
use warnings;

# check the input parameters
 if ( $#ARGV<3 )  {   die "\n Please input 4 parameters as following:\nperl SSEalign_provided_list.pl <total_SSE.fasta> <provided_list.txt> <output_table.txt> <output_align.txt> \n\n Other suggestions: \n1.All SSE sequences must be included in the file <total_SSE.fasta> in FASTA format; \n2.The names of the sequences in the file <The total_SSE.fasta> must be only composed of English letters,numbers and underlines; \n3.The file <provided_list.txt> must be TAB separated, if multiple colunms are provided, the first 2 columns are used. \n4.The file <output_table.txt> is the result in table format, while the file <output_align.txt> contain alignment details. \n\n";    }
 
 # check the files
 unless (-e $ARGV[0] ) { die "The file $ARGV[0] is not existent! \n";  }
 unless (-e $ARGV[1] ) { die "The file $ARGV[1] is not existent! \n";  } 
 unless (-d 'tmp' ) {  system ('mkdir tmp');   }
 
 my $QUERY=$ARGV[0];
 my $LIST=$ARGV[1];
 my $OUTTAB=$ARGV[2];
 my $OUTDEL=$ARGV[3];
 
 my @querys;
 my @subjects;
 my %seq;

open (OUTTAB,">$OUTTAB") or die "cannot open $OUTTAB \n";
open (OUTDEL,">$OUTDEL") or die "cannot open $OUTDEL \n";

print OUTTAB "Query\tSubject\tGap number\tGap open number\tPercent\tWiden%\n";

{
local $/="\n>";
open (QUE,"<$QUERY"); 
 while (my $content=<QUE>)  {
  $content=~s/>//g;
  my $loc=index($content,"\n",0);
  my $firstline=substr($content,0,$loc);
  my @parts=split/\s+/,$firstline;
  
  my $sequence=substr($content,$loc+1);
  $sequence=~s/\n//g;
  $sequence=uc($sequence);
  $seq{$parts[0]}=$content;
} 
 close QUE;
}

my $suc_count=0;
my $pro_count=0;
open (LIS,"<$LIST") or die "cannot open $LIST \n";
while ( my $line=<LIS> ) {
   $pro_count++;
   if ( $pro_count%50==0 )  {   print "......working with $pro_count......\n";   } 
   chomp($line);
   my @data=split/\t/,$line;
   my $flag="NO";
   if ( exists($seq{$data[0]}) and exists($seq{$data[1]}) )   {
       my $query=$data[0];
       my $subject=$data[1];
       my $name="$query\_$subject";
          $name=~s/\W//g;
       my $align_seq="";
       my $queseq="";;
       my $hitseq="";
       my $percent=0;
       my $Widen=0;
       my $CC=0;
       my $HC=0;
       my $EC=0;
       my $HH=0;
       my $HE=0;
       my $EE=0;
       my $gaps=0;
       my $opennum=0;
       my @quechars;
       my @hitchars;
 
       open (TEMP,">tmp/$name.fasta") or die "cannot open $name.fasta \n";
          print TEMP ">$seq{$query}\n"; 
          print TEMP ">$seq{$subject}\n"; 
       close TEMP;

       $align_seq=`./pkg/seq-align/bin/needleman_wunsch --file tmp/$name.fasta --substitution_matrix ./pkg/seq-align/scoring/SSEalign_Matrix.txt --gapopen -2737 --gapextend -1157`;
   
       if ( $align_seq ne "" )  {
         $suc_count++;
         $flag="YES";
         my @seqs=split/\n/,$align_seq;
         $queseq=$seqs[0];
         $hitseq=$seqs[1];

         my @quegap=split/\-+/,$queseq;
         my @hitgap=split/\-+/,$hitseq;
         $opennum=$#quegap+$#hitgap;
         if ( substr($queseq,0,1) eq "-" )  { $opennum=$opennum-1;   }
         if ( substr($hitseq,0,1) eq "-" )  { $opennum=$opennum-1;   }   
       
         @quechars=split(//, $queseq);
         @hitchars=split(//, $hitseq);  
         my $i;
         for ($i=0; $i<$#quechars; $i++ ) {
 
           if (($quechars[$i] eq "-") or ($hitchars[$i] eq "-") ) {  $gaps++;   }

           if (($quechars[$i] eq "C") and ($hitchars[$i] eq "H") ) {  $HC++;   }
           if (($quechars[$i] eq "C") and ($hitchars[$i] eq "E") ) {  $EC++;   }
           if (($quechars[$i] eq "C") and ($hitchars[$i] eq "C") ) {  $CC++;   }

           if (($quechars[$i] eq "H") and ($hitchars[$i] eq "H") ) {  $HH++;   }
           if (($quechars[$i] eq "H") and ($hitchars[$i] eq "E") ) {  $HE++;   }
           if (($quechars[$i] eq "H") and ($hitchars[$i] eq "C") ) {  $HC++;   }

           if (($quechars[$i] eq "E") and ($hitchars[$i] eq "H") ) {  $HE++;   }
           if (($quechars[$i] eq "E") and ($hitchars[$i] eq "E") ) {  $EE++;   }
           if (($quechars[$i] eq "E") and ($hitchars[$i] eq "C") ) {  $EC++;   }   
         }

          $percent=int(($HH+$EE+$CC)/($HH+$EE+$CC+$HE+$HC+$EC+$gaps)*10000+0.5)/10000;
          $Widen=int((1.796*$HH+2.374*$EE+0.351*$CC)/(1.796*$HH+2.374*$EE+0.351*$CC+2.478*$HE+0.557*$HC+0.460*$EC+1.157*($gaps-$opennum)+2.737*$opennum)*10000+0.5)/100;

           print OUTTAB "$query\t$subject\t$gaps\t$opennum\t$percent\t$Widen\n";   
           my $description="Query=$query  Subject=$subject  Gapnum=$gaps  Gapopen=$opennum  percent=$percent  Widen=$Widen"; 
           print OUTDEL "$description\n\n$align_seq\n";   
        }  
     }   
     if ( $flag eq "NO" )   {  print OUTTAB "$line\tNA\n";   }
}
close LIS;

print "\n $suc_count out of total $pro_count SSE sequence alignments had been completed!\n";

close OUTTAB;
close OUTDEL;
