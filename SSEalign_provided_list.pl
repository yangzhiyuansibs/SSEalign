#!/usr/bin/perl
# Description: This script is for Secondary Structure Element Alignment for provided list of protein-protein pair
# usage example: SSEalign_provided_list.pl <total_SSE.fasta> <provided_list.txt> <output_table.txt> <output_align.txt>
# Author: YANG Zhiyuan ( yangzhiyuan@link.cuhk.edu.hk )
# Institute: The Chinese University of Hong Kong
# Date: 2017-03-23

use strict;
use warnings;

# check the input parameter
 if ( $#ARGV<3 )  {   die "\n Please input 4 parameter as following:\nperl SSEalign_provided_list.pl <total_SSE.fasta> <provided_list.txt> <output_table.txt> <output_align.txt> \n\n Other suggestions: \n1.<The total_SSE.fasta> must include all SSE sequence used subsequently in FASTA format; \n2.The name contained NON-word character such as sp|P0ABI8|CYOB_ECOLI is not allowed in <The total_SSE.fasta>; \n3.<provided_list.txt> file must be TAB separated, if multiple colunm is provided, the first 2 columns are used. \n4.<output_table.txt> is the result in table format and <output_align.txt> contain detailed alignment information. \n\n";    }
 
 # check the file
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

{
local $/="\n>";
open (QUE,"<$QUERY"); 
 while (my $content=<QUE>)  {
  $content=~s/>//g;
  my $loc=index($content,"\n",0);
  my $firstline=substr($content,0,$loc-1);
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
       my $name="$query.$subject";
       my $queseq='';;
       my $hitseq='';
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
 
       open (TEMP,">tmp/$query.fasta") or die "cannot open $query.fasta \n";
          print TEMP ">Query.$seq{$query}"; 
       close TEMP;
       open (TEMP2,">tmp/$subject.fasta") or die "cannot open $subject.fasta \n";
          print TEMP2 ">Subje.$seq{$subject}"; 
       close TEMP2;
       `stretcher -asequence tmp/$query.fasta -bsequence tmp/$subject.fasta -outfile tmp/$name.stre -auto Y;`;
   
       if (-e "tmp/$name.stre" )  {
       $suc_count++;
       $flag="YES";
      open (CAL,"<tmp/$name.stre"); 
       while (my $line=<CAL>) {
        if ( $line=~/^Query\.\s+(\S+)\n/ ) {  $queseq=$queseq.$1;   }
        if ( $line=~/^Subje\.\s+(\S+)\n/ ) {  $hitseq=$hitseq.$1;   } 
       }
      close CAL;
    
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

           print OUTTAB "$query\t$subject\t$HH\t$EE\t$CC\t$HE\t$HC\t$EC\t$gaps\t$opennum\t$percent\t$Widen\n";   
           my $description="Query=$query  Subject=$subject  Gapnum=$gaps  Gapopen=$opennum  percent=$percent  Widen=$Widen";  
           `sed -i \"32i $description\" tmp/$name.stre`;
           `cat tmp/$name.stre >> $OUTDEL`; 
  
        }   else  {

        }
    }
     if ( $flag eq "NO" )   {  print OUTTAB "$line\tNA\n";   }
}
close LIS;

print "\n $suc_count out of total $pro_count sequence alignments had been done!\n";

close OUTTAB;
close OUTDEL;



 