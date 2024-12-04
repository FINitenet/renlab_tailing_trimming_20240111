#!usr/bin/perl -w

$sample=$ARGV[0];

%end=();
open(IN,"/bios-store1/chenyc/scripts/Tailing_Trimming/miRNA_start_sequence_length_new.txt")||die "Can not open file: /bios-store1/chenyc/scripts/Tailing_Trimming/miRNA_start_sequence_length_new.txt\n";
while($raw=<IN>)
{
	@inf=split(/\s+/,$raw);
	@inf2=split(/-|\./,$inf[0]);
	$hairpin=$inf2[0]."-".$inf2[1];
	$end{$hairpin."-$inf[2]"}=$inf[2]+length($inf[1])-1;
	# print $hairpin."-$inf[2]\t",$inf[2]+length($inf[1])-1,"\n";				
}

close(IN);


%lengthRead=();
%length5GCM=();
$total=0;
%readName=();
open(IN,$sample)||die "Can not open file: $sample\n";
while($raw=<IN>)
{
	if(!($raw=~/^@/))
	{
		@inf=split(/\s+/,$raw); # Sample-7865227-AAACGCAAAGAAATGG--3	0	ath-MIR156a	7	255	16M	*	0	0	AAACGCAAAGAAATGG	IIIIIIIIIIIIIIII	XA:i:0	MD:Z:13	NM:i:0
		@inf2=split(/-+/,$inf[0]); # Sample 7865227 AAACGCAAAGAAATGG 3
		$lenRead=length($inf2[2]);
		$end5GCM=$inf[3]+$lenRead-$inf2[3]-1;
		$len5GCM=+$lenRead-$inf2[3];
		if(exists($end{$inf[2]."-".$inf[3]})&&(!($readName{$inf2[0]."-".$inf2[1]})))
		{
			if($inf2[3]==0&&$end5GCM>$end{$inf[2]."-".$inf[3]})
			{
				$len5GCM=length($inf2[2]);
			}		
			$lengthRead{$lenRead}++;
			$length5GCM{$len5GCM}++;
			$total++;
		}
		$readName{$inf2[0]."-".$inf2[1]}+=1;
		
	}
	if($total%10000==0&&$total>0)
	{
		#print "$total\n";
	}
}
close(IN);

print "total reads mapped to miRNA: ",$total,"\n";

print "\ndistribution of whole reads\n";
foreach $key (sort(keys(%lengthRead)))
{
	print "$key\t",$lengthRead{$key},"\t";
	printf("%.4f\n",$lengthRead{$key}/$total); 
}
print "\ndistribution of 5GMC\n";
foreach $key (sort(keys(%length5GCM)))
{
	print "$key\t",$length5GCM{$key},"\t";
	printf("%.4f\n",$length5GCM{$key}/$total); 
}
=cut
