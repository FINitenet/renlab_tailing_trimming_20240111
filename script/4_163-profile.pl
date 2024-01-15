#!usr/bin/perl -w

my $metaFile = $ARGV[0];


open(IN,$metaFile)||die "$metaFile\n";

while($raw=<IN>)
{
	$raw=~s/\n|\r//g;
	@inf3=split(/\s+/,$raw);
	$miRNA=$inf3[0];
	$start=$inf3[2];
	$leng=$inf3[3];
	$end=$start+$leng-1;
	@inf=split(/-|\./,$miRNA);
	$hairpin=$inf[0]."-".$inf[1];
		
	$samFile=$ARGV[1];
	$outFile="tmp.txt";
	read_profile($samFile,$miRNA,$start,$leng,$end,$hairpin,$outFile);
}
close(IN);

sub read_profile
{
	my ($samFile,$miRNA,$start,$leng,$end,$hairpin,$outFile)=@_;
	my @inf=();
	my @infName=();
	my @infNt=();
	my $end5GCM=0;
	my $key=0;
	my $raw;
	my $i=0;
	my $total=0;
	my $total1=0;
	my %profile=();
	
	open(SAM,$samFile)||die "Can not open file: $samFile\n";
	while($raw=<SAM>)
	{
		if((!($raw=~/^@/))&&$raw=~/$hairpin/)
		{
			@inf=split(/\s+/,$raw);
			@infName=split(/-+/,$inf[0]);
			$end5GCM=$inf[3]+(length($infName[2])-$infName[3])-1;
			if($end5GCM>$end)
			{
				$end5GCM=$end;     
				$infName[3]=$inf[3]+length($infName[2])-$end-1;	
			}
			if(abs($inf[3]-$start)==0&&$infName[3]<=10&&$end-$end5GCM<=10) #####5GCM+-2 and tailling <=7 and tailed <=7
			{
				$total++;
				$profile{($infName[3])."-".($end-$end5GCM)}+=1;
			}
		}		
	}
	close(SAM);

	for($i=10;$i>=0;$i--)
	{
		for($j=10;$j>=0;$j--)
		{
			if(exists($profile{($i)."-".($j)}))
			{
				print $profile{($i)."-".($j)},"\t";
			}
			else
			{
				print "0\t";
			}
		}
		print "\n";
	}	
	print "$miRNA\t$total\t$total1\n";	
	open(OUT,">".$outFile);
	{
		print OUT "$miRNA\t$total\t$total1\n";
		for($i=10;$i>=0;$i--)
		{
			for($j=10;$j>=0;$j--)
			{
				if(exists($profile{($i)."-".($j)}))
				{
					print OUT  $profile{($i)."-".($j)},"\t";
				}
				else
				{
					print OUT  "0\t";
				}
			}
			print OUT  "\n";
		}			
	}
	close(OUT);
}

