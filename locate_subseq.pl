
# Author: Hengyi Jiang <hengyi.jiang@gmail.com>
# purpose to calculate a subsequence described by regexp from a huge genome sequence
# written by JHY, 2018-7-12
# output format like below:
# >Genome Name_start_end
# matched_subseq

# use Digest::MD5 qw(md5 md5_hex md5_base64);

 sub sub_regexp_pos{
 	my $name = shift; # seq name
 	my $text = shift; # reference to string, which is the sequence
 	my $tg = shift;   # target sequence
 	my $e, $s;        # $e start of next search, $s start of the match
	my @r = ();       # return value

 	return @r if $name eq "" or $$text eq "" or $tg eq "" ; # empty array

	while($$text =~ m/$tg/ig){
		$e = pos($$text); # start of next search
		$s= $e - length($&)+1; 
		$m =$&;
		push(@r,[substr($name,1), $s,$e,$m]);
	}
	return @r;
	# @r is  array reference, this array include :
	#  start_of_target, end_of_target, matched_target
 }

 sub format_print{
 	my $input =shift; 
 	# input is a reference, because passing array as parameter only do shallow copy
 	foreach $h(@$input){
 		if ($h != undef) {
 			print ">";
			print $h->[0]."_";
			print $h->[1]."_";
			print $h->[2]."\n";
			print $h->[3]."\n";
		}
	}
 }

my $file_name = shift;
my $target = shift;
$usage = "perl locate_subseq.pl path_of_file TARGET_SEQ_PATTERN\n";
die "$usage" unless -f $file_name;
die "$usage" unless $target;

open(FASTA,$file_name) or die "open $file_name failed\n";

$started = 0;
$seq="";
$name="";


while(<FASTA>){
	chomp;
	if(/^\s*>.+/){
		if($started){
			@hits = sub_regexp_pos($name,\$seq,$target);
			format_print(\@hits);
			$name = $_;
			$seq = "";
		}else{
			$name = $_;
			$started = 1;
		}
		$n++;
	}else{
		next unless m/\w+/;
		next unless $started;
		$seq .=$_;
	}
}

if($seq){
	@hits = sub_regexp_pos($name,\$seq,$target);
	format_print(\@hits);
	$seq = "";
}