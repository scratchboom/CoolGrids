#!/usr/bin/perl
use strict;
my $amline_root='/serxa/private/www/amline_1.6.2.1/amline';
my $amline_winroot = `cygpath -w $amline_root`;
chomp $amline_winroot;
my $link_root = 'file:///' . $amline_winroot;
$link_root =~ s/\\/\//g;

sub usage() {
    print "plot amline chart\n";
}

sub add_xid() {
    $_[0]->[0]{$_[1]} = $_[2];
}

if (@ARGV < 1) {
    usage();
    exit(1);
}

my $outfile = $ARGV[0];
my $title = 'Title';

# Parse format string
#my $COL_RE='[0-9]+';
#my $TITLE_RE='[^\s]+';
#my $format_str = defined($ARGV[1])? $ARGV[1]: "1:2 $title";
#my %fmt = ();
#foreach (split(';',$format_str)) {
#    m/^($COL_RE):($COL_RE)\s+($TITLE_RE)\s*$/;
#    $fmt{$1} = [$2, $3];
#}

#my $data = [{},{}];
#my $xid = 0;
#while(<>) {
#    chomp;
#    my @f = split /\t/;
#    my $i = 1;
#    foreach my $xcol (sort {$a <=> $b} keys %fmt) {
#	add_xid($data, $xid, $f[$xcol]);
#	$xid++;
#    }
#}

my $csv_data='';
while(my $line = <STDIN>) {
    chomp $line;
    $line =~ s/\t/;/g;
    $line =~ s/;$//g;
    $line =~ s/^;//g;
    $csv_data .= $line . '\n';
}

open OUT, "> $outfile";
print OUT <<HTML;
<html>
<head>
<meta http-equiv="Content-Type" content="text/html; charset=UTF-8" />
<title>$title</title>
</head>
<body>
  <script type="text/javascript" src="$link_root/swfobject.js"></script>
    <div id="flashcontent"><strong>You need to upgrade your Flash Player</strong></div>
      <script type="text/javascript">
	// <![CDATA[
	var so = new SWFObject("$link_root/amline.swf", "amline", "1430", "770", "8", "#FFFFFF");
//	so.addVariable("path", "$link_root");
//	so.addVariable("settings_file", encodeURIComponent("amline/amline_settings.xml"));
//	so.addVariable("data_file", encodeURIComponent("amline/amline_data.xml"));
	so.addVariable("chart_settings", encodeURIComponent("<settings><data_type>csv</data_type></settings>"));
	so.addVariable("chart_data", encodeURIComponent("$csv_data"));
//	so.addVariable("additional_chart_settings", encodeURIComponent("<settings>...</settings>"));
//      so.addVariable("loading_settings", "LOADING SETTINGS");
//      so.addVariable("loading_data", "LOADING DATA");
//	so.addVariable("preloader_color", "#999999");
	so.write("flashcontent");
	// ]]>
      </script>
</body>
</html>
HTML

close OUT
