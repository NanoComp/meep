#!/usr/bin/perl -w

my $amcode = 0;

while (<>) {
  if ($amcode && /^(\\end\{comment\}|\\end\{verbatim\})$/) {
    $amcode = 0
  } elsif ($amcode) {
    print $_;
  } elsif (/^(\\begin\{comment\}|\\begin\{verbatim\})$/) {
    $amcode = 1;
  }
}
