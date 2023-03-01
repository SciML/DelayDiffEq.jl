$pdf_mode = 1;
$bibtex_use = 2;

sub build_header {
  system("ruby ./prep.rb")
}

build_header()
