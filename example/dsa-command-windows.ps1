$dir = ".\"
$stems = "S16R", "S17R"
foreach ($stem in $stems) { `
    & 'C:\Program Files\ccb\dsa\dsa-win64.exe' `
    -g 2 `
    -c horizontal `
    -f TACnnnnnnAGTnnnnnnCGNNNNNNNNNNNGTCCTCTCT `
    -r CCCnnnnnnnnnACGACAAN `
    --template_dna=AAGAAAGTGGTGTTGGCCAAAAAAGGCGATACCGTGGAGCTGACCTGCACCGCAAGCCAGAAGAAGAACATCCAGTTCCACTGGAAGAACTCCAACCAGATCAAGATCCTGGGCAACCAGGGCAGCTTCTTGACCAAGGGACCTAGCAAGCTGAATGACAGAGTGGACTCTCGGAGGAGCCTGTGGGATCAAGGCAACTTCCCTCTGATCATAAAAAACCTGAAGATAGAGGATAGTGATACCTACATCTGTGAAGTGGAAGATCAGAAGGAGGAAGTGCAGCTGCTGGTGTTCGGGCTGACTGCTAACTCCGATACCCATCTgCTgCAGGGGCAGAGCCTAACACTGACACTGGAGAGCCCTCCTGGCAGCAGCCCAAGCCTGCAaTGCCGCAGCCCtGGAGGCAAGAACATCCAAGGTGGCAAAACCCTTTCTGTCAGCCAGCTGGAACTGCAGGATTCTGGAACCTGGACATGTACAGTGCTGCAGGATCAGAAAACCTTGGAGTTCAAGATTGAT `
  $([string]::concat($dir, "", $stem, "2.fastq")) `
  $([string]::concat($dir, "", $stem, "1.fastq")) `
  | Out-File -FilePath $([string]::concat($stem, "_dsa_output.csv")) -Encoding utf8 `
}