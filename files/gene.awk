{
  if(NR==1) print "SNPID CHR POS Z F N ens_coding_exons ens_noncoding_exons tssdist syn nonsyn SEGNUMBER"
  if($2==$9) print $2,$3,$1,$4,$5,$6,$10,$11,$12,$13,$14,$7
}
