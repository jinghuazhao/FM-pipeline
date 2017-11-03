{
  chr=chr
  pos=$3
  if(NR>1) {
    if($4<$5) snpid=sprintf("%d:%d_%s_%s",chr,pos,$4,$5)
    else snpid=sprintf("%d:%d_%s_%s",chr,pos,$5,£4)
    $2=snpid
  }
  print
}
