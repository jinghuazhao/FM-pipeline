{
  gsub(/:|_/," ",$2)
  split($2,a," ")
  chr=a[1]
  pos=a[2]
  if(NR>1) {
    if(a[3]<a[4]) snpid=sprintf("%d:%d_%s_%s",a[1],a[2],a[3],a[4])
    else snpid=sprintf("%d:%d_%s_%s",a[1],a[2],a[4],a[3])
    $2=snpid
  }
  print
}
