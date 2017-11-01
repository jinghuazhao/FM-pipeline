{
   SNP=$1;rsid=$2;pos=$3;a1=$4;a2=$5
   gsub(/:|_/," ",rsid)
   split(rsid,r," ")
   chr=r[1]
   pos=r[2]
   if(r[3]!=a1||r[4]!=a2) print "error"
   else if(r[3]>r[4]) {
      snpid=chr ":" pos "_" a2 "_" a1
      a1=r[4];a2=r[3]
   } else snpid=chr ":" pos "_" a1 "_" a2
   printf SNP " " snpid " " pos " " a1 " " a2
   M=(NF-5)/3
   for(i=0;i<M;i++) {
     i3=i*3
     if(r[3]<a[4]) printf " " $(i3+6) " " $(i3+7) " " $(i3+8)
     else printf " " $(i3+8) " " $(i3+7) " " $(i3+6)
   }
   printf "\n"
}
