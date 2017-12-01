{
   rsid=$2
   pos=$3
   a1=$4
   a2=$5
   if(a1>a2) snpid=chr ":" pos "_" a2 "_" a1
   else snpid=chr ":" pos "_" a1 "_" a2
   printf snpid " " $2 " " pos " " a1 " " a2
   M=(NF-5)/3
   for(i=0;i<M;i++) {
     i3=i*3
     if(a1]<a2) printf " " $(i3+6) " " $(i3+7) " " $(i3+8)
     else printf " " $(i3+8) " " $(i3+7) " " $(i3+6)
   }
   printf "\n"
}
