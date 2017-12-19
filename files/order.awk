{
   rsid=$2
   pos=$3
   a1=$4
   a2=$5
   if (a1>a2) {
      snpid=chr ":" pos "_" a2 "_" a1 
      printf rsid " " snpid " " pos " " a2 " " a1
      for (i=0;i<(NF-5)/3;i++)
      {
          i3=i*3
          printf " " $(i3+8) " " $(i3+7) " " $(i3+6)
      }
   } else {
      snpid=chr ":" pos "_" a1 "_" a2
      printf rsid " " snpid " " pos " " a1 " " a2
      for (i=0;i<(NF-5)/3;i++)
      {
          i3=i*3
          printf " " $(i3+6) " " $(i3+7) " " $(i3+8)
      }
   }
   printf "\n"
}
