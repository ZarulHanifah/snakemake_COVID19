awk '{
    arr[$1]+=$2
   }
   END {
     for (key in arr) printf("%s\t%s\n", key, arr[key])
   }'
