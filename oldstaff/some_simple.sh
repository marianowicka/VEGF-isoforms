for entry in b*
do
  echo tail -n 7617 $entry | bash > v$entry
done

cat v* > wyn_b.txt
