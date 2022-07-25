for i in {1..9}; do
  site1=$(sed -n "${i}p;${i}q" tmp.lst | cut -f 5)
  site2=$(sed -n "${i}p;${i}q" tmp.lst | cut -f 6)
  site3=$(sed -n "${i}p;${i}q" tmp.lst | cut -f 7)
  site4=$(sed -n "${i}p;${i}q" tmp.lst | cut -f 8)
  sed -n "${i}p;${i}q" tmp.lst | perl snp4fasta.pl -f rawseq.fa --stdin >seq"${i}".fa
  RNAfold <seq"${i}".fa >seq"${i}".fold
  RNAplot --pre "${site1} ${site2} 8 0.33 1 0.33 omark ${site3} ${site4} 8 0 1 1 omark" <seq"${i}".fold
  rm seq"${i}".fold seq"${i}".fa
  ls ./*_ss.ps | xargs -n1 -i{} mv {} ps/"${i}"_{}
done
