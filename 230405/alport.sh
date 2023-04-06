mkdir alportrna
cd alportrna || exit
mkdir NJU8126
ln -sf /home/ivan/cater_data/NJU8126/NJU8126_R1.fq.gz \
  /home/ivan/fat/alportrna/NJU8126/R1.fq.gz
ln -sf /home/ivan/cater_data/NJU8126/NJU8126_R2.fq.gz \
  /home/ivan/fat/alportrna/NJU8126/R2.fq.gz

cd NJU8126 || exit
fastqc -t 4 ./*.fq.gz

trim_galore --nextera \
  --length 35 --cores 8 \
  --paired R1.fq.gz R2.fq.gz

fastqc -t 8 R*_val_*.fq.gz