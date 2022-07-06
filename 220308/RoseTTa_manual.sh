for i in mcr mgr-d mgr-pl mgr-s mgr-ul mgr mhr mpr uur mgr-p; do
  mkdir ${i}
done

for i in mcr mgr-d mgr-pl mgr-s mgr-ul mgr mhr mpr uur mgr-p; do
  bsub -n 24 -J "${i}_1" "bash ../step1.sh ${i}.fa ${i}"
done

for i in mcr mgr-d mgr-pl mgr-s mgr-ul mgr mhr mpr uur mgr-p; do
  bsub -n 24 -J "${i}_2" "bash ../step2.sh ${i}.fa ${i}"
done

for i in mcr mgr-d mgr-pl mgr-s mgr-ul mgr mhr mpr uur mgr-p; do
  bsub -q gpu_v100 -gpu "num=1" -J "predict_pyRosetta" "bash ../step3.sh ${i}.fa ${i}"
done

for i in mcr mgr-d mgr-pl mgr-s mgr-ul mgr mhr mpr uur mgr-p; do
  bash ../step4.sh ${i}
done

for i in mcr mgr-d mgr-pl mgr-s mgr-ul mgr mhr mpr uur mgr-p; do
  bsub -n 24 -J "${i}_5" "bash ../step5.sh ${i}.fa ${i}"
done

for i in mcr mgr-d mgr-pl mgr-s mgr-ul mgr mhr mpr uur mgr-p; do
  bsub -n 24 -J "${i}_6" "bash ../step6.sh ${i}.fa ${i}"
done

for i in mcr mgr-d mgr-pl mgr-s mgr-ul mgr mhr mpr uur mgr-p; do
  bsub -n 24 -J "${i}_7" "bash ../step7.sh ${i}.fa ${i}"
done

for i in gde00 gde02 yggl02 yggl05; do
  bsub -n 24 -J "${i}_1" "bash ../step1.sh ${i}.fa ${i}"
done

for i in gde00 gde02 yggl02 yggl05; do
  bsub -n 24 -J "${i}_2" "bash ../step2.sh ${i}.fa ${i}"
done

for i in gde00 gde02 yggl02 yggl05; do
  bsub -q gpu_v100 -gpu "num=1" -J "predict_pyRosetta_${i}" "bash ../step3.sh ${i}.fa ${i}"
done

for i in gde00 gde02 yggl02 yggl05; do
  bash ../step4.sh ${i}
done

for i in gde00 gde02 yggl02 yggl05; do
  bsub -n 24 -J "${i}_5" "bash ../step5.sh ${i}.fa ${i}"
done

for i in gde00 gde02 yggl02 yggl05; do
  bsub -n 24 -J "${i}_6" "bash ../step6.sh ${i}.fa ${i}"
done

for i in gde00 gde02 yggl02 yggl05; do
  bsub -n 24 -J "${i}_7" "bash ../step7.sh ${i}.fa ${i}"
done

for i in brab braz; do
  bsub -n 24 -J "${i}_1" "bash ../step1.sh ${i}.fa ${i}"
done

for i in brab braz; do
  bsub -n 24 -J "${i}_2" "bash ../step2.sh ${i}.fa ${i}"
done

for i in brab braz; do
  bsub -q gpu_v100 -gpu "num=1" -J "predict_pyRosetta_${i}" "bash ../step3.sh ${i}.fa ${i}"
done

for i in brab braz; do
  bash ../step4.sh ${i}
done

for i in brab braz; do
  bsub -n 24 -J "${i}_5" "bash ../step5.sh ${i}.fa ${i}"
done

for i in brab braz; do
  bsub -n 24 -J "${i}_6" "bash ../step6.sh ${i}.fa ${i}"
done

for i in brab braz; do
  bsub -n 24 -J "${i}_7" "bash ../step7.sh ${i}.fa ${i}"
done

for PREFIX in 284 327 352 360 374 384 470 479; do
  bsub -n 24 -J "LC${PREFIX}" "
    ../bismark/bismark --bowtie2 \
      --bam -o LC${PREFIX} \
      --genome hg38_chrom_bismark/ \
      -1 LC${PREFIX}/*_R1_val_1.fq.gz \
      -2 LC${PREFIX}/*_R2_val_2.fq.gz \
      -p 8 --nucleotide_coverage \
      --temp_dir LC${PREFIX}/ \
      1>LC${PREFIX}/bismark.log 2>&1
"
done


