for i in mcr mgr-d mgr-pl mgr-s mgr-ul mgr mhr mpr uur mgr-p; do
  bsub -n 24 -J "${i}_2" "bash ../step1.sh ${i}.fa ${i}"
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

bsub -n 24 -J "RosettaTR" '
source /share/home/wangq/miniconda3/etc/profile.d/conda.sh
conda activate folding
parallel -j 16 \
< /share/home/wangq/wyf/RoseTTAFold/example/MgR_mut/parallel.fold.list \
> /share/home/wangq/wyf/RoseTTAFold/example/MgR_mut/log/folding.stdout \
2> /share/home/wangq/wyf/RoseTTAFold/example/MgR_mut/log/folding.stderr
'

for i in mcr mgr-d mgr-pl mgr-s mgr-ul mgr mhr mpr uur mgr-p; do
  bsub -n 24 -J "${i}_6" 'bash ../step6.sh ${i}.fa ${i}'
done

for i in mcr mgr-d mgr-pl mgr-s mgr-ul mgr mhr mpr uur mgr-p; do
  bsub -n 24 -J "${i}_7" 'bash ../step7.sh ${i}.fa ${i}'
done
