for i in mgr-p
do
bsub -n 24 -J "${i}_2" "bash ../step1.sh ${i}.fa ${i}"
done

for i in mgr-p
do
bsub -n 24 -J "${i}_2" "bash ../step2.sh ${i}.fa ${i}"
done

bsub -q gpu_v100 -gpu "num=1" -J "predict_pyRosetta" 'bash ../step3.sh MgR_mut.fa MgR_mut'

mkdir -p MgR_mut/pdb-3track WDIR=`realpath -s MgR_mut`
PIPEDIR=`dirname /share/home/wangq/wyf/RoseTTAFold/step3.sh`
IN=MgR_mut.fa
for m in 0 1 2
do
for p in 0.05 0.15 0.25 0.35 0.45
do
for ((i=0;i<1;i++))
do
if [ ! -f $WDIR/pdb-3track/model${i}_${m}_${p}.pdb ]; then
echo "python -u $PIPEDIR/folding/RosettaTR.py --roll -r 3 -pd $p -m $m -sg 7,3 $WDIR/t000_.3track.npz $IN $WDIR/pdb-3track/model${i}_${m}_${p}.pdb"
fi
done
done
done > $WDIR/parallel.fold.list

bsub -n 24 -J "RosettaTR" '
source /share/home/wangq/miniconda3/etc/profile.d/conda.sh
conda activate folding
parallel -j 16 \
< /share/home/wangq/wyf/RoseTTAFold/example/MgR_mut/parallel.fold.list \
> /share/home/wangq/wyf/RoseTTAFold/example/MgR_mut/log/folding.stdout \
2> /share/home/wangq/wyf/RoseTTAFold/example/MgR_mut/log/folding.stderr
'

bsub -n 24 -J "DAN-msa" 'bash ../step6.sh MgR_mut.fa MgR_mut'

bsub -n 24 -J "pick_final_models" 'bash ../step7.sh MgR_mut.fa MgR_mut'