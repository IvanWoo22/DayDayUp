bsub -n 24 -q largemem -J "msa" '
    export PIPEDIR=/share/home/wangq/wyf/RoseTTAFold/
    source /share/home/wangq/miniconda3/etc/profile.d/conda.sh
    conda activate RoseTTAFold
    /share/home/wangq/wyf/RoseTTAFold/input_prep/make_msa.sh \
        /share/home/wangq/wyf/RoseTTAFold/example/MgR_mut.fa \
        /share/home/wangq/wyf/RoseTTAFold/example/MgR_mut \
        22 250 > /share/home/wangq/wyf/RoseTTAFold/example/MgR_mut/log/make_msa.stdout \
        2> /share/home/wangq/wyf/RoseTTAFold/example/MgR_mut/log/make_msa.stderr
'

bsub -n 24 -J "ss2" '
    source /share/home/wangq/miniconda3/etc/profile.d/conda.sh conda activate RoseTTAFold
    BLASTMAT=/share/home/wangq/wyf/RoseTTAFold/example/blast-2.2.26/data/
    PIPEDIR=/share/home/wangq/wyf/RoseTTAFold/
    /share/home/wangq/wyf/RoseTTAFold/input_prep/make_ss.sh \
        /share/home/wangq/wyf/RoseTTAFold/example/MgR_mut/t000_.msa0.a3m \
        /share/home/wangq/wyf/RoseTTAFold/example/MgR_mut/t000_.ss2 \
        > /share/home/wangq/wyf/RoseTTAFold/example/MgR_mut/log/make_ss.stdout \
        2> /share/home/wangq/wyf/RoseTTAFold/example/MgR_mut/log/make_ss.stderr
'

bsub -n 24 -J "hhsearch" 'bash ../step2.sh MgR_mut.fa MgR_mut'

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