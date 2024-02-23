#!/usr/bin/env bash
source /share/home/wangq/miniconda3/etc/profile.d/conda.sh
conda activate RF2
export PATH=/share/gpu-apps/cuda/11.7/bin:$PATH
export LD_LIBRARY_PATH=/share/gpu-apps/cuda/11.7/lib64:$LD_LIBRARY_PATH
export LIBRARY_PATH=/share/gpu-apps/cuda/11.7/lib64:$LIBRARY_PATH

echo "Running RoseTTAFold2 to predict structures"

python /share/home/wangq/wyf/RoseTTAFold2/network/predict.py \
	-inputs 7UGF/rcsb_pdb_7UGF_1.msa0.a3m -prefix 7UGF/models/model \
	-model /share/home/wangq/wyf/RoseTTAFold2/network/weights/RF2_apr23.pt \
	-db /share/home/wangq/wyf/RoseTTAFold2/pdb100_2021Mar03/pdb100_2021Mar03 -symm C1
python /share/home/wangq/wyf/RoseTTAFold2/network/predict.py \
	-inputs "${1}"/"${1}"_1.msa0.a3m \
	-prefix "${1}"/models/"${1}" \
	-model /share/home/wangq/wyf/RoseTTAFold2/network/weights/RF2_apr23.pt \
	-db /share/home/wangq/wyf/RoseTTAFold2/pdb100_2021Mar03/pdb100_2021Mar03 -symm C1

echo "Done"

python /share/home/wangq/wyf/RoseTTAFold2/network/predict.py \
	-inputs "${1}"/"${1}"_1.msa0.a3m:"${1}"/"${1}"_1.hhr:"${1}"/"${1}"_1.atab \
	-prefix "${1}"/models/"${1}" \
	-model /share/home/wangq/wyf/RoseTTAFold2/network/weights/RF2_apr23.pt \
	-db /share/home/wangq/wyf/RoseTTAFold2/pdb100_2021Mar03/pdb100_2021Mar03 -symm C1
