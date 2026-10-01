conda activate TLS
python ./code/feature_selection/RNA_mrmr.py --config ./config_rna.json > "./results/PAAD/logs/RNA_mrmr_42.log"
python ./code/feature_selection/DNAmeth_mrmr.py --config ./config_meth.json > "./results/PAAD/logs/DNAmeth_mrmr_42.log"
python ./code/feature_selection/DNAmut_mrmr.py --config ./config_mut.json > "./results/PAAD/logs/DNAmut_mrmr_42.log"
python ./code/feature_selection/CNV_mrmr.py --config ./config_cnv.json > "./results/PAAD/logs/CNV_mrmr_42.log"
python ./code/feature_selection/ensemble_importance.py --config ./config_ensemble.json > "./results/PAAD/logs/ensemble_importance_mrmr_42.log"
python ./code/feature_selection/consensus_selection.py --config ./config_ensemble.json > "./results/PAAD/logs/consensus_selection_mrmr_42.log"

python ./code/modeling/final_model.py --config ./config_final_model.json > ./results/PAAD/logs/final_model.log
python ./code/modeling/final_model.py --config ./config_final_model.json > ./results/all/logs/final_model_consensus.log


sh run_seed.sh ./config_rna.json 0 9 rna all
sh run_seed.sh ./config_meth.json 0 9 meth all
sh run_seed.sh ./config_mut.json 0 9 mut all
sh run_seed.sh ./config_cnv.json 0 9 cnv all
sh run_seed.sh ./config_ensemble.json 0 9 ensemble all



