# DREAM_2024_ANTS

Source code for the Placental Methylation Clock – DREAM 2024 Challenge. Wiki: https://www.synapse.org/Synapse:syn61846522/wiki/629109

This repository contains R scripts used to develop and run the final models for the DREAM 2024 placental methylation clock challenge.

    SC1_synapse_model.R – Final model for the SC1 sub-challenge, incorporating correlation-based feature selection and k-means-based feature clustering.

    Run_SC1_synapse_model.R – Runs the SC1 model on a new dataset for gestational age prediction.

    SC2_paper_model.R – Final model for the SC2 sub-challenge, incorporating correlation-based and EWAS-based feature selection.

    Run_SC2_paper_model.R – Runs the SC2 model on a new dataset for gestational age prediction.

Input format

    A .csv file where: Rows = CpG probe names (cgxxxxxxxx); Columns = Sample IDs

Output result

    Predicted gestational ages (in months) saved to predictions.csv.
