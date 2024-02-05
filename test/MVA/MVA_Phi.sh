#!/bin/bash
python train_trees.py phi
python create_BDT_plots.py phi
python BDT_significance.py phi
