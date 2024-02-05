#!/bin/bash
./Phi_preselection.sh
cd MVA/
./MVA_Phi.sh 
cd ..
./Phi_BDT.sh
./fit_phi.sh
