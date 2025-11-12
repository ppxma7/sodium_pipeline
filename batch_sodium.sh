#!/bin/bash

mypath="/Users/ppzma/Documents/MATLAB/mycode/"
cd "$mypath" || exit
source heph_venv/bin/activate

# sh $mypath/sodium_pipeline/run_sodium_pipeline_nascar.sh Subject1 radial site1
# echo "⏸️ Pausing 2 seconds..."
# sleep 2
# sh $mypath/sodium_pipeline/run_sodium_pipeline_nascar.sh Subject1 floret site2
# echo "⏸️ Pausing 2 seconds..."
# sleep 2
# sh $mypath/sodium_pipeline/run_sodium_pipeline_nascar.sh Subject1 seiffert site3
# echo "⏸️ Pausing 2 seconds..."
# sleep 2

# sh $mypath/sodium_pipeline/run_sodium_pipeline_nascar.sh Subject2 seiffert site1
# echo "⏸️ Pausing 2 seconds..."
# sleep 2
# sh $mypath/sodium_pipeline/run_sodium_pipeline_nascar.sh Subject2 radial site2
# echo "⏸️ Pausing 2 seconds..."
# sleep 2
# sh $mypath/sodium_pipeline/run_sodium_pipeline_nascar.sh Subject2 floret site3
# echo "⏸️ Pausing 2 seconds..."
# sleep 2

# sh $mypath/sodium_pipeline/run_sodium_pipeline_nascar.sh Subject3 floret site1
# echo "⏸️ Pausing 2 seconds..."
# sleep 2
# sh $mypath/sodium_pipeline/run_sodium_pipeline_nascar.sh Subject3 seiffert site2
# echo "⏸️ Pausing 2 seconds..."
# sleep 2
# sh $mypath/sodium_pipeline/run_sodium_pipeline_nascar.sh Subject3 radial site3
# echo "⏸️ Pausing 2 seconds..."
# sleep 2

# sh $mypath/sodium_pipeline/run_sodium_pipeline_nascar.sh Subject4 radial site1
# echo "⏸️ Pausing 2 seconds..."
# sleep 2
sh $mypath/sodium_pipeline/run_sodium_pipeline_nascar.sh Subject4 floret site2
echo "⏸️ Pausing 2 seconds..."
sleep 2
# sh $mypath/sodium_pipeline/run_sodium_pipeline_nascar.sh Subject4 seiffert site3
# echo "⏸️ Pausing 2 seconds..."
# sleep 2

# echo "FINISHED"
