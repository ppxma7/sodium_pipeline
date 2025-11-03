# Sodium pipeline
## Prereq

- You will need a bunch of Python libraries and FSL


- Run `batch_sodium.sh` which runs `run_sodium_pipeline_nascar.sh`
- This does 3 python scripts: 

```
run_sodium_pipeline_nascar.py
run_sodium_pipeline_nascar_to_mni.py 
run_sodium_pipeline_atlasread.py 
```

- The above registers sodium files together, moves everything to MNI space, then reads in an atlas. It also moves things to native space.

- Two final scripts operate on the group results:

```
run_sodium_pipeline_nascar_calcmean.py
run_sodium_pipeline_nascar_fastMNImeans.py
```

- These calculate the across subjects mean and stdev, and also the across sites mean and stdev, and also gets the values from the FAST outputs.
