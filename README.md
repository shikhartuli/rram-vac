# RRAM-VAC: A Variability-Aware Controller for RRAM-based Memory Architectures
 
This is the set of codes to run and test the RRAM-VAC (MMU) on MATLAB. Following is the list of codes and test data and about what they do:

## The different test data (memory traces):
1. reportConv
2. reportCS_full - Complete CS test
3. reportCS_full_2 - Complete CS test with a register size of 2
4. reportCS_full_4 - Complete CS test with a register size of 4
5. reportCS_full_8 - Complete CS test with a register size of 8
6. reportCS_full_16 - Complete CS test with a register size of 16
7. reportCS_full_32 - Complete CS test with a register size of 32
8. reportDT - Decision Tree
9. reportDT_cached - Decision Tree with a cache of 8 words
10. reportDT_cached_2 - Decision Tree with a cache of 16 words
11. reportFE - Feature Extraction
12. reportFE+DT - Complete ML algo with FE and DT
13. reportMM - Matrix Multiplication 16x20 vs 20x16
14. reportMM_new - Matrix Multiplication 30x30 vs 30x30

## Main codes:
1. MMU_rw_parallel_wEnergy.m - outputs stall times, total process time and read/write energies
2. MMU_rw_parallel_wEnergy_winst_new - outputs transient times and energies for every batch as well (for the transients curve in he paper)

## Testing codes:
1. memTrace_to_rwArray.m - Extracts performance and energy gains figure for all applications (Fig. 8 in ASP-DAC 2020)
2. MMU_rw_test_bounded_waitBuffer_vs_writeBuffer - gets the stall times contour plots to determine optimal waitBuffer and bach sizes (Fig. 6 in ASP-DAC 2020)
3. MMU_rw_test_perfGains_vs_rwBurstSize - Extracts the performance improvement with increasing read/write burst size (Fig. 9 in ASP-DAC 2020)

## Debugging and hacking RRAM-VAC

Other testing codes including the MMU codes for serial write, are available in the folder "Other Codes (old)"


## Developer

[Shikhar Tuli](https://www.github.com/shikhartuli) (shikhartuli98@gmail.com)

## Cite this work
If you use our static model, please cite:
```
@INPROCEEDINGS{tuli2020_rram-vac,
  author={S. {Tuli} and M. {Rios} and A. {Levisse} and D. A. {ESL}},
  booktitle={2020 25th Asia and South Pacific Design Automation Conference (ASP-DAC)}, 
  title={RRAM-VAC: A Variability-Aware Controller for RRAM-based Memory Architectures}, 
  year={2020},
  volume={},
  number={},
  pages={181-186},
  doi={10.1109/ASP-DAC47756.2020.9045220}}
```

## Reference

* **Shikhar Tuli, Marco Antonio Rios, Alexandre Sébastien Julien Levisse, David Atienza Alonso, [RRAM-VAC: A Variability-Aware Controller for RRAM-based Memory Architectures.](https://ieeexplore.ieee.org/abstract/document/9045220) 2020 25th Asia and South Pacific Design Automation Conference (ASP-DAC), Beijing, China, 2020, pp. 181-186, doi: 10.1109/ASP-DAC47756.2020.9045220.** 
