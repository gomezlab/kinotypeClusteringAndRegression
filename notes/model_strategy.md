### Some ideas/stuff to work on

Note: GRinf < 0 -> apoptosis/cell-death  
Note: GRinf == 0 -> cytostasis

### Model Outputs
|  | no50 | GR50 | GRmax | GRinf | hGR | GRaoc | GRr2 |
|:---- || :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| **dtype** | boolean | float | float | float | float | float | float |
| **range** | 0/1 | > 0 to inf | -1 to 1 | -1 to 1 | 0 to 5 | 0 (inclusive) to large |  -1 to 1 |
| **activation** | sigmoid | softplus | tanh | tanh | relu6 | relu | tanh |
| **loss** | binary cross- entropy | 0 if no50 & not inf else abs(1 - logfold) | mse | mse | mse | mse | mse |
| **weighting** | 5 | 1 | 1 | 1 | 0.1 | 1 | 1 |
