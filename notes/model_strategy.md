### Some ideas/stuff to work on

Note: GRinf < 0 -> apoptosis/cell-death  
Note: GRinf == 0 -> cytostasis

### Model Outputs
|  | no50 (3) | GR50 | GRmax | GRinf | hGR | GRaoc | GRr2 |
|:---- || :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| **dtype** | ternary | float | float | float | float | float | float |
| **range** | -1/0/1 | -inf & > 0 to inf | -1 to 1 | -1 to 1 | 0 to 5 | < 0 to about 2 |  -1 to 1 |
| **activation** | softmax (over 3) | softplus | tanh | tanh | relu6 (5) | celu | tanh |
| **loss** | cross entropy | 0 if no50 & not inf else abs(1 - logfold) | mse | mse | mse | mse | mse |
| **weighting** | 5 | 1 | 1 | 1 | 0.1 | 1 | 1 |
