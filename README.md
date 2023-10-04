# MassRunning.jl
This package solves the RG equations to perform the running down to the $\overline{\mathrm{MS}}$ scheme.
### Installation
The package is not registered in the general Julia registry. `MassRunning.jl` depends on the `ADerrors.jl` package https://igit.ific.uv.es/alramos/aderrors.jl which should be installed beforehand. Having fulfilled this requirement, `MassRunning.jl` can be installed with the Julia package manager Pkg:
```
julia > ] # open package manager
(@v1.7) pkg> add https://github.com/Bunniies/MassRunning.jl
```

### Features
- running from SF at scale $\mu_{had}$ to $\overline{\mathrm{MS}}$ scheme supported.
- running from RGI mass to $\overline{\mathrm{MS}}$ supported.
- Covariance between $\Lambda$ and $\frac{\bar m_{SF}(\mu_0/2)}{\bar m_{SF}(\mu_{had})}$ properly taken into account in the error propagation.
- Decoupling between $N_f=3$ and $N_f=4$ flavours supported.

### Usage and theory
- In the  folder example  there is a file showing how to perform the running for the charm quark mass in the SF scheme down to the $\overline{\mathrm{MS}}$  scheme.
- In the folder notes there is a theoretical explanation of the running by Carlos Pena on which the code is based.


### Acknowledgments 
- We thanks Gregorio Herdo\'{i}za for debugging the code and implementing the decoupling between $N_f=3$ and N_f=4$ flavours.




