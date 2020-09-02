# keomod

This is the module that contains the subroutine that generates the sinc discrete variable representation kinetic energy operator

## Subroutines

```f90
  !returns sysparams%T with sinc DVR KEO
  subroutine buildkeo(sysparams)
    type(systemparameters) ,intent(inout) :: sysparams
```
