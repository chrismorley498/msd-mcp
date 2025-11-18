# MPC practice library

## Deps
- Eigen
- boostodeint
- NLOpt

## Runtime deps
Need to point your program to the nlopt shared lib
`export LD_LIBRARY_PATH=/usr/local/lib`

## Notes
You need to make sure that pinocchio is included before other boost files to set the default size of boost::variant. Specifically `pinocchio/fwd.hpp`
https://github.com/stack-of-tasks/pinocchio/issues/837