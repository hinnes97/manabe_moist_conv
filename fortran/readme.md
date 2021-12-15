### Instructions for compiling/running:
- Run ./build to compile the source code
- Run ./main to execute code
- Run `python python_comp.py` to compare with python results

### Notes
- Fsolve (More 1980) is an external library that can be found [here](https://people.sc.fsu.edu/~jburkardt/f_src/fsolve/fsolve.html).
It requires double-precision variables, whereas FMS uses single precision, so I have done some hacky type conversions to interface the two.
I don't know if there is a better way of doing this.

