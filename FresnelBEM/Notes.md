## Out of memory error
If the computers runs out of memory, it can be caused by the interpolation table for the Green functions.
The required relative accuracy of the Green function can be set by

```bash
export SCUFF_INTERPOLATION_TOLERANCE=1E-2
```

The default value is 1E-6.
