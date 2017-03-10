## data

```bash

├── metadata
│   └── sci_ns_29923.tbl
├── reduced
│   └── S68N_NIRSPEC.hdf5
└── spectra
    └── S68N_2014Jun_noemis.txt
```



The three directories are:
- `metadata`, contains a file `sci_ns_29923.tbl` a csv file exported from the Keck Archive query of the target containing info all its Keck observations, mostly with NIRSPEC
- `reduced` contains the `.hdf5` file-format required by Starfish at the moment
- `spectra` contains the original reduced spectrum provided by Tom Greene.
