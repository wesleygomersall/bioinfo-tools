# Bioinfo_module

As part of your Python training, you will create and submit your own module containing useful Bioinformatics functions. You should continually update this module as you progress throughout the program.

Find a starter template [here](bioinfo.py). You are not required to use it, but may find it helpful.

Upload your final `bioinfo.py` module to this repo.

## Your module should include AT LEAST the following:

- Functions:
  - `convert_phred()`
  - `qual_score()`
  - `validate_base_seq()`
  - `gc_content()`
  - `calc_median()`
  - `oneline_fasta()` (or similar, from PS7)
- Constants:
  - `DNA_bases`
  - `RNA_bases`

## Remember to include your unit tests:
```
if __name__ == "__main__":
  tests to run
```

## Also be sure that you can successfully: 
1. execute `bioinfo.py` from the command line
   ```
   $ ./bioinfo.py
   ```
3. import and use your module in a separate script
   ```python
   #!/usr/bin/env python
   import bioinfo

   print(bioinfo.validate_base_seq("AATGCGA"))
   ```
Use
   ```python
   #!/usr/bin/env python
   
    import sys
    sys.path.append("/home/wesley/PythonModules/bioinfo-module-wesleygomersall")
    import bioinfo

   print(bioinfo.validate_base_seq("AATGCGA"))
   ```



