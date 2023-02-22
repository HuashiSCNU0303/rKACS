Here is the source code of paper *Top-r keyword-based community search in attributed graphs (ICDE'23)*

### Example of compiling and running (Linux)

Before running, you should:
+ Prepare the dataset and query files according to the templates
+ Check file names and file paths with the macros defined in ./src/Def.h
+ Install required libraries (SNAP, phmap)
+ Specify the paths in the provided Makefile

Then run the following commands,

```shell
$ make
$ ./rKACS_Query Test t ckr 5 30
# For running queries, five parameters should be specified:
# "Test" is the name of dataset
# "t" means that the query time is measured
# "ckr" means that we are running queries on varying |Q|, k, r. Similarly, "ck" means that the queries on varying |Q|, k are conducted
# "5" is the default k value
# "30" is the default r value
```

You can see the result in the corresponding output files (set in ./src/Def.h. You need to create the corresponding directory in advance).
