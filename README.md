# CombLabel

Implemented in Python 3, no extra packages required

Install Python:
https://www.python.org/downloads/


## Contributors

1. Mikhail Yu. Myshkin <mikhail.myshkin@phystech.edu>
2. Maxim A. Dubinnyi <dumkas@yandex.ru>
3. Zakhar O. Shenkarev <zakhar-shenkarev@yandex.ru>

Group of Structural Biology of the Ion Channels, IBCh RAS, Russia,
Moscow <http://www.ibch.ru/en/structure/groups/sbic>

---
## Licence

The program is distributed under AGPL-3.0 license and is
available at https://github.com/mmjmike/CombLabel


---
## Program overview

CombLabel is an application that calculates sequence-specific
combinatorial selective
labeling (CSL) schemes [1,2] to label proteins with stable 13C and 15N
isotopes and perform fast NMR signal assignment. The application takes
protein sequence, laboratory stock of isotopically labeled amino acids
as input and outputs optimal labeling scheme and scheme statistics.
The application is a Python 3 script and is run from the command line.


---

## Package data folders

* **Stocks_and_prices** - here are collected amino acid stocks used in
 calculations in the `examples` folder.
* **examples** - collection of all examples of **CombLabel** use.
    Among the examples are all proteins described in the main paper [1].

---

## Command line commands

The only one required positional argument is the config-file
`[name].task`. The description of the file is presented in section ???.
Example:
````
>python CombLabel.py task_name.task
````
Or if you have specified `(“#!/usr/bin/python3”)`, you can just run
program like follows:
````
>CombLabel.py task_name.task
````
The program by default starts to search for solution starting with one
sample. You can specify the number of samples to start with using `–s
(--samples)` argument as follows:
````
>python CombLabel.py task_name.task –s 5
````
To run the program in check-solution mode you need to add `–c (--check)`
argument with the parameter of the file with existing solution.
````
>python CombLabel.py task_name.task -c existing_solution.txt
````
If you have several solutions in solution file (e.g.
`[job_name]_all_solutions.txt` from the output), marked by `[solution]`
expression, then you can also use `–n (--number)` argument with parameter
 that specifies the number of the solution of interest for which the
 stats will be calculated:
````
>python CombLabel.py task_name.task -c jobname_all_solutions.txt –n 27
````

---


## Input files

### Config-file

Template for config-file is presented below:
```
job_name task_name
sequence_file protein.seq
stock_file LabelStock.txt
optimize_price Yes
prices_file BEST_PRICE.csv
start_samples 1

#Please, specify the use of all 8 labeling types by '1' or '0'.
#Order as follows: "X", "N", "C", "D", "A", "S", "T", "F"
#You can see the description in ...
#Use 1 1 1 0 0 0 0 0 as default
labeling_types 1 1 1 1 0 0 0 0

#Please, specify the use of all 7 spectra by '1' or '0'.
#Order as follows: HSQC, HNCO, HNCA, HNCOCA, COf-HNCA, DQ-HNCA, HNCACO
#You can see the description in ...
#Use 1 1 0 0 0 0 0 as default
spectra_types 1 1 1 0 0 0 0
```

In the first part of the file the user needs to specify the job name,
paths to sequence and labeling stock files. If flag `“optimize_price”`
takes values `“yes”/“y”/“1”`(case insensitive), then the prices file
is required and program continues search for price-optimal solution
after finding the first solution; otherwise it stops after finding
first solution. `“start_samples”` is an option that specifies the
number of samples for the algorithm to use as a starting point
(duplicate of `–s` command line argument; in case the number of
samples is specified by both methods, the command line argument
overrides config value); if the line is missed then the algorithm
uses 1 sample as default value. It’s recommended to use
`“start_samples”` option if program is stuck with finding first
solution in particular number of samples (S); then it’s worth
starting with S+1 samples.

In the second and third parts of config-file user needs to specify the
labeling types that should be used in the calculation of the labeling
scheme. And the same goes to spectra types. The description of the
labeling types can be found in [1] and the descriptions of the
CO-filtered HNCA, DQ-HNCA spectra can be found in [3,4,5], other
(HSQC, HNCO, HNCA, HNCOCA , HNCACO) are standard pulse sequences.

---


### Protein sequence file

Protein sequence in one-letter code is read from the sequence file.
The program allows using `fasta` format and any number of spaces and
blank lines.
````
>VSD-Kv2.1
MAEKRKKLWDLLEKPNSSVAAKILAIISIMFIVLSTIALSLNTLPELQSLDEFGQSTDNPQLAHVEAVSIA
WFTMEYLLRFLSSPKKWKFFKGPLNAIDLLAILPYYVTIFLTESNKSVLQFQNVRRVVQIFRIMRILRILK
LARHSTGLQSLGFTLRRSGSHHHHHH
````

---

### Labeled amino acids stock file

Stock file contains a CSV table with the header of any number of
labeling types and corresponding presence or absence of particular
labeling types for amino acids in rows marked by “1” or “0”.

Example:
```
Res,X,N,C,D
Ala,1,1,1,1
Cys,1,0,0,0
Asp,1,1,1,0
Glu,1,1,1,0
Phe,1,1,1,0
Gly,1,1,1,1
His,1,0,0,0
Ile,1,1,1,1
Lys,1,1,1,0
Leu,1,1,1,1
Met,1,1,1,0
Asn,1,1,0,0
Pro,1,0,1,0
Gln,1,1,0,1
Arg,1,1,1,0
Ser,1,1,1,0
Thr,1,1,0,0
Val,1,1,1,0
Trp,1,1,1,1
Tyr,1,1,1,0
```

---
### Prices file

Prices file contains a table of prices corresponding to different amino
acids and their labeling types. The file is comma-separated. The default
prices are calculated for 1L of 1mM of feeding mixture for cell-free
protein expression system based on the cheapest offer from
[CIL](https://isotope.com),
[SigmaAldrich](https://www.sigmaaldrich.com/) and
[CortecNet](http://www.cortecnet.com/)
companies [2]. The commercially unavailable
labeling types are marked with “-1” value in the table. The full table
and calculation of prices is done in [Google Docs Spreadsheet](https://docs.google.com/spreadsheets/d/1pKwB0ANa-DyLOV3NbCwcikVvmkJOZyj1-1r4LDTbZLo/edit?usp=sharing).
For more up-to-date prices send quotes to the above mentioned companies
or use actual prices your lab payed for particularly labeled amino
acids.

Example:

```
Aa ,       X,          C,          A,          F,          N,          S,          T,          D
Ala, 0.04388,   29.66697,  105.13298,   67.19840,   12.47260,  246.95748,  282.88358,   55.68125
Arg, 0.04172,  526.65000,   -1      ,  126.39600,  231.72600,  588.79470,   -1      ,  153.78180
Asn, 0.05549,   -1      ,   -1      ,  375.35000,  118.90800,  792.73920,   -1      ,  132.12000
Asp, 0.01233,  179.68500,  783.72136,  139.75500,   36.60250,  461.05840,  656.98160,   99.82500
Cys, 0.03235,  216.85850,   -1      , 1023.71750,  445.83200,   -1      ,   -1      ,  768.57560
Gln, 0.06445,  194.08720,  591.90750,  182.68750,   92.80525,   -1      ,   -1      ,  175.38000
Glu, 0.00730,  186.11945, 2098.07380,  117.70400,   16.18430,  543.79248,   -1      ,   95.63450
Gly, 0.00616,    8.55798,   17.58139,   14.26330,    7.50700,   68.49612,   38.35783,   24.77310
His, 0.09170,   -1      ,   -1      ,  838.56000,  318.07800,   -1      ,   -1      ,  411.17400
Ile, 0.10429,  167.64804,   -1      ,  683.71016,   77.39620,   -1      ,   -1      ,  118.06200
Leu, 0.05890,   18.26026,  115.96312,  100.22152,   34.89388,   39.74754,  314.83200,  104.94400
Lys, 0.41591,  612.63156, 3637.22600,  155.25250,   73.06000,   -1      ,   -1      ,  155.25250
Met, 0.07244,   45.38122,   -1      ,  990.75440,   70.27791,  787.82880,   -1      ,  981.80180
Phe, 0.07681,   54.26075,  559.66372,  165.19000,   28.08230,   -1      ,   -1      ,  107.37350
Pro, 0.05377,  306.24580,   -1      ,  115.13000,   -1      ,   -1      ,   -1      ,   -1
Ser, 0.06232,  167.72364,  612.67470,  583.24950,   63.05400,  388.41264,   -1      ,   99.83550
Thr, 0.10447,  458.61200,   -1      ,  572.25248,  216.08368,   -1      ,   -1      ,   83.38400
Trp, 0.12540, 1037.48840,  776.07400, 2501.81750,  265.49900,   -1      ,   -1      , 2782.76049
Tyr, 0.10763,  142.59653,  387.74660,  344.26100,   49.82725,  543.57000,   -1      ,  163.07100
Val, 0.04897,   33.97350,  562.32000,  117.15000,   22.84425,  424.55160,   -1      ,   67.36125
```

---

## Output files

### Log-file

Log-file (`[job_name]_logfile.log`) contains information about all
current events happening during calculations. The program writes out
new solution reports and every 100000 iteration reports, showing current
state of the solution, its depth and price and number of samples.

```
Search for Kv21 solution started
Start time: 27-10-2017 08:38:51

Elapsed time:   0 days 00:00:00

Solutions found: 1
Best price: 6896.08
Iterations: 1164
Elapsed time:   0 days 00:00:00

Solutions found: 2
Best price: 6846.36
Iterations: 1172
Elapsed time:   0 days 00:00:00

...

iteration: 100000
time elapsed: 19.08
iterations/min: 314401.0
curr best price: 6673.28444
Current soution:
Job name: Kv21
Current soution price: 2675.08239
Number of sumples:5
Depth:8
L,NNNNC
S,NNNCN
I,NNDNN
F,NCNCN
A,CNNXN
V,DCDDX
Q,DXXDD
N,DDDND
R,XXXXX
P,XXXXX
H,XXXXX

Price: 2675.08239
Codes calc:1778547
```

---

### All solutions file

This file (`[job_name]_all_solutions.txt`) contains all currently found
solutions for CSL schemes and is updated each time a new solution is
found. In the rows are amino acids and in columns are the samples.

```
[solution]
% Solution number = 1
% Solution price  = 6896.08
Res, S_1, S_2, S_3, S_4, S_5
Leu,   N,   N,   N,   N,   C
Ser,   N,   N,   N,   C,   N
Ile,   N,   N,   D,   N,   N
Phe,   N,   C,   N,   N,   N
Ala,   N,   C,   N,   C,   N
Val,   C,   N,   D,   N,   D
Gln,   D,   D,   X,   X,   D
Asn,   D,   X,   X,   D,   D
Gly,   C,   N,   X,   N,   N
Lys,   N,   N,   X,   N,   N
Trp,   X,   X,   X,   D,   N
Thr,   N,   N,   X,   N,   X
Asp,   X,   X,   N,   N,   N
Glu,   X,   N,   N,   X,   N
Met,   N,   N,   N,   X,   X
Tyr,   N,   N,   X,   X,   X
Arg,   X,   X,   X,   X,   X
Pro,   X,   X,   X,   X,   X
His,   X,   X,   X,   X,   X

[solution]
% Solution number = 2
% Solution price  = 6846.36
Res, S_1, S_2, S_3, S_4, S_5
Leu,   N,   N,   N,   N,   C
Ser,   N,   N,   N,   C,   N
Ile,   N,   N,   D,   N,   N
Phe,   N,   C,   N,   N,   N
Ala,   N,   C,   N,   C,   N
Val,   C,   N,   D,   N,   D
Gln,   D,   D,   X,   X,   D
Asn,   D,   X,   X,   D,   D
Gly,   C,   N,   X,   N,   N
Lys,   N,   N,   X,   N,   N
Trp,   X,   X,   X,   D,   N
Thr,   N,   N,   X,   N,   X
Asp,   X,   X,   N,   N,   N
Glu,   X,   N,   N,   X,   N
Met,   N,   N,   N,   X,   X
Tyr,   N,   X,   X,   X,   X
Arg,   X,   X,   X,   X,   X
Pro,   X,   X,   X,   X,   X
His,   X,   X,   X,   X,   X

...
```

---

### Solution file

Solution file (`[job_name]_solution.txt`) contains all statistical
information about the best found solution and the solution itself. All
the tables and stats have their description in [1].

Example:

---

### Dictionary file

Assignment dictionary (`[job_name]_code_dictionary.txt`) is created for
the last best solution found (if program hasn’t complete computation the
dictionary for the last solution will be available). It shows possible
dipeptide pairs for each code. The code corresponds to NMR response in
the second residue in dipeptide pair. During the process of NMR
assignment you look at the superposition of all NMR spectra used for
each sample and write down the codes for spin systems in 1H-15N plane
according to `[codes]` table in the solution file. Then you look up the
code in the dictionary and if the code corresponds to a single dipeptide
fragment, then you can assign the spin system to the second residue in
the dipeptide. If the code corresponds to multiple dipeptides, then
you can mark it and then resolve the ambiguity by triple-resonance
spectra (if available).

Example of codes dictionary:

```
...

00042: L8 - W9
00100: T74 - M75
00300: I134 - M135
00300: I29 - M30
00343: Q123 - N124
00434: S115 - N116
00443: L41 - N42
00443: L95 - N96
00444: D58 - N59
00444: P15 - N16
01001: D51 - E52

...
```

---

## Recommendations for users
1) We recommend using price optimization. Even if the program doesn’t
finish calculation in reasonable time it will significantly (usually
by 30-40%) optimize price compared to the first solution found. Users
can specify their own prices for all labeling types. After the
calculation was aborted you can calculate stats for the best solution
using solution check mode of the application (see Command Line Commands). Note
that the first solution in the file marked by `[solution]` will be used
for stats calculation.
2) If the program is working for a considerable amount of time and the
`max_depth` value in standard output doesn’t change then we recommend
to start finding the solution with higher number of samples (`–s
(--samples)` command line argument or `start_samples` parameter in config-file).
3) If you decide to use NCD2 or NCD4 coding systems [2] (X, N, C, D
labeling types and HSQC+HNCO spectra without or with HNCA spectrum) we
recommend calculating solutions for both. If calculation for NCD4
system doesn’t result in solution with fewer samples then use the
solution for NCD2 and calculate stats with HNCA spectrum included so
that you get some extra redundancy in the information from the spectra.


---

## Future improvements
1.	As the program is written in Python, it lacks speed. Some of the
code can be rewritten in C++ or partially precompiled using Numba or
Cython. (Up to 100x speed)
2.	It is potentially parallelable as it has iterative tree-search core
algorithm. (up to 100x speed)
3.	Web-application can be created.
4.	Add GUI for creating input files and results analysis. 

---

## References
[1] M.Yu. Myshkin, M.A. Dubinnyi, D.S. Kulbatskii, E.N. Lyukmanova, Z.O. Shenkarev. CombLabel: rational design of optimized combinatorial labeling schemes. Application to backbone assignment of helical membrane proteins with limited stability. J Biomol NMR, 2018 (currently unpublished)

[2] M.A. Dubinnyi, M.Yu. Myshkin, Z.O. Shenkarev. Universal unambiguous combinatorial labeling schemes. J Biomol NMR, 2018 (currently unpublished)

[3] Löhr F., Reckel S., Karbyshev M., Connolly P.J., Abdul-Manan N., Bernhard F., Moore J.M., Dötsch V. (2012) Combinatorial triple-selective labeling as a tool to assist membrane protein backbone resonance assignment. J Biomol NMR 52(3):197-210

[4] Löhr F., Laguerre A., Bock C., Reckel S., Connolly P.J., Abdul-Manan N., Tumulka F., Abele R., Moore J.M., Dötsch V. (2014) Time-shared experiments for efficient assignment of triple-selectively labeled proteins. J Magn Reson 248:81-95

[5] Löhr F., Tumulka F., Bock C., Abele R., Dötsch V. (2015) An extended combinatorial 15N, 13Cα, and 13C' labeling approach to protein backbone resonance assignment. J Biomol NMR 62(3):263-79


