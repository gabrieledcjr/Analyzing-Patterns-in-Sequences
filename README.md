Analyzing-Patterns-in-Sequences
===============================

##Software Requirements
* Tested only on Python 2.7
* BioPython --- `pip install biopython`
* xhtml2dpf python module --- `pip install xhtml2pdf`
* MUST create a `Results` folder at the same level as the `Library` folder. All output of the program are saved in this folder.

##How to use program for a single sample file
```
python seq_identifier.py -p \<patterns.fa\> -s \<sample.fa\>
python seq_identifier.py --min-len=1 --max-gap=3 -p \<patternFile\> -s \<sampleFile\>
python seq_identifier.py [Other Options] [Required Arguments]
```

##How to use program for multiple sample files (Tested on Ubuntu platform only)
Use the `run` bash script to get results for multiple files. Open this file using your favorite text editor (e.g. gedit, vim, etc). In this file, replace `line 3` with the location of the folder relative to the `seq_identifier.py` path. Currently it is set to find sample files inside the `Samples` folder at the same level as the `Library` folder. Also specify the file extension whether the samples uses `*.fa` or `*.fasta`.

At `line 6`, you can change the values of the arguments or add more arguments as required by your experiment.

Once you are satisfied with the setup, you can run this command in the terminal by typing `./run`.

##Required Arguments
| Argument           | Description            |
| ------------------ | ---------------------- |
| -p, --pattern=FILE | Pattern file location  |
| -s, --sample=FILE  | Sample file location"  |

##Optional Arguments
| Argument           | Description                                  |
| ------------------ | -------------------------------------------- |
| -h, --help         | Help                                         |
| -d, --debug        | Turn on debug mode                           |
| --min-len=NUM      | Minimum length for pattern match [Default 4] |
| --min-gap=NUM      | Minimum gap [Default 1]                      |
| --max-gap=NUM      | Maximum gap [Default 4]                      |
| --out-pdf=[0 or 1] | Output pdf file [Default 1 (true)]           |
| --st-anchor=STRING | Starting anchor [Default KWG]                |
| --en-anchor=STRING | Ending anchor [Default GMA]                  |
