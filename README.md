# vcf_processing_analysis
VCF file download, parsing, and SNP analysis. Works on OS X and Linux operating systems.


## Installation
Easiest way to install is via command line:

```Bash
$ cd /path/to/desired/directory
$ git clone https://github.com/ecrum19/vcf_processing_analysis.git
$ cd vcf_processing_analysis/
```

Create a Virtual Environment:
```Bash
$ python3 -m venv venv
$ python3 source ~/path/to/curr/dir/venv/bin/activate
```

Install dependencies:
```Bash
$ pip install rdflib
```

Run a test to ensure dependencies are installed:
```Bash
$ python3 Challenge.py -t
```

## Usage
```Bash
$ python3 Challenge.py --sample sample_ID --query SNP_ID 
```
