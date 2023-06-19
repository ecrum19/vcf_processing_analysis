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
$ source venv/bin/activate
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

# The commands below are equivalent for this challenge because the script will default to using:
# sample huF7A4DE's vcf file (--sample NB72462M.vcf) and SNP ID rs762551 (--query rs762551)
$ python3 Challenge.py
$ python3 Challenge.py -s NB72462M -q rs762551
```
Note that this script can unzip and process any *.vcf.gz file but is only designed to fetch files from https://my.pgp-hms.org/public_genetic_data/ . Please make sure to specify the file name present in FILENAME.vcf.gz within the --sample field for proper fetching and unzipping -- NOT the participant ID.
