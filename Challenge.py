import os, argparse
from pathlib import Path
from rdflib import Graph, URIRef, Literal, Namespace
working_directory = Path(os.getcwd())

''' 
Python script for downloading the 'huF7A4DE' VFC file, converting VCF data to RDF format (with splitting for decentralization),
and retrieving 'rs762551 SNP' information via a local SPARQL query.

Dependencies:
    Linux OS terminal (this script throws multiple errors in windows command prompt)
    RDFLib - Python Package
'''

# Commandline functionality
run = argparse.ArgumentParser(description='Download, parse, and/or query VCF files.')
run.add_argument('-s', '--sample', nargs='?', default='NB72462M',
                 help='-d flag specifies a Personal Genome Project sample ID you would like to download ' +
                      'and/or parse to RDF')
run.add_argument('-q', '--query', nargs='?', default='rs762551',
                 help='-q flag specifies a SNP ID and will return information about it')
run.add_argument('-t', '--test', action='store_true',
                 help='runs a quick test to ensure proper instillation')
args = run.parse_args()

sample = vars(args)['sample']
query = vars(args)['query']
test = vars(args)['test']



def fix_ontology(filename):
    """Alters ontology document to reflect current memory location.

        Args:
            input file name (included in github repo)

        Returns:
            A personalized ontology file
        """
    starting_ont = open(working_directory / filename, "r").readlines()
    better_ont = open('vcf_ontology.ttl', 'w')

    first_prefix = "@prefix vcf: <{}> .".format(str(working_directory / "vcf_ontology.ttl" / "vcf" / '_')[:-1])

    better_ont.write(first_prefix + '\n')
    for line in starting_ont[1:]:
        better_ont.write(line)


def preprocess_vcf(reference_id):
    """Downloads and Unzips VCF file from online source.

    Args:
        reference_id (str): The identity of the file being investigated (sample reference ID)

    Returns:
        A downloaded, unzipped VCF file: a .vcf file that can be parsed
    """

    curr_files = os.listdir(working_directory)
    # Downloads .vcf file if not already in working directory
    if (reference_id+'.vcf' or reference_id+'.vcf.gz') not in curr_files:
        os.system(
            "wget --mirror --no-parent --no-host --cut-dirs=1 " +
            "https://531155966bc06bca5de62439c00ce64b-282.collections.ac2it.arvadosapi.com/_/" +
            "{}".format(reference_id+'.vcf.gz')
        )
        os.system("gunzip {}".format(reference_id+'.vcf.gz'))
    # unzips .vcf file if it is in current directory but not unzipped
    elif reference_id+'.vcf.gz' in curr_files:
        os.system("gunzip {}".format(reference_id+'.vcf.gz'))

    print("\n\n\tPreprocessing completed. Continuing to parsing and RDF conversion.")


def parse_to_rdf(rid):
    """Parses VCF file into RDF triples and stores them in different places based on chromosome. Works by
        iterating over each line of an input sample VCF file. Each line contains data for 1 SNP.
        Only nucleotide location (Location), Reference nucleotide (Reference), and Variant nucleotide (Variant) were
        added to RDF data stores (more semantic info could be added, but for this challenge it was unnecessary).

        Args:
            Sample VCF file ID

        Returns:
            VCF sample information in n# of RDF.ttl files (where n = number of chromosomes in VCF file)
    """
    # basic ontology defined here (just the classes relevant for the challenge - see vcf_ontology.ttl)
    prefix = Namespace(str(working_directory / "vcf_ontology.ttl" / "vcf" / "Info" / "_")[:-1])
    in_data = open(working_directory / (rid + '.vcf'), "r").readlines()     # read in vcf data
    o = 0
    while in_data[o][0:2] == "##":
        ontology_desc.append(in_data[o])      # for skipping over the ## lines of the vcf file (unimportant for current challenge)
        o += 1

    prev_chr = 'chr1'       # for dynamic programming and memory efficiency
    curr_chromosome = 'chr1'     # to split data into multiple locations
    g = Graph()

    # iterates over each line of the VCF file and converts data to RDF
    for i in in_data[o+1:]:
        curr_variant = i.split('\t')

        # because we are only working with known SNPs for this challenge, this statement excludes those without SNP IDs
        if curr_variant[2] != '.':

            # split data by chromosome for decentralization when chromosome number changes
            if curr_variant[0] != prev_chr:
                g.serialize(destination=working_directory / "RDF_Data" / "{}.ttl".format(str(curr_chromosome)))
                g = Graph()
                curr_chromosome = curr_variant[0]

            # Important part -- splits each VCF row (SNP) into RDF data and adds to a chromosome-specific graph
            c_id = URIRef(str(working_directory / "vcf_ontology.ttl" / "vcf" / "ID" / "{}".format(curr_variant[2])))
            g.add((c_id, prefix.Location, Literal(curr_variant[1])))  # ==> position
            g.add((c_id, prefix.Reference, Literal(curr_variant[3])))  # ==> reference base
            g.add((c_id, prefix.Variant, Literal(curr_variant[4])))  # ==> variant base

            # for the last SNP in the VCF file (when the chromosome number does not change)
            if i == in_data[-1]:
                g.serialize(destination=working_directory / "RDF_Data" / "{}.ttl".format(str(curr_chromosome)))

            # tracks current chromosome for data separation
            prev_chr = curr_variant[0]

    print("\n\tVCF to RDF Conversion completed. Now beginning SNP information query process.\n")


def snp_query(snp_id):
    """Searches parsed SNP RDF Triples via SPARQL for one specific SNP by its SNP ID. Process is iterative over the 
        number of chromosome .ttl files produced by the previous parse_to_rdf() method (should be 1-22, X, and potentially Y for human data).

        Args:
            Desired ID of SNP of interest

        Returns:
            SNP_Information.txt file that contains the Variant allele information for designated SNP
    """
    # prefix variables to clean-up SPARQL query
    info_prefix = "file://{}".format(str(working_directory / "vcf_ontology.ttl" / "vcf" / "Info" / "_")[:-1])
    id_prefix = "file://{}".format(str(working_directory / "vcf_ontology.ttl" / "vcf" / "ID" / "_")[:-1])

    # SPARQL Query to retrieve designated SNP variant allele nucleotide
    q = """
    PREFIX inf: <%s>
    PREFIX id: <%s>
    SELECT ?o
    WHERE {
        id:%s inf:Variant ?o .
    }""" % (info_prefix, id_prefix, snp_id)

    # perform SPARQL query on each chromosome's RDF graph until we get the result, break out of loop when result is found (SNP IDs are unique)
    done = False
    output = open(working_directory / 'SNP_Information.txt', 'w')
    # outer loop for chromosome_N.ttl files
    for kg in os.listdir(working_directory / 'RDF_Data'):
        g2 = Graph()
        g2.parse(working_directory / "RDF_Data" / kg)

        # actual SPARQL search of chromosome RDF Triples
        hits = g2.query(q)
      
        # for the case that the SNP is found (the object returned by REFLib is a bit complex to extract literal data from)
        for row in hits:
            for r in row:
                if r is not None:
                    output.write('The SNP %s is found on %s and exhibits the variant "%s" allele'
                                 % (snp_id, kg[:-4], str(r)))
                    # SNP found, now break out of parent loops
                    done = True
                    break
        if done:
            output.close()
            print('\nSuccess! The SNP information you are looking for can be found in the "SNP_Information.txt" file.\n')
            # end iterating over chromosomeN.ttl files because the SNP has been found
            break


# to personalize ontology and ensure there are no file path errors
fix_ontology('u_vcf_ontology.ttl')
if 'RDF_data' not in os.listdir(working_directory):
    os.mkdir(str(working_directory / 'RDF_data'))

# test case
if test:
    print("You selected to run the test version. If there are any errors they will be displayed below.")
    preprocess_vcf('test')
    parse_to_rdf('test')
    snp_query('rs62637819')

# Challenge case
else:
    preprocess_vcf(sample)
    parse_to_rdf(sample)
    snp_query(query)


