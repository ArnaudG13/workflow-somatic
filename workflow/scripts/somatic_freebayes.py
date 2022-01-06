#!/usr/bin/env python

import sys

#Filtering based on bcbio-nextgen see https://github.com/bcbio/bcbio-nextgen/blob/master/bcbio/variation/freebayes.py#L118

def _check_lods(parts, tumor_thresh, normal_thresh, indexes):
    """Ensure likelihoods for tumor and normal pass thresholds.
    Skipped if no FreeBayes GL annotations available.
    """
    try:
        gl_index = parts[8].split(":").index("GL")
    except ValueError:
        return True
    try:
        tumor_gls = [float(x) for x in parts[indexes["tumor"]].strip().split(":")[gl_index].split(",") if x != "."]
        if tumor_gls:
            tumor_lod = max(tumor_gls[i] - tumor_gls[0] for i in range(1, len(tumor_gls)))
        else:
            tumor_lod = -1.0
    # No GL information, no tumor call (so fail it)
    except IndexError:
        tumor_lod = -1.0
    try:
        normal_gls = [float(x) for x in parts[indexes["normal"]].strip().split(":")[gl_index].split(",") if x != "."]
        if normal_gls:
            normal_lod = min(normal_gls[0] - normal_gls[i] for i in range(1, len(normal_gls)))
        else:
            normal_lod = normal_thresh
    # No GL inofmration, no normal call (so pass it)
    except IndexError:
        normal_lod = normal_thresh
    return normal_lod >= normal_thresh and tumor_lod >= tumor_thresh

def _check_freqs(parts, indexes):
    """Ensure frequency of tumor to normal passes a reasonable threshold.
    Avoids calling low frequency tumors also present at low frequency in normals,
    which indicates a contamination or persistent error.
    """
    thresh_ratio = 2.7
    try:  # FreeBayes
        ao_index = parts[8].split(":").index("AO")
        ro_index = parts[8].split(":").index("RO")
    except ValueError:
        ao_index, ro_index = None, None
    try:  # VarDict
        af_index = parts[8].split(":").index("AF")
    except ValueError:
        af_index = None
    if af_index is None and ao_index is None:
        # okay to skip if a gVCF record
        if parts[4].find("<*>") == -1:
            raise NotImplementedError("Unexpected format annotations: %s" % parts[8])
    def _calc_freq(item):
        try:
            if ao_index is not None and ro_index is not None:
                ao = sum([int(x) for x in item.split(":")[ao_index].split(",")])
                ro = int(item.split(":")[ro_index])
                freq = ao / float(ao + ro)
            elif af_index is not None:
                freq = float(item.split(":")[af_index])
            else:
                freq = 0.0
        except (IndexError, ValueError, ZeroDivisionError):
            freq = 0.0
        return freq
    tumor_freq, normal_freq = _calc_freq(parts[indexes["tumor"]]), _calc_freq(parts[indexes["normal"]])
    return normal_freq <= 0.001 or normal_freq <= tumor_freq / thresh_ratio

def remove_missingalt(line):
    """Remove lines that are missing an alternative allele.
    During cleanup of extra alleles, bcftools has an issue in complicated cases
    with duplicate alleles and will end up stripping all alternative alleles.
    This removes those lines to avoid issues downstream.
    """
    if not line.startswith("#"):
        parts = line.split("\t")
        if parts[4] == ".":
            return None
    return line

def call_somatic(vcf, tumor_name, normal_name):
    """Call SOMATIC variants from tumor/normal calls, adding REJECT filters and SOMATIC flag.
    Works from stdin and writes to stdout, finding positions of tumor and normal samples.
    Uses MuTect like somatic filter based on implementation in speedseq:
    https://github.com/cc2qe/speedseq/blob/e6729aa2589eca4e3a946f398c1a2bdc15a7300d/bin/speedseq#L62
    Extracts the genotype likelihoods (GLs) from FreeBayes, which are like phred scores
    except not multiplied by 10.0 (https://en.wikipedia.org/wiki/Phred_quality_score).
    For tumors, we retrieve the best likelihood to not be reference (the first GL) and
    for normal, the best likelhood to be reference.
    After calculating the likelihoods, we compare these to thresholds to pass variants
    at tuned sensitivity/precision. Tuning done on DREAM synthetic 3 dataset evaluations.
    We also check that the frequency of the tumor exceeds the frequency of the normal by
    a threshold to avoid calls that are low frequency in both tumor and normal. This supports
    both FreeBayes and VarDict output frequencies.
    """
    # Thresholds are like phred scores, so 3.5 = phred35
    tumor_thresh, normal_thresh = 3.5, 3.5
    new_headers = ['##INFO=<ID=SOMATIC,Number=0,Type=Flag,Description="Somatic event">\n',
                   ('##FILTER=<ID=REJECT,Description="Not somatic due to normal call frequency '
                    'or phred likelihoods: tumor: %s, normal %s.">\n')
                   % (int(tumor_thresh * 10), int(normal_thresh * 10))]
    def _output_filter_line(line, indexes):
        parts = line.split("\t")
        if _check_lods(parts, tumor_thresh, normal_thresh, indexes) and _check_freqs(parts, indexes):
            parts[7] = parts[7] + ";SOMATIC"
        else:
            if parts[6] in set([".", "PASS"]):
                parts[6] = "REJECT"
            else:
                parts[6] += ";REJECT"
        line = "\t".join(parts)
        sys.stdout.write(line)
    def _write_header(header):
        for hline in header[:-1] + new_headers + [header[-1]]:
            sys.stdout.write(hline)
    header = []
    indexes = None
    for line in open(vcf, 'r'):
        if not indexes:
            if line.startswith("#"):
                header.append(line)
            else:
                parts = header[-1].rstrip().split("\t")
                indexes = {"tumor": parts.index(tumor_name), "normal": parts.index(normal_name)}
                _write_header(header)
                _output_filter_line(line, indexes)
        else:
            _output_filter_line(line, indexes)
    # no calls, only output the header
    if not indexes:
        _write_header(header)

def _clean_freebayes_output(line):
    """Clean FreeBayes output to make post-processing with GATK happy.
    XXX Not applied on recent versions which fix issues to be more compatible
    with bgzip output, but retained in case of need.
    - Remove lines from FreeBayes outputs where REF/ALT are identical:
      2       22816178        .       G       G       0.0339196
      or there are multiple duplicate alleles:
      4       60594753        .       TGAAA   T,T
    - Remove Type=Int specifications which are not valid VCF and GATK chokes
      on.
    """
    if line.startswith("#"):
        line = line.replace("Type=Int,D", "Type=Integer,D")
        return line
    else:
        parts = line.split("\t")
        alleles = [x.strip() for x in parts[4].split(",")] + [parts[3].strip()]
        if len(alleles) == len(set(alleles)):
            return line
    return None

def clean_vcf_output(orig_file, clean_fn, config, name="clean"):
    """Provide framework to clean a file in-place, with the specified clean
    function.
    """
    base, ext = utils.splitext_plus(orig_file)
    out_file = "{0}-{1}{2}".format(base, name, ext)
    if not utils.file_exists(out_file):
        with open(orig_file) as in_handle:
            with file_transaction(config, out_file) as tx_out_file:
                with open(tx_out_file, "w") as out_handle:
                    for line in in_handle:
                        update_line = clean_fn(line)
                        if update_line:
                            out_handle.write(update_line)
        move_vcf(orig_file, "{0}.orig".format(orig_file))
        move_vcf(out_file, orig_file)
        with open(out_file, "w") as out_handle:
            out_handle.write("Moved to {0}".format(orig_file))


input_vcf = sys.argv[1]
tumor = sys.argv[2]
normal = sys.argv[3]
call_somatic(input_vcf, tumor, normal)
