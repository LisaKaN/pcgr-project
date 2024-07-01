def helpMessage() {
  log.info"""
  -----------------------------------------------------------------------
                          PCGR nextflow PIPELINE
  -----------------------------------------------------------------------
  Usage:

  nextflow run main.nf

  Mandatory arguments:

    --samplesheet   [file]      .csv file with columns patient, status, sample, tissue, vcf

    --pcgr_dir      [file]      Directory where PCGR reference data bundle was downloaded and unpacked

    --output_dir    [file]      Output directory
  General Optional Arguments:

    --assembly      [str]       Either grch37 (default) or grch38. Human genome assembly build

    --assay         [str]       WES, WGS or TARGETED (default). Type of DNA sequencing assay performed for input data (VCF)
   
    --tumor_dp_min  [float]     minimum required depth for inclusion in report (default: 0)

    --tumor_af_min  [float]     minimum required AF for inclusion in report (default: 0)
    
    --tmb_algorithm [str]      {all_coding,nonsyn}, Method for calculation of TMB, all coding variants (Chalmers et al., Genome Medicine, 2017), or non-synonymous variants only, default: all_coding
    """.stripIndent()
}

if (params.help) exit 0, helpMessage()


// Check mandatory parameters
if(params.pcgr_dir == null){
  exit 1, "Please include --pcgr_dir <path>"
}

if(params.output_dir == null){
  exit 1, "Please include --output_dir <path>"
}

if(params.samplesheet == null){
  exit 1, "Please include --samplesheet <your_samplesheet>"
}


process makeDirs {
    input:
    tuple val(patient), val(status), val(sample), val(tissue), val(vcf)
    val output_dir

    output:
    tuple val(patient), val(status), val(sample), val(tissue), val(vcf)

    script:
    """
    mkdir -p ${output_dir}/${patient}/${sample}/pcgr
    """
}

process reformatVCF {
    input:
    tuple val(patient), val(status), val(sample), val(tissue), val(vcf)
    val output_dir

    output:
    tuple val(patient), val(status), val(sample), val(tissue), val(vcf)

    script:
    """
    python /Users/lisak/PCGR_project/nf-pcgr-new/reformat_vcf.py ${vcf} ${output_dir}/${patient}/${sample}/${sample}.reformat.vcf
    """
}

process pcgr {
    conda '/Users/lisak/PCGR_project/nf-pcgr-new/PCGR/conda/env/pcgr'

    input:
    tuple val(patient), val(status), val(sample), val(tissue), val(vcf)
    val pcgr_dir
    val output_dir
    val dp_tag
    val af_tag
    val dp_min
    val af_min
    val assembly
    val assay

    output:
    tuple val(patient), val(status), val(sample), val(tissue), val(vcf)

    script:
    """
    pcgr \\
    --pcgr_dir ${pcgr_dir} \\
    --output_dir ${output_dir}/${patient}/${sample}/pcgr \\
    --sample_id ${sample} \\
    --tumor_dp_tag ${dp_tag} \\
    --tumor_af_tag ${af_tag} \\
    --tumor_dp_min ${dp_min} \\
    --tumor_af_min ${af_min} \\
    --genome_assembly ${assembly} \\
    --input_vcf ${output_dir}/${patient}/${sample}/${sample}.reformat.vcf \\
    --include_trials \\
    --assay ${assay} \\
    --tumor_site ${tissue} \\
    --tumor_only \\
    --estimate_tmb \\
    --estimate_msi \\
    --estimate_signatures \\
    --exclude_dbsnp_nonsomatic \\
    --force_overwrite
    """
}

process parsePcgr {
    input:
    tuple val(patient), val(status), val(sample), val(tissue), val(vcf)
    val output_dir

    output:
    val true

    script:
    """
    Rscript /Users/lisak/PCGR_project/nf-pcgr-new/parse_pcgr.R ${output_dir}/${patient}/${sample}/pcgr ${output_dir}/${patient}/${sample} ${sample}
    """
}

process collectVariants {
    input:
    val ready
    val output_dir

    script:
    """
    python /Users/lisak/PCGR_project/nf-pcgr-new/collect_variants.py ${output_dir} ${output_dir}/all_variants.xlsx
    """
}

workflow {
    def pcgr_ch = Channel.fromPath(params.pcgr_dir)
    def output_ch = Channel.fromPath(params.output_dir)
  
    def samples_channel = Channel.fromPath(params.samplesheet)
          .splitCsv(header: true)
          .map{ row -> tuple(row.patient, row.status, row.sample, row.tissue, row.vcf) }
          

    def makeDirs_output = makeDirs(samples_channel, params.output_dir)
    
    def reformatVCF_output = reformatVCF(makeDirs_output, params.output_dir)
    
    def pcgr_output = pcgr(reformatVCF_output, params.pcgr_dir, params.output_dir, 
                           params.tumor_dp_tag, params.tumor_af_tag, 
                           params.dp_min, params.af_min, 
                           params.assembly, 
                           params.assay)
    
    def parsePcgr_output = parsePcgr(pcgr_output, params.output_dir)

    parsePcgr_output
        .collect()

    collectVariants(parsePcgr_output.collect().map { true }, params.output_dir)



}
