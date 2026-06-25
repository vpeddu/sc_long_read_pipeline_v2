process runIsoSeQL {
    tag "C02_isoSeQL"

    cpus 8
    memory '64 GB'

    input:
    tuple val(sample), 
          path(sqanti_gtf), 
          path(sqanti_fasta), 
          path(sqanti_classification), 
          path(tx_xref), 
          path(sqanti_genePred),
          path(tx_dups_xref)

  //   output:
  //  tuple val(sample), path("sqanti_corrected.gtf"), path("sqanti_corrected.fasta"), path("sqanti_RulesFilter_classification.txt"), path("tx_dups_xref.tsv"), path("sqanti_corrected.genePred")

    script:
    """
    source /opt/miniconda3/bin/activate sqanti3

    python ${projectDir}/original_scripts/summarize_classif.py


    # Generate sample config file
    printf "sample_name,tissue,disease,age,sex\\n${sample},unknown,control,NA,NA\\n" > sample_config.txt

    # Generate exp config file
    printf "date,RIN,platform,method,vMap,vReference,vAnnot,vLima,vCCS,vIsoseq3,vCupcake,vSQANTI,exp_name\\n" > exp_config.txt
    printf "\$(date +%Y-%m-%d),NA,PacBio,IsoSeq,minimap2,GRCh38,GENCODE_v38,lima_2.0,ccs_6.0,isoseq3_3.4,cupcake_12.0,SQANTI3_4.2,${sample}_exp\\n" >> exp_config.txt

  git clone https://github.com/christine-liu/isoSeQL

    /opt/miniconda3/envs/sqanti3/bin/python3  isoSeQL/isoSeQL_run.py \
           --classif sqanti_classif_summ.txt \
           --genePred ${sqanti_genePred} \
           --sampleConfig sample_config.txt \
           --expConfig exp_config.txt \
           --db isoforms.db
    """
}
