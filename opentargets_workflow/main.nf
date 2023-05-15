#!/usr/bin/env nextflow
nextflow.enable.dsl=2

workflow {
	main :
	filter_parquets(numlist,parquetpath)
        filter_parquets.out.filtered_files.collect().set{filtered_parquets}
        filter_index(numlist2,indexpath)
        filter_index.out.filtered_index.collect().set{filtered_indexs}
        merge_findings(filtered_parquets,filtered_indexs)
}

numlist = Channel.of(1..1200)
Channel
       	.fromPath(params.parquetlist)
        .splitText()
        .map{it -> it.trim()}
        .set{parquetpath}

process filter_parquets {
  conda 'my-env.yaml'

  input:
  val x
  val parquetfile

  output:
  path('*.tsv'), emit: filtered_files

  script:
  """
  mkdir -p ${params.workdir}/Output/${params.project}/

  python3 ${params.workdir}/filter_parquet.py ${params.parquetfolder}$parquetfile output_$x ${params.workdir$
  """

}

numlist2 = Channel.of(1..200)
Channel
       	.fromPath(params.indexlist)
        .splitText()
        .map{it -> it.trim()}
        .set{indexpath}

process filter_index {
  conda 'my-env.yaml'

  input:
  val x
  val indexfile

  output:
  path('*.tsv'), emit: filtered_index

  script:
  """
  mkdir -p ${params.workdir}/Output/${params.project}/

  python3 ${params.workdir}/filter_parquet.py ${params.indexfolder}$indexfile output_index_$x ${params.workd$
  """

}

process merge_findings {
  publishDir "${params.workdir}/Output/${params.project}"

  input:
  file(filtered_parquets)
  file(filtered_indexs)

  output:
  file('combined_output.txt')
  file('combined_index.txt')
  path('*.tsv')

  module 'lang/R/4.1.0-foss-2019b'

  script:
  """
  awk 'FNR==1 && NR!=1{next;}{print}' ${filtered_parquets} > combined_output.txt
  awk 'FNR==1 && NR!=1{next;}{print}' ${filtered_indexs} > combined_index.txt

  Rscript --vanilla ${params.workdir}/process_final_outputs.R --args combined_output.txt $params.qtlthreshol$
  """
}
