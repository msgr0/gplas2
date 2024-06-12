params.gfa = ""

process GPLAS {

  input:
  tuple val(meta), path(graph), path(pred)

  output:
  tuple val(meta), path(res_path), emit: res
  tuple val(meta), path(bin_path), emit: bins

  script:
  res_path = "results/${meta.id}_results.tab"
  bin_path = "results/${meta.id}_bins.tab"
  """
  gplas -c predict -i ${graph} -P ${pred} -n ${meta.id}
  """
}


workflow BINNING {
  take:
  prediction

  main:
  GPLAS ( prediction )

  emit:
  results = GPLAS.out.res
  bins = GPLAS.out.bins
}


process EXTRACT {
  input:
  tuple val(meta), path(graph)

  output:
  tuple val(meta), path(fasta)

  script:
  fasta = "gplas_input/${meta.id}_contigs.fasta"
  """
  gplas -c extract -i ${graph} -n ${meta.id}
  """
}
process MLPLASMIDS {
  conda 'envs/r_packages.yaml'
  
  input:
  tuple val(meta), path(fasta)

  output:
  tuple val(meta), path(pred), emit: pred

  script:
  pred = "${meta.id}.ml.pred"
  """
  #!/usr/bin/Rscript

  Rscript ${projectDir}/gplas/scripts/run_mlplasmids.R ${fasta} ${pred} 0.5 '${meta.species}'
  """
}

workflow PREDICTION {
  take:
  graph

  main:
  MLPLASMIDS ( EXTRACT ( graph ) )

  emit:
  prediction = MLPLASMIDS.out.pred
}


workflow {
  sample = Channel.fromPath("$params.gfa", type: 'file')
  input = sample.map{file -> tuple ([id: file.getName(), species: "Klebsiella Pneumoniae"], file)}
  
  PREDICTION ( input )
  BINNING ( PREDICTION.out.prediction )
  BINNING.out.results | view

}
