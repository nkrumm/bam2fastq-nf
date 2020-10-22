import groovy.json.JsonOutput

Channel.fromPath("s3://ncgl-prod.sample-bucket/${params.sample}/${params.library_type}/analyses/**.bam")
        .toSortedList({a, b -> a.lastModified() <=> b.lastModified()}).flatten().last() // take the most recent bam
        .set { bam_to_fastqs_ch }
    
process bam_to_fastqs {
    label 'bamutils'
    echo true
    input: 
        file(bam) from bam_to_fastqs_ch
    output:
        path("output/**.fastq.gz") into fastq_group_ch
        file("log.txt")
    shell:
    '''
    mkdir output
    bam bam2FastQ --in !{bam} --splitRG --outBase split 2> log.txt

    READGROUPS=$(bam dumpHeader !{bam} | grep "^@RG" | cut -f2 | cut -f2 -d:)

    for RG in $READGROUPS; do
        mkdir -p output/${RG}
        cat $(ls split.${RG}*_1.fastq) | gzip > output/${RG}/1.fastq.gz
        cat $(ls split.${RG}*_2.fastq) | gzip > output/${RG}/2.fastq.gz
    done
    '''
}

fastq_group_ch.flatten().map { path -> 
        def (filename, readgroup_id, rest) = path.toString().tokenize('/').reverse() // tokenize path
        return [readgroup_id, file(path)]}
    .groupTuple()
    .map { readgroup_id, fastqs -> 
        def (fcid, lane, barcodes) = readgroup_id.tokenize(".")
        def config = [
            sample_id: params.sample,
            library_id: params.sample,
            readgroup_id: readgroup_id,
            fcid: fcid, lane: lane, barcodes: barcodes,
            library_type: params.library_type
        ]
        return tuple(readgroup_id, config, fastqs)
    }
    .view()
    .set { finalize_libraries_ch }

process finalize_libraries {
    label 'finalize'
    echo true
    input:
        tuple val(readgroup_id), val(config), path(fastqs) from finalize_libraries_ch

    script:
        def library_id = config.sample_id
        def meta_json_data = [
          "sample_id": config.sample_id,
          "barcode": config.barcodes,
          "library_id": library_id,
          "flowcell_id": config.fcid,
          "lane_id": Integer.parseInt(config.lane),
          "patient_test_id": config.sample_id,
          "project_id": "Clinical"
        ]
        def readgroup = "${config.fcid}.${config.lane}.${config.barcodes}"

        def meta_json = JsonOutput.prettyPrint(JsonOutput.toJson(meta_json_data))
   
        def library_path = "${params.out_base}/${config.sample_id}/${config.library_type}/libraries/${library_id}/${readgroup}"
        """
        echo '${meta_json}' > meta.json
        
        aws s3 cp --recursive --only-show-errors --acl bucket-owner-full-control --exclude="*" --include="*.fastq.gz" --include="*.json" . ${library_path}/
        """
}

