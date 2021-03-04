#! /usr/bin/env bash
set -e
set -x

say_usage(){  ##todo: edit this usage function, in case you forget what this function does
    echo "Usage:"
}

readonly NUM_ARG=4   ## todo: change arg numbers here

log_it(){
    local log_type=$1
    local log_text=$2
    echo "[$(date +"%Y-%m-%d %T")] [${log_type}] ${log_text}"
}

check_args_num(){
    if [ $# -ne ${NUM_ARG} ];  ##todo:change the way of args checking
    then
        log_it "ERROR" "script take ${NUM_ARG} argment, caught $# argment(s)."
        say_usage
        exit 1
    fi
}


#### your code starts here

prepare_gff(){
    
    local gff_in=$1
    local path_list=$2
    local gff_out=$3
    
    for line in $(cat ${path_list});
    do
        grep ${line} ${gff_in} | awk '{if($3~"exon") print $0}' >> ${gff_out}
    done
}


main(){  # wrap your code here
    
    local path_list=$1
    local gff_in=$2
    local fa_in=$3
    
    local out_prefix=$4
    local gff_out=${out_prefix}.gff
    local fa_out=${out_prefix}.fa
    
    prepare_gff ${gff_in} ${path_list} ${gff_out}
    log_it "gff is ok , now extract fa"
    
}

#### start running
check_args_num $@
main $@
