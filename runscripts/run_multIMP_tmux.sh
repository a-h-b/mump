#! /bin/bash -i

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
VARCONFIG=$DIR/VARIABLE_CONFIG

while read var val; do unset $var ; declare $var="$val" ; done < $VARCONFIG

if [ -z "$MAX_THREADS" ]; then
    MAX_THREADS=50
fi

if [ -z "$BIGMEM_CORES" ]; then
    BIGMEM_CORES=5
fi

usage() {
    echo "Usage: $0 [-u|d|c|f|l|i] [-g goTo] [-t number] [-r] [-n name] /absolute_path/to/config_file " 1>&2
    echo "       -n <name for main job>, only works with -c and -f" 1>&2
    echo "       -r if set, a report is generated (it's recommended to run -c, -f and -l with -r)" 1>&2
    echo "       -d if set, a dryrun is performed" 1>&2
    echo "       -c if set, the whole thing is submitted to the cluster" 1>&2
    echo "       -i if set, only the conda environments will be installed, if they don't exist" 1>&2
    echo "       -u if set, the working directory will be unlocked (only necessary for crash/kill recovery)" 1>&2
    echo "       -f if set, the main snakemake thread and indivdual rules are run in a tmux session" 1>&2
    echo "       -l if set, the main snakemake thread and indivdual rules are run in the current terminal session" 1>&2
    echo "       -t <max_threads> maximum number of cpus to use for all rules at a time. Defaults to $MAX_THREADS for -c, and to 1 for -l and -f. No effect on -r, -d or -u only." 1>&2
    echo "       -g <goTo> alternative target of the workflow." 1>&2
}

while getopts n:t:g:udflcrhi flag
do
    case $flag in
        i)
            INITIAL=true;;
        u)
            UNLOCK=true;;
        d)
            DRYRUN=true;;
        c)
            CLUSTER=true;;
        n)
            JNAME=$OPTARG;;
        r)
            REPORT=true;;
        f)
            FRONTEND=true;;
        l)
            LAPTOP=true;;
        t)
            THREADS=$OPTARG;;
        g)
            TARGET=$OPTARG;;
        h)
            usage
            exit;;
        *)  
            echo "Unimplemented option: -$OPTARG" >&2 
            usage
            exit 1;;
        :) 
            echo "Missing option argument for -$OPTARG" >&2 
            usage
            exit 1;;
        ?)
            usage
            exit
             ;;
    esac
done

shift $((OPTIND-1))

if [ -z "$1" ]; then
    echo "missing input"
    usage
    exit 1
else
    CONFIGFILE=$1
fi


#if the file cannot be found
if [[ !  -e "$1" ]]; then
   echo "Configfile "$1" was not found."
   echo "Provide full path."
   exit 1
fi

if [ "$SNAKEMAKE_VIA_CONDA" = true ]; then
   CONDA_START="conda activate $DIR/conda/snakemake_env"
   CONDA_END="conda deactivate"
   CONDA_END_t="conda deactivate;"
else
   CONDA_START=""
   CONDA_END=""
   CONDA_END_t=""
fi

START_TIME=`date +%s`
NAMEHASH=`echo $START_TIME| cksum | awk '{print $1}'`
if [ -z "$JNAME" ]; then
    JNAME="IMP3_${NAMEHASH}"
else
    JNAME="${JNAME}_${NAMEHASH}"
fi

if [ "$UNLOCK" = true ]; then
    echo "Unlocking working directory."
    eval $LOADING_MODULES
    eval $CONDA_START
    snakemake --cores 1 -s $DIR/multi_Snakefile --unlock --configfile $CONFIGFILE
    eval $CONDA_END
elif [ "$DRYRUN" = true ]; then
    echo "Dryrun."
    eval $LOADING_MODULES
    eval $CONDA_START
    if [ "$LOCK_SETTINGS" = true ]; then
        snakemake --cores 1 -s $DIR/multi_Snakefile --configfile $CONFIGFILE --config sessionName=$JNAME locked_normalMem=$NORMAL_MEM_EACH locked_bigMem=$BIGMEM_MEM_EACH locked_bigCores=$BIGMEM_CORES locked_bigTotal=$BIGMEM_TOTAL settingsLocked=true sessionKind=dryrun --dryrun $TARGET
    else
        snakemake --cores 1 -s $DIR/multi_Snakefile --config sessionName=$JNAME sessionKind=dryrun --configfile $CONFIGFILE --dryrun $TARGET
    fi
    eval $CONDA_END
elif [ "$INITIAL" = true ]; then
    echo "Initializing conda environments."
    eval $LOADING_MODULES
    eval $CONDA_START
    if [ "$LOCK_SETTINGS" = true ]; then
        snakemake --cores 1 -s $DIR/multi_Snakefile --config locked_normalMem=$NORMAL_MEM_EACH locked_bigMem=$BIGMEM_MEM_EACH locked_bigCores=$BIGMEM_CORES locked_bigTotal=$BIGMEM_TOTAL settingsLocked=true sessionKind=dryrun --conda-create-envs-only --use-conda --conda-prefix $DIR/conda --configfile $CONFIGFILE --local-cores 1 $TARGET
    else
        snakemake --cores 1 -s $DIR/multi_Snakefile --conda-create-envs-only --use-conda --conda-prefix $DIR/conda --configfile $CONFIGFILE --local-cores 1 $TARGET 
    fi
elif [ "$CLUSTER" = true ]; then
    if [ -z "$THREADS" ]; then
        THREADS=$MAX_THREADS
    fi
    echo "Submitting workflow to cluster - snakemake stays on the frontend, via tmux."
    tmux new -s $JNAME -d
    tmux send-keys -t $JNAME "$LOADING_MODULES >> $JNAME.stdout 2>> $JNAME.stderr" C-m
    tmux send-keys -t $JNAME "$CONDA_START >> $JNAME.stdout 2>> $JNAME.stderr" C-m
    if [ "$REPORT" = true ]; then
        if [ "$LOCK_SETTINGS" = true ]; then
            tmux send-keys -t $JNAME "snakemake --cores $THREADS -s $DIR/multi_Snakefile --keep-going --local-cores 1 --cluster-config $DIR/config/$SCHEDULER.config.yaml --cluster \"{cluster.call} {cluster.runtime}{resources.runtime} {cluster.mem_per_cpu}{resources.mem} {cluster.threads}{threads} {cluster.partition} {cluster.stdout}\" --configfile $CONFIGFILE --config sessionName=$JNAME locked_normalMem=$NORMAL_MEM_EACH locked_bigMem=$BIGMEM_MEM_EACH locked_bigCores=$BIGMEM_CORES locked_bigTotal=$BIGMEM_TOTAL settingsLocked=true sessionKind=cluster --use-conda --conda-prefix $DIR/conda $TARGET >> $JNAME.stdout 2>> $JNAME.stderr; snakemake --cores 1 -s $DIR/multi_Snakefile --configfile $CONFIGFILE --use-conda --conda-prefix $DIR/conda --report report.html $TARGET >> $JNAME.stdout 2>> $JNAME.stderr; $CONDA_END_t tmux kill-session" C-m
        else
            tmux send-keys -t $JNAME "snakemake --cores $THREADS -s $DIR/multi_Snakefile --keep-going --local-cores 1 --cluster-config $DIR/config/$SCHEDULER.config.yaml --cluster \"{cluster.call} {cluster.runtime}{resources.runtime} {cluster.mem_per_cpu}{resources.mem} {cluster.threads}{threads} {cluster.partition} {cluster.stdout}\" --configfile $CONFIGFILE --config sessionName=$JNAME sessionKind=cluster --use-conda --conda-prefix $DIR/conda $TARGET >> $JNAME.stdout 2>> $JNAME.stderr; snakemake --cores 1 -s $DIR/multi_Snakefile --configfile $CONFIGFILE --use-conda --conda-prefix $DIR/conda --report report.html $TARGET >> $JNAME.stdout 2>> $JNAME.stderr; $CONDA_END_t tmux kill-session" C-m
        fi
    else
        if [ "$LOCK_SETTINGS" = true ]; then
            tmux send-keys -t $JNAME "snakemake --cores $THREADS -s $DIR/multi_Snakefile --keep-going --local-cores 1 --cluster-config $DIR/config/$SCHEDULER.config.yaml --cluster \"{cluster.call} {cluster.runtime}{resources.runtime} {cluster.mem_per_cpu}{resources.mem} {cluster.threads}{threads} {cluster.partition} {cluster.stdout}\" --configfile $CONFIGFILE --config sessionName=$JNAME locked_normalMem=$NORMAL_MEM_EACH locked_bigMem=$BIGMEM_MEM_EACH locked_bigCores=$BIGMEM_CORES locked_bigTotal=$BIGMEM_TOTAL settingsLocked=true sessionKind=cluster --use-conda --conda-prefix $DIR/conda $TARGET >> $JNAME.stdout 2>> $JNAME.stderr; $CONDA_END_t tmux kill-session" C-m
        else
            tmux send-keys -t $JNAME "snakemake --cores $THREADS -s $DIR/multi_Snakefile --keep-going --local-cores 1 --cluster-config $DIR/config/$SCHEDULER.config.yaml --cluster \"{cluster.call} {cluster.runtime}{resources.runtime} {cluster.mem_per_cpu}{resources.mem} {cluster.threads}{threads} {cluster.partition} {cluster.stdout}\" --configfile $CONFIGFILE --config sessionName=$JNAME sessionKind=cluster --use-conda --conda-prefix $DIR/conda $TARGET >> $JNAME.stdout 2>> $JNAME.stderr; $CONDA_END_t tmux kill-session" C-m
        fi
    fi
elif [ "$FRONTEND" = true ]; then
    echo "Running workflow on frontend - don't use this setting except with small datasets and with no more than one run at a time."
    if [ -z "$THREADS" ]; then
        THREADS=1
    fi
    tmux new -s $JNAME -d
    tmux send-keys -t $JNAME "$LOADING_MODULES >> $JNAME.stdout 2>> $JNAME.stderr" C-m
    tmux send-keys -t $JNAME "$CONDA_START >> $JNAME.stdout 2>> $JNAME.stderr" C-m
    if [ "$REPORT" = true ]; then
        tmux send-keys -t $JNAME "snakemake --cores $THREADS -s $DIR/multi_Snakefile --keep-going --configfile $CONFIGFILE --config sessionName=$JNAME sessionKind=tmux --use-conda --conda-prefix $DIR/conda $TARGET >> $JNAME.stdout 2>> $JNAME.stderr; snakemake --cores 1 -s $DIR/multi_Snakefile --configfile $CONFIGFILE --use-conda --conda-prefix $DIR/conda --report report.html $TARGET >> $JNAME.stdout 2>> $JNAME.stderr; $CONDA_END_t tmux kill-session" C-m
    else
        tmux send-keys -t $JNAME "snakemake --cores $THREADS -s $DIR/multi_Snakefile --keep-going --configfile $CONFIGFILE --config sessionName=$JNAME sessionKind=tmux --use-conda --conda-prefix $DIR/conda $TARGET >> $JNAME.stdout 2>> $JNAME.stderr; $CONDA_END_t tmux kill-session" C-m
    fi
elif [ "$LAPTOP" = true ]; then
    echo "Running workflow in current session - don't use this setting except with small datasets and databases."
    JNAME=${JNAME//./_}
    if [ -z "$THREADS" ]; then
        THREADS=1
    fi
    eval $LOADING_MODULES
    eval $CONDA_START
    if [ "$REPORT" = true ]; then
        snakemake --cores $THREADS -s $DIR/multi_Snakefile --keep-going --configfile $CONFIGFILE --config sessionName=$JNAME sessionKind=laptop --use-conda --conda-prefix $DIR/conda $TARGET 
        snakemake --cores $THREADS -s $DIR/multi_Snakefile --configfile $CONFIGFILE --use-conda --conda-prefix $DIR/conda --report report.html $TARGET 
        eval $CONDA_END
    else
        snakemake --cores $THREADS -s $DIR/multi_Snakefile --keep-going --configfile $CONFIGFILE --config sessionName=$JNAME sessionKind=laptop --use-conda --conda-prefix $DIR/conda $TARGET 
        eval $CONDA_END
    fi    
elif [ "$REPORT" = true ]; then
    echo "Writing report."
    eval $LOADING_MODULES
    eval $CONDA_START
    snakemake --cores 1 -s $DIR/multi_Snakefile --report report.html --configfile $CONFIGFILE --use-conda --conda-prefix $DIR/conda $TARGET
    eval $CONDA_END
else
    echo "Nothing was done, please give -u, -d, -r, -c, -f, -i, or -l to start anything."
fi


