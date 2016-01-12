#!/bin/bash

# TODO docs
# Required parameters
# * PROJECT - project base directory
# * INTERVALS - target intervals
# * SAMPLES - list of prefixes, format: <sample_id>.<run_id>

### Settings
FASTQ_EXT=fastq.gz
FASTQ_CAT=zcat

# Processing
DATA=${PROJECT}
RESULTS=${DATA}/$(date +%F)-${REFERENCE}
if [[ $1 == '--results' ]]; then
	shift
	RESULTS="$1"
	shift
fi

# TODO detectCores
CPUS=24
JAVA_OPTS=-Xmx4g
CLUSTER_NODES=1

# Resources
RESOURCES=/resources
HG=${RESOURCES}/${REFERENCE}/${REFERENCE_NAME}
ADAPTORS=${RESOURCES}/qc/adaptors.fa
SCRIPTS=$(dirname $0)/../scripts

function log {
	echo '***' $@
}

function kill_all_jobs {
	# TODO is there some cleaner way?!?
	echo 'Killing all jobs on cluster'
	squeue -h | awk '{ print $1; }' | xargs scancel
}

function main {
	trap kill_all_jobs SIGINT

	if [[ $* ]]; then
		temp=/tmp/pipeline-$(date +%F-%H%M%S)-$$
		[[ -d ${temp} ]] || mkdir -p ${temp}
		[[ -d ${RESULTS} ]] || mkdir -p ${RESULTS}
		$*
		[[ -d ${temp} ]] && rm -rf ${temp}
	else
		echo "Usage: $0 [--results <results directory>] [command] [args]"
		echo
		echo "Commands:"
		for cmd in $(declare | grep '^[a-zA-Z_-]\+_usage=' | cut -d '=' -f 1 | rev | cut -d '_' -f 2- | rev); do
			echo -n " * ${cmd} "; eval echo \$${cmd}_usage
		done
	fi
}

function batch {
	name=$1
	shift
	sbatch -J${name} -c ${CPUS} -- $0 --results ${RESULTS} $@
}

function wait_for_jobs {
	while [[ $(squeue -h) ]]; do sleep 1; done
}
