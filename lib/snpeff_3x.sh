#!/bin/bash

SNPEFF_HOME=${RESOURCES}/software/snpEff-3.6

function snpeff {
	java ${JAVA_OPTS} -jar ${SNPEFF_HOME}/snpEff.jar "$@"
}

function snpsift {
	java ${JAVA_OPTS} -jar ${SNPEFF_HOME}/SnpSift.jar "$@"
}
