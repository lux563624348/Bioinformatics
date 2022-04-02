#!/bin/sh

#
# Downloads sequence for the GRCh38 release 84 version of H. sapiens (human) from
# Ensembl.
#
# Note that Ensembl's GRCh38 build has three categories of compressed fasta
# files:
#
# The base files, named ??.fa.gz
#
# By default, this script builds and index for just the base files,
# since alignments to those sequences are the most useful.  To change
# which categories are built by this script, edit the CHRS_TO_INDEX
# variable below.
#

ENSEMBL_RELEASE=84
ENSEMBL_GRCh38_BASE=ftp://ftp.ensembl.org/pub/release-${ENSEMBL_RELEASE}/fasta/homo_sapiens/dna

get() {
	file=$1
	if ! wget --version >/dev/null 2>/dev/null ; then
		if ! curl --version >/dev/null 2>/dev/null ; then
			echo "Please install wget or curl somewhere in your PATH"
			exit 1
		fi
		curl -o `basename $1` $1
		return $?
	else
		wget $1
		return $?
	fi
}

HISAT2_BUILD_EXE=./hisat2-build
if [ ! -x "$HISAT2_BUILD_EXE" ] ; then
	if ! which hisat2-build ; then
		echo "Could not find hisat2-build in current directory or in PATH"
		exit 1
	else
		HISAT2_BUILD_EXE=`which hisat2-build`
	fi
fi

HISAT2_REPEAT_EXE=./hisat2-repeat
if [ ! -x "$HISAT2_REPEAT_EXE" ]; then
	if ! which hisat2-repeat ; then
		echo "Could not find hisat2-repeat in current directory or in PATH"
		exit 1
	else
		HISAT2_REPEAT_EXE=`which hisat2-repeat`
	fi
fi

rm -f genome.fa
F=Homo_sapiens.GRCh38.dna.primary_assembly.fa
if [ ! -f $F ] ; then
	get ${ENSEMBL_GRCh38_BASE}/$F.gz || (echo "Error getting $F" && exit 1)
	gunzip $F.gz || (echo "Error unzipping $F" && exit 1)
	mv $F genome.fa
fi

# Build repeat
CMD="${HISAT2_REPEAT_EXE} -p 4 --repeat-count 5 --repeat-length 51-300,76-300,100-300,101-300,151-300 genome.fa genome_rep"
echo Running $CMD
$CMD 
[[ "$?" -eq 0 ]] || { echo "Index building failed; see error message";  exit 1; };

CMD="${HISAT2_BUILD_EXE} -p 4 --repeat-ref genome_rep.rep.fa --repeat-info genome_rep.rep.info genome.fa genome_rep"
echo Running $CMD
if $CMD ; then
	echo "genome index built; you may remove fasta files"
else
	echo "Index building failed; see error message"
fi
