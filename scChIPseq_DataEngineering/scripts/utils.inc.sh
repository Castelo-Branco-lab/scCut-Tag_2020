## utils.inc.sh
##
## Copyright (c) 2018 Institut Curie                               
## Author(s): Eric Nicolas Servant
## Contact: nicolas.servant@curie.fr
## This software is distributed without any guarantee under the terms of the CECILL licence.
## See the LICENCE file for details

###########################
## trap handler
###########################

function trap_error()
{
    exec >&3 2>&4
    echo -e "Error: $(basename $0) exit with status: $3" >&2
    echo -e "Command (line $1): $2" >&2
    echo -e
    echo -e "Please, look at the logs folder for more details." >&2
    exit 1
}

function trap_exit()
{
    exec >&3 2>&4
    ##Since bash-4.0 $LINENO is reset to 1 when the trap is triggered
    if [ "$?" != "0" ]; then
	echo "Error: exit status detected." >&2
    fi

    if [ -e ${ODIR}/mapping/tmp ]; then 
	echo -e "Cleaning temporary folders ..." >&2
	/bin/rm -rf ${ODIR}/mapping/tmp; 
    fi
}

exec 3>&1 4>&2
trap 'trap_error "$LINENO" "$BASH_COMMAND" "$?"' ERR
trap 'trap_exit' 0 1 2 3

set -E ## export trap to functions
set -o pipefail  ## trace ERR through pipes         
#set -o errexit   ## set -e : exit the script if any statement returns a non-true return value


###########################
## Subroutine for pipelines
###########################

die() 
{ 
    echo "Exit: $@" >&2
    exit 1
}

exec_cmd()
{
    echo $*
    if [ -z "${DRY_RUN+x}" ]; then
	eval "$@" #|| die 'Error during exec !'
    fi
}

exec_ret()
{
    if [ -z "${DRY_RUN+x}" ]; then
	eval "$@" #|| die 'Error during exec !'
    fi
}

abspath() 
{
    [[ $1 = /* ]] && echo "$1" || echo "$PWD/${1#./}"
}

read_config()
{
    eval "$(sed -e '/^$/d' -e '/^#/d' -e 's/ =/=/' -e 's/= /=/' $1 | \
awk -F"=" '{printf("%s=\"%s\"; export %s;\n", $1, $2, $1)} $1~"PATH"{printf("export PATH=%s:$PATH;\n", $2)}')"
}

is_in_path()
{
    type -P $1 > /dev/null && echo 1 || echo 0
}


#sort_big TMP 	infile outfile kstring
sort_big()
{

	
	LARGETMP=${1}
	infile=$2
	outfile=$3
	$kstring=$4

	random_string=$(head /dev/urandom | tr -dc A-Za-z0-9 | head -c 13 ; echo '')
	mkdir -p ${LARGETMP}/${random_string}
	N_LINES=1000000 # Adjust when to still too large or too small
	split --lines=${N_LINES} ${infile} ${LARGETMP}/${random_string}/splitted_
	for small in ${LARGETMP}/${random_string}/splitted_*; do
	   sort  -T ${LARGETMP} ${small} > ${LARGETMP}/${random_string}/sorted_${small}
	   rm ${small}
	done
	echo "Done with sorting the splitted files, now concate the stuff"
	sort -m -T ${LARGETMP} ${LARGETMP}/${random_string}/sorted_* > ${infile}.sorted
}

file_exists()
{

    if [ ! -e $1 ]; then
 	echo -e "Error: The file ${1} was not found. Exit" >&2
 	echo
 	exit 1
    elif [ ! -s $1 ]; then
 	echo -e "Error: The file ${1} was found but is empty. Exit" >&2
 	echo
 	exit 1
    fi
}

get_fastq_prefix()
{
    basename $1 | sed -e 's/[\._]R*[12].*.fastq\(.gz\)*//'
}
