Bootstrap:docker
From:chrismoth/image_phase9

%runscript 
	/bin/bash

%environment
    CONTAINER_TYPE=Singularity
    export CONTAINER_TYPE
    LD_LIBRARY_PATH=/psbadmin/rosetta3.7/main/source/build/src/release/linux/3.10/64/x86/gcc/4.8/default:\
/psbadmin/rosetta3.7/main/source/build/external/release/linux/3.10/64/x86/gcc/4.8/default:$LD_LIBRARY_PATH
    export LD_LIBRARY_PATH

%post -c bin/bash
    CUSTOM_PS1_SCRIPT=/.singularity.d/env/pipeline_ps1_setup.sh
    cat > $CUSTOM_PS1_SCRIPT << EOF
PS1='\r\n\[$(/usr/bin/tput setaf 195; /usr/bin/tput setab 033)\]   PSB Pipeline    host:\h  user:\u      Singularity Container\[$(/usr/bin/tput sgr0)\]\r\n\w\$ '
#export PS1
EOF
    chmod 755 $CUSTOM_PS1_SCRIPT
    mkdir /scratch /data /gpfs22 /gpfs23

