# Building a VUStruct Pipeline Container for Production Deployment
Running the VUStruct pipeline from a container is the preferred installaton and execution mechanism for production.

The distribution Singularity image is build from a Docker image which is built up in stages.  Deploying a container replaces the legacy concept of downloading
pipeline components and installing them in a local filesystem.  This container build process, described below, as well as the applications which support containerization, 
is "in flux", and you must to approach it with a sense of flexibility.

Comparing the final container's performance to previous "known good" Singularity images, is a good way to quickly sanity-test new container builds.

## Overview
You will download and build container components using a user id on a local machine.  Then, you will log in
as root to that same machine build the Docker image.  The singularity .simg image file, needed to run VUStruct in most HPC environments, 
is created from that Docker image.  The .simg, today about 3.7GB in size, can only be built by "root", as well.

I have achieved the end to end process on both a personal PC that boots to CentOS 7 and a similarly configured 
CentOS 7 booting VirtualBox Virtual Machine (VM) configured with at least 200GB of virtual hard drive storage, 4 CPU cores, and 8GB of RAM (12+ helps)   
(VM software downloadable from oracle.com)  [CentOS 7](https://www.centos.org/download/) reaches end of life November 2024.  It should not be hard to move the VUStruct pipeline to an alternate linux distribution, but that has not been implemented as yet.

Part of the build-up process includes downloading and building Rosetta 3.7 and 3.13, which requires a [Rosetta License]  (https://www.rosettacommons.org/software/license-and-download)

## Steps in Container Creation
1. Install [CentOS 7.latest](https://www.centos.org/download/) workstation/full development configuration from the full ISO, downloaded from any of the CentOS partner sites.
+ This process includes creating a root user, and a "yourlocaluserid"
+ Also , dring either the initial OS install process, or in followup "yum install"s, full g++ development and python environments for the Rosetta build to come.
+ As root, install docker and singularity following current online guidance from vendors.  Recently, the commands added some repos and ended with:

```# yum install docker-ce-20.10.9```  
and  
```# yum install singularity```


2. Boot the machine created above and login as yourlocaluserid.  Download Rosetta 3.7 and 3.13 full source to /home/yourlocaluserid/rosetta{3.7,3.13}
+ Build each set of binaries and libraries following the [Rosetta Build Documentation](https://www.rosettacommons.org/docs/latest/build_documentation/Build-Documentation) ->  ```./scons.py -j4 mode=release bin```  
+ You may wish to run ddg_monomer (3.7) and cartesian (3.13) binary files to check for any problems with missing libraries.  For example:  
```$ ~/rosetta3.13/complex_path_to/cartesian_ddg.default.linuxgccrelease --help```  
and  
```$ ~/rosetta3.7/complex_path_to/ddg_monomer.default.linuxgccrelease --help```  
should give command options and NOT error about missing libraries.

Note: It is quite reasonable to expect that Rosetta's static pre-built binaries would download easily and function well on the pipeline.  That avenue may be pursued in future.

3. Continuing as yourlocaluserid, check out the psbadmin (pipeline) and pdbmap github projects into your local machine

For example:

    1. Checkout CapraLab/psbadmin from github to ~/psbadmin.  (You may use other local directories, as you like)

```$ cd ~
$ git clone https://github.com/CapraLab/psbadmin/
$ cd psbadmin
$ git status # << Check that you are on the branch, typically master, you expect
```

    2. Checkout CapraLab/pdbmap to ~/psbadmin/pdbmap (pdbmap/ must be under the ~/psbadmin directory created in the first clone above.

```$ cd ~/psbadmin # << you may already be there.  This should be the same psbadmin directory from the above step.
$ git clone https://github.com/CapraLab/pdbmap/
$ cd pdbmap
$ git status # << Check that you are on the branch, typically master, you expect
```

3) Copy [headless chimera](https://www.cgl.ucsf.edu/chimera/download.html) file chimera-1.15-linux_x86_64_osmesa.bin from ucsf chimera legacy site to ~/psbadmin  Do NOT "run the file" to install it.
4) The rosetta builds are quite large, and scripts are provided here to transfer only needed excerpts.  Edit ONLY the ROSETTA_MAIN_SOURCE= variable in the scripts ~/psbadmin/docker/xfer_rosetta{3.7,3.13}/.bash to copy these rosetta excerpts to ~/psbadmin/rosetta{3.7,3.13}/main  Run both these scripts, and watch for any warnings over libraries not found.
5) Logout and login as root.  You are now ready to build the docker container:
+ Start the docker daemon and run the script to create the image layers:
```$ systemctl start docker
$ cd /home/yourlocalusername/psbadmin
$ source build_docker_images.bash
```

In the docker instance in which you are working, are many other things that you _might_ do.  Notably, if you get a build up of lots of dead images clogging local storage, you can do ```docker image prune -a``` to clear disk space.

7) You must push your docker image to the dockerhub repo before it can be integrated into singularity
```
docker image push chrismoth/image_phase9
```

8) Finally, to build the singularity image from the docker image:  
```source build_singularity.bash```

Now place image_phase9.simg in the UDN case root directory.  Typically,  
``` $ singularity shell ./image_phase9.simg ```  
is sufficient to start the VUstruct interactive pipeline component.  Though often, volume mapping and other parameter settings may be required for the $ inside the container to see the UDN subdirectory from whence it was started, and to create subdirectories for new cases.


## Helpful resources
[Vanderbilt ACCRE page on Singularity construction](https://www.vanderbilt.edu/accre/documentation/singularity/)  
O'Reilly's [Docker Up and Running](https://www.amazon.com/Docker-Shipping-Reliable-Containers-Production/dp/1492036730/)

