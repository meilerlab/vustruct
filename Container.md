Running the pipeline from a container is the preferred mechanism for production applications.

The distribution Docker image is built up in stages, which replace the legacy concept of downloading
pipeline components and installing them in a local filesystem.

Importantly, you can only build the Docker image on a machine where you have root access.
The Singularity image, needed to run in most HPC environments, is created from the final Docker
image.

Begin the process by checking out the psbadmin (pipeline) and pdbmap github projects into your local machine, using
you usual non-root user id.

For example:

1) Checkout CapraLab/psbadmin to ~/psbadmin.  (You may use other local directories, as you like)

```$ cd ~
$ git clone https://github.com/CapraLab/psbadmin/
$ cd psbadmin
$ git status # << Check that you are on the branch, typically master, you expect
```

2) Checkout CapraLab/pdbmap to ~/psbadmin/pdbmap (pdbmap/ must be under the directory target of 1) above)

```$ cd ~/psbadmin # << you may already be there.  This should be the same psbadmin directory from the above step.
$ git clone https://github.com/CapraLab/pdbmap/
$ cd pdbmap
$ git status # << Check that you are on the branch, typically master, you expect
```

3) You must download and install Rosetta versions 3.7 and 3.13 to support ddg_monomer and ddg_cartesian calculations, respectively 


