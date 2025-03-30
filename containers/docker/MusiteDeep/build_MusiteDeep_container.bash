# Build chrismoth/musitedeep in one shot
docker build --progress=plain --rm -f Dockerfile.0.centos_miniconda -t chrismoth/musitedeep . 2>&1 | tee image_musitedeep.stderr_out
