# Dockerfile to add PERL modules for Ensembl VEP

FROM chrismoth/image_phase0

RUN \
  yum -y install mysql && \
  yum -y install which && \
  yum -y install perl-App-cpanminus && \
  cpanm DBI && \
  cpanm DBD::mysql && \
  cpanm autodie && \
  cpanm parent && \
  cpanm Try && \
  yum -y install git 

CMD [ "/bin/bash" ]
