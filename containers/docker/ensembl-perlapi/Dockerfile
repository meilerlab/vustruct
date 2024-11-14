FROM rockylinux:9

# ENSEMBL has created a container for VEP, but strangely, not for the general PERL API
# I have create this one from non-container install instructions at
#  https://useast.ensembl.org/info/docs/api/api_installation.html

ENV PATH="/opt/vustruct_perl:/opt/ensembl_perlapi/src:$PATH"
ENV PERL5LIB=/opt/ensembl_perlapi/src/bioperl-1.6.924
ENV PERL5LIB=${PERL5LIB}:/opt/ensembl_perlapi/src/ensembl/modules
ENV PERL5LIB=${PERL5LIB}:/opt/ensembl_perlapi/src/ensembl-compara/modules
ENV PERL5LIB=${PERL5LIB}:/opt/ensembl_perlapi/src/ensembl-variation/modules
ENV PERL5LIB=${PERL5LIB}:/opt/ensembl_perlapi/src/ensembl-funcgen/modules

COPY ./vustruct_perl /opt/vustruct_perl

# I am installing Rocky9 intrinsic perl
# and the database interfaces may go better when
# this is done a dnf level and not via cpanminus

RUN echo 'Building perl5 LIB' && \
     chmod a+x /opt/vustruct_perl/* && \
     dnf -y update && \
     dnf -y install zip && \
     dnf -y install wget && \
     dnf -y install perl && \
     dnf -y install mysql && \
     dnf -y install perl-DBI && \
     dnf -y install perl-DBD-mysql && \
     dnf -y install perl-autodie && \
     dnf -y install which && \
     dnf -y install perl-App-cpanminus && \
     cpanm parent && \
     cpanm Try && \
     dnf -y install git && \
     mkdir /opt/ensembl_perlapi && \
     cd /opt/vustruct_perl && \
     git clone https://github.com/Ensembl/ensembl-git-tools.git && \
     PATH=$PWD/ensembl-git-tools/bin:$PATH && \
     git ensembl --clone api && \
     git ensembl --checkout --branch release/111 api && \
     wget https://github.com/bioperl/bioperl-live/archive/release-1-6-924.zip && \
     unzip release-1-6-924.zip && \
     rm release-1-6-924.zip && \
     mv bioperl-live-release-1-6-924 bioperl-1.6.924
