FROM rocker/r-ver
ENV CRAN_MIRROR http://cran.rstudio.com

RUN apt-get update && apt-get install -y gnupg2 \
      python3-pip \
      libboost-all-dev \
      libgsl0-dev \
      libeigen3-dev \
      cmake \
      wget \
      default-jre-headless

#Install gsutil
RUN apt-get update --fix-missing \
	 && apt-get install -y \
            apt-transport-https \
            curl \
            gnupg \
            lsb-release \
            openssh-client 

ENV CLOUD_SDK_VERSION=226.0.0
RUN export CLOUD_SDK_REPO="cloud-sdk-$(lsb_release -c -s)" && \
    echo "deb https://packages.cloud.google.com/apt $CLOUD_SDK_REPO main" > /etc/apt/sources.list.d/google-cloud-sdk.list && \
    curl https://packages.cloud.google.com/apt/doc/apt-key.gpg | apt-key add - && \
    apt-get update && apt-get install -y google-cloud-sdk=${CLOUD_SDK_VERSION}-0 $INSTALL_COMPONENTS && \
    gcloud config set core/disable_usage_reporting true && \
    gcloud config set component_manager/disable_update_check true && \
    gcloud config set metrics/environment github_docker_image && \
    gcloud --version
VOLUME ["/root/.config"]

#Install fuse
RUN export GCSFUSE_REPO="gcsfuse-$(lsb_release -c -s)" && \
    echo "deb http://packages.cloud.google.com/apt $GCSFUSE_REPO main" > /etc/apt/sources.list.d/gcsfuse.list && \
    curl https://packages.cloud.google.com/apt/doc/apt-key.gpg | apt-key add - && \
    apt-get update && apt-get install -y gcsfuse

RUN install2.r --repos ${CRAN_MIRROR}\
    Matrix\
    future\
    future.apply\
    rARPACK\
    Rcpp\
    RcppAnnoy\
    irlba\
    nmslibR\
    data.table

RUN pip3 install --upgrade setuptools \
    scipy\
    pybind11\
    numpy\
    nmslib      

RUN mkdir gephi scsvatools

#Get Gephi Toolkit
RUN wget https://github.com/gephi/gephi-toolkit/releases/download/v0.9.2/gephi-toolkit-0.9.2-all.jar
RUN mv gephi-toolkit-0.9.2-all.jar /gephi/ 
COPY force-atlas-3d.jar /gephi/
COPY R/*.R /scsvatools/

