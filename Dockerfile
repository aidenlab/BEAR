FROM debian:latest as builder

RUN apt-get update && apt-get install -y \
    openssh-client \
    git \
    && rm -rf /var/lib/apt/lists/*
COPY github_key .
RUN eval $(ssh-agent) && \
    ssh-add github_key && \
    ssh-keyscan -H github.com >> /etc/ssh/ssh_known_hosts && \
    git clone git@github.com:aidenlab/Polar.git /opt/Polar

FROM conda/miniconda3:latest    
LABEL maintainer="weisz@bcm.edu"

RUN useradd -u 999 --user-group --system --create-home --no-log-init --shell /bin/bash polar
COPY --from=builder /opt/Polar /opt/Polar
RUN conda config --set always_yes yes --set changeps1 no && conda update -q conda && conda init bash
RUN conda env create -n Polar_cond_env -f /opt/Polar/Polar_conda_env.yml
RUN echo "source /usr/local/etc/profile.d/conda.sh\nPATH=/opt/Polar:\$PATH\nconda activate Polar_cond_env" >>/etc/profile.d/polar.sh
RUN conda clean --all -f -y


ENV BASH_ENV "/etc/profile.d/polar.sh"

SHELL ["/bin/bash", "--login", "-c"]

ENV PATH=/opt/Polar:$PATH
WORKDIR /fastq
ENTRYPOINT ["align_serial.sh"]