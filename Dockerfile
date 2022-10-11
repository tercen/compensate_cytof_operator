FROM tercen/runtime-flowsuite:3.15-2

RUN R -e "BiocManager::install('CATALYST')"
RUN R -e "remotes::install_github('tercen/tim', ref = '0.0.20')"

COPY . /operator
WORKDIR /operator

ENV TERCEN_SERVICE_URI https://tercen.com

ENTRYPOINT ["R", "--no-save", "--no-restore", "--no-environ", "--slave", "-f", "main.R", "--args"]
CMD ["--taskId", "someid", "--serviceUri", "https://tercen.com", "--token", "sometoken"]