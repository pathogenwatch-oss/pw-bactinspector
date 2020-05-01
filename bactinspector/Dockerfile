FROM python:3

RUN pip3 install pandas

RUN curl -L https://github.com/marbl/Mash/releases/download/v2.2/mash-Linux64-v2.2.tar | tar xv \
    && mv mash-Linux64-v2.2/mash /usr/bin/ \
    && rm -rf mash-Linux64-v2.2

RUN mkdir -p /bactinspector

COPY . /bactinspector/

RUN chmod +x /bactinspector/entrypoint.sh

WORKDIR /bactinspector

ENTRYPOINT ["/bactinspector/entrypoint.sh"]
