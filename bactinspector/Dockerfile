FROM python:3

RUN pip3 install pandas pyarrow profilehooks

RUN curl -L https://github.com/marbl/Mash/releases/download/v2.2/mash-Linux64-v2.2.tar | tar xv \
    && mv mash-Linux64-v2.2/mash /usr/bin/ \
    && rm -rf mash-Linux64-v2.2

RUN mkdir -p /bactinspector

COPY data /bactinspector/data

COPY libraries /bactinspector/libraries

COPY *.py /bactinspector/

COPY entrypoint.sh /bactinspector/

RUN chmod +x /bactinspector/entrypoint.sh

WORKDIR /bactinspector

ENTRYPOINT ["/bactinspector/entrypoint.sh"]
