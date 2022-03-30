# Latest version of Python
FROM python

RUN apt update -y && apt install vim -y && pip3 install biopython

# Change directory ('cd ' like)
WORKDIR home

RUN git clone https://github.com/FranGASTRIN/RNA_assessment_V2.git

VOLUME /home/RNA_assessment_V2/data_DI

WORKDIR RNA_assessment_V2

RUN python setup.py install
