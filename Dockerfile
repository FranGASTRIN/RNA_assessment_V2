# Latest version of Python
FROM python

RUN apt update && pip3 install biopython

# Change directory ('cd ' like)
WORKDIR home

RUN git clone https://github.com/FranGASTRIN/RNA_assessment_V2.git

#WORKDIR RNA_assessment

#RUN python setup.py install
