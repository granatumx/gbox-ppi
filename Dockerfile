FROM granatumx/gbox-py-sdk:1.0.0

RUN wget https://downloads.thebiogrid.org/Download/BioGRID/Release-Archive/BIOGRID-4.1.190/BIOGRID-ALL-4.1.190.tab3.zip

RUN pip install pycairo
RUN pip install python-igraph==0.8.2

COPY . .

RUN ./GBOXtranslateVERinYAMLS.sh
RUN tar zcvf /gbox.tgz package.yaml yamls/*.yaml
