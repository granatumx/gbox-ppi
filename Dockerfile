FROM granatumx/gbox-py-sdk:1.0.0

RUN wget https://downloads.thebiogrid.org/Download/BioGRID/Release-Archive/BIOGRID-4.2.193/BIOGRID-ALL-4.2.193.tab3.zip -O BIOGRID.zip

RUN apt-get remove -y python3-igraph
RUN pip install pycairo
RUN pip install python-igraph==0.7.1.post6

COPY . .

RUN ./GBOXtranslateVERinYAMLS.sh
RUN tar zcvf /gbox.tgz package.yaml yamls/*.yaml
