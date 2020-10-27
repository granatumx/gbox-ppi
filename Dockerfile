FROM granatumx/gbox-py-sdk:1.0.0

RUN wget https://downloads.thebiogrid.org/File/BioGRID/Release-Archive/BIOGRID-4.1.190/BIOGRID-ALL-4.1.190.tab3.zip

COPY . .

RUN ./GBOXtranslateVERinYAMLS.sh
RUN tar zcvf /gbox.tgz package.yaml yamls/*.yaml
