#!/bin/bash

echo "Check that the granatum path is set correctly:"
echo $GRANATUM_SWD
echo ""
echo "I am currently in directory..."
res=`pwd`
echo "$res"
echo ""
echo "Running zgsea"
python3 ./ppi.py --top_scoring_genes=100 --ppi_table ./BIOGRID-ALL-4.1.190.tab3.zip
