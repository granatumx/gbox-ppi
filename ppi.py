#!/usr/bin/env python3 -tt
"""
20180731 tkwolf

This warning message may be ignored:
/usr/lib/python3.6/importlib/_bootstrap.py:219: RuntimeWarning: numpy.dtype size changed, may indicate binary incompatibility.
"""

# Imports
# Standard library imports
import sys
#import os
import argparse as ap
# Related third party imports
import pandas as pd
import igraph as ig
# Local application/library specific imports
import granatum_sdk as gsdk
import matplotlib as mpl
import matplotlib.pyplot as plt
import cairo
from matplotlib.artist import Artist
from matplotlib.backends.backend_cairo import RendererCairo

# Global variables

# Class declarations

class GraphArtist(Artist):
    """Matplotlib artist class that draws igraph graphs.

    Only Cairo-based backends are supported.
    """

    def __init__(self, graph, bbox, palette=None, *args, **kwds):
        """Constructs a graph artist that draws the given graph within
        the given bounding box.

        `graph` must be an instance of `igraph.Graph`.
        `bbox` must either be an instance of `igraph.drawing.BoundingBox`
        or a 4-tuple (`left`, `top`, `width`, `height`). The tuple
        will be passed on to the constructor of `BoundingBox`.
        `palette` is an igraph palette that is used to transform
        numeric color IDs to RGB values. If `None`, a default grayscale
        palette is used from igraph.

        All the remaining positional and keyword arguments are passed
        on intact to `igraph.Graph.__plot__`.
        """
        Artist.__init__(self)

        if not isinstance(graph, ig.Graph):
            raise TypeError("expected igraph.Graph, got %r" % type(graph))

        self.graph = graph
        self.palette = palette or ig.palettes["rainbow"]
        self.bbox = ig.BoundingBox(bbox)
        self.args = args
        self.kwds = kwds

    def draw(self, renderer):
        if not isinstance(renderer, RendererCairo):
            raise TypeError("graph plotting is supported only on Cairo backends")
        self.graph.__plot__(renderer.gc.ctx, self.bbox, self.palette)


# Function declarations
def get_igraph_ppi():
    return 0

def main():
    args = sys.argv[1:]
    if not args:
        print('usage: [-h, --help]')
        sys.exit(1)
    parser = ap.ArgumentParser()
    parser.add_argument(
        '--test_mode',
        action='store_true',
        help='Generate test graphs instead of using the Granatum SDK.'
    )
    parser.add_argument(
        '--ppi_table',
        type=str,
        metavar='TAB_DELIMITED_TABLE',
        help='BIOGRID table with "Official Symbol Interactor A" and '+
            'Official Symbol Interactor B" named columns containing '+
            'the gene symbol pairs.'
    )
    parser.add_argument(
        '--top_scoring_genes',
        type=int,
        default=100,
        metavar='INTEGER',
        help='Number of top absolute value scoring genes to include.'
    )
    # parsed_args, unparsed = parser.parse_known_args()
    parsed_args = parser.parse_args()
    gene_pairs = pd.read_table(
        parsed_args.ppi_table, 
        low_memory=False
        # compression='bz2',
        # header=1
    )
    number_of_top_genes = parsed_args.top_scoring_genes
    print(gene_pairs.head()) # TEST - can be removed later
    # Populate species and gene-scores depending on gbox_mode
    species = None
    gene_score_dict = None
    if (parsed_args.test_mode == False):
        print("Gbox mode")
        gn = gsdk.Granatum()
        species = gn.get_arg('species')
        gene_score_dict = gn.get_import('genesAndScores')
    else:
        print("Test mode")
        # Select vertices according to 
        species = 'Hs' # Keeping Hs or Mm convention from Granatum (v1)
        gene_score_dict = dict()
        if species == 'Hs':
            gene_score_dict['BRCA1'] = 2
            gene_score_dict['PRRC2A'] = 3
            gene_score_dict['HNRNPM'] = 1
            gene_score_dict['QARS'] = -1
            gene_score_dict['MSH2'] = -0.5
            gene_score_dict['HAND2'] = 0.25
            gene_score_dict['EIF3K'] = 1
            gene_score_dict['TCF3'] = -0.15
            # if 'TCF3' in gene_score_dict:
            #     print("OK "+'TCF3 '+str(gene_score_dict['TCF3'])) # TEST
    # Subset to top absolute values (otherwise graph will get too busy)
    gene_score_top_genes = sorted(gene_score_dict, reverse=True, key=lambda key: abs(gene_score_dict[key]))[0:number_of_top_genes] 
    print("TEST:") # TEST
    print(gene_score_top_genes) # TEST
    # Set species code
    species_code = None
    if species == 'Hs':
        species_code = 9606
    elif species == 'Mm':
        species_code = 10090
    else:
        # Error
        raise ValueError('Could not understand species arguement.')
    # Get min/max of gene scores for scaling to adjust colors in graph
    score_min = min(gene_score_dict.values())
    score_max = max(gene_score_dict.values())
    score_scale = mpl.colors.Normalize(vmin=score_min, vmax=score_max)
    # Colors: https://matplotlib.org/users/colormaps.html
    #color_mapper = mpl.cm.ScalarMappable(norm=score_scale, cmap=mpl.cm.Greys) # Grey scale
    # If check for +/-, etc. may go with colors and scaled +/- the highest absolute value
    color_mapper = mpl.cm.ScalarMappable(norm=score_scale, cmap=mpl.cm.RdBu) # Red to blue
    print("min: "+str(score_min)) # TEST
    print("max: "+str(score_max)) # TEST
    # Load graph, filtering on input
    g = ig.Graph()
    # Without further filtering, multiple edges may be drawn for each
    # record, e.g., the gene-gene pair may be referenced in multiple
    # publications. Another idea could be to deduplicate entirely or
    # incrase the edge size for each additional pair. The same gene
    # may also reference its self, and currently will show an edge that
    # gets drawn back on its self, which may be reasonable to show.
    # To speed up processing, the human and mouse datasets could be separated,
    # however, to make it easy to update the PPI database, one file is used.
    #for index, row in gene_pairs.head(200).iterrows(): # TEST
    for index, row in gene_pairs.iterrows():
        gene_A = row['Official Symbol Interactor A']
        gene_B = row['Official Symbol Interactor B']
        organism_A = row['Organism Interactor A']
        organism_B = row['Organism Interactor B']
        if ((organism_A == species_code) and
          (organism_A == organism_B)):
            # if ((gene_A in gene_score_dict) and
            #   (gene_B in gene_score_dict)):
            if ((gene_A in gene_score_top_genes) and
              (gene_B in gene_score_top_genes)):
                # Keep the "label" attribute, since it is recognized by the plot function
                # Add colors according to score and min/max scaling
                score_A = gene_score_dict[gene_A]
                score_B = gene_score_dict[gene_B]
                color_A = color_mapper.to_rgba(score_A)
                color_B = color_mapper.to_rgba(score_B)
                g.add_vertex(gene_A, label=gene_A, color=color_A, label_size=14, label_angle=1.5708, label_dist=1)
                g.add_vertex(gene_B, label=gene_B, color=color_B, label_size=14, label_angle=1.5708, label_dist=1)
                g.add_edge(gene_A, gene_B)
    print(g) # TEST
    # Subset if want to load entire graph then subset
    # g_subgraph = g.induced_subgraph(vertices=gene_score_dict.keys())
    # Not subsetting
    g_subgraph = g
    print(g_subgraph) # TEST
    print(g_subgraph.degree()) # TEST
    # Use vertices with at least one edge (degree > 0)
    g_seq = g_subgraph.vs.select(_degree_gt=0)
    g_subgraph_degree_gt_0 = g_subgraph.induced_subgraph(g_seq)
    # kk_layout = g_subgraph.layout("kamada_kawai")
    # plot_of_all = ig.plot(g_subgraph, "ppi_kk_plot_deg_all.png", layout=kk_layout, bbox=(1000,1000), legend=(1,2,3))
    kk_layout = g_subgraph_degree_gt_0.layout("kamada_kawai")
    plot = ig.plot(g_subgraph_degree_gt_0, "ppi_kk_plot_deg_gt_0.png", layout=kk_layout, bbox=(1000,1000))
    context = cairo.Context(plot.surface) # Get plot surface for manipulation
    ## Add title
    #context.set_font_size(24)
    #drawer = ig.drawing.text.TextDrawer(context, "PPI network", halign=ig.drawing.text.TextDrawer.CENTER)
    #drawer.draw_at(300, 50, width=400)
    # Add legend
    # This part could use improvement: consider switching to put the igraph
    # into matplotlib, scale to absolute min/max, and include tick marks
    pattern = cairo.LinearGradient(0.0, 100.0, 0.0, 200.0) # cairo.LinearGradient(x0, y0, x1, y1)
    pattern.add_color_stop_rgba(0.0, 0, 0, 0, 0.9) 
    pattern.add_color_stop_rgba(1.0, 1, 1, 1, 0.9) 
    context.rectangle(900, 100, 20, 100)  # Rectangle(x0, y0, x1, y1)
    context.set_source(pattern)
    context.fill() 
    context = cairo.Context(plot.surface) # Get plot surface for manipulation
    context.set_font_size(20)
    drawer = ig.drawing.text.TextDrawer(context,
      score_max,
      valign=ig.drawing.text.TextDrawer.BOTTOM,
      halign=ig.drawing.text.TextDrawer.CENTER)
    drawer.draw_at(900, 95, width=20)
    drawer = ig.drawing.text.TextDrawer(context,
      score_min,
      valign=ig.drawing.text.TextDrawer.TOP,
      halign=ig.drawing.text.TextDrawer.CENTER)
    drawer.draw_at(900, 205, width=20)

    
    mpl.use("cairo")
    # Create the figure
    fig = plt.figure()

    # Create a basic plot
    axes = fig.add_subplot(111)
    graph_artist = GraphArtist(g_subgraph_degree_gt_0, (1000, 1000, 150, 150), layout=kk_layout)
    graph_artist.set_zorder(float('inf'))
    axes.artists.append(graph_artist)
    axes.axis('off')

    #plot.show()
    
    gn.add_current_figure_to_results(description="PPI-plot")
    gn.commit()

# Main body
if __name__ == '__main__':
    main()
