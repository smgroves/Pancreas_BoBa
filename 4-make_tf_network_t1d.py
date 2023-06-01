
from pythontools import enrichr
import networkx as nx
import time


def prune(G_orig):
    G = G_orig.copy()
    n = len(G.nodes())
    nold = n + 1

    while (n != nold):
        nold = n
        for tf in list(G.nodes()):
            if G.in_degree(tf) == 0 or G.out_degree(tf)==0:
                if tf == 'INS' or tf == 'GCG':
                    continue
                else: G.remove_node(tf)
            # if G.in_degree(tf) == 0: G.remove_node(tf)

        n = len(G.nodes())
    return G


def prune_info(G_orig, prune_self_loops=True):
    G = G_orig.copy()
    for tf in list(G.nodes()):
        edges = G.adj[tf]
        for target in list(edges.keys()):
            if tf == target and prune_self_loops:
                G.remove_edge(tf, target)
                continue
            if 'db' not in edges[target]:
                G.remove_edge(tf, target)
            elif len(edges[target]['db']) < 2:
                G.remove_edge(tf, target)
    return prune(G)


def prune_to_chea(G_orig, prune_self_loops=True):
    G = G_orig.copy()
    for tf in list(G.nodes()):
        edges = G.adj[tf]
        for target in list(edges.keys()):
            if tf == target and prune_self_loops:
                G.remove_edge(tf, target)
                continue
            if 'db' in edges[target]:
                if not True in ['ChEA' in i for i in edges[target]['db']]: G.remove_edge(tf, target)
    #                if len(edges[target]['db']) < 2: G.remove_edge(tf, target)
    return prune(G)


# Pseudotime top 100
# tfs = ['ABL1', 'ARID3B', 'ASB9', 'BACH2', 'BEX1', 'BMP5', 'BTG2', 'CALR', 'CBFA2T2', 'CBX6', 'CD36', 'CERS6', 'CHD7',
#        'CNOT1', 'CRIP2', 'CSDE1', 'DSP', 'DUSP1', 'EEF1A1', 'ELAVL2', 'ENO1', 'EPAS1', 'FEV', 'FLNA', 'FOS', 'FOXA2',
#        'FOXP1', 'H1F0', 'HELZ2', 'HES6', 'HMGB3', 'HOPX', 'ID1', 'ID3', 'IER2', 'INSM1', 'IRX2', 'JARID2', 'JUN',
#        'LMO1', 'LMO2', 'LMO4', 'LMO7', 'MAF', 'MAFB', 'MBNL1', 'MET', 'MLLT11', 'MNX1', 'NDN', 'NEUROD1', 'NFKB1',
#        'NKX2-2', 'NKX6-1', 'NPM1', 'NR2F1', 'NR3C1', 'ONECUT2', 'PAM', 'PAX6', 'PBX1', 'PBX3', 'PLXNA2', 'POU2F2',
#        'PRDM16', 'PTGER3', 'RAN', 'RBP1', 'RBP4', 'RFX3', 'RFX6', 'RGS9', 'RPL7', 'RPL7A', 'RPS3', 'RUFY3', 'RYBP',
#        'SCG2', 'SIX2', 'SLC4A10', 'SMARCD1', 'SOX1', 'SOX4', 'SOX6', 'SQSTM1', 'STAT1', 'TBX3', 'TCF4', 'TCF7L2',
#        'TEAD2', 'TERF2IP', 'TGIF2', 'TLE1', 'TLE3', 'TLE4', 'TXNIP', 'VTN', 'ZDHHC9', 'ZNF503', 'ZNF703', 'INS', 'GCG']
# Pseudotime top 250
# tfs = ['ABL1', 'ABT1', 'ANKS1A', 'ANP32A', 'ARID3B', 'ARIH2', 'ASB9', 'ATF4', 'ATXN7', 'BACH2', 'BARX2', 'BCL11A',
#         'BEX1', 'BEX2', 'BMP5', 'BTBD11', 'BTBD3', 'BTG1', 'BTG2', 'CALR', 'CASZ1', 'CBFA2T2', 'CBFB', 'CBX2', 'CBX6',
#        'CCT4', 'CD36', 'CDKN1C', 'CERS6', 'CHD7', 'CHD9', 'CLOCK', 'CNOT1', 'CNOT6', 'CREB3L2', 'CRIP2', 'CSDE1',
#        'CSRNP3', 'CSRP2', 'DDB1', 'DDIT3', 'DDX5', 'DSP', 'DUSP1', 'E2F5', 'EEF1A1', 'EGR1', 'EID1', 'ELAVL2', 'ELL2',
#        'ENC1', 'ENO1', 'EPAS1', 'ESRRG', 'ETV1', 'ETV4', 'EZH2', 'FEV', 'FLNA', 'FMNL2', 'FOS', 'FOXA2', 'FOXO1',
#        'FOXP1', 'GLIS3', 'GREB1', 'H1F0', 'H1FX', 'HDAC9', 'HELZ2', 'HES6', 'HIC2', 'HMGA1', 'HMGB3', 'HMGN2', 'HMGN3',
#        'HNRNPK', 'HOPX', 'HOXB2', 'HOXC4', 'HPCA', 'HSF4', 'ID1', 'ID3', 'IER2', 'INSM1', 'IRX1', 'IRX2', 'ISL1',
#        'JARID2', 'JUN', 'JUNB', 'KDM4B', 'KDM5B', 'KHDRBS3', 'KIAA1549', 'KLF10', 'KLF11', 'KLF9', 'KLHL3', 'KLHL8',
#        'LARP1', 'LDB2', 'LGR4', 'LIN28A', 'LIN28B', 'LMO1', 'LMO2', 'LMO4', 'LMO7', 'LMX1B', 'LRRFIP1', 'LSR', 'LZTS1',
#        'MACF1', 'MAF', 'MAFA', 'MAFB', 'MAML3', 'MAPK3', 'MBNL1', 'MBNL2', 'MED24', 'MEN1', 'MET', 'MLLT11', 'MLXIPL',
#        'MNX1', 'MXD1', 'MXI1', 'MYCL', 'NAP1L3', 'NDN', 'NEUROD1', 'NFATC2', 'NFATC2IP', 'NFKB1', 'NKX2-2', 'NKX6-1',
#        'NOTCH3', 'NPM1', 'NR0B1', 'NR2F1', 'NR3C1', 'NR6A1', 'NRIP1', 'ONECUT2', 'PAM', 'PAX6', 'PBX1', 'PBX3', 'PHF21B',
#        'PLAGL1', 'PLXNA2', 'PLXNB1', 'PLXNB2', 'PLXNC1', 'POLR1D', 'POU2F2', 'PRDM16', 'PROX1', 'PTF1A', 'PTGER3',
#        'PURA', 'RAN', 'RBM15B', 'RBP1', 'RBP4', 'RBPMS', 'REXO2', 'RFX3', 'RFX6', 'RGS9', 'RNF8', 'ROCK1', 'ROR1',
#        'RPA1', 'RPL7', 'RPL7A', 'RPS3', 'RUFY3', 'RUNX1T1', 'RUVBL1', 'RYBP', 'SALL4', 'SCG2', 'SETBP1', 'SIX2', 'SKOR2',
#        'SLC30A9', 'SLC39A10', 'SLC4A10', 'SMARCD1', 'SNW1', 'SOX1', 'SOX13', 'SOX2', 'SOX4', 'SOX5', 'SOX6', 'SQSTM1',
#        'SRSF9', 'SSRP1', 'STAT1', 'TAX1BP1', 'TBL1X', 'TBX3', 'TCEA3', 'TCF3', 'TCF4', 'TCF7L1', 'TCF7L2', 'TDG', 'TEAD2',
#        'TERF2IP', 'TGIF2', 'THRB', 'TLE1', 'TLE3', 'TLE4', 'TOP2B', 'TRAC', 'TRIM13', 'TRIOBP', 'TRIP10', 'TRIP6', 'TRPS1',
#        'TSC22D1', 'TSC22D3', 'TSHZ2', 'TXNIP', 'VTN', 'WBP11', 'YBX1', 'YWHAB', 'ZBTB18', 'ZDHHC9', 'ZFP36L2', 'ZHX2',
#        'ZMYND10', 'ZNF385A', 'ZNF428', 'ZNF431', 'ZNF467', 'ZNF503', 'ZNF585A', 'ZNF697', 'ZNF703', 'ZNF827', 'ZYX']

# Filtered DiffEx using scanpy's filter_rank_gene_groups on 200 genes per cell type
# tfs = ['ASCL1', 'CBX1', 'CBX6', 'CENPF', 'CERS6', 'CRABP2', 'DEK', 'DUSP1', 'EID1', 'ELL2', 'EPAS1', 'ETS2', 'ETV1',
# 'FLNA', 'FOS', 'FOXA1', 'FOXJ1', 'GATA4', 'HES1', 'HES6', 'HIF1A', 'HMGA2', 'HMGB2', 'HMGN2', 'HOXA1', 'ID2', 'ILF3',
# 'INSM1', 'IRX2', 'ISL1', 'LCOR', 'LCORL', 'LIN28A', 'LMO4', 'MAFB', 'MAPK3', 'MCM7', 'MNX1', 'MYCN', 'NEUROD1',
# 'NEUROG3', 'NR3C1', 'NR6A1', 'ONECUT3', 'PAM', 'PCNA', 'PDLIM1', 'PDX1', 'PLXNA2', 'POLR2I', 'PTTG1', 'RBP1',
# 'RBPJ', 'RUFY3', 'RUNX1T1', 'SALL4', 'SMARCA2', 'SND1', 'SOX9', 'STAT1', 'TCF7L2', 'TEAD2', 'TERF2IP', 'TOP2B',
# 'TSHZ1', 'VTN', 'ZFP36L1', 'ZFP36L2', 'ZNF503', 'ZNF706', 'INS','GCG']

# Unfiltered DiffEx
tfs = ['ABL1', 'AES', 'AFF1', 'AFF3', 'AMOT', 'ANP32A', 'APC2', 'APEX1', 'ARNT2', 'ASCC3', 'ASCL1', 'ATOH7', 'BACH2',
      'BAZ2B', 'BEX1', 'BEX2', 'BMP7', 'BTBD11', 'BTBD3', 'BTG1', 'BTG2', 'CARM1', 'CASZ1', 'CBFA2T2', 'CBFA2T3',
      'CBX1', 'CBX2', 'CBX3', 'CBX5', 'CBX6', 'CCNA2', 'CDCA7', 'CDKN1C', 'CDX2', 'CELF3', 'CENPF', 'CENPK', 'CERS6',
      'CHD1L', 'CHD7', 'CRABP2', 'CRIP2', 'CSDE1', 'CSRP2', 'CTBP2', 'CTDSPL', 'CTNND1', 'DAB2', 'DACH1', 'DDX5', 'DEK',
      'DEPTOR', 'DISP3', 'DSP', 'DUSP1', 'DUT', 'DYNLL1', 'EEF1A1', 'EID1', 'ELF3', 'ELL2', 'ENO1', 'EPAS1', 'ESPL1',
      'ESRRG', 'ETS1', 'ETS2', 'ETV1', 'ETV4', 'EZH2', 'FLNA', 'FOS', 'FOSL2', 'FOXA1', 'FOXA2', 'FOXA3', 'FOXJ1',
      'FOXM1', 'FOXO1', 'FOXP2', 'FOXQ1', 'FUS', 'GADD45A', 'GATA4', 'GATA6', 'GLI3', 'GTF2I', 'H1FX', 'HDAC2',
      'HES1', 'HES4', 'HES6', 'HEY1', 'HHEX', 'HIF1A', 'HIPK2', 'HIPK3', 'HIST2H2BE', 'HMGA1', 'HMGA2', 'HMGB1',
      'HMGB2', 'HMGN1', 'HMGN2', 'HNRNPA3', 'HNRNPAB', 'HNRNPC', 'HNRNPD', 'HNRNPR', 'HOXA1', 'HOXA3', 'ID1', 'ID2',
      'ID3', 'ID4', 'ILF2', 'ILF3', 'INSM1', 'IRX1', 'IRX2', 'ISL1', 'JARID2', 'JUN', 'JUP', 'KDM1A', 'KDM5B',
      'KHDRBS1', 'KLF5', 'LCOR', 'LCORL', 'LDB2', 'LGR4', 'LHX1', 'LIN28A', 'LMO2', 'LMO3', 'LMO4', 'LMX1B', 'LRRFIP1',
      'LZTS1', 'MAFB', 'MAPK3', 'MAZ', 'MBNL2', 'MCM2', 'MCM3', 'MCM4', 'MCM6', 'MCM7', 'MDFI', 'MEIS1', 'MIS18BP1',
      'MLLT1', 'MLLT11', 'MLXIPL', 'MNX1', 'MXD4', 'MYBL2', 'MYCL', 'MYCN', 'MYOCD', 'MYT1', 'NEUROD1', 'NEUROG3',
      'NFATC2', 'NFATC4', 'NFE2L2', 'NFE2L3', 'NONO', 'NOTCH1', 'NOTCH2', 'NOTCH3', 'NPAS3', 'NPM1', 'NR2F1', 'NR2F2',
      'NR3C1', 'NR6A1', 'NSD2', 'NSMCE3', 'ONECUT2', 'ONECUT3', 'PAM', 'PBX1', 'PCBD1', 'PCBP4', 'PCNA', 'PDLIM1',
      'PDX1', 'PHF21B', 'PKNOX2', 'PLXNA2', 'PLXNB2', 'PLXNC1', 'PNRC1', 'POLR2I', 'POLR2L', 'POU2F2', 'POU3F4', 'PROX1',
      'PSIP1', 'PTBP1', 'PTMA', 'PTTG1', 'RAD21', 'RAD54L', 'RARA', 'RARB', 'RBFOX2', 'RBMX', 'RBP1', 'RBPJ', 'RBPMS',
      'RCOR2', 'RFX3', 'RFX6', 'RGS7', 'RGS9', 'RNF24', 'ROR1', 'RPA1', 'RPL7', 'RPL7A', 'RPS3', 'RUFY3', 'RUNX1T1',
      'SALL1', 'SALL4', 'SAP30', 'SFPQ', 'SIX2', 'SMARCA1', 'SMARCA2', 'SMARCD1', 'SMYD3', 'SND1', 'SNRPB', 'SOX11',
      'SOX4', 'SOX6', 'SOX9', 'SRSF2', 'SRSF3', 'SRSF5', 'SSBP2', 'SSBP3', 'SSBP4', 'STAT1', 'SUPT16H', 'TBX3', 'TCF3',
      'TCF7L2', 'TEAD2', 'TERF2IP', 'TFDP2', 'THRB', 'TIMELESS', 'TOP2B', 'TOX', 'TRIM9', 'TSC22D1', 'TSHZ1', 'UBE2I',
      'UHRF1', 'VTN', 'YAP1', 'YBX1', 'ZBTB18', 'ZFP36L1', 'ZFP36L2', 'ZNF428', 'ZNF503', 'ZNF608', 'ZNF704',
       'ZNF706', 'ZYX', 'INS','GCG']


G = nx.DiGraph()
prelim_G = nx.DiGraph()
with open(
        "/Users/sarahmaddox/Dropbox (Vanderbilt)/Parthenon_Pancreas_Project/BooleaBayes/Network/TF-Lit-Network/tf-lit-network.csv") as  infile:
    for line in infile:
        line = line.strip().split(',')
        prelim_G.add_edge(line[0], line[1])

for tf in tfs: G.add_node(tf)

for tf in tfs:
    enrichr.build_tf_network(G, tf, tfs)
    time.sleep(1)

for edge in prelim_G.edges():
    if edge[0] in tfs and edge[1] in tfs:
        G.add_edge(edge[0], edge[1])

outfile = open("/Users/sarahmaddox/Dropbox (Vanderbilt)/Parthenon_Pancreas_Project/BooleaBayes/Network/DiffEx"
               "-Network/unfiltered/_0_network_updated.csv", "w")
for edge in G.edges(): outfile.write("%s,%s\n" % (edge[0], edge[1]))
outfile.close()

Gp = prune(G)

outfile = open("/Users/sarahmaddox/Dropbox (Vanderbilt)/Parthenon_Pancreas_Project/BooleaBayes/Network/DiffEx"
               "-Network/unfiltered/_1_network_updated.csv", "w")
for edge in G.edges(): outfile.write("%s,%s\n" % (edge[0], edge[1]))
outfile.close()

Gpp = prune_info(Gp)
Gpp = prune_to_chea(Gp)

outfile = open("/Users/sarahmaddox/Dropbox (Vanderbilt)/Parthenon_Pancreas_Project/BooleaBayes/Network/DiffEx"
               "-Network/unfiltered/_2_network_updated.csv", "w")
for edge in Gpp.edges(): outfile.write("%s,%s\n" % (edge[0], edge[1]))
outfile.close()

