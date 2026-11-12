#!/usr/bin/env python3
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg') # Headless safe
import matplotlib.patches as patches
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import matplotlib.ticker as ticker
import argparse
import os
import sys
import json
import gzip

# Hardcoded offline gene structures
GENES_STRUCTURE_38 = {
  "MLH1": {
    "gene_start": 36992834,
    "gene_end": 37051505,
    "strand": 1,
    "chrom": "3",
    "transcript_id": "ENST00000231790",
    "exons": [
      {
        "start": 36993518,
        "end": 36993663
      },
      {
        "start": 36996619,
        "end": 36996709
      },
      {
        "start": 37000955,
        "end": 37001053
      },
      {
        "start": 37004401,
        "end": 37004474
      },
      {
        "start": 37006991,
        "end": 37007063
      },
      {
        "start": 37008814,
        "end": 37008905
      },
      {
        "start": 37011820,
        "end": 37011862
      },
      {
        "start": 37012011,
        "end": 37012099
      },
      {
        "start": 37014432,
        "end": 37014544
      },
      {
        "start": 37017506,
        "end": 37017599
      },
      {
        "start": 37020310,
        "end": 37020463
      },
      {
        "start": 37025637,
        "end": 37026007
      },
      {
        "start": 37028784,
        "end": 37028932
      },
      {
        "start": 37040186,
        "end": 37040294
      },
      {
        "start": 37042268,
        "end": 37042331
      },
      {
        "start": 37047519,
        "end": 37047683
      },
      {
        "start": 37048517,
        "end": 37048609
      },
      {
        "start": 37048904,
        "end": 37049017
      },
      {
        "start": 37050486,
        "end": 37050846
      }
    ]
  },
  "MSH2": {
    "gene_start": 47403067,
    "gene_end": 47663146,
    "strand": 1,
    "chrom": "2",
    "transcript_id": "ENST00000233146",
    "exons": [
      {
        "start": 47403156,
        "end": 47403402
      },
      {
        "start": 47408401,
        "end": 47408555
      },
      {
        "start": 47410094,
        "end": 47410372
      },
      {
        "start": 47412414,
        "end": 47412560
      },
      {
        "start": 47414269,
        "end": 47414418
      },
      {
        "start": 47416296,
        "end": 47416429
      },
      {
        "start": 47429742,
        "end": 47429941
      },
      {
        "start": 47445548,
        "end": 47445657
      },
      {
        "start": 47463031,
        "end": 47463154
      },
      {
        "start": 47466658,
        "end": 47466808
      },
      {
        "start": 47470965,
        "end": 47471062
      },
      {
        "start": 47475025,
        "end": 47475270
      },
      {
        "start": 47476367,
        "end": 47476571
      },
      {
        "start": 47478272,
        "end": 47478519
      },
      {
        "start": 47480696,
        "end": 47480871
      },
      {
        "start": 47482779,
        "end": 47483223
      }
    ]
  },
  "MSH6": {
    "gene_start": 47695530,
    "gene_end": 47810109,
    "strand": 1,
    "chrom": "2",
    "transcript_id": "ENST00000234420",
    "exons": [
      {
        "start": 47783145,
        "end": 47783493
      },
      {
        "start": 47790927,
        "end": 47791123
      },
      {
        "start": 47795894,
        "end": 47796063
      },
      {
        "start": 47798611,
        "end": 47801155
      },
      {
        "start": 47803420,
        "end": 47803685
      },
      {
        "start": 47804910,
        "end": 47805027
      },
      {
        "start": 47805618,
        "end": 47805707
      },
      {
        "start": 47806204,
        "end": 47806358
      },
      {
        "start": 47806452,
        "end": 47806651
      },
      {
        "start": 47806779,
        "end": 47806953
      }
    ]
  },
  "PMS2": {
    "gene_start": 5970925,
    "gene_end": 6009130,
    "strand": -1,
    "chrom": "7",
    "transcript_id": "ENST00000265849",
    "exons": [
      {
        "start": 5970925,
        "end": 5973542
      },
      {
        "start": 5977588,
        "end": 5977757
      },
      {
        "start": 5978596,
        "end": 5978696
      },
      {
        "start": 5982824,
        "end": 5982991
      },
      {
        "start": 5986759,
        "end": 5987620
      },
      {
        "start": 5989800,
        "end": 5989955
      },
      {
        "start": 5991973,
        "end": 5992057
      },
      {
        "start": 5995534,
        "end": 5995633
      },
      {
        "start": 5997326,
        "end": 5997423
      },
      {
        "start": 5999108,
        "end": 5999275
      },
      {
        "start": 6002453,
        "end": 6002636
      },
      {
        "start": 6003690,
        "end": 6003792
      },
      {
        "start": 6003972,
        "end": 6004058
      },
      {
        "start": 6005892,
        "end": 6006031
      },
      {
        "start": 6008997,
        "end": 6009049
      }
    ]
  },
  "EPCAM": {
    "gene_start": 47345158,
    "gene_end": 47387601,
    "strand": 1,
    "chrom": "2",
    "transcript_id": "ENST00000263735",
    "exons": [
      {
        "start": 47369311,
        "end": 47369581
      },
      {
        "start": 47373463,
        "end": 47373570
      },
      {
        "start": 47373808,
        "end": 47374048
      },
      {
        "start": 47375234,
        "end": 47375299
      },
      {
        "start": 47377014,
        "end": 47377077
      },
      {
        "start": 47378953,
        "end": 47379054
      },
      {
        "start": 47379769,
        "end": 47379969
      },
      {
        "start": 47385166,
        "end": 47385210
      },
      {
        "start": 47386572,
        "end": 47387020
      }
    ]
  },
  "MLH3": {
    "gene_start": 75013224,
    "gene_end": 75051532,
    "strand": -1,
    "chrom": "14",
    "transcript_id": "ENST00000355774",
    "exons": [
      {
        "start": 75013775,
        "end": 75017201
      },
      {
        "start": 75018829,
        "end": 75018980
      },
      {
        "start": 75022814,
        "end": 75022892
      },
      {
        "start": 75022995,
        "end": 75023018
      },
      {
        "start": 75030543,
        "end": 75030702
      },
      {
        "start": 75032068,
        "end": 75032179
      },
      {
        "start": 75033419,
        "end": 75033490
      },
      {
        "start": 75038340,
        "end": 75038412
      },
      {
        "start": 75039911,
        "end": 75040015
      },
      {
        "start": 75041615,
        "end": 75041700
      },
      {
        "start": 75042379,
        "end": 75042477
      },
      {
        "start": 75046376,
        "end": 75049718
      },
      {
        "start": 75051380,
        "end": 75051467
      }
    ]
  },
  "POLD1": {
    "gene_start": 50384204,
    "gene_end": 50418018,
    "strand": 1,
    "chrom": "19",
    "transcript_id": "ENST00000440232",
    "exons": [
      {
        "start": 50384347,
        "end": 50384390
      },
      {
        "start": 50398851,
        "end": 50399053
      },
      {
        "start": 50399371,
        "end": 50399484
      },
      {
        "start": 50401778,
        "end": 50401924
      },
      {
        "start": 50401999,
        "end": 50402124
      },
      {
        "start": 50402205,
        "end": 50402373
      },
      {
        "start": 50402454,
        "end": 50402535
      },
      {
        "start": 50402612,
        "end": 50402741
      },
      {
        "start": 50403053,
        "end": 50403219
      },
      {
        "start": 50403493,
        "end": 50403597
      },
      {
        "start": 50406182,
        "end": 50406322
      },
      {
        "start": 50406407,
        "end": 50406517
      },
      {
        "start": 50406983,
        "end": 50407174
      },
      {
        "start": 50407327,
        "end": 50407415
      },
      {
        "start": 50408785,
        "end": 50408901
      },
      {
        "start": 50409122,
        "end": 50409235
      },
      {
        "start": 50409519,
        "end": 50409666
      },
      {
        "start": 50413426,
        "end": 50413521
      },
      {
        "start": 50413742,
        "end": 50413879
      },
      {
        "start": 50414815,
        "end": 50414990
      },
      {
        "start": 50415438,
        "end": 50415590
      },
      {
        "start": 50415724,
        "end": 50415826
      },
      {
        "start": 50416396,
        "end": 50416528
      },
      {
        "start": 50416610,
        "end": 50416723
      },
      {
        "start": 50417045,
        "end": 50417097
      },
      {
        "start": 50417172,
        "end": 50417269
      },
      {
        "start": 50417842,
        "end": 50418014
      }
    ]
  },
  "POLE": {
    "gene_start": 132623753,
    "gene_end": 132687376,
    "strand": -1,
    "chrom": "12",
    "transcript_id": "ENST00000320574",
    "exons": [
      {
        "start": 132623762,
        "end": 132624810
      },
      {
        "start": 132624905,
        "end": 132624994
      },
      {
        "start": 132625645,
        "end": 132625770
      },
      {
        "start": 132626117,
        "end": 132626317
      },
      {
        "start": 132632315,
        "end": 132632508
      },
      {
        "start": 132632664,
        "end": 132632795
      },
      {
        "start": 132634186,
        "end": 132634378
      },
      {
        "start": 132635892,
        "end": 132636024
      },
      {
        "start": 132638014,
        "end": 132638139
      },
      {
        "start": 132639125,
        "end": 132639298
      },
      {
        "start": 132641647,
        "end": 132641851
      },
      {
        "start": 132642177,
        "end": 132642397
      },
      {
        "start": 132642506,
        "end": 132642729
      },
      {
        "start": 132642820,
        "end": 132642996
      },
      {
        "start": 132643224,
        "end": 132643330
      },
      {
        "start": 132643407,
        "end": 132643560
      },
      {
        "start": 132643837,
        "end": 132643977
      },
      {
        "start": 132648929,
        "end": 132649072
      },
      {
        "start": 132649306,
        "end": 132649515
      },
      {
        "start": 132649677,
        "end": 132649889
      },
      {
        "start": 132657136,
        "end": 132657258
      },
      {
        "start": 132657349,
        "end": 132657429
      },
      {
        "start": 132657868,
        "end": 132657970
      },
      {
        "start": 132659295,
        "end": 132659509
      },
      {
        "start": 132660969,
        "end": 132661164
      },
      {
        "start": 132661527,
        "end": 132661684
      },
      {
        "start": 132664004,
        "end": 132664148
      },
      {
        "start": 132664370,
        "end": 132664462
      },
      {
        "start": 132665302,
        "end": 132665450
      },
      {
        "start": 132667503,
        "end": 132667648
      },
      {
        "start": 132668356,
        "end": 132668502
      },
      {
        "start": 132668635,
        "end": 132668737
      },
      {
        "start": 132668811,
        "end": 132668939
      },
      {
        "start": 132672215,
        "end": 132672322
      },
      {
        "start": 132672627,
        "end": 132672839
      },
      {
        "start": 132673164,
        "end": 132673277
      },
      {
        "start": 132673575,
        "end": 132673707
      },
      {
        "start": 132675398,
        "end": 132675517
      },
      {
        "start": 132675735,
        "end": 132675820
      },
      {
        "start": 132676094,
        "end": 132676204
      },
      {
        "start": 132676546,
        "end": 132676653
      },
      {
        "start": 132677363,
        "end": 132677443
      },
      {
        "start": 132677578,
        "end": 132677719
      },
      {
        "start": 132679497,
        "end": 132679651
      },
      {
        "start": 132679954,
        "end": 132680046
      },
      {
        "start": 132680178,
        "end": 132680222
      },
      {
        "start": 132680607,
        "end": 132680687
      },
      {
        "start": 132681138,
        "end": 132681279
      },
      {
        "start": 132687254,
        "end": 132687342
      }
    ]
  },
  "MUTYH": {
    "gene_start": 45329163,
    "gene_end": 45340893,
    "strand": -1,
    "chrom": "1",
    "transcript_id": "ENST00000456914",
    "exons": [
      {
        "start": 45329242,
        "end": 45329437
      },
      {
        "start": 45330516,
        "end": 45330557
      },
      {
        "start": 45331182,
        "end": 45331334
      },
      {
        "start": 45331420,
        "end": 45331556
      },
      {
        "start": 45331661,
        "end": 45331849
      },
      {
        "start": 45332023,
        "end": 45332086
      },
      {
        "start": 45332166,
        "end": 45332310
      },
      {
        "start": 45332391,
        "end": 45332488
      },
      {
        "start": 45332574,
        "end": 45332687
      },
      {
        "start": 45332763,
        "end": 45332834
      },
      {
        "start": 45332918,
        "end": 45332959
      },
      {
        "start": 45333097,
        "end": 45333170
      },
      {
        "start": 45333285,
        "end": 45333324
      },
      {
        "start": 45333413,
        "end": 45333561
      },
      {
        "start": 45334391,
        "end": 45334511
      },
      {
        "start": 45339899,
        "end": 45339970
      }
    ]
  },
  "TP53": {
    "gene_start": 7661779,
    "gene_end": 7687546,
    "strand": -1,
    "chrom": "17",
    "transcript_id": "ENST00000269305",
    "exons": [
      {
        "start": 7668421,
        "end": 7669690
      },
      {
        "start": 7670609,
        "end": 7670715
      },
      {
        "start": 7673535,
        "end": 7673608
      },
      {
        "start": 7673701,
        "end": 7673837
      },
      {
        "start": 7674181,
        "end": 7674290
      },
      {
        "start": 7674859,
        "end": 7674971
      },
      {
        "start": 7675053,
        "end": 7675236
      },
      {
        "start": 7675994,
        "end": 7676272
      },
      {
        "start": 7676382,
        "end": 7676403
      },
      {
        "start": 7676521,
        "end": 7676622
      },
      {
        "start": 7687377,
        "end": 7687490
      }
    ]
  },
  "BRAF": {
    "gene_start": 140719327,
    "gene_end": 140925199,
    "strand": -1,
    "chrom": "7",
    "transcript_id": "ENST00000646891",
    "exons": [
      {
        "start": 140730665,
        "end": 140734770
      },
      {
        "start": 140739812,
        "end": 140739946
      },
      {
        "start": 140749287,
        "end": 140749418
      },
      {
        "start": 140753275,
        "end": 140753393
      },
      {
        "start": 140754187,
        "end": 140754233
      },
      {
        "start": 140776912,
        "end": 140777088
      },
      {
        "start": 140777991,
        "end": 140778075
      },
      {
        "start": 140781576,
        "end": 140781693
      },
      {
        "start": 140783021,
        "end": 140783157
      },
      {
        "start": 140787548,
        "end": 140787584
      },
      {
        "start": 140794308,
        "end": 140794467
      },
      {
        "start": 140800362,
        "end": 140800481
      },
      {
        "start": 140801412,
        "end": 140801560
      },
      {
        "start": 140807960,
        "end": 140808062
      },
      {
        "start": 140808892,
        "end": 140808995
      },
      {
        "start": 140834609,
        "end": 140834872
      },
      {
        "start": 140850111,
        "end": 140850212
      },
      {
        "start": 140924566,
        "end": 140924929
      }
    ]
  }
}

GENES_STRUCTURE_37 = {
  "MLH1": {
    "gene_start": 37034823,
    "gene_end": 37107380,
    "strand": 1,
    "chrom": "3",
    "transcript_id": "ENST00000231790",
    "exons": [
      {
        "start": 37034823,
        "end": 37035154
      },
      {
        "start": 37038110,
        "end": 37038200
      },
      {
        "start": 37042446,
        "end": 37042544
      },
      {
        "start": 37045892,
        "end": 37045965
      },
      {
        "start": 37048482,
        "end": 37048554
      },
      {
        "start": 37050305,
        "end": 37050396
      },
      {
        "start": 37053311,
        "end": 37053353
      },
      {
        "start": 37053502,
        "end": 37053590
      },
      {
        "start": 37055923,
        "end": 37056035
      },
      {
        "start": 37058997,
        "end": 37059090
      },
      {
        "start": 37061801,
        "end": 37061954
      },
      {
        "start": 37067128,
        "end": 37067498
      },
      {
        "start": 37070275,
        "end": 37070423
      },
      {
        "start": 37081677,
        "end": 37081785
      },
      {
        "start": 37083759,
        "end": 37083822
      },
      {
        "start": 37089010,
        "end": 37089174
      },
      {
        "start": 37090008,
        "end": 37090100
      },
      {
        "start": 37090395,
        "end": 37090508
      },
      {
        "start": 37091977,
        "end": 37092409
      }
    ]
  },
  "MSH2": {
    "gene_start": 47630108,
    "gene_end": 47789450,
    "strand": 1,
    "chrom": "2",
    "transcript_id": "ENST00000233146",
    "exons": [
      {
        "start": 47630108,
        "end": 47630541
      },
      {
        "start": 47635540,
        "end": 47635694
      },
      {
        "start": 47637233,
        "end": 47637511
      },
      {
        "start": 47639553,
        "end": 47639699
      },
      {
        "start": 47641408,
        "end": 47641557
      },
      {
        "start": 47643435,
        "end": 47643568
      },
      {
        "start": 47656881,
        "end": 47657080
      },
      {
        "start": 47672687,
        "end": 47672796
      },
      {
        "start": 47690170,
        "end": 47690293
      },
      {
        "start": 47693797,
        "end": 47693947
      },
      {
        "start": 47698104,
        "end": 47698201
      },
      {
        "start": 47702164,
        "end": 47702409
      },
      {
        "start": 47703506,
        "end": 47703710
      },
      {
        "start": 47705411,
        "end": 47705658
      },
      {
        "start": 47707835,
        "end": 47708010
      },
      {
        "start": 47709918,
        "end": 47710367
      }
    ]
  },
  "MSH6": {
    "gene_start": 47922669,
    "gene_end": 48037240,
    "strand": 1,
    "chrom": "2",
    "transcript_id": "ENST00000234420",
    "exons": [
      {
        "start": 48010221,
        "end": 48010632
      },
      {
        "start": 48018066,
        "end": 48018262
      },
      {
        "start": 48023033,
        "end": 48023202
      },
      {
        "start": 48025750,
        "end": 48028294
      },
      {
        "start": 48030559,
        "end": 48030824
      },
      {
        "start": 48032049,
        "end": 48032166
      },
      {
        "start": 48032757,
        "end": 48032846
      },
      {
        "start": 48033343,
        "end": 48033497
      },
      {
        "start": 48033591,
        "end": 48033790
      },
      {
        "start": 48033918,
        "end": 48037240
      }
    ]
  },
  "PMS2": {
    "gene_start": 6012870,
    "gene_end": 6048756,
    "strand": -1,
    "chrom": "7",
    "transcript_id": "ENST00000265849",
    "exons": [
      {
        "start": 6012870,
        "end": 6013173
      },
      {
        "start": 6017219,
        "end": 6017388
      },
      {
        "start": 6018227,
        "end": 6018327
      },
      {
        "start": 6022455,
        "end": 6022622
      },
      {
        "start": 6026390,
        "end": 6027251
      },
      {
        "start": 6029431,
        "end": 6029586
      },
      {
        "start": 6031604,
        "end": 6031688
      },
      {
        "start": 6035165,
        "end": 6035264
      },
      {
        "start": 6036957,
        "end": 6037054
      },
      {
        "start": 6038739,
        "end": 6038906
      },
      {
        "start": 6042084,
        "end": 6042267
      },
      {
        "start": 6043321,
        "end": 6043423
      },
      {
        "start": 6043603,
        "end": 6043689
      },
      {
        "start": 6045523,
        "end": 6045662
      },
      {
        "start": 6048628,
        "end": 6048756
      }
    ]
  },
  "EPCAM": {
    "gene_start": 47572297,
    "gene_end": 47614740,
    "strand": 1,
    "chrom": "2",
    "transcript_id": "ENST00000263735",
    "exons": [
      {
        "start": 47596287,
        "end": 47596720
      },
      {
        "start": 47600602,
        "end": 47600709
      },
      {
        "start": 47600947,
        "end": 47601187
      },
      {
        "start": 47602373,
        "end": 47602438
      },
      {
        "start": 47604153,
        "end": 47604216
      },
      {
        "start": 47606092,
        "end": 47606193
      },
      {
        "start": 47606908,
        "end": 47607108
      },
      {
        "start": 47612305,
        "end": 47612349
      },
      {
        "start": 47613711,
        "end": 47614173
      }
    ]
  },
  "MLH3": {
    "gene_start": 75480467,
    "gene_end": 75518235,
    "strand": -1,
    "chrom": "14",
    "transcript_id": "ENST00000355774",
    "exons": [
      {
        "start": 75480467,
        "end": 75483904
      },
      {
        "start": 75485532,
        "end": 75485683
      },
      {
        "start": 75489517,
        "end": 75489595
      },
      {
        "start": 75489698,
        "end": 75489721
      },
      {
        "start": 75497246,
        "end": 75497405
      },
      {
        "start": 75498771,
        "end": 75498882
      },
      {
        "start": 75500122,
        "end": 75500193
      },
      {
        "start": 75505043,
        "end": 75505115
      },
      {
        "start": 75506614,
        "end": 75506718
      },
      {
        "start": 75508318,
        "end": 75508403
      },
      {
        "start": 75509082,
        "end": 75509180
      },
      {
        "start": 75513079,
        "end": 75516421
      },
      {
        "start": 75518083,
        "end": 75518235
      }
    ]
  },
  "POLD1": {
    "gene_start": 50887461,
    "gene_end": 50921273,
    "strand": 1,
    "chrom": "19",
    "transcript_id": "ENST00000440232",
    "exons": [
      {
        "start": 50887596,
        "end": 50887647
      },
      {
        "start": 50902108,
        "end": 50902310
      },
      {
        "start": 50902628,
        "end": 50902741
      },
      {
        "start": 50905035,
        "end": 50905181
      },
      {
        "start": 50905256,
        "end": 50905381
      },
      {
        "start": 50905462,
        "end": 50905630
      },
      {
        "start": 50905711,
        "end": 50905792
      },
      {
        "start": 50905869,
        "end": 50905998
      },
      {
        "start": 50906310,
        "end": 50906476
      },
      {
        "start": 50906750,
        "end": 50906854
      },
      {
        "start": 50909439,
        "end": 50909579
      },
      {
        "start": 50909664,
        "end": 50909774
      },
      {
        "start": 50910240,
        "end": 50910431
      },
      {
        "start": 50910584,
        "end": 50910672
      },
      {
        "start": 50912042,
        "end": 50912158
      },
      {
        "start": 50912379,
        "end": 50912492
      },
      {
        "start": 50912776,
        "end": 50912923
      },
      {
        "start": 50916683,
        "end": 50916778
      },
      {
        "start": 50916999,
        "end": 50917136
      },
      {
        "start": 50918072,
        "end": 50918247
      },
      {
        "start": 50918695,
        "end": 50918847
      },
      {
        "start": 50918981,
        "end": 50919083
      },
      {
        "start": 50919653,
        "end": 50919785
      },
      {
        "start": 50919867,
        "end": 50919980
      },
      {
        "start": 50920302,
        "end": 50920354
      },
      {
        "start": 50920429,
        "end": 50920526
      },
      {
        "start": 50921099,
        "end": 50921271
      }
    ]
  },
  "POLE": {
    "gene_start": 133200348,
    "gene_end": 133263951,
    "strand": -1,
    "chrom": "12",
    "transcript_id": "ENST00000320574",
    "exons": [
      {
        "start": 133200348,
        "end": 133201396
      },
      {
        "start": 133201491,
        "end": 133201580
      },
      {
        "start": 133202231,
        "end": 133202356
      },
      {
        "start": 133202703,
        "end": 133202903
      },
      {
        "start": 133208901,
        "end": 133209094
      },
      {
        "start": 133209250,
        "end": 133209381
      },
      {
        "start": 133210772,
        "end": 133210964
      },
      {
        "start": 133212478,
        "end": 133212610
      },
      {
        "start": 133214600,
        "end": 133214725
      },
      {
        "start": 133215711,
        "end": 133215884
      },
      {
        "start": 133218233,
        "end": 133218437
      },
      {
        "start": 133218763,
        "end": 133218983
      },
      {
        "start": 133219092,
        "end": 133219315
      },
      {
        "start": 133219406,
        "end": 133219582
      },
      {
        "start": 133219810,
        "end": 133219916
      },
      {
        "start": 133219993,
        "end": 133220146
      },
      {
        "start": 133220423,
        "end": 133220563
      },
      {
        "start": 133225515,
        "end": 133225658
      },
      {
        "start": 133225892,
        "end": 133226101
      },
      {
        "start": 133226263,
        "end": 133226475
      },
      {
        "start": 133233722,
        "end": 133233844
      },
      {
        "start": 133233935,
        "end": 133234015
      },
      {
        "start": 133234454,
        "end": 133234556
      },
      {
        "start": 133235881,
        "end": 133236095
      },
      {
        "start": 133237555,
        "end": 133237750
      },
      {
        "start": 133238113,
        "end": 133238270
      },
      {
        "start": 133240590,
        "end": 133240734
      },
      {
        "start": 133240956,
        "end": 133241048
      },
      {
        "start": 133241888,
        "end": 133242036
      },
      {
        "start": 133244089,
        "end": 133244234
      },
      {
        "start": 133244942,
        "end": 133245088
      },
      {
        "start": 133245221,
        "end": 133245323
      },
      {
        "start": 133245397,
        "end": 133245525
      },
      {
        "start": 133248801,
        "end": 133248908
      },
      {
        "start": 133249213,
        "end": 133249425
      },
      {
        "start": 133249750,
        "end": 133249863
      },
      {
        "start": 133250161,
        "end": 133250293
      },
      {
        "start": 133251984,
        "end": 133252103
      },
      {
        "start": 133252321,
        "end": 133252406
      },
      {
        "start": 133252680,
        "end": 133252790
      },
      {
        "start": 133253132,
        "end": 133253239
      },
      {
        "start": 133253949,
        "end": 133254029
      },
      {
        "start": 133254164,
        "end": 133254305
      },
      {
        "start": 133256083,
        "end": 133256237
      },
      {
        "start": 133256540,
        "end": 133256632
      },
      {
        "start": 133256764,
        "end": 133256808
      },
      {
        "start": 133257193,
        "end": 133257273
      },
      {
        "start": 133257724,
        "end": 133257865
      },
      {
        "start": 133263840,
        "end": 133263945
      }
    ]
  },
  "MUTYH": {
    "gene_start": 45794835,
    "gene_end": 45806142,
    "strand": -1,
    "chrom": "1",
    "transcript_id": "ENST00000372098",
    "exons": [
      {
        "start": 45794914,
        "end": 45795109
      },
      {
        "start": 45796188,
        "end": 45796229
      },
      {
        "start": 45796854,
        "end": 45797006
      },
      {
        "start": 45797092,
        "end": 45797228
      },
      {
        "start": 45797333,
        "end": 45797521
      },
      {
        "start": 45797695,
        "end": 45797758
      },
      {
        "start": 45797838,
        "end": 45797982
      },
      {
        "start": 45798063,
        "end": 45798160
      },
      {
        "start": 45798246,
        "end": 45798359
      },
      {
        "start": 45798435,
        "end": 45798506
      },
      {
        "start": 45798590,
        "end": 45798631
      },
      {
        "start": 45798769,
        "end": 45798842
      },
      {
        "start": 45798957,
        "end": 45798996
      },
      {
        "start": 45799085,
        "end": 45799266
      },
      {
        "start": 45800063,
        "end": 45800183
      },
      {
        "start": 45805891,
        "end": 45806060
      }
    ]
  },
  "TP53": {
    "gene_start": 7565097,
    "gene_end": 7590856,
    "strand": -1,
    "chrom": "17",
    "transcript_id": "ENST00000269305",
    "exons": [
      {
        "start": 7571720,
        "end": 7573008
      },
      {
        "start": 7573927,
        "end": 7574033
      },
      {
        "start": 7576853,
        "end": 7576926
      },
      {
        "start": 7577019,
        "end": 7577155
      },
      {
        "start": 7577499,
        "end": 7577608
      },
      {
        "start": 7578177,
        "end": 7578289
      },
      {
        "start": 7578371,
        "end": 7578554
      },
      {
        "start": 7579312,
        "end": 7579590
      },
      {
        "start": 7579700,
        "end": 7579721
      },
      {
        "start": 7579839,
        "end": 7579940
      },
      {
        "start": 7590695,
        "end": 7590856
      }
    ]
  },
  "BRAF": {
    "gene_start": 140419127,
    "gene_end": 140624564,
    "strand": -1,
    "chrom": "7",
    "transcript_id": "ENST00000288602",
    "exons": [
      {
        "start": 140434279,
        "end": 140434570
      },
      {
        "start": 140439612,
        "end": 140439746
      },
      {
        "start": 140449087,
        "end": 140449218
      },
      {
        "start": 140453075,
        "end": 140453193
      },
      {
        "start": 140453987,
        "end": 140454033
      },
      {
        "start": 140476712,
        "end": 140476888
      },
      {
        "start": 140477791,
        "end": 140477875
      },
      {
        "start": 140481376,
        "end": 140481493
      },
      {
        "start": 140482821,
        "end": 140482957
      },
      {
        "start": 140487348,
        "end": 140487384
      },
      {
        "start": 140494108,
        "end": 140494267
      },
      {
        "start": 140500162,
        "end": 140500281
      },
      {
        "start": 140501212,
        "end": 140501360
      },
      {
        "start": 140507760,
        "end": 140507862
      },
      {
        "start": 140508692,
        "end": 140508795
      },
      {
        "start": 140534409,
        "end": 140534672
      },
      {
        "start": 140549911,
        "end": 140550012
      },
      {
        "start": 140624366,
        "end": 140624564
      }
    ]
  }
}


import os
import sys
import json
import gzip
import re
import csv
import argparse
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.ticker as ticker
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D

# Import gene structures from the original script

# Chromosome lengths and centromere definitions (GRCh38 & GRCh37)
CHROM_INFO_38 = {
    "1":  {"len": 248956422, "cen": (121535434, 126135345)},
    "2":  {"len": 242193529, "cen": (92326171, 95326171)},
    "3":  {"len": 198295559, "cen": (90500000, 93700000)},
    "4":  {"len": 190214555, "cen": (49660117, 52660117)},
    "5":  {"len": 181538259, "cen": (46405641, 49405641)},
    "6":  {"len": 170805979, "cen": (58830166, 61830166)},
    "7":  {"len": 159345973, "cen": (58054331, 61054331)},
    "8":  {"len": 145138636, "cen": (43838887, 46838887)},
    "9":  {"len": 138394717, "cen": (47367679, 50367679)},
    "10": {"len": 133797422, "cen": (39254935, 42254935)},
    "11": {"len": 135086622, "cen": (51644205, 54644205)},
    "12": {"len": 133275309, "cen": (34856694, 37856694)},
    "13": {"len": 114364328, "cen": (16000000, 19000000)},
    "14": {"len": 107043718, "cen": (16000000, 19000000)},
    "15": {"len": 101991189, "cen": (17000000, 20000000)},
    "16": {"len": 90338345,  "cen": (35335801, 38335801)},
    "17": {"len": 83257441,  "cen": (22263006, 25263006)},
    "18": {"len": 80373285,  "cen": (15460898, 18460898)},
    "19": {"len": 58617616,  "cen": (24681782, 27681782)},
    "20": {"len": 64444167,  "cen": (26369569, 29369569)},
    "21": {"len": 46709983,  "cen": (11288129, 14288129)},
    "22": {"len": 50818468,  "cen": (13000000, 16000000)},
    "X":  {"len": 156040895, "cen": (58605579, 61605579)},
    "Y":  {"len": 57227415,  "cen": (10104553, 13104553)}
}

CHROM_INFO_37 = {
    "1":  {"len": 249250621, "cen": (121535434, 126135345)},
    "2":  {"len": 243199373, "cen": (92326171, 95326171)},
    "3":  {"len": 198022430, "cen": (90500000, 93700000)},
    "4":  {"len": 191154276, "cen": (49660117, 52660117)},
    "5":  {"len": 180915260, "cen": (46405641, 49405641)},
    "6":  {"len": 171115067, "cen": (58830166, 61830166)},
    "7":  {"len": 159138663, "cen": (58054331, 61054331)},
    "8":  {"len": 146364022, "cen": (43838887, 46838887)},
    "9":  {"len": 141213431, "cen": (47367679, 50367679)},
    "10": {"len": 135534747, "cen": (39254935, 42254935)},
    "11": {"len": 135006516, "cen": (51644205, 54644205)},
    "12": {"len": 133851895, "cen": (34856694, 37856694)},
    "13": {"len": 115169878, "cen": (16000000, 19000000)},
    "14": {"len": 107349540, "cen": (16000000, 19000000)},
    "15": {"len": 102531392, "cen": (17000000, 20000000)},
    "16": {"len": 90354753,  "cen": (35335801, 38335801)},
    "17": {"len": 81195210,  "cen": (22263006, 25263006)},
    "18": {"len": 78077248,  "cen": (15460898, 18460898)},
    "19": {"len": 59128983,  "cen": (24681782, 27681782)},
    "20": {"len": 63025520,  "cen": (26369569, 29369569)},
    "21": {"len": 48129895,  "cen": (11288129, 14288129)},
    "22": {"len": 51304566,  "cen": (13000000, 16000000)},
    "X":  {"len": 155270560, "cen": (58605579, 61605579)},
    "Y":  {"len": 59373566,  "cen": (10104553, 13104553)}
}

LANDMARK_GENES = {
    "1": [{"name": "ARID1A", "pos": 26.7e6}, {"name": "MUTYH", "pos": 45.3e6}],
    "2": [{"name": "MSH2/6", "pos": 47.6e6}, {"name": "EPCAM", "pos": 47.6e6}, {"name": "IDH1", "pos": 208.2e6}],
    "3": [{"name": "SETD2/VHL", "pos": 10.15e6}, {"name": "MLH1", "pos": 37.02e6}, {"name": "BAP1", "pos": 52.44e6}, {"name": "MITF", "pos": 69.85e6}, {"name": "PIK3CA", "pos": 186.3e6}],
    "4": [{"name": "FGFR3", "pos": 1.8e6}, {"name": "KIT/PDGFRA", "pos": 54.7e6}],
    "5": [{"name": "TERT", "pos": 1.2e6}, {"name": "APC", "pos": 112.7e6}],
    "7": [{"name": "PMS2", "pos": 5.97e6}, {"name": "EGFR", "pos": 55.1e6}, {"name": "MET", "pos": 116.4e6}, {"name": "BRAF", "pos": 140.7e6}],
    "8": [{"name": "FGFR1", "pos": 38.3e6}, {"name": "MYC", "pos": 127.7e6}],
    "9": [{"name": "JAK2", "pos": 5.0e6}, {"name": "CDKN2A", "pos": 21.9e6}, {"name": "NOTCH1", "pos": 136.5e6}],
    "10": [{"name": "RET", "pos": 43.1e6}, {"name": "PTEN", "pos": 87.8e6}],
    "11": [{"name": "WT1", "pos": 32.4e6}, {"name": "CCND1", "pos": 69.4e6}, {"name": "ATM", "pos": 108.2e6}],
    "12": [{"name": "KRAS", "pos": 25.2e6}, {"name": "CDK4", "pos": 58.1e6}, {"name": "MDM2", "pos": 68.8e6}, {"name": "POLE", "pos": 132.6e6}],
    "13": [{"name": "BRCA2", "pos": 32.3e6}, {"name": "RB1", "pos": 48.3e6}],
    "14": [{"name": "MLH3", "pos": 75.0e6}, {"name": "AKT1", "pos": 104.8e6}],
    "15": [{"name": "SMAD4", "pos": 48.6e6}],
    "16": [{"name": "CDH1", "pos": 68.8e6}],
    "17": [{"name": "TP53", "pos": 7.67e6}, {"name": "NF1", "pos": 31.2e6}, {"name": "ERBB2", "pos": 39.7e6}, {"name": "BRCA1", "pos": 43.1e6}],
    "18": [{"name": "SMAD4", "pos": 51.0e6}],
    "19": [{"name": "STK11", "pos": 1.2e6}, {"name": "POLD1", "pos": 50.4e6}],
    "20": [{"name": "SRC", "pos": 37.4e6}, {"name": "GNAS", "pos": 58.9e6}],
    "21": [{"name": "RUNX1", "pos": 34.8e6}],
    "22": [{"name": "CHEK2", "pos": 28.7e6}, {"name": "NF2", "pos": 29.6e6}],
    "X":  [{"name": "AR", "pos": 67.5e6}],
}


def compute_rolling_cn_stats(x_eval, x_obs, y_obs, bandwidth=3.5e6, max_gap=7.0e6):
    """
    Computes kernel-weighted rolling median, MAD, and Standard Error for Copy Number probes.
    Breaks across deserts where distance to nearest probe > max_gap.
    """
    res_center = []
    res_mad = []
    res_se = []
    res_neff = []
    for x in x_eval:
        dists = np.abs(x_obs - x)
        if len(dists) == 0 or np.min(dists) > max_gap:
            res_center.append(np.nan)
            res_mad.append(np.nan)
            res_se.append(np.nan)
            res_neff.append(0.0)
            continue
        weights = np.exp(-0.5 * (dists / bandwidth) ** 2)
        idx = np.where(weights > 0.02)[0]
        if len(idx) == 0:
            res_center.append(np.nan)
            res_mad.append(np.nan)
            res_se.append(np.nan)
            res_neff.append(0.0)
            continue
        w_sub = weights[idx]
        y_sub = y_obs[idx]
        neff = np.sum(w_sub)

        # Weighted median
        sorter = np.argsort(y_sub)
        y_s, w_s = y_sub[sorter], w_sub[sorter]
        cw = np.cumsum(w_s)
        med_val = y_s[min(np.searchsorted(cw, 0.5 * cw[-1]), len(y_s) - 1)]

        # Weighted MAD
        devs = np.abs(y_sub - med_val)
        sorter_dev = np.argsort(devs)
        dev_s, w_dev_s = devs[sorter_dev], w_sub[sorter_dev]
        cw_dev = np.cumsum(w_dev_s)
        mad_val = dev_s[min(np.searchsorted(cw_dev, 0.5 * cw_dev[-1]), len(dev_s) - 1)]

        # Standard error of median: 1.4826 * (MAD + reg) / sqrt(neff)
        se_val = 1.4826 * (mad_val + 0.03) / np.sqrt(max(1.0, neff))

        res_center.append(med_val)
        res_mad.append(mad_val)
        res_se.append(se_val)
        res_neff.append(neff)
    return np.array(res_center), np.array(res_mad), np.array(res_se), np.array(res_neff)


def compute_rolling_baf_stats(x_eval, x_obs, y_obs, bandwidth=6.0e6, max_gap=12.0e6):
    """
    Computes kernel-weighted rolling median and MAD for BAF on folded BAF (|BAF - 0.5| + 0.5).
    Produces major allele center and minor allele center with corresponding confidence ribbons.
    """
    y_folded = np.abs(y_obs - 0.5) + 0.5
    res_center = []
    res_mad = []
    res_se = []
    res_neff = []
    for x in x_eval:
        dists = np.abs(x_obs - x)
        if len(dists) == 0 or np.min(dists) > max_gap:
            res_center.append(np.nan)
            res_mad.append(np.nan)
            res_se.append(np.nan)
            res_neff.append(0.0)
            continue
        weights = np.exp(-0.5 * (dists / bandwidth) ** 2)
        idx = np.where(weights > 0.02)[0]
        if len(idx) == 0:
            res_center.append(np.nan)
            res_mad.append(np.nan)
            res_se.append(np.nan)
            res_neff.append(0.0)
            continue
        w_sub = weights[idx]
        y_sub = y_folded[idx]
        neff = np.sum(w_sub)

        sorter = np.argsort(y_sub)
        y_s, w_s = y_sub[sorter], w_sub[sorter]
        cw = np.cumsum(w_s)
        med_val = y_s[min(np.searchsorted(cw, 0.5 * cw[-1]), len(y_s) - 1)]

        devs = np.abs(y_sub - med_val)
        sorter_dev = np.argsort(devs)
        dev_s, w_dev_s = devs[sorter_dev], w_sub[sorter_dev]
        cw_dev = np.cumsum(w_dev_s)
        mad_val = dev_s[min(np.searchsorted(cw_dev, 0.5 * cw_dev[-1]), len(dev_s) - 1)]

        se_val = 1.4826 * (mad_val + 0.02) / np.sqrt(max(1.0, neff))

        res_center.append(med_val)
        res_mad.append(mad_val)
        res_se.append(se_val)
        res_neff.append(neff)
    return np.array(res_center), np.array(res_mad), np.array(res_se), np.array(res_neff)


def plot_ideogram_bar(ax, chrom, chr_len, centromere, gene_start, gene_end, gene_name, is_grch37=False):
    """Draws a clean, publication-ready cytoband ideogram for the chromosome."""
    ax.set_xlim(0, chr_len)
    ax.set_ylim(-0.2, 1.95)
    ax.axis('off')

    cen_s, cen_e = centromere

    # Draw p-arm
    rect_p = patches.FancyBboxPatch((0, 0.0), cen_s, 0.65,
                                    boxstyle="round,pad=0.02,rounding_size=2e6",
                                    facecolor='#f1f5f9', edgecolor='#475569', linewidth=1.2)
    ax.add_patch(rect_p)
    ax.text(cen_s * 0.45, 0.32, f"{chrom}p arm (0–{cen_s/1e6:.1f} Mb)", ha='center', va='center',
            fontsize=8.5, fontweight='bold', color='#334155')

    # Draw Centromere
    rect_cen = patches.Rectangle((cen_s, 0.0), cen_e - cen_s, 0.65,
                                 facecolor='#e2e8f0', edgecolor='#475569', linewidth=1.2, hatch='///')
    ax.add_patch(rect_cen)
    ax.text((cen_s + cen_e) / 2.0, 0.95, f"Cen ({cen_s/1e6:.1f}–{cen_e/1e6:.1f}M)", ha='center', va='bottom',
            fontsize=7.5, fontweight='bold', color='#64748b')

    # Draw q-arm
    rect_q = patches.FancyBboxPatch((cen_e, 0.0), chr_len - cen_e, 0.65,
                                    boxstyle="round,pad=0.02,rounding_size=2e6",
                                    facecolor='#f1f5f9', edgecolor='#475569', linewidth=1.2)
    ax.add_patch(rect_q)
    ax.text((cen_e + chr_len) / 2.0, 0.32, f"{chrom}q arm ({cen_e/1e6:.1f}–{chr_len/1e6:.1f} Mb)", ha='center', va='center',
            fontsize=8.5, fontweight='bold', color='#334155')

    # Draw panel landmarks
    raw_landmarks = LANDMARK_GENES.get(chrom, [])
    landmarks = []
    has_target = False
    for lm in raw_landmarks:
        if lm['name'].startswith(gene_name) or lm['name'] == gene_name:
            has_target = True
            landmarks.append({"name": f"{gene_name} (Target)", "pos": (gene_start + gene_end)/2.0, "is_target": True})
        else:
            landmarks.append(lm)
    if not has_target:
        landmarks.append({"name": f"{gene_name} (Target)", "pos": (gene_start + gene_end)/2.0, "is_target": True})

    for g in landmarks:
        pos = g['pos']
        is_t = g.get('is_target', False)
        color = '#dc2626' if is_t else '#0284c7'
        y_tick = 1.15 if is_t else 0.85
        y_text = 1.25 if is_t else 0.95
        ax.plot([pos, pos], [0.65, y_tick], color=color, linewidth=1.8 if is_t else 1.0)
        ax.text(pos, y_text, g['name'], ha='center', va='bottom', fontsize=8.2 if is_t else 7.2,
                fontweight='bold', color=color,
                bbox=dict(boxstyle='round,pad=0.2', facecolor='#ffffff', edgecolor=color, linewidth=1.2 if is_t else 0.8, alpha=0.95))

    # Highlight star for target gene on ideogram bar
    g_mid = (gene_start + gene_end) / 2.0
    ax.scatter([g_mid], [0.32], c='#dc2626', marker='*', s=160, zorder=5, edgecolor='#881337', linewidth=1.2)


def render_gene_baf_v10(
    gene,
    gene_info,
    amber_baf_file,
    purple_segments=None,
    target_cn_file=None,
    variants=None,
    methylation_data=None,
    purity=1.0,
    sample_id="Sample",
    genome="GRCh38_hmf",
    output_png="plot.png",
    output_html=None
):
    chrom = str(gene_info['chrom']).upper().replace('CHR', '')
    g_start = gene_info['gene_start']
    g_end = gene_info['gene_end']
    strand = gene_info['strand']
    transcript_id = gene_info.get('transcript_id', '')

    is_grch37 = "37" in str(genome)
    chrom_dict = CHROM_INFO_37 if is_grch37 else CHROM_INFO_38
    chr_len = chrom_dict.get(chrom, {}).get('len', 200e6)
    centromere = chrom_dict.get(chrom, {}).get('cen', (chr_len * 0.45, chr_len * 0.50))
    
    # Process Exons
    raw_exons = gene_info['exons']
    num_exons = len(raw_exons)
    exons = []
    for idx, ex in enumerate(raw_exons):
        e_num = ex.get('exon_num')
        if e_num is None:
            e_num = (idx + 1) if strand == 1 else (num_exons - idx)
        exons.append({
            'start': ex['start'],
            'end': ex['end'],
            'exon_num': int(e_num)
        })

    # Promoter Coordinates
    if strand == 1:
        promoter_start = gene_info.get('promoter_start', g_start - 3500)
        promoter_end = gene_info.get('promoter_end', g_start + 500)
    else:
        promoter_start = gene_info.get('promoter_start', g_end - 500)
        promoter_end = gene_info.get('promoter_end', g_end + 3500)

    # 1. Load Segments
    df_seg = pd.DataFrame()
    if purple_segments and os.path.exists(purple_segments):
        try:
            df_seg = pd.read_csv(purple_segments, sep='\t')
            df_seg['norm_chrom'] = df_seg['chromosome'].astype(str).str.replace('chr', '', case=False)
            df_seg = df_seg[df_seg['norm_chrom'] == chrom].copy()
        except Exception as e:
            print(f"Warning reading segments: {e}", file=sys.stderr)

    # 2. Load Probes
    df_raw = pd.DataFrame()
    if target_cn_file and os.path.exists(target_cn_file):
        try:
            df_raw = pd.read_csv(target_cn_file, sep='\t')
            df_raw['chromosome'] = df_raw['chromosome'].astype(str).str.replace('chr', '', case=False)
        except Exception as e:
            print(f"Warning reading target cn: {e}", file=sys.stderr)
    elif purple_segments and os.path.exists(purple_segments):
        cand = purple_segments.replace('.segment.tsv', '.target_region_cn.tsv')
        if os.path.exists(cand):
            try:
                df_raw = pd.read_csv(cand, sep='\t')
                df_raw['chromosome'] = df_raw['chromosome'].astype(str).str.replace('chr', '', case=False)
            except Exception:
                pass
        else:
            p_dir = os.path.dirname(purple_segments)
            for fname in os.listdir(p_dir):
                if fname.endswith('.target_region_cn.tsv'):
                    try:
                        df_raw = pd.read_csv(os.path.join(p_dir, fname), sep='\t')
                        df_raw['chromosome'] = df_raw['chromosome'].astype(str).str.replace('chr', '', case=False)
                        break
                    except Exception:
                        pass

    # Extract chromosome-wide probes and target gene probes
    chr_probes = pd.DataFrame()
    gene_probes = pd.DataFrame()
    seg_cn = 2.0
    minor_cn = 1.0
    if not df_raw.empty:
        chr_probes = df_raw[
            (df_raw['chromosome'] == chrom) &
            (df_raw['windowTumorRatio'] > 0)
        ].copy()
        if not chr_probes.empty:
            chr_probes['probe_mid'] = (chr_probes['windowStart'] + chr_probes['windowEnd']) / 2.0
            chr_probes['norm_cn'] = chr_probes['windowTumorRatio'] * 2.0
            gene_probes = chr_probes[
                (chr_probes['windowStart'] >= g_start - 10000) &
                (chr_probes['windowEnd'] <= g_end + 10000)
            ].copy()
            if not gene_probes.empty:
                if 'copyNumber' in gene_probes.columns:
                    seg_cn = float(gene_probes['copyNumber'].iloc[0])
                if 'minorAlleleCopyNumber' in gene_probes.columns:
                    minor_cn = float(gene_probes['minorAlleleCopyNumber'].iloc[0])

    if seg_cn == 2.0 and not df_seg.empty:
        gene_segs = df_seg[
            (df_seg['end'] >= g_start) &
            (df_seg['start'] <= g_end)
        ]
        if not gene_segs.empty:
            s_row = gene_segs.iloc[0]
            if 'tumorCopyNumber' in s_row:
                seg_cn = float(s_row['tumorCopyNumber'])
            elif 'copyNumber' in s_row:
                seg_cn = float(s_row['copyNumber'])
            if 'minorAlleleCopyNumber' in s_row:
                minor_cn = float(s_row['minorAlleleCopyNumber'])
            elif 'minorAllelePloidy' in s_row:
                minor_cn = float(s_row['minorAllelePloidy'])

    major_cn = max(0.0, seg_cn - minor_cn)
    frac_minor = (minor_cn / seg_cn) if seg_cn > 0 else 0.5
    frac_major = (major_cn / seg_cn) if seg_cn > 0 else 0.5
    is_balanced_diploid = (1.7 <= seg_cn <= 2.3 and minor_cn >= 0.70)

    # 3. Load AMBER BAF data
    df_baf = pd.DataFrame()
    if amber_baf_file and os.path.exists(amber_baf_file):
        try:
            df_baf = pd.read_csv(amber_baf_file, sep='\t', compression='gzip' if amber_baf_file.endswith('.gz') else None)
            df_baf['norm_chrom'] = df_baf['chromosome'].astype(str).str.replace('chr', '', case=False)
            df_baf = df_baf[df_baf['norm_chrom'] == chrom].copy()
            df_baf = df_baf[(df_baf['tumorBAF'] > 0.05) & (df_baf['tumorBAF'] < 0.95)]
        except Exception as e:
            print(f"Warning reading AMBER BAF: {e}", file=sys.stderr)

    # Adjust chr_len if data extends beyond nominal boundaries
    max_data_pos = max(
        g_end + 2e6,
        df_seg['end'].max() if not df_seg.empty else 0,
        chr_probes['windowEnd'].max() if not chr_probes.empty else 0,
        df_baf['position'].max() if not df_baf.empty else 0
    )
    chr_len = max(chr_len, max_data_pos)

    # 4. Process Variants
    if variants is None:
        variants = []
    gene_vars = [v for v in variants if v.get('gene') == gene or v.get('gene_name') == gene]
    if not gene_vars and variants:
        gene_vars = variants

    parsed_vars = []
    purity_clamped = max(0.01, min(1.0, float(purity)))
    for v in gene_vars:
        # Skip epigenetic promoter methylation entries bundled in variants list
        v_type = str(v.get('var_type', '')).lower()
        v_phasing = str(v.get('phasing', '')).lower()
        v_desc_check = str(v.get('change') or v.get('desc') or v.get('label') or '').lower()
        if 'methylation' in v_type or 'epigenetic' in v_phasing or 'methylation' in v_desc_check:
            continue
            
        v_pos = v.get('pos')
        if not v_pos:
            continue
        v_pos = int(v_pos)
        v_desc = v.get('change') or v.get('desc')
        if not v_desc and (v.get('label') or v.get('variant_desc')):
            lbl = v.get('label') or v.get('variant_desc')
            if ' (' in lbl:
                lbl = lbl.split(' (')[0]
            for pfx in [f"Somatic {gene} ", f"Germline {gene} ", f"{gene} "]:
                if lbl.startswith(pfx):
                    lbl = lbl[len(pfx):]
            v_desc = lbl
        if not v_desc:
            v_desc = f"chr{chrom}:{v_pos}"
        v_exon = v.get('exon')
        if v_exon is None:
            for ex in exons:
                if ex['start'] <= v_pos <= ex['end']:
                    v_exon = ex['exon_num']
                    break
        v_af = v.get('sample_vaf') if v.get('sample_vaf') is not None else (v.get('tumor_vaf') if v.get('tumor_vaf') is not None else v.get('vaf'))
        try:
            v_af = float(v_af) if v_af is not None else 0.0
        except ValueError:
            v_af = 0.0

        p_af = v.get('purple_af')
        if p_af is None:
            if purity_clamped * seg_cn > 0:
                p_af = (v_af * (purity_clamped * seg_cn + 2.0 * (1.0 - purity_clamped))) / (purity_clamped * seg_cn)
            else:
                p_af = v_af
            p_af = min(1.0, max(0.0, float(p_af)))
        else:
            p_af = float(p_af)

        v_vcn = v.get('purple_vcn')
        if v_vcn is None:
            v_vcn = p_af * seg_cn
        else:
            v_vcn = float(v_vcn)

        # Refined Allele Phasing Logic
        if is_balanced_diploid:
            if v_vcn >= 1.6:
                phased_allele_str = "Homozygous / Duplicated on Diploid (2 copies)"
                phasing_call_str = f"Homozygous somatic variant ({v_vcn:.2f} copies on balanced 1:1 diploid background)"
            else:
                phased_allele_str = "Heterozygous Diploid (1 copy)"
                phasing_call_str = f"Monoallelic mutation ({v_vcn:.2f} copies on balanced 1:1 diploid background)"
        else:
            # Asymmetric copy number / LOH background
            if major_cn >= 1.4 and v_vcn >= major_cn * 0.75:
                phased_allele_str = "Major Allele (Duplicated / Both Copies)"
                phasing_call_str = f"Clonal on Major Allele (duplicated across {v_vcn:.2f} copies; pre-duplication hit)"
            elif major_cn >= 1.4 and v_vcn > (minor_cn + 0.35):
                phased_allele_str = "Major Allele (Single Copy)"
                phasing_call_str = f"Mutation on 1 copy of Major Allele ({v_vcn:.2f} of {major_cn:.2f} major copies; post-duplication hit)"
            elif minor_cn < 0.35 and v_vcn <= 0.35:
                phased_allele_str = "Minor Allele (Lost / Residual)"
                phasing_call_str = f"Residual / subclonal hit on deleted Minor Allele ({v_vcn:.2f} copies)"
            elif v_vcn <= minor_cn * 1.25:
                phased_allele_str = "Minor Allele"
                phasing_call_str = f"Clonal mutation on Minor Allele ({v_vcn:.2f} of {minor_cn:.2f} minor copies)"
            elif v_vcn > minor_cn:
                phased_allele_str = "Major Allele"
                phasing_call_str = f"Mutation on Major Allele ({v_vcn:.2f} of {major_cn:.2f} copies)"
            else:
                phased_allele_str = "Minor Allele"
                phasing_call_str = f"Mutation on Minor Allele ({v_vcn:.2f} of {minor_cn:.2f} copies)"

        parsed_vars.append({
            'pos': v_pos,
            'desc': v_desc,
            'exon': v_exon,
            'sample_vaf': v_af,
            'tumor_af': p_af,
            'vcn': v_vcn,
            'phased_allele': phased_allele_str,
            'phasing_call': phasing_call_str,
            'is_somatic': (v.get('type', 'somatic').lower() == 'somatic')
        })

    # 5. Process Methylation Data
    has_meth = False
    tumor_promoter_beta = None
    normal_promoter_beta = None
    methylation_status = "NOT ASSAYED"
    cpg_list = []
    if methylation_data:
        has_meth = True
        tumor_promoter_beta = float(methylation_data.get('tumor_beta', 0.0) or 0.0)
        norm_b = methylation_data.get('normal_beta')
        normal_promoter_beta = float(norm_b) if norm_b is not None else None
        methylation_status = str(methylation_data.get('status', 'unmethylated')).upper()
        cpg_list = methylation_data.get('cpgs', [])

    is_unmeth = (methylation_status.lower() == 'unmethylated' or (tumor_promoter_beta is not None and tumor_promoter_beta < 0.20))

    # Focal padding
    pad_focal = max(9000, int((g_end - g_start) * 0.18))
    x_min_f = g_start - pad_focal
    x_max_f = g_end + pad_focal

    # -------------------------------------------------------------
    # SETUP FIGURE: OPTION A STACKED MACRO-TO-MICRO LAYOUT
    # -------------------------------------------------------------
    fig = plt.figure(figsize=(19, 15.0), dpi=220)
    outer_gs = GridSpec(2, 1, height_ratios=[4.8, 5.4], hspace=0.30,
                        left=0.07, right=0.96, top=0.93, bottom=0.05)

    # Top Macro Grid (Ideogram + Macro CN + Macro BAF)
    macro_gs = GridSpec(3, 1, height_ratios=[0.5, 2.5, 2.3], hspace=0.18)
    macro_gs.update(top=outer_gs[0].get_position(fig).y1, bottom=outer_gs[0].get_position(fig).y0,
                    left=outer_gs[0].get_position(fig).x0, right=outer_gs[0].get_position(fig).x1)

    ax_ideo = fig.add_subplot(macro_gs[0])
    ax_macro_cn = fig.add_subplot(macro_gs[1])
    ax_macro_baf = fig.add_subplot(macro_gs[2], sharex=ax_macro_cn)

    # Bottom Micro Grid (Micro CN + Micro BAF + Micro Gene/Methylation)
    micro_gs = GridSpec(3, 1, height_ratios=[2.5, 2.2, 2.5], hspace=0.18)
    micro_gs.update(top=outer_gs[1].get_position(fig).y1, bottom=outer_gs[1].get_position(fig).y0,
                    left=outer_gs[1].get_position(fig).x0, right=outer_gs[1].get_position(fig).x1)

    ax_micro_cn = fig.add_subplot(micro_gs[0])
    ax_micro_baf = fig.add_subplot(micro_gs[1], sharex=ax_micro_cn)
    ax_micro_gene = fig.add_subplot(micro_gs[2], sharex=ax_micro_cn)

    # -------------------------------------------------------------
    # TOP SECTION: WHOLE CHROMOSOME MACRO TRACKS (0–chr_len)
    # -------------------------------------------------------------
    plot_ideogram_bar(ax_ideo, chrom, chr_len, centromere, g_start, g_end, gene, is_grch37=is_grch37)

    # Rolling stats across chromosome
    x_grid = np.linspace(0, chr_len, 600)
    c_cn, m_cn, se_cn, n_cn = np.full(600, np.nan), np.full(600, np.nan), np.full(600, np.nan), np.zeros(600)
    if not chr_probes.empty:
        c_cn, m_cn, se_cn, n_cn = compute_rolling_cn_stats(x_grid, chr_probes['probe_mid'].values, chr_probes['norm_cn'].values)

    c_baf, m_baf, se_baf, n_baf = np.full(600, np.nan), np.full(600, np.nan), np.full(600, np.nan), np.zeros(600)
    if not df_baf.empty:
        c_baf, m_baf, se_baf, n_baf = compute_rolling_baf_stats(x_grid, df_baf['position'].values, df_baf['tumorBAF'].values)

    # 1. Macro Copy Number
    ax_macro_cn.set_facecolor('#ffffff')
    ax_macro_cn.grid(True, linestyle='--', alpha=0.35, color='#cbd5e1')
    ax_macro_cn.set_xlim(0, chr_len)
    ax_macro_cn.set_ylim(-0.25, 4.4)
    ax_macro_cn.axhline(2.0, color='#64748b', linestyle='--', linewidth=1.1, alpha=0.85)
    ax_macro_cn.axhline(1.0, color='#f59e0b', linestyle=':', linewidth=0.9, alpha=0.7)
    ax_macro_cn.axhline(3.0, color='#10b981', linestyle=':', linewidth=0.9, alpha=0.7)
    ax_macro_cn.axhspan(1.75, 2.25, color='#f8fafc', alpha=0.9, zorder=0)

    # Centromere shaded gap
    ax_macro_cn.axvspan(centromere[0], centromere[1], color='#f1f5f9', hatch='///', alpha=0.8, zorder=1)
    ax_macro_baf.axvspan(centromere[0], centromere[1], color='#f1f5f9', hatch='///', alpha=0.8, zorder=1)

    # Shaded vertical focus strip for target gene
    pad_macro_strip = max(1.5e6, (g_end - g_start) * 2.5)
    ax_macro_cn.axvspan(g_start - pad_macro_strip, g_end + pad_macro_strip, color='#fee2e2', alpha=0.6, zorder=1)
    ax_macro_baf.axvspan(g_start - pad_macro_strip, g_end + pad_macro_strip, color='#fee2e2', alpha=0.6, zorder=1)

    # Plot Rolling Confidence Band for CNV
    cn_upper = c_cn + 1.96 * se_cn
    cn_lower = c_cn - 1.96 * se_cn
    ax_macro_cn.fill_between(x_grid, cn_lower, cn_upper, color='#38bdf8', alpha=0.38, zorder=2,
                             label='Rolling Confidence Band (±1.96 SE, Local MAD / √Neff)')
    ax_macro_cn.plot(x_grid, c_cn, color='#0284c7', linewidth=1.8, linestyle='-', zorder=3, label='Rolling Robust Median CN')

    # Plot PURPLE Segments
    if not df_seg.empty:
        for _, srow in df_seg.iterrows():
            s_s, s_e = srow['start'], srow['end']
            t_cn = srow.get('tumorCopyNumber', srow.get('copyNumber', 0.0))
            if srow.get('depthWindowCount', 1) > 0 and t_cn > 0:
                ax_macro_cn.plot([s_s, s_e], [t_cn, t_cn], color='#0f172a', linewidth=3.2, zorder=4)

    # Plot Probe points
    if not chr_probes.empty:
        p_colors = ['#dc2626' if v < 1.3 else '#f59e0b' if v < 1.7 else '#10b981' if v > 2.5 else '#0284c7' for v in chr_probes['norm_cn']]
        ax_macro_cn.scatter(chr_probes['probe_mid'], chr_probes['norm_cn'], c=p_colors, s=28, edgecolor='#0f172a', linewidth=0.5, zorder=5,
                            label=f'Capture Probes (N={len(chr_probes)})')

    ax_macro_cn.set_ylabel('Copy Number (CN)', fontsize=9.5, fontweight='bold', color='#0f172a')
    plt.setp(ax_macro_cn.get_xticklabels(), visible=False)
    ax_macro_cn.tick_params(axis='x', length=0)
    ax_macro_cn.set_title(f"Whole Chromosome {chrom} Profile to Scale (0–{chr_len/1e6:.1f} Mb) — Multi-Scale Genomic Architecture with Rolling MAD Confidence",
                          fontsize=11.5, fontweight='bold', pad=8, color='#0f172a')
    ax_macro_cn.legend(loc='upper right', frameon=True, facecolor='#ffffff', edgecolor='#cbd5e1', fontsize=7.2, ncol=4)

    # 2. Macro BAF Track
    ax_macro_baf.set_facecolor('#ffffff')
    ax_macro_baf.grid(True, linestyle='--', alpha=0.35, color='#cbd5e1')
    ax_macro_baf.set_ylim(-0.05, 1.05)
    ax_macro_baf.set_yticks([0.0, 0.25, 0.5, 0.75, 1.0])
    ax_macro_baf.axhline(0.5, color='#64748b', linestyle='--', linewidth=1.0, alpha=0.7)
    ax_macro_baf.axhspan(0.0, 0.15, color='#fee2e2', alpha=0.35)
    ax_macro_baf.axhspan(0.85, 1.0, color='#fee2e2', alpha=0.35)

    maj_upper = np.clip(c_baf + 1.96 * se_baf, 0.5, 1.0)
    maj_lower = np.clip(c_baf - 1.96 * se_baf, 0.5, 1.0)
    min_upper = np.clip((1.0 - c_baf) + 1.96 * se_baf, 0.0, 0.5)
    min_lower = np.clip((1.0 - c_baf) - 1.96 * se_baf, 0.0, 0.5)

    ax_macro_baf.fill_between(x_grid, maj_lower, maj_upper, color='#93c5fd', alpha=0.45, zorder=2,
                              label='Major Allele Confidence Band')
    ax_macro_baf.plot(x_grid, c_baf, color='#2563eb', linewidth=1.6, linestyle='-', zorder=3)

    ax_macro_baf.fill_between(x_grid, min_lower, min_upper, color='#d8b4fe', alpha=0.45, zorder=2,
                              label='Minor Allele Confidence Band')
    ax_macro_baf.plot(x_grid, 1.0 - c_baf, color='#9333ea', linewidth=1.6, linestyle='-', zorder=3)

    # Plot PURPLE Segment BAFs
    if not df_seg.empty:
        for _, srow in df_seg.iterrows():
            s_s, s_e = srow['start'], srow['end']
            s_baf = srow.get('tumorBAF', 0.5)
            if srow.get('bafCount', 0) > 0 and s_baf > 0:
                ax_macro_baf.plot([s_s, s_e], [s_baf, s_baf], color='#1e3a8a', linewidth=2.8, zorder=4)
                if abs(s_baf - 0.5) > 0.04:
                    ax_macro_baf.plot([s_s, s_e], [1.0 - s_baf, 1.0 - s_baf], color='#6d28d9', linewidth=2.6, zorder=4)

    # Plot Het SNPs
    if not df_baf.empty:
        maj_snps = df_baf[df_baf['tumorBAF'] >= 0.50]
        min_snps = df_baf[df_baf['tumorBAF'] < 0.50]
        if not maj_snps.empty:
            ax_macro_baf.scatter(maj_snps['position'], maj_snps['tumorBAF'], c='#2563eb', s=32, edgecolor='#1e3a8a',
                                 linewidth=0.5, zorder=5, label=f'Het SNPs (Major, N={len(maj_snps)})')
        if not min_snps.empty:
            ax_macro_baf.scatter(min_snps['position'], min_snps['tumorBAF'], c='#9333ea', s=32, edgecolor='#581c87',
                                 linewidth=0.5, zorder=5, label=f'Het SNPs (Minor, N={len(min_snps)})')

    ax_macro_baf.set_ylabel('Het SNP AF (BAF)', fontsize=9.5, fontweight='bold', color='#0f172a')
    ax_macro_baf.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, pos: f"{x/1e6:.0f} Mb"))
    ax_macro_baf.tick_params(axis='x', labelsize=8.5)
    ax_macro_baf.legend(loc='lower right', frameon=True, facecolor='#ffffff', edgecolor='#cbd5e1', fontsize=7.2, ncol=4)

    # -------------------------------------------------------------
    # BOTTOM SECTION: ZOOMED FOCAL GENE ARCHITECTURE
    # -------------------------------------------------------------
    ax_micro_cn.set_facecolor('#fffdfa')
    ax_micro_cn.grid(True, linestyle='--', alpha=0.35, color='#cbd5e1')
    ax_micro_cn.set_xlim(x_min_f, x_max_f)
    ax_micro_cn.set_ylim(-0.25, 4.4)
    ax_micro_cn.axhline(2.0, color='#64748b', linestyle='--', linewidth=1.1, alpha=0.85)
    ax_micro_cn.axhline(1.0, color='#f59e0b', linestyle=':', linewidth=1.0, alpha=0.8)
    ax_micro_cn.axhline(0.0, color='#ef4444', linestyle=':', linewidth=1.0, alpha=0.8)
    ax_micro_cn.axhspan(1.75, 2.25, color='#f8fafc', alpha=0.9, zorder=0)

    ax_micro_cn.plot([x_min_f, x_max_f], [seg_cn, seg_cn], color='#0f172a', linewidth=3.2, zorder=3,
                     label=f'Focal Segment Total CN ({seg_cn:.2f})')
    ax_micro_cn.plot([x_min_f, x_max_f], [minor_cn, minor_cn], color='#6366f1', linewidth=2.4, linestyle='-.', zorder=3,
                     label=f'Minor Allele CN ({minor_cn:.2f})')

    # Exon shading and exon medians
    for idx_ex, ex in enumerate(exons):
        e_s = ex["start"]
        e_e = ex["end"]
        e_num = ex["exon_num"]
        is_mut = any(pv['exon'] == e_num for pv in parsed_vars)
        vspan = ax_micro_cn.axvspan(e_s, e_e, color='#fee2e2' if is_mut else '#eff6ff', alpha=0.65, zorder=1)
        vspan.set_gid(f"exon_span_{e_num}")

        if not gene_probes.empty and 'probe_mid' in gene_probes.columns:
            sub_p = gene_probes[(gene_probes['probe_mid'] >= e_s - 1200) & (gene_probes['probe_mid'] <= e_e + 1200)]
            if not sub_p.empty:
                e_median = sub_p['norm_cn'].median()
                ax_micro_cn.plot([e_s, e_e], [e_median, e_median], color='#dc2626' if is_mut else '#0284c7', linewidth=2.8, zorder=4)

        y_text = 3.80 if (idx_ex % 2 == 0) else 3.50
        ax_micro_cn.text((e_s + e_e)/2.0, y_text, f"E{e_num}", ha='center', va='bottom',
                         fontsize=7.2, fontweight='bold', color='#dc2626' if is_mut else '#334155',
                         bbox=dict(boxstyle='square,pad=0.1', facecolor='#ffffff', edgecolor='none', alpha=0.75))

    # Focal probes scatter
    if not gene_probes.empty:
        c_colors = ['#dc2626' if v < 1.3 else '#f59e0b' if v < 1.7 else '#10b981' if v > 2.5 else '#0284c7' for v in gene_probes['norm_cn']]
        sc_probes = ax_micro_cn.scatter(gene_probes['probe_mid'], gene_probes['norm_cn'], c=c_colors, s=46, edgecolor='#0f172a', linewidth=0.6, zorder=5)
        sc_probes.set_gid("target_probes")
        ax_micro_cn.plot(gene_probes['probe_mid'], gene_probes['norm_cn'], color='#94a3b8', alpha=0.45, linewidth=1.0, zorder=2)

    # Mutation ploidy (VCN) point and drop line
    for p_idx, pv in enumerate(parsed_vars):
        if pv.get('vcn') is not None:
            mut_vcn = max(0.05, float(pv['vcn']))
            vl = ax_micro_cn.axvline(pv['pos'], color='#dc2626', linestyle=':', linewidth=2.0, alpha=0.95)
            vl.set_gid(f"mut_line_c_{p_idx}")
            ax_micro_cn.plot([pv['pos'], pv['pos']], [seg_cn, mut_vcn], color='#dc2626', linewidth=2.0, zorder=5)
            mut_sc = ax_micro_cn.scatter([pv['pos']], [mut_vcn], c='#dc2626', marker='*', s=260, edgecolor='#881337', linewidth=1.3, zorder=6)
            mut_sc.set_gid(f"mut_point_c_{p_idx}")
            ax_micro_cn.text(pv['pos'], min(4.1, mut_vcn + 0.22),
                             f"{pv['desc']}\n({mut_vcn:.2f} copies, VAF {pv['sample_vaf']:.1%})",
                             ha='center', va='bottom', fontsize=7.2, fontweight='bold', color='#991b1b',
                             bbox=dict(boxstyle='round,pad=0.2', facecolor='#ffffff', edgecolor='#fca5a5', alpha=0.9), zorder=7)

    legend_elements_micro_cn = [
        Line2D([0], [0], color='#0f172a', lw=3.0, label=f'Segment Total CN ({seg_cn:.2f})'),
        Line2D([0], [0], color='#6366f1', lw=2.2, linestyle='-.', label=f'Minor Allele CN ({minor_cn:.2f})'),
        Line2D([0], [0], marker='o', color='w', markerfacecolor='#0284c7', markeredgecolor='#0f172a', markersize=7.5, label='Neutral Probe (1.7–2.5)'),
        Line2D([0], [0], marker='o', color='w', markerfacecolor='#f59e0b', markeredgecolor='#0f172a', markersize=7.5, label='Shallow Loss (1.3–1.7)'),
        Line2D([0], [0], marker='o', color='w', markerfacecolor='#dc2626', markeredgecolor='#0f172a', markersize=7.5, label='Deep Loss (<1.3)'),
        Line2D([0], [0], marker='o', color='w', markerfacecolor='#10b981', markeredgecolor='#0f172a', markersize=7.5, label='Gain (>2.5)'),
        patches.Patch(facecolor='#eff6ff', edgecolor='#93c5fd', label=f'Exon (E1–E{num_exons})')
    ]
    if parsed_vars:
        legend_elements_micro_cn.append(patches.Patch(facecolor='#fee2e2', edgecolor='#fca5a5', label='Mutated Exon'))
        legend_elements_micro_cn.append(Line2D([0], [0], marker='*', color='w', markerfacecolor='#dc2626', markeredgecolor='#881337', markersize=11, label=f"Mutation Ploidy (VCN: {parsed_vars[0]['vcn']:.2f} copies)"))
    ax_micro_cn.legend(handles=legend_elements_micro_cn, loc='upper right', frameon=True, facecolor='#ffffff', edgecolor='#cbd5e1', fontsize=7.2, ncol=4)

    # 2. Focal BAF Track
    ax_micro_baf.set_facecolor('#fffdfa')
    ax_micro_baf.grid(True, linestyle='--', alpha=0.35, color='#cbd5e1')
    ax_micro_baf.set_xlim(x_min_f, x_max_f)
    ax_micro_baf.set_ylim(-0.05, 1.05)
    ax_micro_baf.set_yticks([0.0, 0.25, 0.5, 0.75, 1.0])
    ax_micro_baf.axhline(0.5, color='#64748b', linestyle='--', linewidth=1.0, alpha=0.7)
    ax_micro_baf.axhspan(0.0, 0.15, color='#fee2e2', alpha=0.45)
    ax_micro_baf.axhspan(0.85, 1.0, color='#fee2e2', alpha=0.45)

    focal_baf = pd.DataFrame()
    if not df_baf.empty:
        focal_baf = df_baf[(df_baf['position'] >= x_min_f) & (df_baf['position'] <= x_max_f)].sort_values('position')

    if not focal_baf.empty:
        maj = focal_baf[focal_baf['tumorBAF'] >= 0.50]
        min_b = focal_baf[focal_baf['tumorBAF'] < 0.50]
        if not maj.empty:
            if len(maj) > 1:
                ax_micro_baf.plot(maj['position'], maj['tumorBAF'], color='#93c5fd', alpha=0.70, linewidth=1.2, zorder=2)
            ax_micro_baf.scatter(maj['position'], maj['tumorBAF'], c='#2563eb', s=46, edgecolor='#1e3a8a', linewidth=0.6, zorder=3,
                                 label='Major Allele (AF ≥ 0.50)')
        if not min_b.empty:
            if len(min_b) > 1:
                ax_micro_baf.plot(min_b['position'], min_b['tumorBAF'], color='#d8b4fe', alpha=0.70, linewidth=1.2, zorder=2)
            ax_micro_baf.scatter(min_b['position'], min_b['tumorBAF'], c='#9333ea', s=46, edgecolor='#581c87', linewidth=0.6, zorder=3,
                                 label='Minor Allele (AF < 0.50)')

    # Segment allele fit lines
    if abs(frac_major - 0.5) > 0.03:
        ax_micro_baf.plot([x_min_f, x_max_f], [frac_major, frac_major], color='#6366f1', linewidth=2.8, label=f'Major Phase ({frac_major:.0%})')
        ax_micro_baf.plot([x_min_f, x_max_f], [frac_minor, frac_minor], color='#a855f7', linewidth=2.8, label=f'Minor Phase ({frac_minor:.0%})')
    else:
        ax_micro_baf.plot([x_min_f, x_max_f], [0.5, 0.5], color='#1d4ed8', linewidth=2.8, label='Diploid Balance (0.50)')

    # Somatic mutation VAF points on BAF track
    for p_idx, pv in enumerate(parsed_vars):
        v_sample_af = pv['sample_vaf']
        v_tumor_af = pv['tumor_af']
        ax_micro_baf.plot([pv['pos'], pv['pos']], [v_sample_af, v_tumor_af], color='#dc2626', linestyle='--', linewidth=1.4, zorder=5)
        dia_sc = ax_micro_baf.scatter([pv['pos']], [v_sample_af], c='#f59e0b', marker='D', s=160, edgecolor='#78350f', linewidth=1.3, zorder=6,
                                     label=f"Sample VAF ({v_sample_af:.1%})" if p_idx == 0 else "")
        dia_sc.set_gid(f"somatic_sample_vaf_diamond_{p_idx}")
        star_sc = ax_micro_baf.scatter([pv['pos']], [v_tumor_af], c='#dc2626', marker='*', s=280, edgecolor='#881337', linewidth=1.3, zorder=7,
                                      label=f"Tumor AF ({v_tumor_af:.1%})" if p_idx == 0 else "")
        star_sc.set_gid(f"somatic_tumor_af_star_{p_idx}")

    # Exon shading on BAF track
    for ex in exons:
        e_s, e_e = ex["start"], ex["end"]
        is_mut = any(pv['exon'] == ex['exon_num'] for pv in parsed_vars)
        ax_micro_baf.axvspan(e_s, e_e, color='#fee2e2' if is_mut else '#eff6ff', alpha=0.55, zorder=0)

    # 3. Focal Gene Structure & Quantitative Methylation Track
    ax_micro_gene.set_facecolor('#ffffff')
    ax_micro_gene.set_xlim(x_min_f, x_max_f)
    ax_micro_gene.set_ylim(-0.95, 1.30)
    ax_micro_gene.spines['top'].set_visible(False)
    ax_micro_gene.spines['right'].set_visible(False)
    ax_micro_gene.set_yticks([0.0, 0.5, 1.0])
    ax_micro_gene.set_yticklabels(['0.0', '0.5', '1.0'], fontsize=8.5, fontweight='bold', color='#334155')
    ax_micro_gene.grid(axis='y', linestyle=':', color='#cbd5e1', alpha=0.6)

    for ex in exons:
        e_s, e_e = ex["start"], ex["end"]
        is_mut = any(pv['exon'] == ex['exon_num'] for pv in parsed_vars)
        ax_micro_gene.axvspan(e_s, e_e, color='#fee2e2' if is_mut else '#eff6ff', alpha=0.55, zorder=0)

    # Intron backbone
    ax_micro_gene.plot([g_start, g_end], [-0.35, -0.35], color='#475569', linewidth=2.4, zorder=2)
    arrow_xs = np.linspace(g_start + (g_end - g_start)*0.08, g_end - (g_end - g_start)*0.08, 7)
    for ax_pos in arrow_xs:
        ax_micro_gene.text(ax_pos, -0.35, "▶" if strand == 1 else "◀", ha='center', va='center',
                           fontsize=7.5, color='#64748b', zorder=3)

    # Exon boxes
    for ex in exons:
        e_s, e_e = ex["start"], ex["end"]
        e_num = ex["exon_num"]
        w = max(e_e - e_s, (x_max_f - x_min_f) * 0.0075)
        is_mut = any(pv['exon'] == e_num for pv in parsed_vars)
        rect = patches.Rectangle((e_s, -0.55), w, 0.40, facecolor='#dc2626' if is_mut else '#2563eb',
                                 edgecolor='#991b1b' if is_mut else '#1e3a8a', linewidth=1.1, zorder=4)
        ax_micro_gene.add_patch(rect)
        y_lbl = -0.66 if (e_num % 2 == 1) else -0.78
        ax_micro_gene.text(e_s + w/2.0, y_lbl, f"E{e_num}", ha='center', va='top', fontsize=7.2,
                           fontweight='bold', color='#991b1b' if is_mut else '#334155', zorder=5)

    p_w = max(500, promoter_end - promoter_start)
    prom_bg = '#fee2e2' if not is_unmeth else '#dcfce7'
    prom_border = '#dc2626' if not is_unmeth else '#16a34a'
    rect_prom = patches.Rectangle((promoter_start, -0.58), p_w, 0.46, facecolor=prom_bg, edgecolor=prom_border,
                                  linewidth=1.6, linestyle='--', zorder=4)
    rect_prom.set_gid("promoter_box")
    ax_micro_gene.add_patch(rect_prom)
    prom_beta_str = f"β={tumor_promoter_beta:.2f}" if tumor_promoter_beta is not None else ""
    prom_lbl = f"Promoter\n({methylation_status}\n{prom_beta_str})" if not is_unmeth and tumor_promoter_beta is not None else f"Promoter\n({methylation_status})"
    ax_micro_gene.text(promoter_start + p_w/2.0, -0.35, prom_lbl, ha='center', va='center',
                       fontsize=7.0, fontweight='bold', color=prom_border, zorder=6)

    ax_micro_gene.text(promoter_start - 800, -0.35, "5'", ha='right', va='center', fontsize=8.5, fontweight='bold', color='#334155')
    ax_micro_gene.text(g_end + 800, -0.35, f"3' ({transcript_id})", ha='left', va='center', fontsize=8.2, fontweight='bold', color='#334155')

    # CpG methylation points and mean line
    ax_micro_gene.axhline(0.20, color='#94a3b8', linestyle=':', linewidth=0.8, alpha=0.5, zorder=1)
    ax_micro_gene.axhline(0.50, color='#94a3b8', linestyle='--', linewidth=0.8, alpha=0.5, zorder=1)
    if tumor_promoter_beta is not None:
        ax_micro_gene.plot([promoter_start, promoter_end], [tumor_promoter_beta, tumor_promoter_beta],
                           color='#dc2626' if not is_unmeth else '#059669', linewidth=3.2, zorder=7)

    if cpg_list:
        c_pos = [c['pos'] for c in cpg_list]
        c_tumor = [c['tumor'] for c in cpg_list]
        cpg_sc = ax_micro_gene.scatter(c_pos, c_tumor, c='#10b981', s=45, edgecolor='#047857', linewidth=0.6, zorder=8,
                                       label=f"Tumour CpG (Mean β={tumor_promoter_beta:.3f})" if tumor_promoter_beta is not None else "Tumour CpG")
        cpg_sc.set_gid("tumor_prom_cpgs")
        if any(c.get('normal') is not None for c in cpg_list):
            c_norm = [c['normal'] for c in cpg_list if c.get('normal') is not None]
            c_npos = [c['pos'] for c in cpg_list if c.get('normal') is not None]
            norm_sc = ax_micro_gene.scatter(c_npos, c_norm, c='#94a3b8', s=32, edgecolor='#475569', linewidth=0.5, alpha=0.85, zorder=7,
                                            label="Normal Baseline CpG")
            norm_sc.set_gid("normal_prom_cpgs")

    # Somatic mutation lollipop pin
    for p_idx, pv in enumerate(parsed_vars):
        ax_micro_gene.plot([pv['pos'], pv['pos']], [-0.15, 0.45], color='#dc2626', linewidth=2.2, zorder=9)
        pin_sc = ax_micro_gene.scatter([pv['pos']], [0.45], c='#dc2626', marker='o', s=90, edgecolor='#881337', linewidth=1.3, zorder=10)
        pin_sc.set_gid(f"mut_pin_{p_idx}")
        ax_micro_gene.text(pv['pos'], 0.52, f"{pv['desc']}\nExon {pv['exon']} (VAF: {pv['sample_vaf']:.1%})",
                           ha='center', va='bottom', fontsize=7.5, fontweight='bold', color='#991b1b',
                           bbox=dict(boxstyle='round,pad=0.25', facecolor='#fff1f2', edgecolor='#f43f5e', alpha=0.95), zorder=11)

    # Format x-axis in Mb
    for ax in [ax_micro_cn, ax_micro_baf]:
        plt.setp(ax.get_xticklabels(), visible=False)
        ax.tick_params(axis='x', length=0)
    ax_micro_gene.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, pos: f"{x/1e6:.3f}M"))
    ax_micro_gene.tick_params(axis='x', labelsize=8, rotation=25)

    ax_micro_cn.set_ylabel('Copy Number (CN)', fontsize=9.5, fontweight='bold', color='#0f172a')
    ax_micro_baf.set_ylabel('Het SNP AF (BAF)', fontsize=9.5, fontweight='bold', color='#0f172a')
    ax_micro_gene.set_ylabel('Methylation (β)', fontsize=9.5, fontweight='bold', color='#0f172a')

    ax_micro_cn.set_title(f"Zoomed Focal Architecture: {gene} (chr{chrom}:{x_min_f/1e6:.3f}–{x_max_f/1e6:.3f} Mb) — Exon Structure, Allelic Trajectory & Quantitative Methylation",
                          fontsize=11.5, fontweight='bold', pad=8, color='#991b1b')

    # Zoom connector polygon linking macro tracks to micro tracks
    fig.canvas.draw()
    trans_macro = ax_macro_baf.transData
    trans_fig = fig.transFigure.inverted()

    p_top_l = trans_fig.transform(trans_macro.transform((g_start - pad_macro_strip, -0.05)))
    p_top_r = trans_fig.transform(trans_macro.transform((g_end + pad_macro_strip, -0.05)))

    bbox_micro = ax_micro_cn.get_position()
    p_bot_l = (bbox_micro.x0, bbox_micro.y1 + 0.015)
    p_bot_r = (bbox_micro.x1, bbox_micro.y1 + 0.015)

    poly = patches.Polygon([p_top_l, p_top_r, p_bot_r, p_bot_l],
                           closed=True, transform=fig.transFigure,
                           facecolor='#fee2e2', edgecolor='#fca5a5', linestyle='--', linewidth=1.2, alpha=0.45, zorder=0)
    fig.patches.append(poly)

    fig.suptitle(f"Multi-Scale Genomic, Allelic & Epigenetic Architecture: {gene} on Chromosome {chrom} ({sample_id})",
                 fontsize=15.0, fontweight='bold', y=0.99)

    # Save PNG and SVG via local temp files
    import tempfile
    import shutil

    os.makedirs(os.path.dirname(os.path.abspath(output_png)), exist_ok=True)
    if output_html is None:
        output_html = output_png.rsplit('.', 1)[0] + '.html'
    os.makedirs(os.path.dirname(os.path.abspath(output_html)), exist_ok=True)

    with tempfile.NamedTemporaryFile(suffix='.png', delete=False) as tf_png, \
         tempfile.NamedTemporaryFile(suffix='.svg', delete=False) as tf_svg:
        temp_png_path = tf_png.name
        temp_svg_path = tf_svg.name

    try:
        plt.savefig(temp_png_path, dpi=220, bbox_inches='tight')
        shutil.copyfile(temp_png_path, output_png)
        print(f"Saved clean PNG: {output_png}")

        plt.savefig(temp_svg_path, format='svg', bbox_inches='tight')
        plt.close(fig)

        with open(temp_svg_path, 'r', encoding='utf-8') as f:
            svg_content = f.read()
    finally:
        if os.path.exists(temp_png_path):
            os.remove(temp_png_path)
        if os.path.exists(temp_svg_path):
            os.remove(temp_svg_path)

    mut_tooltips_json = {}
    for p_idx, pv in enumerate(parsed_vars):
        if is_balanced_diploid:
            phasing_html = f"""
            <div style='margin-top:6px; padding:4px 6px; background:#f0fdf4; border-left:3px solid #16a34a; border-radius:3px;'>
                <strong>Allelic State:</strong> <strong>Heterozygous (1 copy on Balanced Diploid Background)</strong><br>
                • Segment Copy Numbers: Total = {seg_cn:.2f} | Minor = {minor_cn:.2f} (Diploid 1:1 baseline)<br>
                • State: Monoallelic somatic mutation (no copy asymmetry; distinguishing major/minor allele is non-applicable)
            </div>
            """
        else:
            badge_color = "#dc2626" if "Major" in pv['phased_allele'] else "#d97706"
            bg_color = "#fef2f2" if "Major" in pv['phased_allele'] else "#fffbeb"
            phasing_html = f"""
            <div style='margin-top:6px; padding:4px 6px; background:{bg_color}; border-left:3px solid {badge_color}; border-radius:3px;'>
                <strong>Allele Phasing:</strong> Mutation is phased to the <strong>{pv['phased_allele'].upper()}</strong><br>
                • Segment Copy Numbers: Total = {seg_cn:.2f} | Major = {major_cn:.2f} | Minor = {minor_cn:.2f}<br>
                • Phasing Call: {pv['phasing_call']}
            </div>
            """

        mut_tooltips_json[f"mut_{p_idx}"] = f"""
        <div style='font-size:12.5px; line-height:1.4;'>
            <div style='font-weight:700; color:#991b1b; font-size:13.5px; border-bottom:1px solid #fecdd3; padding-bottom:4px; margin-bottom:6px;'>
                ⭐ Somatic Mutation: {gene} {pv['desc']}
            </div>
            <div><strong>Locus:</strong> chr{chrom}:{pv['pos']:,} (Exon {pv['exon']})</div>
            <div style='margin-top:4px;'><strong>Measured Sample VAF:</strong> <span style='color:#b45309; font-weight:700;'>{pv['sample_vaf']:.1%}</span> (Tumor Purity = {purity_clamped:.0%})</div>
            <div><strong>Purity-Adjusted Tumor AF:</strong> <span style='color:#dc2626; font-weight:700;'>{pv['tumor_af']:.1%}</span></div>
            <div><strong>Variant Ploidy (VCN):</strong> {pv['vcn']:.2f} copies</div>
            {phasing_html}
        </div>
        """

    if has_meth and tumor_promoter_beta is not None:
        norm_str = f"{normal_promoter_beta:.3f} ({normal_promoter_beta*100:.1f}%)" if normal_promoter_beta is not None else "N/A (No Normal EM-Seq)"
        delta_str = f"{tumor_promoter_beta - (normal_promoter_beta or 0.0):+.3f}" if normal_promoter_beta is not None else "N/A"
        meth_tooltip_html = f"""
        <div style='font-size:12.5px; line-height:1.4;'>
            <div style='font-weight:700; color:#047857; font-size:13.5px; border-bottom:1px solid #a7f3d0; padding-bottom:4px; margin-bottom:6px;'>
                🌿 Promoter CpG Island: {gene}
            </div>
            <div><strong>Promoter Coordinates:</strong> chr{chrom}:{min(promoter_start, promoter_end):,} – {max(promoter_start, promoter_end):,}</div>
            <div style='margin-top:4px;'><strong>Status:</strong> <span style='font-weight:700; color:{"#15803d" if is_unmeth else "#b91c1c"};'>{methylation_status}</span></div>
            <div><strong>Tumour Promoter Mean Beta (β):</strong> <span style='color:#059669; font-weight:700;'>{tumor_promoter_beta:.3f} ({tumor_promoter_beta*100:.1f}%)</span></div>
            <div><strong>Normal Baseline Mean Beta (β):</strong> <span style='color:#475569; font-weight:700;'>{norm_str}</span></div>
            <div><strong>Epigenetic Delta (Δβ):</strong> {delta_str}</div>
            <div style='margin-top:5px; font-size:11px; color:#475569;'>
                • Closed circles: Promoter CpGs | Open circles: Gene body CpGs<br>
                • Green: Tumour sample | Grey: Normal reference baseline
            </div>
        </div>
        """
    else:
        meth_tooltip_html = f"""
        <div style='font-size:12.5px; line-height:1.4;'>
            <div style='font-weight:700; color:#475569; font-size:13.5px; border-bottom:1px solid #cbd5e1; padding-bottom:4px; margin-bottom:6px;'>
                Promoter CpG Island: {gene}
            </div>
            <div><strong>Promoter Coordinates:</strong> chr{chrom}:{min(promoter_start, promoter_end):,} – {max(promoter_start, promoter_end):,}</div>
            <div style='margin-top:4px;'><strong>Status:</strong> <span style='font-weight:700; color:#64748b;'>NOT ASSAYED / NO DATA</span></div>
            <div style='margin-top:4px; color:#64748b; font-style:italic;'>EM-Seq DNA methylation sequencing was not performed or failed QC for this sample.</div>
        </div>
        """

    probe_tooltip_html = f"""
    <div style='font-size:12px;'>
        <div style='font-weight:700; color:#0369a1;'>Target Capture Probes ({gene})</div>
        <div>Obs/Exp Normalized Ratio × 2.0</div>
        <div>Blue: Neutral (1.7–2.5) | Amber: Loss (1.3–1.7) | Red: Deep Loss (&lt;1.3) | Green: Gain (&gt;2.5)</div>
    </div>
    """

    header_mut_badge = ""
    if parsed_vars:
        header_mut_badge = f"<span class='badge badge-mut'>Somatic: {parsed_vars[0]['desc']}</span>"

    html_page = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>{gene} ({sample_id}) - Genomic, Allelic & Methylation Architecture</title>
    <style>
        body {{
            font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, Helvetica, Arial, sans-serif;
            margin: 0;
            padding: 16px 24px;
            background: #f8fafc;
            color: #0f172a;
        }}
        .header-bar {{
            display: flex;
            justify-content: space-between;
            align-items: center;
            background: #ffffff;
            padding: 12px 20px;
            border-radius: 8px;
            box-shadow: 0 1px 3px rgba(0,0,0,0.08);
            margin-bottom: 16px;
        }}
        .header-title {{
            font-size: 18px;
            font-weight: 700;
            color: #0f172a;
        }}
        .header-badges {{
            display: flex;
            gap: 10px;
        }}
        .badge {{
            padding: 4px 10px;
            border-radius: 6px;
            font-size: 12px;
            font-weight: 600;
        }}
        .badge-sample {{ background: #eff6ff; color: #1d4ed8; border: 1px solid #bfdbfe; }}
        .badge-mut {{ background: #fef2f2; color: #b91c1c; border: 1px solid #fecdd3; }}
        .badge-meth {{ background: {"#f0fdf4; color: #15803d; border: 1px solid #bbf7d0;" if is_unmeth else "#fef2f2; color: #b91c1c; border: 1px solid #fecdd3;"} }}
        .plot-container {{
            background: #ffffff;
            border-radius: 8px;
            box-shadow: 0 1px 4px rgba(0,0,0,0.08);
            padding: 16px;
            position: relative;
            overflow: auto;
        }}
        svg {{
            width: 100%;
            height: auto;
            min-width: 1100px;
            display: block;
        }}
        #interactive-tooltip {{
            position: absolute;
            display: none;
            background: #ffffff;
            color: #1e293b;
            padding: 12px 14px;
            border-radius: 8px;
            box-shadow: 0 10px 25px -5px rgba(0,0,0,0.18), 0 8px 10px -6px rgba(0,0,0,0.12);
            border: 1px solid #cbd5e1;
            pointer-events: none;
            z-index: 10000;
            max-width: 360px;
            transition: opacity 0.12s ease-out;
        }}
        [id^="somatic_tumor_af_star"], [id^="somatic_sample_vaf_diamond"], [id^="mut_point_c"], [id^="mut_line_c"], #promoter_box, #tumor_prom_cpgs, #normal_prom_cpgs, #target_probes {{
            cursor: pointer;
            transition: all 0.15s ease;
        }}
        [id^="somatic_tumor_af_star"]:hover, [id^="mut_point_c"]:hover {{
            filter: drop-shadow(0 0 6px #dc2626);
        }}
        #promoter_box:hover {{
            filter: brightness(0.95) drop-shadow(0 0 4px #10b981);
        }}
        .instruction-note {{
            font-size: 12px;
            color: #64748b;
            text-align: right;
            margin-top: 8px;
        }}
    </style>
</head>
<body>
    <div class="header-bar">
        <div class="header-title">{gene} Genomic, Allelic & Methylation Profile</div>
        <div class="header-badges">
            <span class="badge badge-sample">Sample: {sample_id}</span>
            {header_mut_badge}
            <span class="badge badge-meth">Promoter: {methylation_status}</span>
        </div>
    </div>

    <div class="plot-container" id="plot-wrapper">
        <div id="interactive-tooltip"></div>
        {svg_content}
    </div>
    <div class="instruction-note">💡 Hover over the somatic mutation star/diamond, promoter CpG island, or probe points to inspect detailed metrics.</div>

    <script>
        const tooltip = document.getElementById('interactive-tooltip');
        const wrapper = document.getElementById('plot-wrapper');

        const mutPayloads = {json.dumps(mut_tooltips_json)};
        const methPayload = `{meth_tooltip_html}`;
        const probePayload = `{probe_tooltip_html}`;

        function showTooltip(e, html) {{
            tooltip.innerHTML = html;
            tooltip.style.display = 'block';
            updateTooltipPos(e);
        }}

        function hideTooltip() {{
            tooltip.style.display = 'none';
        }}

        function updateTooltipPos(e) {{
            const rect = wrapper.getBoundingClientRect();
            let x = e.clientX - rect.left + 16;
            let y = e.clientY - rect.top + 16;
            if (x + 360 > rect.width) {{
                x = e.clientX - rect.left - 370;
            }}
            tooltip.style.left = x + 'px';
            tooltip.style.top = y + 'px';
        }}

        wrapper.addEventListener('mousemove', (e) => {{
            if (tooltip.style.display === 'block') {{
                updateTooltipPos(e);
            }}
        }});

        function attachListeners() {{
            Object.keys(mutPayloads).forEach(k => {{
                const p_idx = k.replace('mut_', '');
                const sel = `[id$="_star_${{p_idx}}"], [id$="_diamond_${{p_idx}}"], [id$="_point_c_${{p_idx}}"], [id$="_line_c_${{p_idx}}"], [id$="_pin_${{p_idx}}"]`;
                document.querySelectorAll(sel).forEach(el => {{
                    el.addEventListener('mouseenter', (e) => showTooltip(e, mutPayloads[k]));
                    el.addEventListener('mouseleave', hideTooltip);
                }});
            }});

            document.querySelectorAll('#promoter_box, #tumor_prom_cpgs, #normal_prom_cpgs, #tumor_body_cpgs, #normal_body_cpgs').forEach(el => {{
                el.addEventListener('mouseenter', (e) => showTooltip(e, methPayload));
                el.addEventListener('mouseleave', hideTooltip);
            }});

            document.querySelectorAll('#target_probes').forEach(el => {{
                el.addEventListener('mouseenter', (e) => showTooltip(e, probePayload));
                el.addEventListener('mouseleave', hideTooltip);
            }});
        }}

        attachListeners();
    </script>
</body>
</html>
"""

    # Write HTML via local temp file and atomic copy
    with tempfile.NamedTemporaryFile(mode='w', suffix='.html', delete=False, encoding='utf-8') as tf_html:
        temp_html_path = tf_html.name
        tf_html.write(html_page)

    try:
        shutil.copyfile(temp_html_path, output_html)
        print(f"Saved interactive HTML: {output_html}")
    finally:
        if os.path.exists(temp_html_path):
            os.remove(temp_html_path)


def main():
    parser = argparse.ArgumentParser(description='Generate a gene-focused BAF & Copy Number plot from AMBER BAF')
    parser.add_argument('--amber-baf', required=True, help='Path to AMBER BAF tsv.gz file')
    parser.add_argument('--purple-segments', help='Path to purple.segment.tsv file')
    parser.add_argument('--target-region-cn', help='Path to purple.target_region_cn.tsv file')
    parser.add_argument('--ichor-dir', help='Path to ichorCNA output directory')
    parser.add_argument('--sample-id', help='Sample ID to match')
    parser.add_argument('--patient-id', help='Patient ID to match')
    parser.add_argument('--variants-json', help='JSON string or file path containing variants to annotate')
    parser.add_argument('--gene', required=True, help='Gene of interest')
    parser.add_argument('--output', required=True, help='Path to output PNG plot')
    parser.add_argument('--output-html', help='Path to output HTML plot (optional)')
    parser.add_argument('--genome', default='GRCh38_hmf', help='Genome build (GRCh38_hmf or GRCh37_hmf)')
    parser.add_argument('--min-depth', type=int, default=100, help='Minimum depth for SNPs')
    parser.add_argument('--purity', type=float, default=1.0, help='Tumor purity (0.0 to 1.0)')
    parser.add_argument('--methylation-json', help='JSON string containing CpG-level promoter methylation data')
    
    args = parser.parse_args()

    gene = args.gene
    is_grch37 = "37" in args.genome
    structures = GENES_STRUCTURE_37 if is_grch37 else GENES_STRUCTURE_38

    if gene not in structures:
        print(f"Error: Gene {gene} is not in the offline database.", file=sys.stderr)
        sys.exit(1)

    gene_info = structures[gene]

    variants = []
    if args.variants_json:
        try:
            if os.path.isfile(args.variants_json):
                with open(args.variants_json, 'r', encoding='utf-8') as f:
                    variants = json.load(f)
            else:
                variants = json.loads(args.variants_json)
        except Exception as e:
            print(f"Warning: Failed to parse variants JSON: {e}", file=sys.stderr)

    methylation_data = None
    if args.methylation_json:
        try:
            if os.path.isfile(args.methylation_json):
                with open(args.methylation_json, 'r', encoding='utf-8') as f:
                    methylation_data = json.load(f)
            else:
                methylation_data = json.loads(args.methylation_json)
        except Exception as e:
            print(f"Warning: Failed to parse methylation JSON: {e}", file=sys.stderr)

    sample_id = args.sample_id or "Sample"
    if not args.sample_id and args.output:
        base = os.path.basename(args.output)
        if '_' in base:
            sample_id = base.split('_')[0]

    render_gene_baf_v10(
        gene=gene,
        gene_info=gene_info,
        amber_baf_file=args.amber_baf,
        purple_segments=args.purple_segments,
        target_cn_file=args.target_region_cn,
        variants=variants,
        methylation_data=methylation_data,
        purity=args.purity,
        sample_id=sample_id,
        genome=args.genome,
        output_png=args.output,
        output_html=args.output_html
    )

if __name__ == '__main__':
    main()
