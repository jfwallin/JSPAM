var stateInfos = [['-1.01385,-2.25736,-1.07197,-0.09418,-0.41519,-1.10167,0.38946,0.60577,0.69777,0.52425,142.48649,65.29434,258.15299,238.65171,0.3,0.3,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0,0,1.35749,-2.81999,0.66757,0.0,0.0,1.0,0.0,0.0,0.0','-6.696014176130678','4.068901123074229','SDSS 1237678620102623480'],
['-9.93853,-4.5805,15.43348,-1.54025,-2.97899,4.69234,34.41671,29.49607,2.97913,4.2501,301.0,35.5,310.80745,321.98757,0.3,0.3,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0,0,0.24008,-2.41999,5.62715,0.0,1.0,0.0,1.0,0.0,0.0','-6.904942981293307','24.945037352311584','SDSS 587722984435351614'],
['1.08119,-1.66798,-2.06826,-0.18438,-0.39693,-0.02919,0.51654,0.36192,0.94043,0.43877,12.97415,135.11902,97.54035,305.01225,0.3,0.3,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0,0,0.68139,-4.67999,1.02298,0.0,0.0,1.0,0.0,0.0,0.0','-9.639999999999883','2.8158787054319507','SDSS 587724234257137777'],
['-0.49344,-1.43805,0.05985,0.07349,-1.02296,0.89921,1.09707,1.54176,1.27761,0.33429,48.78911,43.8006,50.06943,135.05795,0.3,0.3,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0,0,5.08621,-0.9,0.43029,0.0,0.0,0.0,0.0,0.0,0.0','-5.0','3.6092302720033733','SDSS 587726033843585146'],
['0.93974,0.17991,2.66529,0.0,-0.42418,1.16041,0.36689,1.02587,0.65803,0.61394,6.2096,95.0,200.0,289.23076,0.3,0.3,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0,0,0.73231,-1.98,1.00208,0.0,1.0,1.0,0.0,0.0,0.0','-5.454348810387462','2.244202687779696','SDSS 587727177926508595'],
['-0.40793,-1.92518,1.95772,-1.14918,-1.02439,0.77963,1.3042,1.25831,0.94075,0.54373,97.78523,60.52239,144.60117,216.50117,0.3,0.3,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0,0,0.77689,-1.44,1.17399,0.0,0.0,1.0,0.0,0.0,0.0','-5.307844817363995','3.6642833366547523','SDSS 587727178988388373'],
['0.66665,-2.12701,1.61528,-0.18579,-0.23147,0.20647,0.19195,0.1555,0.81482,0.75573,15.0,87.5,203.47826,22.36024,0.3,0.3,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0,0,0.35895,-5.65999,1.10062,0.0,1.0,1.0,1.0,0.0,0.0','-5.0','3.362365819667262','SDSS 587727222471131318'],
['0.05958,-0.64055,-0.62252,0.19343,-0.15746,-0.31657,0.87551,0.47456,0.09603,0.05894,78.08622,87.64218,36.69664,347.76946,0.3,0.3,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0,0,93.92675,-0.9,0.07628,0.0,0.0,1.0,0.0,0.0,0.0','-4.159999999999999','0.9533878782563739','SDSS 587728676861051075'],
['-0.12623,0.35028,3.35743,-0.12979,0.55971,4.02072,0.78258,0.98072,0.30541,0.41196,152.77952,195.62801,147.22542,138.42519,0.3,0.3,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0,0,45.47487,-0.49999,0.07399,0.0,0.0,0.0,0.0,0.0,0.0','-1.3485536705527896','1.770658089034511','SDSS 587729227151704160'],
['-1.1109,1.42179,1.69164,-0.49869,0.50953,0.08717,0.016,0.3198,0.33392,0.38043,51.18749,-35.0,295.15527,249.31677,0.3,0.3,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0,0,0.31636,-2.83999,1.01393,0.0,1.0,0.0,0.0,0.0,0.0','-6.187313530679593','2.1223294285207097','SDSS 587731913110650988'],
['-2.87284,0.90246,0.0622,-0.33826,0.18874,-0.12684,0.11758,0.10218,0.53431,0.4208,135.70152,354.86088,140.82002,138.13649,0.3,0.3,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0,0,0.6753,-5.73999,0.62153,0.0,0.0,0.0,0.0,0.0,0.0','-5.0','3.8505110069557897','SDSS 587732136993882121'],
['-0.10526,0.14243,-0.20682,-2.00321,-0.17917,0.6021,0.79131,1.11639,0.47156,0.1715,117.875,147.04046,326.45962,65.75701,0.3,0.3,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0,0,12.38433,0.0,0.26993,0.0,1.0,1.0,0.0,0.0,0.0','-3.1400000000000023','4.1739074186920515','SDSS 587732772130652231'],
['0.08788,-0.63517,1.09776,-0.66456,-0.47795,0.39667,0.20617,1.16277,0.27916,0.33541,14.93183,68.71883,103.69424,239.11166,0.3,0.3,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0,0,0.89859,-2.83999,1.16114,0.0,1.0,0.0,1.0,0.0,0.0','-5.0','1.6523463700624257','SDSS 587733080814583863'],
['-1.23326,1.43243,1.17391,-0.24002,-0.52906,0.74803,2.30689,0.43804,1.01348,0.78231,140.375,-4.875,317.51552,221.36645,0.3,0.3,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0,0,1.42528,-3.03999,0.99198,0.0,1.0,0.0,2.0,0.0,0.0','-5.0','3.921951431940543','SDSS 587734862680752822'],
['-0.48646,-1.45004,-0.62585,0.40688,-0.53238,0.04629,1.48573,0.74812,0.39338,0.19293,26.64819,49.48416,342.60021,232.84096,0.3,0.3,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0,0,5.02108,-1.52,0.42493,0.0,0.0,0.0,0.0,0.0,0.0','-9.77999999999988','2.394917346758852','SDSS 587735043609329845'],
['-0.8126,0.26807,1.42974,-0.61636,0.10013,0.34708,0.92793,0.18718,0.49658,0.14763,58.81685,94.3125,235.38311,335.40372,0.3,0.3,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0,0,14.6183,-1.28,0.18054,0.0,1.0,0.0,0.0,0.0,0.0','-10.23999999999987','1.429735488107589','SDSS 587735665840881790'],
['1.29969,-0.31467,3.70378,1.33823,0.91009,-0.88972,0.55893,0.17255,1.18348,0.84617,126.88209,51.3787,296.58787,232.55479,0.3,0.3,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0,0,0.01868,0.0,3.9662,0.0,0.0,1.0,0.0,0.0,0.0','-6.409059757572907','3.3279442892451963','SDSS 587736523764334706'],
['-0.27341,0.18994,0.74436,1.28828,1.04149,0.49963,0.17637,0.11239,0.30307,0.3086,96.61117,49.83149,152.65312,109.12104,0.3,0.3,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0,0,0.25949,0.0,0.8198,0.0,1.0,1.0,1.0,0.0,0.0','-1.8091978836128806','0.8753428381605564','SDSS 587736941981466667'],
['-4.62932,-3.48666,8.53665,-1.02385,-0.94978,2.5994,1.40808,0.97714,1.1802,0.65139,176.93944,71.08069,337.85444,335.44237,0.3,0.3,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0,0,1.4927,-3.55999,0.65576,0.0,1.0,0.0,0.0,0.0,0.0','-6.622481217085435','6.9853132699771905','SDSS 587738569246376675'],
['-0.82909,2.18149,-1.22409,0.08683,-0.22468,-0.92517,0.91339,2.31295,0.73445,1.10831,153.56774,43.83877,168.63224,326.12543,0.3,0.3,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0,0,0.63048,-3.03999,1.67104,0.0,1.0,0.0,0.0,0.0,0.0','-4.839999999999985','4.6856489298212605','SDSS 587738569249390718'],
['0.47569,-0.71987,-2.73633,0.15406,-0.31429,-0.25394,0.74212,0.35877,0.47819,0.22201,55.45,45.65,319.41176,335.29411,0.3,0.3,0.0,0.0,0.0,0.0,0.0,0.0,2.0,0,0,1.6323,-4.87999,0.64674,0.0,1.0,1.0,1.0,0.0,0.0','-8.019999999999918','1.6236927448351341','SDSS 587739153356095531'],
['-0.30275,-0.54877,-0.03638,-0.4393,0.47019,0.0,0.15827,0.23379,0.24804,0.34462,96.07505,107.19568,227.47677,229.43557,0.3,0.3,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0,0,2.04922,-3.81999,0.49257,0.0,1.0,0.0,1.0,0.0,0.0','-5.0','2.179981536023529','SDSS 587739407868690486'],
['-0.89489,0.7566,-2.61146,-0.83865,-0.36327,-1.79472,0.51513,0.27306,0.30949,0.7403,7.91644,42.80034,115.71619,304.2272,0.3,0.3,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0,0,0.26036,-1.32,1.18627,0.0,0.0,0.0,0.0,0.0,0.0','-2.678171154437894','1.8378469459476634','SDSS 587739505541578866'],
['-1.9621,2.22832,-0.63225,0.00118,0.11473,0.27316,0.97107,0.71787,1.09937,0.43672,82.03236,271.02679,191.07553,296.26974,0.3,0.3,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0,0,2.21449,-5.23999,0.60942,0.0,0.0,0.0,0.0,0.0,0.0','-7.719999999999923','5.048231960279983','SDSS 587739646743412797'],
['-0.30712,0.91123,-4.13962,-0.72369,-0.3113,-0.14825,0.42912,0.66235,0.57938,0.5465,70.38163,178.66037,59.54056,109.24137,0.3,0.3,0.0,0.0,0.0,0.0,0.0,0.0,2.0,0,0,0.0781,-1.4,4.15674,0.0,1.0,1.0,1.0,0.0,0.0','-8.563936785041877','1.7279746764071704','SDSS 587739647284805725'],
['-0.4413,-1.22532,-2.06561,-0.25131,-0.19632,-0.30839,1.82612,1.30459,0.58709,0.53904,213.73768,67.3073,203.49577,241.01371,0.3,0.3,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0,0,4.1093,-3.71999,0.51102,0.0,0.0,1.0,0.0,0.0,0.0','-3.5000000000000027','5.154879003376974','SDSS 587739707420967061'],
['2.37024,-1.36863,0.25617,0.43873,-0.34901,-0.15217,0.54922,0.13328,0.73772,0.42516,148.30063,49.41001,231.49951,261.44764,0.3,0.3,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0,0,2.70624,-3.31999,0.39999,0.0,0.0,0.0,0.0,0.0,0.0','-5.0','3.8355954659543916','SDSS 587739720308818095'],
['0.18431,-1.70765,4.96754,-1.15744,-1.80231,0.57509,0.5629,0.20912,0.94419,0.63482,6.239,36.32914,283.48382,345.33778,0.3,0.3,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0,0,0.01389,-0.59999,5.08991,0.0,0.0,1.0,1.0,0.0,0.0','-5.418063128971349','3.2598635430193523','SDSS 587739721376202860'],
['-0.43119,-0.92661,-1.79015,0.02476,0.3569,-0.01856,0.3307,0.55382,0.71235,0.39054,7.45496,57.96958,85.87695,212.74889,0.3,0.3,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0,0,0.86861,-5.33999,0.8924,0.0,0.0,0.0,0.0,0.0,0.0','-7.39999999999993','4.697271459647266','SDSS 587739721900163101'],
['1.11881,-0.58701,1.7904,0.04865,-0.14019,1.01095,0.13149,0.09535,0.61225,0.17248,6.80196,90.49367,224.25723,228.89454,0.3,0.3,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0,0,2.01801,-4.69999,0.34829,0.0,0.0,0.0,0.0,0.0,0.0','-7.203187541653982','1.487879649503476','SDSS 587739810496708646'],
['0.96821,-2.63044,-2.22943,0.16074,-0.78512,-0.84691,0.51707,0.59847,0.43977,0.2968,45.10324,358.74357,62.33573,207.08113,0.3,0.3,0.0,0.0,0.0,0.0,0.0,0.0,3.0,0,0,2.07794,-2.39999,0.49411,0.0,0.0,1.0,2.0,0.0,0.0','-4.335352069594705','4.143636596149168','SDSS 587739845393580192'],
['-1.92706,-1.6735,-1.02462,0.11886,-0.93364,0.06848,2.11334,0.35556,1.29023,0.1899,173.28798,45.57512,44.62654,230.34231,0.3,0.3,0.0,0.0,0.0,0.0,0.0,0.0,2.0,0,0,0.73028,-2.95999,1.42974,0.0,0.0,1.0,1.0,0.0,0.0','-7.659999999999925','2.8849602237178806','SDSS 587741391565422775'],
['-0.10214,-2.12459,5.33525,0.21877,-0.02317,1.56857,0.59529,2.82108,0.61541,0.79173,-34.375,141.75,306.3354,0.0,0.3,0.3,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0,0,0.6252,-3.15999,1.4541,0.0,1.0,1.0,0.0,0.0,0.0','-5.6214129945742055','4.312806123287474','SDSS 587741532784361481'],
['-0.03825,1.80871,0.30038,-0.2417,0.18815,0.26631,0.09611,0.12302,0.54242,0.24137,53.17979,146.64794,46.33769,24.75554,0.3,0.3,0.0,0.0,0.0,0.0,0.0,0.0,2.0,0,0,0.17788,-2.55999,1.41533,0.0,0.0,1.0,1.0,0.0,0.0','-5.0','2.797773003905253','SDSS 587741534400217110'],
['2.96778,3.61902,-0.12375,0.30078,0.67945,0.07129,0.33782,0.18678,0.70035,0.46849,45.8125,166.875,286.21118,234.7826,0.3,0.3,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0,0,0.60809,-5.45999,0.82469,0.0,1.0,1.0,0.0,0.0,0.0','-11.009531000819964','4.763339101259693','SDSS 587741602030026825'],
['-0.42567,-0.83668,-1.93196,-0.4365,0.00717,-1.09404,1.33245,0.1993,0.59839,0.1963,108.92933,39.06386,319.92442,160.84304,0.3,0.3,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0,0,2.98555,-2.77999,0.50916,0.0,0.0,0.0,0.0,0.0,0.0','-4.046843729837199','2.1198043043269137','SDSS 587741722819493915'],
['1.06338,1.04559,1.63578,1.00283,0.31386,1.55365,0.56343,0.23131,0.9187,0.50148,108.04444,30.31226,327.91147,146.78275,0.3,0.3,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0,0,1.10371,-1.1,0.56561,0.0,1.0,1.0,1.0,0.0,0.0','-5.1015176494969365','3.8707124709669136','SDSS 587741817851674654'],
['-1.15174,-0.32733,-1.2717,-1.85991,0.14873,-0.01861,1.13467,1.49592,0.68248,1.11228,9.1875,78.0,160.99378,252.6708,0.3,0.3,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0,0,1.04768,-0.78,1.01893,0.0,1.0,0.0,2.0,0.0,0.0','-3.665174450750971','2.788528271038801','SDSS 587741829658181698'],
['-0.98742,-0.65926,-2.46557,0.12428,-0.48102,-0.12119,0.51217,0.30955,0.80325,0.31048,149.1875,110.6875,78.26086,70.43478,0.3,0.3,0.0,0.0,0.0,0.0,0.0,0.0,3.0,0,0,0.30397,-4.97999,1.64905,0.0,1.0,1.0,1.0,0.0,0.0','-9.099999999999895','1.7004384780019914','SDSS 587742010583941189'],
['-0.23056,1.79463,-2.58944,-0.24992,0.26734,-0.62973,0.1458,0.16789,0.53863,0.36106,116.68822,65.78378,32.04906,249.07924,0.3,0.3,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0,0,0.40524,-3.37999,0.8248,0.0,0.0,1.0,0.0,0.0,0.0','-8.749211157819898','3.190470818792097','SDSS 587742014353702970'],
['0.60209,-0.36672,-3.1872,-0.01698,0.20511,-1.60752,0.73935,1.22722,0.49155,0.54816,101.875,53.0625,2.23602,216.8944,0.3,0.3,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0,0,1.466,-1.36,0.67055,0.0,1.0,0.0,1.0,0.0,0.0','-3.533949208942654','2.8024626673950612','SDSS 587742571610243080'],
['-0.54962,0.12782,0.12509,0.50274,-0.83732,0.28701,1.08587,0.45008,0.52331,0.34502,86.86455,55.83273,348.72428,148.42085,0.3,0.3,0.0,0.0,0.0,0.0,0.0,0.0,3.0,0,0,11.06724,-1.64,0.26106,0.0,1.0,1.0,1.0,0.0,0.0','-5.0','3.272215785934329','SDSS 587745402001817662'],
['-0.20706,-0.12072,-0.51456,-1.17676,-0.37393,-1.19254,0.48543,0.86206,0.10418,0.11542,96.0,65.6875,127.45341,221.36645,0.3,0.3,0.0,0.0,0.0,0.0,0.0,0.0,2.0,0,0,25.66484,-0.26,0.14193,0.0,1.0,1.0,3.0,0.0,0.0','-5.0','0.4093291686888595','SDSS 587746029596311590'],
['1.20338,0.67312,0.64081,0.43046,0.50439,0.43521,0.50371,0.60483,0.58896,0.86039,72.5,48.8125,111.80124,67.08074,0.3,0.3,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0,0,19.59198,-1.22,0.15549,0.0,1.0,0.0,1.0,0.0,0.0','-4.1000000000000005','2.9133436375228223','SDSS 587747120521216156'],
['-0.26053,-1.07287,1.61671,0.31334,-0.12918,0.49763,0.41938,0.14077,0.48995,0.51658,45.22657,32.37807,236.54081,61.01804,0.3,0.3,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0,0,0.79193,-2.67999,0.80776,0.0,0.0,0.0,0.0,0.0,0.0','-11.519999999999843','2.139856664761071','SDSS 588007005230530750'],
['0.96515,0.50909,0.06211,1.04124,-1.16798,0.25165,1.53689,1.4,0.48536,0.50568,117.0625,116.8125,131.92546,228.07453,0.3,0.3,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0,0,1.93691,-0.72,0.89473,0.0,1.0,0.0,0.0,0.0,0.0','-1.8400000000000012','2.715167447594645','SDSS 588011124116422756'],
['0.57769,-0.37316,-1.10901,1.16477,0.09758,-1.35063,1.24805,0.59457,0.19797,0.23287,85.43445,97.5088,317.60473,4.75135,0.3,0.3,0.0,0.0,0.0,0.0,0.0,0.0,3.0,0,0,5.00857,-0.65999,0.37311,0.0,1.0,1.0,1.0,0.0,0.0','-1.2917952266792656','1.83713078600622','SDSS 588013383816904792'],
['0.77006,-0.18571,-2.39132,0.22907,-0.3688,0.12392,0.69896,0.30098,0.71739,0.34935,51.4375,136.0625,199.00621,243.7267,0.3,0.3,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0,0,0.55906,-6.86,1.26156,0.0,1.0,1.0,0.0,0.0,0.0','-9.81999999999988','1.6804039671144984','SDSS 588017604696408086'],
['-1.65653,-2.3295,-1.43985,-0.3518,-0.50926,-0.03246,0.4397,0.2497,0.30783,0.3846,141.0625,110.25,138.63354,26.83229,0.3,0.3,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0,0,5.17158,-3.51999,0.27189,0.0,1.0,0.0,2.0,0.0,0.0','-6.043090755868784','3.7863714576827427','SDSS 588017604696408195'],
['0.77539,-1.76333,-3.66139,0.4448,0.03362,-0.82359,10.87707,0.36014,0.46316,0.45767,83.53667,329.68415,310.01058,317.60661,0.3,0.3,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0,0,18.78104,-2.05999,0.29848,0.0,0.0,0.0,0.0,0.0,0.0','-5.0','4.728701363248437','SDSS 588017702948962343'],
['0.50703,1.54971,0.15061,-0.14252,0.26477,0.15418,0.46797,0.11446,1.14919,0.27229,11.0625,65.5625,111.80124,290.68322,0.3,0.3,0.0,0.0,0.0,0.0,0.0,0.0,2.0,0,0,6.16516,-2.63999,0.25323,0.0,1.0,1.0,1.0,0.0,0.0','-6.359999999999952','2.7058056331443594','SDSS 588017978901528612'],
['-0.10914,-1.06691,-1.77066,0.25931,0.18152,0.03363,0.38702,0.4096,0.56671,0.46368,73.18749,-63.6875,346.58385,187.82608,0.3,0.3,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0,0,0.78828,-5.21999,0.92797,0.0,1.0,1.0,0.0,0.0,0.0','-7.979999999999918','3.8654346341427477','SDSS 588018055130710322'],
['-0.63148,0.63778,0.26845,-1.00831,-0.19458,0.42192,0.90726,0.95484,0.29162,0.27274,97.25,-4.75,194.53416,35.77639,0.3,0.3,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0,0,6.14846,-0.72,0.3755,0.0,1.0,1.0,3.0,0.0,0.0','-2.740000000000002','2.176974099437279','SDSS 758874299603222717'],
['0.61603,-2.03504,-1.99154,0.24821,-0.75856,-0.00241,0.60656,0.47041,0.91437,0.5201,149.9,136.1875,50.4,38.01242,0.3,0.3,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0,0,0.44147,-3.17999,1.35853,0.0,1.0,0.0,0.0,0.0,0.0','-5.0','3.9137884481791456','SDSS 758877153600208945'],
['0.82889,0.04034,-0.62191,1.71911,0.608,-0.20669,0.04024,0.10375,0.66818,0.88701,84.24978,24.27873,79.47009,227.46727,0.3,0.3,0.0,0.0,0.0,0.0,0.0,0.0,0.81481,0,0,0.17757,-0.42,0.66069,0.0,0.0,1.0,0.0,0.0,0.0','-4.011681038327637','2.347310086225794','HST Arp148'],
['1.46799,-1.90751,-0.55072,0.07554,-0.38851,-0.60193,0.61882,0.16625,0.77494,0.88675,72.6875,118.0,116.27329,245.96273,0.3,0.3,0.0,0.0,0.0,0.0,0.0,0.0,0.75,0,0,0.30059,-2.59999,1.56147,0.0,1.0,0.0,0.0,0.0,0.0','-5.0','4.500688042506512','HST CGCG'],
['-0.07728,-0.80372,-0.13938,-0.31746,-0.75224,-0.08412,0.33607,0.44353,0.91692,0.89357,79.55898,35.13127,141.02645,306.89238,0.3,0.3,0.0,0.0,0.0,0.0,0.0,0.0,0.92857,0,0,20.24723,-0.67999,0.14396,0.0,0.0,0.0,0.0,0.0,0.0','-5.0','5.45603392190466','HST Arp244'],
['2.69455,-0.04067,0.02039,0.44033,-0.03982,0.20649,0.22587,0.12288,2.52756,1.21528,339.26402,62.125,320.75111,181.11801,0.3,0.3,0.0,0.0,0.0,0.0,0.0,0.0,0.78571,0,0,56.95151,-3.85999,0.06188,0.0,1.0,0.0,2.0,0.0,0.0','-5.0','5.206083269292581','HST Arp272'],
['-0.10257,-4.89429,-0.45646,-0.25191,-0.99373,-0.80271,0.70144,0.35002,2.86951,2.36963,95.61553,-46.0,233.87594,82.8,0.3,0.3,0.0,0.0,0.0,0.0,0.0,0.0,0.96296,0,0,0.44674,-4.03999,1.20341,0.0,1.0,1.0,0.0,0.0,0.0','-5.0','7.414712606387326','HST Arp273'],
['0.14477,-0.08121,-0.12185,0.44786,0.37667,0.34975,0.25506,0.19451,0.0943,0.07058,83.375,-86.6875,309.68944,259.37888,0.3,0.3,0.0,0.0,0.0,0.0,0.0,0.0,0.71428,0,0,12.34862,-1.98,0.1894,0.0,1.0,0.0,0.0,0.0,0.0','-5.0','0.9039706159318648','HST ESO77-14'],
['0.55851,1.6516,1.33066,0.75369,0.91969,0.75846,0.35894,0.22324,1.62444,0.91594,118.16774,13.18068,108.74559,124.37474,0.3,0.3,0.0,0.0,0.0,0.0,0.0,0.0,0.65517,0,0,0.54217,-1.56,0.83259,0.0,0.0,1.0,1.0,0.0,0.0','-10.196992911457555','5.968129701154468','HST NGC5331'],
['-1.76016,-0.40892,-4.48686,-0.2624,-0.29988,-0.53499,0.6195,0.25405,0.75273,0.45032,88.9,91.3,334.8,0.0,0.3,0.3,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0,0,0.94579,-5.23999,0.76974,0.0,1.0,1.0,0.0,0.0,0.0','-11.160599280811024','3.034361684208167','HST NGC6786']];
