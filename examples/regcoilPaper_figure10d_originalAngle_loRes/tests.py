#!/usr/bin/env python

# This python script checks the output file for an example to 
# see if the results are close to expected values.  This script may be
# run directly, and it is also called when "make test" is run from the
# main REGCOIL directory.

execfile('../testsCommon.py')
absoluteTolerance = 1e-100

numFailures = 0

f = readOutputFile()

variableName = 'lambda'
data = f.variables[variableName][()]
relativeTolerance = 1e-12
numFailures += arrayShouldBe(data, [0, 1e-15, 1.33352143216332e-15, 1.77827941003892e-15,\
                                    2.37137370566166e-15, 3.16227766016838e-15, 4.21696503428582e-15,\
                                    5.62341325190349e-15, 7.49894209332456e-15, 1e-14],relativeTolerance,absoluteTolerance)


variableName = 'chi2_B'
data = f.variables[variableName][()]
relativeTolerance = 0.03
# Skip the lambda=0 case, which is pathological, and make sure the other ones are within a few % of the high-res results:
numFailures += arrayShouldBe(data[1:], [0.174519878306313, 0.233148551834137, \
    0.310446213587783, 0.411565852326815, 0.542616084841207, \
    0.710568315783369, 0.922890554110523, 1.18673864681126, 1.50770781013225],relativeTolerance,absoluteTolerance)

variableName = 'chi2_K'
data = f.variables[variableName][()]
relativeTolerance = 0.01
numFailures += arrayShouldBe(data[1:], [1.74957088182873e+15, 1.6989658516062e+15,\
    1.64892520707378e+15, 1.59982527129183e+15, 1.55209554042516e+15,\
    1.50621123410416e+15, 1.46269693947801e+15, 1.42212837760326e+15,\
    1.38509938973785e+15], relativeTolerance,absoluteTolerance)

variableName = 'max_Bnormal'
data = f.variables[variableName][()]
relativeTolerance = 0.03
numFailures += arrayShouldBe(data[1:], [0.117555544291034, 0.134561224876993,\
    0.153923792979921, 0.176347144447289, 0.201603899107085,\
    0.229977435728317, 0.261839981656698, 0.296625161466877, 0.333988791941713], relativeTolerance,absoluteTolerance)

variableName = 'max_K'
data = f.variables[variableName][()]
relativeTolerance = 0.06
numFailures += arrayShouldBe(data[1:], [ 8824481.60102707, 8355450.23653953,\
    7884212.39105115, 7412969.17310808, 6945281.57578389, 6485916.75531679,\
    6040435.92861007, 5614638.75966927, 5214018.48872153], relativeTolerance,absoluteTolerance)

variableName = 'single_valued_current_potential_mn'
data = f.variables[variableName][()]
#print data.shape
relativeTolerance = 1e-4
numFailures += arrayShouldBe(data[0,:], [-636481.60416242, -102835.333542226, -39124.0186294881, -6816.08113940242,
   542.692475698608, -2686.1108530638, -601.765971600708, 26.2180240803913,
    -150.847194134617, -58.2551976061389, 60.0317419077943, -15.997364372303,
    -12.3528609642142, -28.8468306445688, 120.907353313668, -90.444227136533,
    -254.610877079356, 117.750619995718, -586.859801429981,
    -3240.79714512101, -345.706420298433, -8000.93888813538,
    -33571.3281108956, -146463.897042633, 261525.559042215, 875162.411027093,
    223424.734276641, 48679.3704112507, 24653.0844662479, 3207.55972139988,
    -971.300457622961, 1735.36333359508, 392.289700069504, 58.3403278403076,
    69.2088816044866, 28.574277764148, -32.6049569184015, -3.19468817585156,
    18.9509927995784, 40.9417395029467, -191.329948990297, 102.505974883377,
    476.918019850139, -313.219021953551, 228.539785025719, 3535.38858987503,
    1994.33005106025, -2723.95384286104, 25378.8241944312, 211577.172631438,
    -509987.705552534, -305882.370067181, -64140.5967994057, 5526.489569249,
    -21143.6598694639, -1312.30942574156, 1256.40679971238,
    -866.753793091458, -323.218082280049, -62.4712285310379,
    -22.3887872985142, 11.572463519096, 18.8670999860767, 6.09601174794048,
    -39.626829330563, -58.1376718144676, 258.637725082145, -93.2264049514627,
    -690.771352763538, 137.091910375363, 618.390992293105, -2858.58032050993,
    -8887.84722613869, 18273.2756559291, 31254.1463773911, -183820.187349618,
    317263.663711709, 136546.037104119, 1714.17033452626, -6965.82579245848,
    19700.2927688041, 279.360825551664, -585.394321551508, 457.56013983205,
    247.414984963762, 36.063011937723, -2.88433594557794, -14.657242495135,
    -29.7277484867182, -5.15155109726196, 79.9609755678001, 96.2382594231763,
    -307.126611682774, -217.681183192655, 1164.78374124385, 1189.74972998963,
    -2200.7785485338, -901.659579742522, 18143.0098542936, -12823.6723296253,
    -18793.9927181867, 191936.052698242, -134572.41665727, -112479.974965958,
    -19248.1496384529, 8028.67003637305, -18917.5018103108, -663.00856123082,
    -130.269218689209, -359.598076210987, -121.339449458628,
    -42.6774334572464, 17.9474379704981, 23.4796365098314, 39.903936366997,
    6.71424592396383, -163.619373671154, -164.958633554062, 457.9492928255,
    929.023995244022, -2170.03690293255, -4109.54155806639, 4423.02245595864,
    11848.7759951353, -19183.4447881742, 20698.488290232, 98376.4092161681,
    -246652.945255868, 36177.2383322875, 120007.953328824, 11484.7657544551,
    -4712.15687254718, 17600.9544860059, -109.759176704224, 858.005478267692,
    250.801901773195, 112.595088965925, -29.0802679750806, -41.9900107268828,
    -32.446943536193, -50.3151682883662, -5.27400182515745, 314.549242682154,
    268.959151540374, -921.339292653633, -2110.06565879334, 4471.04983495344,
    8373.18809218046, -8063.05480853702, -24107.6728123965, 13542.8083194455,
    -49083.7639169206, -130918.200746624, 265974.857549389, 19277.1970752313,
    -94636.1096069281, -12455.6660070208, 205.614622738596,
    -13352.0580607509, 88.0140555713081, -751.856883151767, -231.62025244395,
    24.2787971083689, 79.1321103557964, 68.1790577884637, 33.0664705300337,
    62.1791798540125, -8.96468981509356, -551.338416161408   ],relativeTolerance,absoluteTolerance,requireSameLength=False)

numFailures += arrayShouldBe(data[1,:], [-636845.139192448, -93726.100985817, -12808.2679337137, -6436.11331882515, 
    -1681.81086973519, -485.628613896459, -51.0035530084606, 
    2.57965700650486, -3.08359685930898, -0.95839259239159, 
    -0.0542860537675785, 0.0589544960724215, 0.0184659230256529, 
    0.0769552670319529, -0.262669749613766, -1.85762870238826, 
    -4.48294147405658, 9.77123536044076, -28.4339720480356, 
    -450.405437429645, -2048.68587258934, -9542.84769319848, 
    -11669.3480566478, -126993.914047055, 200134.524174578, 780678.271122884, 
    221364.046785469, 48874.9084250086, 1732.76398357889, 2614.62551563965, 
    799.077029570748, 219.002026009396, 29.0188384278004, -3.1024899442833, 
    0.314721570012263, 0.0893881351479318, -0.0812851357532177, 
    0.00290851130928144, -0.0550572574676785, -0.147780529776928, 
    0.38376209245763, 2.36663288136853, 4.74859467921797, -19.3144909584734, 
    -76.8150486658404, -48.6780930746745, 956.722007151166, 3224.97803831513, 
    14791.0744397381, 157700.163165199, -487877.544270743, -165831.030675654, 
    -52754.7805491044, -3106.22754847172, 473.627446170829, 
    -561.542644462128, -117.003129564884, -25.9990546377287, 
    -0.621612316021711, 4.52433375939335, 0.78366288299259, 
    0.111023309228341, -0.00758576015460735, -0.00643956505715865, 
    0.154816558534507, 0.401147312418584, 0.0352085229246653, 
    -4.24674093374629, -8.32654759807846, 49.1679592736964, 177.523100278182, 
    418.431282522965, -1919.9181912971, 8164.06507384419, 12742.9035731632, 
    -125313.597804208, 251605.800949885, -26875.6370303012, 
    -6551.74001875138, 152.613961229345, -606.027206935648, 
    -222.877823486467, -62.1808469036173, -18.5690428992164, 
    -7.12477659577854, -2.97106388609079, -0.382819851517509, 
    -0.00503060124001372, 0.0161359845510441, -0.00998783802528895, 
    -0.368609249523516, -1.11292554363902, -2.01207279148648, 
    7.33366678816116, 14.9415662122255, -66.4299343505126, -65.7382483985067, 
    -499.849037166508, 5598.6008557347, -3111.17077838504, 8389.66498705812, 
    87851.854449707, -34120.1030577796, 38812.9747478082, -2556.15730402831, 
    -430.042936935821, 319.725970749222, 137.093492800773, 33.7420385155216, 
    10.2930579648087, 3.30062400243539, 0.861064214735285, 
    0.00367925130425598, 0.0139232158257035, -0.0222371785709498, 
    -0.0300175983684826, 0.645683067183313, 1.90363938028169, 
    6.55063199568604, -0.0326354551966688, -16.2130112965381, 
    -50.275102393623, -123.525433819699, 3257.47329534058, -8305.30343341391, 
    -7971.01752571554, 39625.692610642, -73762.1713333272, -48092.6536727011, 
    -11355.3822208888, -510.262304707652, 479.718899459551, 43.5788931587508, 
    -16.5903117683302, -7.70867912158858, 1.36885139892951, 
    0.0990160316681461, 0.00267920733865999, -0.00685727034461888, 
    -0.0283406794140528, 0.0356207864623813, 0.174783727405189, 
    -0.748386349705818, -0.799894357903513, -12.6643380901238, 
    -19.6013780623049, -14.5281290410233, 221.471071053636, 770.289884864596, 
    -6076.24785297526, -1618.96547832442, 5993.3173828301, -44505.6636867058,    ],relativeTolerance,absoluteTolerance,requireSameLength=False)


del data
f.close()
print "numFailures:",numFailures
exit(numFailures > 0)
